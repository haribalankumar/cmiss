      SUBROUTINE CALC_SPARSE_SOLVE_1DTREE(ISCMAX,ISRMAX,ISC,ISR,M,N,NBJ,
     &  NEELEM,NENP,NHST,NPNE,NPNODE,nr,nx,NXI,NYNE,NYNP,NZMAX,NZTOT,
     &  WORK_ISC,WORK_ISR,ERROR,*) 

C#### Subroutine: CALC_SPARSE_SOLVE_1DTREE
C###  Description:
C###    CALC_SPARSE_SOLVE_1DTREE calculates sparsity patterns for the
C###    matrices used in 1D tree problems.
C***    Only uses compressed row format.  For multiple nh, nh=1 will be
C***    stored in the first half of nzz, and nh=2 in second etc.
C***  Created by : Merryn Howatson Tawhai, October 1997

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
!     Parameter List
      INTEGER ISCMAX,ISRMAX,ISC(ISCMAX),ISR(ISRMAX),M,N,NBJ(NJM,NEM),
     &  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NHST,
     &  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,nx,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NZMAX,NZTOT,WORK_ISC(M*8),
     &  WORK_ISR(M*8)
      CHARACTER ERROR*(*)
C      INTEGER*1 WORK_ARRAY(N,M)
!     Local Variables
      INTEGER i,IDATA(10),na,nb,ndb,NDG,ne,ne2,nh,nhs,nhx,nn,noelem,
     '  noelem2,nonode,np,NPLIST2(0:NPM),
     '  np1,np2,ny,ny1,ny2,ny2c,ny3,nz,nzz,nzz_row,
     '  nz_count,NZ_MAX,nz_work
      LOGICAL FLOW_BALANCED
      CHARACTER CHAR*1
!     EXTERNAL FUNCTION
      INTEGER IDIGITS

      CALL ENTERS('CALC_SPARSE_SOLVE_1DTREE',*9999)
      CALL ASSERT(USE_SPARSE.EQ.1,
     '  '>>Set USE_SPARSE to 1 to use sparsity arrays',ERROR,*9999)
C***  Calculate sparsity pattern
      IF(ISRMAX.LT.(N+1)) THEN
        WRITE(CHAR,'(I1)') IDIGITS(N+1)
        WRITE(ERROR,'(''>>Increase NISR_GKM to '',I'//CHAR//')')
     '    (N+1)
        GO TO 9999
      ENDIF
C      CALL ASSERT((N+1).LE.ISRMAX,'>>Increase ISRMAX',ERROR,*9999)
C***  Initalise sparsity arrays
      DO i=1,M+1
        ISR(i)=1
      ENDDO !i
      DO nz=1,ISCMAX
        ISC(nz)=0
      ENDDO !nz
      nzz=1
      IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5)THEN
        IF(LUNG_SPARSITY)THEN !use direct assembly of ISC and ISR
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            IF(NXI(-1,0,ne).EQ.0)THEN !stem branch: only two in row
              np1=NPNE(1,nb,ne) !start node #
              DO nhx=1,NHST !easiest to loop over nh here?
                nh=NH_LOC(nhx,nx)
                ny1=NYNP(1,1,nh,np1,0,1,nr) !row #
                DO nhs=1,2 !for each row entry
                  IF(nzz.LE.ISCMAX) THEN
                    np2=NPNE(nhs,nb,ne)
                    ny2=NYNP(1,1,nh,np2,0,1,nr) !col #
                    ISC(nzz)=ny2 !store the col #
                  ENDIF
                  nzz=nzz+1
                ENDDO !nhs
                IF(nh.EQ.1.AND.NHST.EQ.2)THEN
                  DO nhs=1,2 !for each row entry
                    IF(nzz.LE.ISCMAX) THEN
                      np2=NPNE(nhs,nb,ne)
                      ny2c=NYNP(1,1,2,np2,0,1,nr) !col #
                      ISC(nzz)=ny2c !store the col #
                    ENDIF
                    nzz=nzz+1
                  ENDDO !nhs
                ENDIF
                DO na=1,NAT(nb)
                  ny3=NYNE(na,nh,1,1,ne)
                  ISC(nzz)=ny3
                  nzz=nzz+1
                ENDDO !na for auxiliary variables i.e. flow
                ISR(ny1+1)=nzz !store the row #
              ENDDO
            ENDIF
            np1=NPNE(2,nb,ne) !end node # of parent
            DO nhx=1,NHST
              nh=NH_LOC(nhx,nx)
              ny1=NYNP(1,1,nh,np1,0,1,nr) !row #
              DO nhs=1,2 !for each row entry
                IF(nzz.LE.ISCMAX) THEN
                  np2=NPNE(nhs,nb,ne)
                  ny2=NYNP(1,1,nh,np2,0,1,nr) !col #
                  ISC(nzz)=ny2 !store the col #
                ENDIF
                nzz=nzz+1
              ENDDO !nhs
              IF(nh.EQ.1.AND.NHST.EQ.2)THEN
                DO nhs=1,2 !for each row entry
                  IF(nzz.LE.ISCMAX) THEN
                    np2=NPNE(nhs,nb,ne)
                    ny2c=NYNP(1,1,2,np2,0,1,nr) !col #
                    ISC(nzz)=ny2c !store the col #
                  ENDIF
                  nzz=nzz+1
                ENDDO !nhs
              ENDIF
              DO ndb=1,NXI(1,0,ne) !for each daughter branch
                ne2=NXI(1,ndb,ne) !daughter element #
                np2=NPNE(2,nb,ne2) !end node #
                IF(nzz.LE.ISCMAX) THEN
                  ny2=NYNP(1,1,nh,np2,0,1,nr) !col #
                  ISC(nzz)=ny2
                ENDIF
                nzz=nzz+1
              ENDDO !ndb
              IF(nh.EQ.1.AND.NHST.EQ.2)THEN
                DO ndb=1,NXI(1,0,ne) !for each row entry
                  ne2=NXI(1,ndb,ne) !daughter element #
                  np2=NPNE(2,nb,ne2) !end node #
                  IF(nzz.LE.ISCMAX) THEN
                    ny2c=NYNP(1,1,2,np2,0,1,nr) !col #
                    ISC(nzz)=ny2c !store the col #
                  ENDIF
                  nzz=nzz+1
                ENDDO !ndb
              ENDIF
              ISR(ny1+1)=nzz !store the row #
            ENDDO !nh
          ENDDO !noelem
        ELSE !use 'work' arrays
          nz_work=0
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nhx=1,NHST
              nh=NH_LOC(nhx,nx)
              ny=NYNP(1,1,nh,np,0,1,nr)
              nz_work=nz_work+1
              WORK_ISC(nz_work)=ny !store diagonal non-zeros
              WORK_ISR(nz_work)=ny
            ENDDO
          ENDDO
          NDG=nz_work
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            DO nhx=1,NHST
              nh=NH_LOC(nhx,nx)
              ny1=NYNP(1,1,nh,np1,0,1,nr)
              ny2=NYNP(1,1,nh,np2,0,1,nr)
              nz_work=nz_work+1
              WORK_ISC(nz_work)=ny1
              WORK_ISR(nz_work)=ny2
              nz_work=nz_work+1
              WORK_ISC(nz_work)=ny2
              WORK_ISR(nz_work)=ny1
            ENDDO
          ENDDO
          NZ_MAX=nz_work
          ISR(1)=1
          nzz=0
          DO nz_work=1,NDG
            nz_count=0
            nz=1
            DO WHILE(nz_count.LT.4.AND.nz.LT.NZ_MAX+1)
              IF(WORK_ISR(nz).EQ.nz_work)THEN
                nz_count=nz_count+1
                IDATA(nz_count)=WORK_ISC(nz)
              ENDIF
              nz=nz+1
            ENDDO
            CALL ISORT(nz_count,IDATA) !reorder into ascending order
            DO nz=1,nz_count
              nzz=nzz+1
              ISC(nzz)=IDATA(nz)
            ENDDO
            ISR(nz_work+1)=nzz+1
          ENDDO
        ENDIF
      ELSE IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.ITYP3(nr,
     &    nx).EQ.6)THEN !capillaries
        nzz_row=1 !ISR index
        NPLIST2(0)=0
        DO nonode=1,NPM
          NPLIST2(nonode)=0 !initialise
        ENDDO !nonode
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          np1=NPNE(1,nb,ne) !start node #
          nh=NH_LOC(1,nx)
          ny1=NYNP(1,1,nh,np1,0,1,nr) !row #
          DO nhs=1,2 !for each row entry i.e. pressures
            IF(nzz.LE.ISCMAX) THEN
              np2=NPNE(nhs,nb,ne)
              ny2=NYNP(1,1,nh,np2,0,1,nr) !col #
              ISC(nzz)=ny2 !store the col #
            ENDIF
            nzz=nzz+1
          ENDDO !nhs
          DO na=1,NAT(nb) !auxiliary variable - flow
            IF(nzz.LE.ISCMAX) THEN
              ny3=NYNE(na,nh,1,1,ne)
              ISC(nzz)=ny3
            ENDIF
            nzz=nzz+1
          ENDDO !na
          nzz_row=nzz_row+1
          CALL ASSERT(nzz_row.LT.ISRMAX,
     '      '>>Increase ISR_GK in .ippara file',ERROR,*9999)
          ISR(nzz_row)=nzz !store the row #
          DO nn=1,NNT(nb)
            FLOW_BALANCED=.FALSE. !initialise
            np=NPNE(nn,nb,ne) !balances each node of element
            DO nonode=1,NPLIST2(0)
              IF(np.EQ.NPLIST2(nonode)) THEN
                FLOW_BALANCED=.TRUE. !flow balance already done at np
              ENDIF
            ENDDO !nonode
            IF((NENP(np,0,nr).GT.1).AND.(.NOT.FLOW_BALANCED))THEN
              DO na=1,NAT(nb) !balance flow at junction
                DO noelem2=1,NENP(np,0,nr) !only for ne with subtended
                  IF(nzz.LE.ISCMAX) THEN
                    ne2=NENP(np,noelem2,nr) !branches, not outlet branch
                    ny3=NYNE(na,nh,1,1,ne2)
                    ISC(nzz)=ny3
                  ENDIF
                  nzz=nzz+1
                ENDDO !noelem2
              ENDDO !na
              NPLIST2(0)=NPLIST2(0)+1 !ensures balance only done
              NPLIST2(NPLIST2(0))=np !once at each junction node
              nzz_row=nzz_row+1
              CALL ASSERT(nzz_row.LT.ISRMAX,
     '          '>>Increase ISR_GK in .ippara file',ERROR,*9999)
              ISR(nzz_row)=nzz !store the row #
            ENDIF !np.EQ.NPLIST2(nonode)
          ENDDO !nn
        ENDDO !noelem
        nzz=nzz-1 !total number of non-zero entries
      ENDIF !(ITYP3)
      NZTOT=nzz !total # of non-zero entries
C      WRITE(OP_STRING,'(/'' ISCMAX,NZMAX= '',I8)') NZTOT
C      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      IF(ISCMAX.LT.NZTOT) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NZTOT)
        WRITE(ERROR,'(''>>Increase NISC_GKM to '',I'//CHAR//')')
     '    NZTOT
        GO TO 9999
      ELSE IF(NZMAX.LT.NZTOT) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NZTOT)
        WRITE(ERROR,'(''>>Increase NZ_GK_M to '',I'//CHAR//')')
     '    NZTOT
        GO TO 9999
      ENDIF
C      CALL ASSERT(NZTOT.LE.ISCMAX,'>>Increase ISCMAX',ERROR,*9999)
C      CALL ASSERT(NZTOT.LE.NZMAX,'>>Increase NZMAX',ERROR,*9999)
      
      CALL EXITS('CALC_SPARSE_SOLVE_1DTREE')
      RETURN
 9999 CALL ERRORS('CALC_SPARSE_SOLVE_1DTREE',ERROR)
      CALL EXITS('CALC_SPARSE_SOLVE_1DTREE')
      RETURN 1
      END



