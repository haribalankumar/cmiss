      SUBROUTINE CALC_SPARSE_GKK_1DTREE(ISCMAX,ISRMAX,ISC_GK,ISC,ISR_GK,
     '  ISR,M,N,NBJ,NDG,NEELEM,NENP,NHST,NPLIST2,NPNE,nr,nx,NXI,NYNE,
     &  NYNP,NYNR,NZTOT,SPARSENESS,FIX,ERROR,*)

C#### Subroutine: CALC_SPARSE_GKK_1DTREE
C###  Description:
C###    CALC_SPARSE_GKK_1DTREE calculates sparsity patterns for global
C###    solution matrix GKK for 1D tree problems.
C***  Created by : Merryn Howatson Tawhai, October 1997


      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER ISCMAX,ISRMAX,ISC_GK(NISC_GKM),ISC(ISCMAX),
     &  ISR_GK(NISR_GKM),ISR(ISRMAX),M,N,NBJ(NJM,NEM),NDG,
     &  NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),NHST,NPLIST2(0:NPM),
     &  NPNE(NNM,NBFM,NEM),nr,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     &  NYNP(NKM,NVM,NHM,NPM),NYNR(0:NY_R_M,0:NRCM,NCM),nx,
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NZTOT,SPARSENESS
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IDATA(50),na,nb,ndb,ne,ne2,nh,nhs,nhx,nn,noelem,noelem2,
     '  nonode,no_ny,no_nynr,np,np1,np2,ny,ny1,ny2,
     '  ny2c,nz,nzz,nzz2,nzz_row,ost1,ost2,ost3
      LOGICAL FLOW_BALANCED
      CHARACTER CHAR*1
!     External functions
      INTEGER IDIGITS

      CALL ENTERS('CALC_SPARSE_GKK_1DTREE',*9999)

      CALL ASSERT(USE_SPARSE.EQ.1,
     '  '>>Set USE_SPARSE to 1 to use sparsity arrays',ERROR,*9999)
C***  Calculate sparsity pattern
      ost1=0
      ost2=0
      NDG=0
      IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.ITYP3(nr,nx).EQ.6)
     &  THEN !pulmonary capillaries or P-R-F
        NPLIST2(0)=0
        DO nonode=1,NPM
          NPLIST2(nonode)=0 !initialise
        ENDDO !nonode
        nzz_row=1 !for position in ISR (for capillaries)
      ENDIF
C***  Initalise sparsity arrays
      IF(SPARSENESS.NE.5)THEN
        IF((N+1).GT.ISRMAX) THEN
          WRITE(CHAR,'(I1)') IDIGITS(N+1)
          WRITE(ERROR,'(''>>Increase NISR_GKKM to '',I'//CHAR//')')
     &      (N+1)
          GO TO 9999
        ENDIF
      ENDIF
      DO i=1,M+1
        ISR(i)=1
      ENDDO !i
      DO nz=1,ISCMAX
        ISC(nz)=0
      ENDDO !nz
      nzz=1
      IF(LUNG_SPARSITY)THEN !use direct assembly of ISC and ISR
        DO noelem=1,NEELEM(0)
          ne=NEELEM(noelem)
          nb=NBJ(1,ne)
          IF(NXI(-1,0,ne).EQ.0) THEN !stem branch: only two in row
            IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5) THEN !.OR.(ITYP3(nr,nx).EQ.3).AND.
              np1=NPNE(1,nb,ne) !start node #
              DO nhx=1,NHST
                nh=NH_LOC(nhx,nx)
                ny1=NYNP(1,1,nh,np1) !row #
                IF(.NOT.FIX(ny1,1))THEN
                  NDG=NDG+1
                  DO nhs=1,2 !for each row entry
                    IF(nzz.LE.ISCMAX) THEN
                      np2=NPNE(nhs,nb,ne)
                      ny2=NYNP(1,1,nh,np2) !col #
                      IF(.NOT.FIX(ny2,1))THEN
                        ost2=0
                        DO ny=ny2-1,1,-1 !loop over previous columns
                          IF(FIX(ny,1))THEN !new line
                            ost2=ost2+1 !count # of fixed BC
                          ENDIF
                        ENDDO !ny
                        ISC(nzz)=ny2-ost2 !store the col #
                        nzz=nzz+1
                        IF(nh.EQ.1.AND.NHST.EQ.2)THEN
                          ny2c=NYNP(1,1,2,np2) !col # for w.v. conc
                          ISC(nzz)=ny2c-ost2 !store the col #
                          nzz=nzz+1
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !nhs
                ELSE
                  ost1=ost1+1
                ENDIF
                ISR(ny1-ost1+1)=nzz
              ENDDO !nh
            ENDIF !ITYP3
          ENDIF !NXI
          np1=NPNE(2,nb,ne)
          DO nhx=1,NHST
            nh=NH_LOC(nhx,nx)
            ny1=NYNP(1,1,nh,np1) !row #
            IF((.NOT.FIX(ny1,1)).OR.(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,
     &        nx).EQ.4.OR.ITYP3(nr,nx).EQ.6).AND.(NENP(np1,0).LE.1))THEN !temporary (must test)
              NDG=NDG+1
              DO nhs=1,2 !for each row entry
                IF(nzz.LE.ISCMAX) THEN
                  np2=NPNE(nhs,nb,ne)
                  ny2=NYNP(1,1,nh,np2) !col #
                  IF(.NOT.FIX(ny2,1))THEN
                    ost2=0
                    DO ny=ny2-1,1,-1 !loop over previous columns
                      IF(FIX(ny,1))THEN !new line
                        ost2=ost2+1 !count # of fixed BC
                      ENDIF
                    ENDDO !ny
                    ISC(nzz)=ny2-ost2 !store the col #
                    nzz=nzz+1
                    IF(nh.EQ.1.AND.NHST.EQ.2)THEN
                      ny2c=NYNP(1,1,2,np2) !col # for w.v. conc
                      ISC(nzz)=ny2c-ost2 !store the col #
                      nzz=nzz+1
                    ENDIF
                  ENDIF
                ELSE
                  WRITE(CHAR,'(I1)') IDIGITS(nzz)
                  WRITE(ERROR,
     &              '(''>>Increase NISC_GKKM to at least '',I'//CHAR
     &              //')') nzz
                  GO TO 9999
                ENDIF
              ENDDO !nhs
              IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.ITYP3(nr,
     &          nx).EQ.6)THEN !capillaries
                DO na=1,NAT(nb) !auxiliary variable - flow
                  IF(nzz.LE.ISCMAX) THEN
                    ny2c=NYNE(na,nh,1,1,ne)
                    IF(.NOT.FIX(ny2c,1)) THEN!.OR.(NENP(np1,0).LE.1)) THEN
                      ost2=0
                      DO ny=ny2c-1,1,-1 !loop over previous columns
                        IF(FIX(ny,1))THEN !new line
                          ost2=ost2+1 !count # of fixed BC
                        ENDIF
                      ENDDO !ny
                      ISC(nzz)=ny2c-ost2 !store column #
                      nzz=nzz+1
                    ENDIF !FIX
                  ELSE !iff nzz.GT.ISCMAX -> ERROR
                    WRITE(CHAR,'(I1)') IDIGITS(nzz)
                    WRITE(ERROR,
     &                '(''>>Increase NISC_GKKM to at least '',I'//CHAR
     &                //')') nzz
                    GO TO 9999
                  ENDIF !nzz
                ENDDO !na
                nzz_row=nzz_row+1
                ISR(nzz_row)=nzz !stores row #
              ENDIF !ITYP3
              IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5) THEN
                DO ndb=1,NXI(1,0,ne) !for each daughter branch
                  ne2=NXI(1,ndb,ne) !daughter element #
                  np2=NPNE(2,nb,ne2) !end node #
                  IF(nzz.LE.ISCMAX) THEN
                    ny2=NYNP(1,1,nh,np2) !col #
                    IF(.NOT.FIX(ny2,1))THEN
                      ost2=0
                      DO ny=ny2-1,1,-1 !loop over previous columns
                        IF(FIX(ny,1))THEN !new line
                          ost2=ost2+1 !count # of fixed BC
                        ENDIF
                      ENDDO !ny
                      ISC(nzz)=ny2-ost2
                      nzz=nzz+1
                      IF(nh.EQ.1.AND.NHST.EQ.2)THEN
                        ny2c=NYNP(1,1,2,np2) !col # for w.v. conc
                        ISC(nzz)=ny2c-ost2 !store the col #
                        nzz=nzz+1
                      ENDIF
                    ENDIF
                  ELSE
                    WRITE(CHAR,'(I1)') IDIGITS(nzz)
                    WRITE(ERROR,
     &                '(''>>Increase NISC_GKKM to at least '',I'//CHAR
     &                //')') nzz
                    GO TO 9999
                  ENDIF
                ENDDO !ndb
              ELSEIF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.ITYP3(nr,
     &            nx).EQ.6) THEN !capillaries
                DO nn=1,NNT(nb) !balances at each node of ne
                  FLOW_BALANCED=.FALSE. !initialise
                  np=NPNE(nn,nb,ne)
                  DO nonode=1,NPLIST2(0) !junctions balanced at already
                    IF(np.EQ.NPLIST2(nonode)) THEN
                      FLOW_BALANCED=.TRUE. !already balanced at np
                    ENDIF
                  ENDDO !nonode
                  IF((NENP(np,0).GT.1).AND.(.NOT.FLOW_BALANCED))THEN
                    DO na=1,NAT(nb) !at an unbalanced junction
                      DO noelem2=1,NENP(np,0) !for elems with
                        ne2=NENP(np,noelem2) !subtended branches only
                        IF(nzz.LT.ISCMAX) THEN
                          ny2=NYNE(na,nh,1,1,ne2)
                          IF(.NOT.FIX(ny2,1))THEN
                            ost2=0
                            DO ny=ny2-1,1,-1 !loop over previous columns
                              IF(FIX(ny,1))THEN !new line
                                ost2=ost2+1 !count # of fixed BC
                              ENDIF
                            ENDDO !ny
                            ISC(nzz)=ny2-ost2
                            nzz=nzz+1
                          ENDIF
                        ELSE
                          WRITE(CHAR,'(I1)') IDIGITS(nzz)
                          WRITE(ERROR,
     &                      '(''>>Increase NISC_GKKM to at least '',I'
     &                      //CHAR//')') nzz
                          GO TO 9999
                        ENDIF
                      ENDDO !noelem2
                    ENDDO !na
                    NPLIST2(0)=NPLIST2(0)+1 !stores nodes where flow
                    NPLIST2(NPLIST2(0))=np !balance been done
                  ENDIF !NENP(np,0).GT.1
                  IF((.NOT.FLOW_BALANCED).AND.(NENP(np,0).GT.1))THEN
                    nzz_row=nzz_row+1
                    ISR(nzz_row)=nzz !store the row #
                    NDG=NDG+1 !# entries on diagonal
                  ENDIF !.NOT.FLOW_BALANCED
                ENDDO !nn
              ENDIF !ITYP3
            ELSE
              ost1=ost1+1
            ENDIF
            IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5) THEN
              ISR(ny1-ost1+1)=nzz !store the row #
            ENDIF !ITYP3
          ENDDO !nh
        ENDDO !noelem
      ELSE
        ost1=0
        no_ny=1
        nzz=1
        ISR(1)=1
        DO ny1=1,NYNR(0,0,1)
          IF(.NOT.FIX(ny1,1))THEN !put entries into ISC,ISR
            ost3=0
            DO no_nynr=ISR_GK(ny1),ISR_GK(ny1+1)-1 !for each row entry
              ny2=ISC_GK(no_nynr) !unreduced ny
              IF(.NOT.FIX(ny2,1))THEN
                ost2=0
                DO ny=ny2-1,1,-1 !loop over previous columns
                  IF(FIX(ny,1))THEN
                    ost2=ost2+1 !count # of fixed bc
                  ENDIF
                ENDDO
                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,NZT(1,nx),
     '            ISC_GK,ISR_GK,6,ERROR,*9999)
                ISC(nzz)=ISC_GK(nz)-ost2
                nzz=nzz+1
              ELSE
                ost3=ost3+1
              ENDIF
            ENDDO !ny2
            no_ny=no_ny+1
            ISR(no_ny)=nzz
          ELSE
            ost1=ost1+(ISR_GK(ny1+1)-ISR_GK(ny1))
          ENDIF
        ENDDO
        NDG=no_ny-1
      ENDIF
      NZTOT=nzz !total # of non-zero entries
      IF(NZTOT.GT.ISCMAX) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NZTOT)
        WRITE(ERROR,'(''>>Increase NISC_GKKM to '',I'//CHAR//')')
     &    NZTOT
        GO TO 9999
      ENDIF
C KSB (26/09/02): ISC still needs to be reordered for capillary problems
      IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.4.OR.ITYP3(nr,nx).EQ.6)
     &  THEN
        nzz=0
        nzz2=0
        DO ost1=1,NDG
          DO ost2=1,(ISR(ost1+1)-ISR(ost1))
            nzz=nzz+1
            IDATA(ost2)=ISC(nzz)
          ENDDO !ost2
          CALL ISORT(ISR(ost1+1)-ISR(ost1),IDATA)
          DO ost2=1,(ISR(ost1+1)-ISR(ost1))
            nzz2=nzz2+1
            ISC(nzz2)=IDATA(ost2) !sorts into increasing value
          ENDDO !ost2
        ENDDO !ost1
      ENDIF !ITYP3

      IF(NZTOT.GT.NZ_GKK_M) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NZTOT)
        WRITE(ERROR,'(''>>Increase NZ_GKK_M to '',I'//CHAR//')')
     &    NZTOT
        GO TO 9999
      ENDIF

      CALL EXITS('CALC_SPARSE_GKK_1DTREE')
      RETURN
 9999 CALL ERRORS('CALC_SPARSE_GKK_1DTREE',ERROR)
      CALL EXITS('CALC_SPARSE_GKK_1DTREE')
      RETURN 1
      END



