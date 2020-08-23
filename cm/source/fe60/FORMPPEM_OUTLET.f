      SUBROUTINE FORMPPEM_OUTLET(FACE_TOT,MAXLOC_FACES,NBJ,NENP,
     '  NFVC_OUTLET,NPLIST,NPNE,NPNODE,nr,NVCB,NVCBNODE_OUTLET,
     '  NVCB_OUTLET,NVCNODE,N_OUTLET,OUTLET_MATRIX,XNFV_OUTLET,XP,ZA,
     '  ERROR,*)

C#### Subroutine: FORMPPEM_OUTLET
C###  Description:
C###    <HTML>
C###    Forms the pressure Poisson outlet matrix.
C###    <PRE>
C###
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER N_OUTLET,FACE_TOT,MAXLOC_FACES,NBJ(NJM,NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NFVC_OUTLET(2,0:MAXLOC_FACES,N_OUTLET),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVCB(-1:3,NVCBM),NVCBNODE_OUTLET(NP_R_M),
     '  NVCB_OUTLET(N_OUTLET),NVCNODE(2,NP_R_M)
      REAL*8 OUTLET_MATRIX(N_OUTLET,N_OUTLET),
     '  XNFV_OUTLET(-1:3,FACE_TOT),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER bnvc,LOC_N_OUTLET,LOC_FACE_TOT,nonode,np,ELEMCOUNT,
     '  noelem,ne,nb,INTCOUNT,SIMPLEXLIST(60),ADJNODE,nn,cnonode,
     '  nfvl,cnp,FTOT,njj,nj,m,noelem2,ne2,nb2,nn2,CNE(2),bnvc_outlet,
     '  cbnvc_outlet,nfv,i,j,BNENP(0:100),noelem3,ne3
      REAL*8 DXP(3),DIAG,TEMP
      LOGICAL ADD,FOUND

      CALL ENTERS('FORMPPEM_OUTLET',*9999)

      LOC_N_OUTLET=0
      LOC_FACE_TOT=0

C     ..Set up NPLIST to get the nonode number of each np
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPLIST(np)=nonode
        NVCBNODE_OUTLET(nonode)=0
      ENDDO

      DO bnvc=1,NVBT

C       ..Outlet nodes only
        IF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN
          nonode=NVCB(1,bnvc)
          LOC_N_OUTLET=LOC_N_OUTLET+1
          NVCB_OUTLET(LOC_N_OUTLET)=bnvc
          NVCBNODE_OUTLET(nonode)=LOC_N_OUTLET
          np=NPNODE(nonode,nr)

C         ..Compute the boundary simplices only
          ELEMCOUNT=0
          DO noelem=1,NENP(np,0,nr)
            ne=NENP(np,noelem,nr)
            nb=NBJ(1,ne)
            INTCOUNT=0
            DO nn=1,NNT(nb)
              ADJNODE=NPNE(nn,nb,ne)
              cnonode=NPLIST(ADJNODE)
C             ..Split up into if..elseif statement in case we get an
C               array out of bounds on NVCB
              IF(NVCNODE(TYPE,cnonode).EQ.INTERNAL.OR.
     '          NVCNODE(TYPE,cnonode).EQ.INTBOUN) THEN
                INTCOUNT=INTCOUNT+1
              ELSEIF(NVCNODE(TYPE,cnonode).EQ.BOUNDARY.AND.
     '            NVCB(BCTYPE,NVCNODE(MAP,cnonode)).EQ.WALL) THEN
                INTCOUNT=INTCOUNT+1
              ENDIF
            ENDDO
            IF(INTCOUNT.LT.NNT(nb)) THEN
              ELEMCOUNT=ELEMCOUNT+1
              SIMPLEXLIST(ELEMCOUNT)=ne
            ENDIF
          ENDDO

          NFVC_OUTLET(1,0,LOC_N_OUTLET)=0
          DO noelem=1,ELEMCOUNT
            ne=SIMPLEXLIST(noelem)

C           ..Loop over the element's local nodes to get the surrounding
C           nodes
            nb=NBJ(1,ne)
            DO nn=1,NNT(nb)
              ADJNODE=NPNE(nn,nb,ne)

              IF(ADJNODE.NE.np) THEN

                cnonode=NPLIST(ADJNODE)

C               ..Only IB-IB connections
                IF(NVCNODE(TYPE,cnonode).EQ.INTBOUN) THEN

C                 ..Make sure that we haven't done this node already
                  ADD=.TRUE.
                  FOUND=.FALSE.
                  nfvl=1
                  DO WHILE(.NOT.FOUND.AND.nfvl.LE.
     '              NFVC_OUTLET(1,0,LOC_N_OUTLET))
                    cnp=NPNODE(NFVC_OUTLET(1,nfvl,LOC_N_OUTLET),nr)
                    IF(ADJNODE.EQ.cnp) THEN
                      ADD=.FALSE.
                      FOUND=.TRUE.
                    ELSE
                      nfvl=nfvl+1
                    ENDIF
                  ENDDO
                ELSE
                  ADD=.FALSE.
                ENDIF

C               ..Increment arrays
                IF(ADD) THEN
                  LOC_FACE_TOT=LOC_FACE_TOT+1
                  FTOT=NFVC_OUTLET(1,0,LOC_N_OUTLET)+1
                  NFVC_OUTLET(1,FTOT,LOC_N_OUTLET)=cnonode
                  NFVC_OUTLET(2,FTOT,LOC_N_OUTLET)=LOC_FACE_TOT
                  NFVC_OUTLET(1,0,LOC_N_OUTLET)=FTOT
                  NFVC_OUTLET(2,0,LOC_N_OUTLET)=FTOT
                  IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN

C                   ..Compute distance between nodes
                    XNFV_OUTLET(IDIST,LOC_FACE_TOT)=0.d0
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                      DXP(njj)=XP(1,1,nj,ADJNODE)-XP(1,1,nj,np)
                      XNFV_OUTLET(IDIST,LOC_FACE_TOT)=
     '                  XNFV_OUTLET(IDIST,LOC_FACE_TOT)+DXP(njj)**2
                    ENDDO

C                   ..Invert
                    XNFV_OUTLET(IDIST,LOC_FACE_TOT)=
     '                1.d0/DSQRT(XNFV_OUTLET(IDIST,LOC_FACE_TOT))

C                   ..Get the unit normal
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      XNFV_OUTLET(njj,LOC_FACE_TOT)=
     '                  DXP(njj)*XNFV_OUTLET(IDIST,LOC_FACE_TOT)
                    ENDDO

C                   ..In 2d the facial area will always be 1
                    XNFV_OUTLET(FAREA,LOC_FACE_TOT)=1.d0

                    IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                    call mp_setlock()
                      WRITE(OP_STRING,'(/,'' Voronoi face:  '',I12)')
     '                  LOC_FACE_TOT
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' Normal:        '',2F12.6)')
     '                  (XNFV_OUTLET(j,LOC_FACE_TOT),j=1,2)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' Dist:          '',F12.6)')
     '                  1.d0/XNFV_OUTLET(IDIST,LOC_FACE_TOT)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                    call mp_unsetlock()
                    ENDIF

                  ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN

C                   ..Compute distance between nodes
                    XNFV_OUTLET(IDIST,LOC_FACE_TOT)=0.d0
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nj=NJ_LOC(NJL_GEOM,njj,nr)
                      DXP(njj)=XP(1,1,nj,ADJNODE)-XP(1,1,nj,np)
                      XNFV_OUTLET(IDIST,LOC_FACE_TOT)=
     '                  XNFV_OUTLET(IDIST,LOC_FACE_TOT)+DXP(njj)**2
                    ENDDO

C                   ..Invert
                    XNFV_OUTLET(IDIST,LOC_FACE_TOT)=
     '                1.d0/DSQRT(XNFV_OUTLET(IDIST,LOC_FACE_TOT))

C                   ..Get the unit normal
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      XNFV_OUTLET(njj,LOC_FACE_TOT)=
     '                  DXP(njj)*XNFV_OUTLET(IDIST,LOC_FACE_TOT)
                    ENDDO

C                   ..Compute the face area

C                   ..First construct the boundary simplices surrounding
C                     node np
                    BNENP(0)=0
                    DO noelem2=1,NENP(np,0,nr)
                      ne2=NENP(np,noelem2,nr)
                      DO noelem3=1,ELEMCOUNT
                        ne3=SIMPLEXLIST(noelem3)
                        IF(ne3.EQ.ne2) THEN
                          BNENP(0)=BNENP(0)+1
                          BNENP(BNENP(0))=ne3
                        ENDIF
                      ENDDO
                    ENDDO

C                   ..Now find the 2 joining elements, and the distance
C                     between their circumcentres that will be the
C                     face area
                    m=1
                    noelem2=1
                    DO WHILE(m.LE.2.AND.noelem2.LE.BNENP(0))
                      ne2=BNENP(noelem2)
                      nb2=NBJ(1,ne2)
                      DO nn2=1,NNT(nb2)
                        IF(NPNE(nn2,nb2,ne2).EQ.ADJNODE) THEN
                          CNE(m)=ne2
                          m=m+1
                        ENDIF
                      ENDDO
                      noelem2=noelem2+1
                    ENDDO
                    XNFV_OUTLET(FAREA,LOC_FACE_TOT)=0.d0
                    DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
                      XNFV_OUTLET(FAREA,LOC_FACE_TOT)=
     '                  XNFV_OUTLET(FAREA,LOC_FACE_TOT)+
     '                  (ZA(1,njj,1,CNE(1))-ZA(1,njj,1,CNE(2)))**2
                    ENDDO
                    XNFV_OUTLET(FAREA,LOC_FACE_TOT)=
     '                DSQRT(XNFV_OUTLET(FAREA,LOC_FACE_TOT))
                  ENDIF
                ENDIF
              ENDIF !intboun
            ENDDO !nn
          ENDDO !noelem
        ENDIF !outlet
      ENDDO !bnvc


      CALL ASSERT(LOC_N_OUTLET.EQ.N_OUTLET,
     '  '>> different number of outlet cells found',ERROR,*9999)
      CALL ASSERT(LOC_FACE_TOT.EQ.FACE_TOT,
     '  '>> different number of outlet faces found',ERROR,*9999)

C     ..Zero all the matrix entries
      DO i=1,N_OUTLET
        DO j=1,N_OUTLET
          OUTLET_MATRIX(i,j)=0.d0
        ENDDO
      ENDDO

C     ..Form the outlet pressure Poisson matrix
      DO bnvc_outlet=1,N_OUTLET
        DIAG=0.d0
        DO nfvl=1,NFVC_OUTLET(1,0,bnvc_outlet)
          nfv=NFVC_OUTLET(2,nfvl,bnvc_outlet)
          cnonode=NFVC_OUTLET(1,nfvl,bnvc_outlet)
          cbnvc_outlet=NVCBNODE_OUTLET(cnonode)
          TEMP=XNFV_OUTLET(FAREA,nfv)*XNFV_OUTLET(IDIST,nfv)
          DIAG=DIAG+TEMP
          OUTLET_MATRIX(bnvc_outlet,cbnvc_outlet)=-TEMP
        ENDDO
        OUTLET_MATRIX(bnvc_outlet,bnvc_outlet)=DIAG
      ENDDO

C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$    call mp_setlock()
      WRITE(OP_STRING,'(/$,'' Outlet Pressure matrix:'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' ########################'')')
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      DO i=1,N_OUTLET
        WRITE(OP_STRING,'(200D12.4)')
     '    (OUTLET_MATRIX(i,j),j=1,N_OUTLET)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDDO
      DO bnvc=1,N_OUTLET
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/$,''Adjacency Molecule:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''==================='')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/$,''bnvc_outlet = '',I4,'' bnvc '//
     '    '= '',I5)') bnvc,NVCB_OUTLET(bnvc)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/,''  Node  Face        Area'//
     '    '        Dist        Norm'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nfvl=1,NFVC_OUTLET(1,0,bnvc)
          WRITE(OP_STRING,
     '      '(I6,I6,D12.4,D12.4,''      '',3D12.4)')
     '      NFVC_OUTLET(1,nfvl,bnvc),NFVC_OUTLET(2,nfvl,bnvc),
     '      XNFV_OUTLET(FAREA,NFVC_OUTLET(2,nfvl,bnvc)),
     '      1.0D0/XNFV_OUTLET(IDIST,NFVC_OUTLET(2,nfvl,bnvc)),
     '      (XNFV_OUTLET(nj,NFVC_OUTLET(2,nfvl,bnvc)),
     '      nj=1,NJ_LOC(NJL_GEOM,0,nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDDO
CC$    call mp_unsetlock()

      CALL EXITS('FORMPPEM_OUTLET')
      RETURN
 9999 CALL ERRORS('FORMPPEM_OUTLET',ERROR)
      CALL EXITS('FORMPPEM_OUTLET')
      RETURN 1
      END


