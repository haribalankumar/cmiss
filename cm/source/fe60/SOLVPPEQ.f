      SUBROUTINE SOLVPPEQ(FACETOT_OUTLET,IPIVOT_OUTLET,ISC_GKK,
     '  ISR_GKK,MAXLOC_OUTLET_FACES,NFVC,NFVC_OUTLET,NODENVC,NODENVCB,
     '  NPLIST,NPNODE,nr,NVCB,NVCBNODE_OUTLET,NVCB_OUTLET,NVCNODE,
     '  nx,NYNP,N_OUTLET,GKK,GRR,OUTLET_MATRIX,OUTLET_RHS,
     '  XNFV,XNFV_OUTLET,XO,YP,FIRST_A,REFACT_OUTLET,
     '  UPDATE_MATRIX,x_INIT,ERROR,*)

C#### Subroutine: SOLVPPEQ
C###  Description:
C###    Solves the pressure Poisson equation.
C###    System represented as Ax=b

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER N_OUTLET,FACETOT_OUTLET,IPIVOT_OUTLET(N_OUTLET),
     '  ISC_GKK(NISC_GKKM),ISR_GKK(NISR_GKKM),MAXLOC_OUTLET_FACES,
     '  NFVC(2,0:NFVCM,NVCM),
     '  NFVC_OUTLET(2,0:MAXLOC_OUTLET_FACES,N_OUTLET),NODENVC(NVCM),
     '  NODENVCB(NVCBM),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVCB(-1:3,NVCBM),NVCBNODE_OUTLET(NP_R_M),
     '  NVCB_OUTLET(N_OUTLET),NVCNODE(2,NP_R_M),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 GKK(NZ_GKK_M),GRR(NOM),OUTLET_MATRIX(N_OUTLET,N_OUTLET),
     '  OUTLET_RHS(N_OUTLET),XNFV(-(NJM+1):NJM,NFVM),
     '  XNFV_OUTLET(-1:3,FACETOT_OUTLET),XO(NOM),YP(NYM,NIYM)
      LOGICAL FIRST_A,REFACT_OUTLET,UPDATE_MATRIX,x_INIT
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER fixnonode,fixnvc,nz,nfvl,cnonode,cnvc,cnonode2,cnvc2,
     '  ibnvc,ibnonode,nvc,nonode,np,ny,nfvl2,bnvc_outlet,INFO,
     '  ibnp,cbnonode,cbnp,nfv,nh,nhx,ibny,cbny,inonode,bnvc,i_outlet,
     '  fixnp,FIXOUTNOD
      REAL*8 DIV,TEMP

      CALL ENTERS('SOLVPPEQ',*9999)

C     ..Compute the list of nonode numbers foreach np
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        NPLIST(np)=nonode
      ENDDO
      IF(CONFINED) THEN

        fixnonode=NPLIST(FIXDNODE)
        fixnvc=NVCNODE(MAP,fixnonode)
        GRR(fixnvc)=0.d0

C       ..Need to rezero the fixed node entry in the matrix if it needs
C         to be refactorized
        IF((FIRST_A.AND.UPDATE_MATRIX).OR.SOLVEROPTION(nx).NE.1) THEN

C         ..Find the corresponding nonode of ADJNODE

C         Fix the pressure at the fixed node. Adjusts the pressure
C         matrix so that row and column fixnvc are zero, except for
C         A(fixnvc,fixnvc), which is 1. The rhs of fixnvc should
C         also equal zero.

C         ..Get the sparse non-zero diagonal entry nz
          CALL SPARSE(fixnvc,fixnvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),
     '      ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)

C         ..Fix the matrix entry to be 1 and the RHS to be 0
          GKK(nz)=1.d0

C         ..Loop over the adjacent nodes of Voronoi cell fixnvc
C         (corresponding to the rows of the Pressure matrix)
          DO nfvl=1,NFVC(1,0,fixnvc)
            cnonode=NFVC(1,nfvl,fixnvc)

C           ..Internal nodes only
            IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
              cnvc=NVCNODE(MAP,cnonode)

C             ..Get the sparse non-zero entry nz and set it to be zero
              CALL SPARSE(fixnvc,cnvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),
     '          ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
              GKK(nz)=0.d0

C             ..Loop over the adjoining nodes of the adjoining nodes
C             (this corresponds to the columns of the pressure matrix)
              DO nfvl2=1,NFVC(1,0,cnvc)
                cnonode2=NFVC(1,nfvl2,cnvc)

                IF(NVCNODE(TYPE,cnonode2).NE.BOUNDARY) THEN
                  cnvc2=NVCNODE(MAP,cnonode2)

C                 ..Fix the column entry
                  IF(cnvc2.EQ.fixnvc) THEN
                    CALL SPARSE(cnvc,cnvc2,NVCT,nz,NZ_GKK_M,NZT(1,nx),
     '                ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
                    GKK(nz)=0.d0
                  ENDIF !cnvc2.eq.fixnvc
                ENDIF !internal
              ENDDO !nfvl2
            ENDIF !internal
          ENDDO !nfvl
        ENDIF !first_a.and.update_matrix

      ELSE !.not.confined

        IF(DUCTFLOW) THEN

C         ..Compute the planar velocity divergence at the outlet
          DO bnvc_outlet=1,N_OUTLET
            ibnonode=NVCB(1,NVCB_OUTLET(bnvc_outlet))
            ibnvc=NVCNODE(MAP,ibnonode)
            ibnp=NPNODE(ibnonode,nr)

            DIV=0.d0
            DO nfvl=1,NFVC_OUTLET(1,0,bnvc_outlet)
              cbnonode=NFVC_OUTLET(1,nfvl,bnvc_outlet)
              cbnp=NPNODE(cbnonode,nr)
              nfv=NFVC_OUTLET(2,nfvl,bnvc_outlet)

              TEMP=0.d0
              DO nhx=1,nh_loc(0,nx)-1
                nh=nh_loc(nhx,nx)
                ibny=NYNP(1,1,nh,ibnp,0,1,nr)
                cbny=NYNP(1,1,nh,cbnp,0,1,nr)

                TEMP=TEMP+(YP(ibny,1)+YP(cbny,1))*XNFV_OUTLET(nhx,nfv)
              ENDDO
              DIV=DIV+TEMP*0.5d0*XNFV_OUTLET(FAREA,nfv)
            ENDDO
            OUTLET_RHS(bnvc_outlet)=-DIV/DT
          ENDDO

C         ..Zero an entry
          DO i_outlet=1,NUM_OUTLETS
            fixnp=OUTLET_FIXDNODES(i_outlet)
            fixnonode=NPLIST(fixnp)
            fixnvc=NVCNODE(MAP,fixnonode)
            FIXOUTNOD=NVCBNODE_OUTLET(NVCB(1,fixnvc))
            IF(FIXOUTNOD.NE.0) THEN
              OUTLET_RHS(FIXOUTNOD)=0.d0
            ELSE
              ERROR='Could not find fixed outlet node'
              GOTO 9999
            ENDIF
          ENDDO

           IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/$,'' Outlet Pressure RHS'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ####################'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nvc=1,N_OUTLET
              nonode=NODENVCB(NVCB_OUTLET(nvc))
              np=NPNODE(nonode,nr)
              WRITE(OP_STRING,'('' Voronoi outlet: '',I6,
     '            '' (np = '',I6,'') RHS: '',D14.6)')
     '          nvc,np,OUTLET_RHS(nvc)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nvc
CC$          call mp_unsetlock()
          ENDIF !dop


C         ..Zero the matrix entry if matrix needs to be refactorized
          INFO=0
          IF(REFACT_OUTLET) THEN

            DO i_outlet=1,NUM_OUTLETS

C             ..Find the fixed node for the duct flow outlet(s)
              fixnp=OUTLET_FIXDNODES(i_outlet)
              fixnonode=NPLIST(fixnp)
              fixnvc=NVCNODE(MAP,fixnonode)
              FIXOUTNOD=NVCBNODE_OUTLET(NVCB(1,fixnvc))
              CALL ASSERT(FIXOUTNOD.NE.0,'Fixoutnod=0',ERROR,*9999)

C             ..Fix the diagonal fixnod matrix entry to be 1
              OUTLET_MATRIX(FIXOUTNOD,FIXOUTNOD)=1.d0

C             ..Loop over the adjacent nodes of Voronoi cell fixnvc
C               (corresponding to the rows of the Pressure matrix)
              DO nfvl=1,NFVC_OUTLET(1,0,FIXOUTNOD)
                cnonode=NFVC_OUTLET(1,nfvl,FIXOUTNOD)
                cnvc=NVCBNODE_OUTLET(cnonode)
                OUTLET_MATRIX(FIXOUTNOD,cnvc)=0.d0

C               ..Loop over the adjoining nodes of the adjoining nodes
C                 (this corresponds to the cols of the pressure matrix)
                DO nfvl2=1,NFVC_OUTLET(1,0,cnvc)
                  cnonode2=NFVC_OUTLET(1,nfvl2,cnvc)
                  cnvc2=NVCBNODE_OUTLET(cnonode2)

C                 ..Fix the column entry
                  IF(cnvc2.EQ.FIXOUTNOD) OUTLET_MATRIX(cnvc,cnvc2)=0.d0
                ENDDO !nfvl2
              ENDDO !nfvl
            ENDDO !i_outlet

            CALL DGETRF(N_OUTLET,N_OUTLET,OUTLET_MATRIX,N_OUTLET,
     '        IPIVOT_OUTLET,INFO)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(/,'' Factorised duct flow matrix'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ###########################'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nvc=1,N_OUTLET
                WRITE(OP_STRING,'(200D12.4)')
     '            (OUTLET_MATRIX(nvc,cnvc),cnvc=1,N_OUTLET)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
CC$            call mp_unsetlock()
            ENDIF !dop
          ENDIF !refact_outlet
          REFACT_OUTLET=.FALSE.

          IF(INFO.EQ.0) THEN
            CALL DGETRS('N',N_OUTLET,1,OUTLET_MATRIX,N_OUTLET,
     '        IPIVOT_OUTLET,OUTLET_RHS,N_OUTLET,INFO)
            IF(INFO.LT.0) THEN
              WRITE(ERROR,'('' >>DGETRS. The '',I1,''th argument '//
     '          'had an error.'')') INFO
              GOTO 9999
            ENDIF
          ELSE
            IF(INFO.LT.0) THEN
              WRITE(ERROR,'('' >>DGETRF. The '',I1,''th argument '//
     '          'had an error.'')') INFO
              GOTO 9999
            ELSE
              WRITE(ERROR,'('' >>DGETRF. OUTLET_MATRIX('',I3,'','',I3,
     '          '') is exactly zero.'')') INFO,INFO
              GOTO 9999
            ENDIF
          ENDIF

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/$,'' Outlet Pressure solution'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ########################'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nvc=1,N_OUTLET
              nonode=NODENVCB(NVCB_OUTLET(nvc))
              np=NPNODE(nonode,nr)
              WRITE(OP_STRING,'('' Voronoi outlet: '',I6,
     '            '' (np = '',I6,'') Pressure: '',D14.6)')
     '          nvc,np,OUTLET_RHS(nvc)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nvc
CC$          call mp_unsetlock()
          ENDIF !dop

          DO bnvc_outlet=1,N_OUTLET
            inonode=NVCB(1,NVCB_OUTLET(bnvc_outlet))
            nvc=NVCNODE(MAP,inonode)

            DO nfvl=1,NFVC(1,0,nvc)
              nonode=NFVC(1,nfvl,nvc)
              bnvc=NVCNODE(MAP,nonode)
              IF(NVCNODE(TYPE,nonode).EQ.BOUNDARY) THEN
                IF(NVCB(BCTYPE,bnvc).EQ.OUTLET.AND.
     '            NVCB(1,bnvc).EQ.inonode) THEN
                  IF(FIRST_A.AND.UPDATE_MATRIX) THEN
                    nfv=NFVC(2,nfvl,nvc)
                    TEMP=XNFV(IDIST,nfv)*XNFV(FAREA,nfv)
                    CALL SPARSE(nvc,nvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),
     '                ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
                    GKK(nz)=GKK(nz)-TEMP
                  ENDIF
                  GRR(nvc)=GRR(nvc)-TEMP*OUTLET_RHS(bnvc_outlet)
                ENDIF
              ENDIF
            ENDDO
          ENDDO

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/$,'' Duct flow RHS listing of PPE:'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' ###################'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nvc=1,NVCT
              WRITE(OP_STRING,
     '          '('' Voronoi cell: '',I7,'' RHS:'',D16.8)')
     '          nvc,GRR(nvc)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nvc
CC$          call mp_unsetlock()
          ENDIF !dop

        ELSE !.not.ductflow

C         ..Zero the RHS pressures
          DO bnvc=1,NVBT
            IF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN
              IF(NVCB(0,bnvc).EQ.1) THEN
                ibnonode=NVCB(1,bnvc)
                ibnvc=NVCNODE(MAP,ibnonode)
              ELSE
                ERROR='>>Should only be one internal node '//
     '            ' connected to an outlet boundary node'
                GOTO 9999
              ENDIF
              GRR(ibnvc)=0.d0
            ENDIF
          ENDDO

C         ..Need to zero the outlet pressures if not using LU or this
C           is the unfactorised matrix
          IF((FIRST_A.AND.UPDATE_MATRIX).OR.SOLVEROPTION(nx).NE.1) THEN
            DO bnvc=1,NVBT

              IF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN

                IF(NVCB(0,bnvc).EQ.1) THEN
                  ibnonode=NVCB(1,bnvc)
                  ibnvc=NVCNODE(MAP,ibnonode)
                ELSE
                  ERROR='>>Should only be one internal node '//
     '              ' connected to an outlet boundary node'
                  GOTO 9999
                ENDIF

C               ..Get the diagonal entry & set it to 1, the rhs to 0
                CALL SPARSE(ibnvc,ibnvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),
     '            ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
                GKK(nz)=1.d0

C               ..Loop over the adjacent nodes of Voronoi cell bnvc
C                 (corresponding to the rows of the Pressure matrix)
                DO nfvl=1,NFVC(1,0,ibnvc)
                  cnonode=NFVC(1,nfvl,ibnvc)

C                 ..Internal nodes only
                  IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
                    cnvc=NVCNODE(MAP,cnonode)

C                 ..Get the sparse entry nz and set it to be zero
                    CALL SPARSE(ibnvc,cnvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),
     '              ISC_GKK,ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
                    GKK(nz)=0.d0

C                   ..Loop over the adjacent nodes of Voronoi cell cnvc
C                     (ie the columns of the pressure matrix)
                    DO nfvl2=1,NFVC(1,0,cnvc)
                      cnonode2=NFVC(1,nfvl2,cnvc)

C                     ..Internal nodes only
                      IF(NVCNODE(TYPE,cnonode2).NE.BOUNDARY) THEN
                        cnvc2=NVCNODE(MAP,cnonode2)

C                       ..Column that matches up with the boundary cell
                        IF(cnvc2.EQ.ibnvc) THEN
                          CALL SPARSE(cnvc,cnvc2,NVCT,nz,NZ_GKK_M,
     '                      NZT(1,nx),ISC_GKK,ISR_GKK,SPARSEGKK(nx),
     '                      ERROR,*9999)
                          GKK(nz)=0.d0
                        ENDIF !cnvc2.eq.ibnvc
                      ENDIF !internal
                    ENDDO !nfvl2
                  ENDIF !internal
                ENDDO !nfvl
              ENDIF !outlet
            ENDDO !bnvc
          ENDIF !(first_a.and.update_matrix).or.solver.ne.LU
        ENDIF !ductflow
      ENDIF !confined

C     ..Solve the system of equations
      CALL SOLVE_SYSTEM(ISC_GKK,ISR_GKK,NVCT,NVCT,NVCT,NZT(1,nx),
     '  IWRIT4(nr,nx),PRECON_CODE(nx),SOLVEROPTION(nx),SPARSEGKK(nx),
     '  GKK,GRR,XO,FIRST_A,UPDATE_MATRIX,x_INIT,nx,ERROR,*9999)

C     ..Now transfer the pressures back into YP
      DO nvc=1,NVCT
        nonode=NODENVC(nvc)
        np=NPNODE(nonode,nr)
        ny=NYNP(1,1,nh_loc(0,nx),np,0,1,nr)
        YP(ny,1)=XO(nvc)
      ENDDO

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/$,'' Pressure solution'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' #################'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nvc=1,NVCT
          nonode=NODENVC(nvc)
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,'('' Voronoi cell: '',I6,'' (np = '',I6,
     '      '') Pressure: '',D14.6)') nvc,np,XO(nvc)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nvc
CC$     call mp_unsetlock()
      ENDIF !dop

      CALL EXITS('SOLVPPEQ')
      RETURN
 9999 CALL ERRORS('SOLVPPEQ',ERROR)
      CALL EXITS('SOLVPPEQ')
      RETURN 1
      END


