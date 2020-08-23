      SUBROUTINE FORMPPEM(ISC_GKK,ISR_GKK,NFVC,nr,NVCNODE,
     '  nx,GKK,XNFV,ERROR,*)

C#### Subroutine: FORMPPEM
C###  Description:
C###    <HTML>
C###    Forms the pressure Poisson matrix.
C###    <PRE>
C###
C###    Laplacian(P) is pressure Poisson matrix
C###
C###    L(pi) = Sumj[(pj - pi)(Aij/dij)]
C###
C###    Where P  = pressure
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER ISC_GKK(NISC_GKKM),ISR_GKK(NISR_GKKM),
     '  NFVC(2,0:NFVCM,NVCM),nr,NVCNODE(2,NP_R_M),nx
      REAL*8 GKK(NZ_GKK_M),XNFV(-(NJM+1):NJM,NFVM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*4 WORK_PTR
      INTEGER VECTMAX
      PARAMETER(VECTMAX=100)
      INTEGER nvc,nfvl,nfv,nz,i,j,ROWNZTOT,COL(VECTMAX),cnonode,cnvc,
     '  NCOL
      REAL*8 DIAG,TEMP,ROWNZ(VECTMAX)
      LOGICAL MEM_INIT

      CALL ENTERS('FORMPPEM',*9999)

C     ..Calculate sparsity pattern if required

C The reason why I am calculating the sparsity pattern here for
C compressed row is that the Image calculation is *TOO* slow.
C It uses an O(n^2) loop. If the mesh is moving, then the matrix
C has to be recalculated at each step, so it takes too long. It is
C only done for compressed row as this is the sparsity pattern that
C the conjugate gradient method uses, which is the only sort of solver
C that you would want to use for a moving mesh anyway.

      IF(SPARSEGKK(nx).EQ.0) THEN ! No sparsity pattern
        DO i=1,NVCT
          DO j=1,NVCT
            nz=i+(j-1)*NVCT
            CALL ASSERT(nz.LE.NZ_GKK_M,'>>Increase NZ_GKK_M',
     '        ERROR,*9999)
            GKK(nz)=0.d0
          ENDDO
        ENDDO
      ELSEIF(SPARSEGKK(nx).EQ.1.AND..NOT.MESHFIXD) THEN ! Compressed row
        CALL ASSERT((NVCT+1).LE.NISR_GKKM,'>>Increase NISR_GKKM',
     '    ERROR,*9999)
        nz=1
        DO nvc=1,NVCT
          ISR_GKK(nvc)=nz

C         ..Get the column entries
          NCOL=1
          DO nfvl=1,NFVC(1,0,nvc)
            cnonode=NFVC(1,nfvl,nvc)
            IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
              cnvc=NVCNODE(MAP,cnonode)
              IF(NCOL.LE.VECTMAX) THEN
                COL(NCOL)=cnvc
              ENDIF
              NCOL=NCOL+1
            ENDIF
          ENDDO

C         ..Don't forget the diagonal entry !
          COL(NCOL)=nvc

          CALL ASSERT(NCOL.LE.VECTMAX,'>>Increase vectmax',ERROR,*9999)

C         ..Sort the column entries
          IF(NCOL.LE.50) THEN !Use shell sort for small N
            CALL ISHELLSORT(NCOL,COL)
          ELSE
            CALL IHEAPSORT(NCOL,COL)
          ENDIF

C         ..Put them into the sparsity structure
          DO j=1,NCOL
            cnvc=COL(j)
            IF(nz.LE.NISC_GKKM) THEN
              ISC_GKK(nz)=cnvc
            ENDIF
            nz=nz+1
          ENDDO
        ENDDO
        ISR_GKK(NVCT+1)=nz
        NZT(1,nx)=nz-1
        NZZT(1,nr,nx)=NZT(1,nx)
        CALL ASSERT(NZT(1,nx).LE.NZ_GKK_M,'>>Increase NZ_GKK_M',
     '    ERROR,*9999)
        CALL ASSERT(NZT(1,nx).LE.NISC_GKKM,'>>Increase NISC_GKKM',
     '    ERROR,*9999)
      ELSE ! All others - use Image matrix
        WORK_PTR=0
        MEM_INIT=.FALSE.
        CALL ALLOCATE_MEMORY(NVCT**2,1,
     '    CHARTYPE,WORK_PTR,MEM_INIT,ERROR,*9999)
        CALL SPRSVORO(ISC_GKK,ISR_GKK,NVCT,NVCT,NFVC,nr,
     '    NVCNODE,nx,%VAL(WORK_PTR),ERROR,*9999)
        CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
      ENDIF

C     ..Compute matrix entries
      DO nvc=1,NVCT
        DIAG=0.d0
        DO nfvl=1,NFVC(1,0,nvc)
          cnonode=NFVC(1,nfvl,nvc)
          IF(NVCNODE(TYPE,cnonode).NE.BOUNDARY) THEN
            cnvc=NVCNODE(MAP,cnonode)
            nfv=NFVC(2,nfvl,nvc)
            CALL SPARSE(nvc,cnvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),ISC_GKK,
     '        ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
            TEMP=XNFV(FAREA,nfv)*XNFV(IDIST,nfv)
            GKK(nz)=TEMP
            DIAG=DIAG-TEMP
          ENDIF
        ENDDO
        CALL SPARSE(nvc,nvc,NVCT,nz,NZ_GKK_M,NZT(1,nx),ISC_GKK,
     '    ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
        GKK(nz)=DIAG
      ENDDO


      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/$,'' Pressure correction matrix:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' ###########################'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO i=1,NVCT
          ROWNZTOT=1
          DO j=1,NVCT
            CALL SPARSE(i,j,NVCT,nz,NZ_GKK_M,NZT(1,nx),ISC_GKK,
     '        ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
            IF(nz.NE.0) THEN
              ROWNZ(ROWNZTOT)=GKK(nz)
              COL(ROWNZTOT)=j
              ROWNZTOT=ROWNZTOT+1
            ENDIF
          ENDDO
          ROWNZTOT=ROWNZTOT-1
          WRITE(OP_STRING,'(''Row:     '',I12)') i
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''Cols:    '',100I12)')
     '      (COL(j),j=1,ROWNZTOT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''Colents: '',100D12.4)')
     '      (ROWNZ(j),j=1,ROWNZTOT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$     call mp_unsetlock()
      ENDIF

      CALL EXITS('FORMPPEM')
      RETURN
 9999 CALL ERRORS('FORMPPEM',ERROR)
      CALL EXITS('FORMPPEM')
      RETURN 1
      END


