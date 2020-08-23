      SUBROUTINE IPSOLU_ITERATVE(NX,SOLVER,NOQUES,FILEIP,ERROR,*)

C#### Subroutine: IPSOLU_ITERATVE
C###  Description:
C###    IPSOLU_ITERATVE gets solution parameters for the iterative solver

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'bmg00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER NX,SOLVER,NOQUES
      CHARACTER ERROR*(*)
      LOGICAL FILEIP
!     Local Variables
      INTEGER INFO
      REAL*8 EPS
!     External functions
      REAL*8 DLAMCH
      EXTERNAL DLAMCH


      CALL ENTERS('IPSOLU_ITERATVE',*9999)

      EPS=DLAMCH('EPS')**2

C     Start reading in the values ....
      IDEFLT(1)=NYT(1,1,nx)/2
      WRITE(FORMAT,'(A,I7,A)') '($,'' The maximum number of iterations '
     '  //'is [',IDEFLT(1),']: '',I7)'
      IF(IOTYPE.EQ.3) IDATA(1)=SOLV_MAXIT(nx)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,999999,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) SOLV_MAXIT(nx)=IDATA(1)

      RDEFLT(1)=1.0d-6
      FORMAT='($,'' Enter the solver tolerance [1.0d-6]: '',D11.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=SOLV_TOL(nx)
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,EPS,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) SOLV_TOL(nx)=RDATA(1)

C     Jacobi/SOR relaxation parameter
      IF(SOLVER.EQ.SOLV_JACOBI .OR. SOLVER.EQ.SOLV_SOR) THEN
        RDEFLT(1)=0.6667D0
        FORMAT='($,'' The relaxation parameter [0.6667D0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=SOLV_OMEGA(nx)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0001D0,2.0D0,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SOLV_OMEGA(nx)=RDATA(1)

C     Krylov space methods
      ELSE IF(SOLVER.EQ.SOLV_CG .OR. SOLVER.EQ.SOLV_BICGSTAB .OR.
     '    SOLVER.EQ.SOLV_GMRES) THEN

        IDEFLT(1)=100 ! arbitrary - any better ideas
        FORMAT='($,'' The number of orthogonalisations before '
     '    //'restarting is [10]: '',I6)'
        IF(IOTYPE.EQ.3) IDATA(1)=SOLV_NRES(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,999999,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SOLV_NRES(nx)=IDATA(1)

        FORMAT='('' Specify type of preconditioning [2]: '''//
     '    '/''   (1) None'''//
     '    '/''   (2) Point Jacobi'''//
     '    '/''   (3) Jacobi Iteration'''//
     '    '/''   (4) Symmetric Succesive Over Relaxation'''//
     '    '/''   (5) Incomplete LU Decomposition(0)'''//
     '    '/''   (6) Incomplete LU Decomposition(1)'''//
     '    '/''   (7) Row scale'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=PRECON_CODE(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,7,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) PRECON_CODE(nx)=IDATA(1)

        IDEFLT(1)=2
        FORMAT='($,'' The number of preconditioner iterations per loop'
     '    //' [2]: '',I6)'
        IF(IOTYPE.EQ.3) IDATA(1)=SOLV_NPRECON(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,20,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SOLV_NPRECON(nx)=IDATA(1)

C       Hardcode the preconditioner relaxation parameter.
        SOLV_OMEGA(nx)=0.667D0

C     Special using Preconditioned CG from BOXMG
      ELSE IF(SOLVER.EQ.SOLV_PCG_BMG) THEN

        FORMAT='('' Specify type of preconditioning [2]: '''//
     '    '/''   (1) None'''//
     '    '/''   (2) Point Jacobi'''//
     '    '/''   (3) BOXMG'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=PRECON_CODE(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) PRECON_CODE(nx)=IDATA(1)

        IDEFLT(1)=1
        FORMAT='($,'' The number of preconditioner iterations per loop'
     '    //' [1]: '',I6)'
        IF(IOTYPE.EQ.3) IDATA(1)=SOLV_NPRECON(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,20,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SOLV_NPRECON(nx)=IDATA(1)

        IF(PRECON_CODE(nx).EQ.3) THEN

           IDEFLT(1)=1 !V-Cycle
           FORMAT='('' Specify type of Multigrid Cycle [1]: '''//
     '          '/''   (1) V-Cycle'''//
     '          '/''   (2) W-Cycle'''//
     '          '/$,''    '',I1)' 
           IF(IOTYPE.EQ.3) IDATA(1)=MG_CYCLE(nx) 
           CALL GINOUT(
     '          IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,2,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
           IF(IOTYPE.NE.3) MG_CYCLE(nx)=IDATA(1)

           IDEFLT(1)=1 ! Pointwise Gauss-Seidel
           FORMAT='('' Specify type of multigrid smoother [1]: '''//  
     '          '/''   (1) Pointwise Gauss-Seidel'''//
     '          '/''   (2) X-Line Relaxation (2D only)'''//
     '          '/''   (3) Y-Line Relaxation (2D only)'''//
     '          '/''   (4) XY-Line Relaxation (2D only)'''//
     '          '/''   (5) Plane Relaxation (3D only)'''//
     '          '/$,''    '',I1)'
           IF(IOTYPE.EQ.3) IDATA(1)=MG_RELAX_TYPE(nx)
           CALL GINOUT(
     '          IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '          1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,5,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
           IF(IOTYPE.NE.3) MG_RELAX_TYPE(nx)=IDATA(1)

        ENDIF

        IDEFLT(1)=0          ! Boxmg Output
        FORMAT='('' Specify Level of Output from PCG in BOXMG '
     '    //'[0]: '''//
     '       '/''   (0) No Output '''//
     '       '/''   (1) Residual norm of PCG iterations'''//
     '       '/''   (2) 1 + BOXMG Post-iterations residual norm'''//
     '       '/''   (3) 2 + BOXMG Post-Relaxation residual norm'''//
     '       '/''   (4) 3 + BOXMG Cycle Times'''//
     '       '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MG_OUTPUT(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '       FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,0,4,
     '       LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MG_OUTPUT(nx)=IDATA(1)

C     Multigrid methods
      ELSE IF(SOLVER.EQ.SOLV_AMG .OR. SOLVER.EQ.SOLV_BMG ) THEN
         
        IDEFLT(1)=1 ! V-cycles
        FORMAT='('' Specify type of Multigrid Cycle [1]: '''//
     '    '/''   (1) V-Cycle'''//
     '    '/''   (2) V*-Cycle (only AMG)'''//
     '    '/''   (3) F-Cycle'''//
     '    '/''   (4) W-Cycle'''//
     '    '/$,''    '',I1)' 
        IF(IOTYPE.EQ.3) IDATA(1)=MG_CYCLE(nx) 
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MG_CYCLE(nx)=IDATA(1)

      ENDIF

      IF(SOLVER.EQ.SOLV_BMG  ) THEN 

        IDEFLT(1)=1 ! Pointwise Gauss-Seidel
        FORMAT='('' Specify type of multigrid smoother [1]: '''//  
     '    '/''   (1) Pointwise Gauss-Seidel'''//
     '    '/''   (2) X-Line Relaxation (2D only)'''//
     '    '/''   (3) Y-Line Relaxation (2D only)'''//
     '    '/''   (4) XY-Line Relaxation (2D only)'''//
     '    '/''   (5) Plane Relaxation (3D only)'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MG_RELAX_TYPE(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MG_RELAX_TYPE(nx)=IDATA(1)

        IDEFLT(1)=1 ! Symmetric or Nonsymmetric relaxation
        FORMAT='('' Specify Symmetry of Cycle [1]: '''//
     '    '/''   (1) Symmetric Relaxation'''//
     '    '/''   (2) Nonsymmetric Relaxation'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MG_RELAX_SYM(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MG_RELAX_SYM(nx)=IDATA(1) 

        IDEFLT(1)=0 ! Boxmg Output
        FORMAT='('' Specify Level of Output from BOXMG [0]: '''//
     '    '/''   (0) No Output '''//
     '    '/''   (1) Post-iterations residual norm'''//
     '    '/''   (2) 1 + Post-Relaxation residual norm'''//
     '    '/''   (3) 2 + Cycle Times'''//
     '    '/''   (4) 3 + Setup Times and Memory Usage'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MG_OUTPUT(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,0,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MG_OUTPUT(nx)=IDATA(1) 

      ENDIF

      CALL EXITS('IPSOLU_ITERATVE')
      RETURN
 9999 CALL ERRORS('IPSOLU_ITERATVE',ERROR)
      CALL EXITS('IPSOLU_ITERATVE')
      RETURN 1
      END


