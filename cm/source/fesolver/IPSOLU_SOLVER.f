      SUBROUTINE IPSOLU_SOLVER(nr,nx,NOQUES,FILEIP,ERROR,*)

C#### Subroutine: IPSOLU_SOLVER
C###  Description:
C###    Defines additional solution parameters for region nr for
C###    solvers in libsolver.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='IPSOLU_SOLVER')
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER NOQUES,nr,nx
      LOGICAL FILEIP
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER INFO,OLD_SOLVEROPTION
      LOGICAL SPARSELU,ITERATIVE
      CHARACTER SOLV_FORM*800


      CALL ENTERS(ROUTINENAME,*9999)

      FILEIP=.FALSE.  
      NOQUES=0

      SOLV_FORM='('' Specify type of linear solution procedure '
     '  //'[1]: '''//
     '  '/''   (1)  LU Decomposition'''//
     '  '/''   (2)  Single Value Decomposition'''//
     '  '/''   (3)  Least Squares'''//
     '  '/''   (4)  Cholesky Decomposition'''//
     '  '/''   (5)  Jacobi Iteration'''//
     '  '/''   (6)  Succesive Over Relaxation'''//
     '  '/''   (7)  Incomplete LU Decomposition(0)'''//
     '  '/''   (8)  Incomplete LU Decomposition(1)'''//
     '  '/''   (9)  Conjugate Gradient'''//
     '  '/''   (10) Biconjugate Gradient Stabilised'''//
     '  '/''   (11) Generalised Minimum Residual'''//
     '  '/''   (12) Black Box Multigrid (BOXMG)'''//
     '  '/''   (13) Preconditioned CG from BOXMG'''//
     '  '/''   (14) Algebraic Multigrid (AMG1R6)'''//
     '  '/$,''    '',I2)'
      OLD_SOLVEROPTION=SOLVEROPTION(nx)
      IF(OLD_SOLVEROPTION.GT.SOLV_SPARSELU) SOLVEROPTION(nx)=SOLV_LU
      IF(IOTYPE.EQ.3) IDATA(1)=SOLVEROPTION(nx)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '  SOLV_FORM,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IONE,1,14,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) SOLVEROPTION(nx)=IDATA(1)
      ITERATIVE=IDATA(1).GE.5

      FORMAT='($,'' Do you want the solution matrix stored as '
     '  //'a sparse matrix [Y]? '',A)'
      IF(ITYP4(nr,nx).EQ.2) THEN !BEM
        CDEFLT(1)='N'
      ELSE !all other cases
        CDEFLT(1)='Y'
      ENDIF
      SPARSELU=.FALSE.
      IF(IOTYPE.EQ.3) THEN
        IF(SPARSEGKK(nx).GT.0) THEN
          ADATA(1)='Y'
        ELSE
          ADATA(1)='N'
        ENDIF
        SPARSELU=(SOLVEROPTION(nx).EQ.1).AND.(SPARSEGKK(nx).GT.0)
      ENDIF
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,AYES,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

      IF(IOTYPE.NE.3) THEN
        IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
          IF(SOLVEROPTION(nx).EQ.SOLV_LU) THEN !LU Factorisation
            SPARSELU=.TRUE.
          ELSE IF(SOLVEROPTION(nx).EQ.SOLV_SVD) THEN !SVD Factorisation
            ERROR='>>Must be non-sparse for SVD Factorisation'
            GOTO 9999
          ELSE IF(SOLVEROPTION(nx).EQ.SOLV_LSQ) THEN !least squares
            ERROR='>>Must be non-sparse for least squares'
            GOTO 9999
          ELSE IF(SOLVEROPTION(nx).EQ.SOLV_CHOLESKY) THEN !least squares
            ERROR='>>Must be non-sparse for Cholesky'
            GOTO 9999
          ENDIF

          SPARSEGKK(nx)=1
          CALL ASSERT(USE_SPARSE.NE.0,
     '      '>>Set USE_SPARSE to 1 to use sparse matrices',
     '      ERROR,*9999)
        ELSE
           IF(SOLVEROPTION(nx).EQ.SOLV_AMG) THEN !AMG1R6
             ERROR='>>Must be sparse for calling AMG1R6.'
             GOTO 9999
          ELSEIF(SOLVEROPTION(nx).EQ.SOLV_BMG) THEN !BOXMG
             ERROR='>>Must be sparse for calling BOXMG'
             GOTO 9999
          ELSEIF(SOLVEROPTION(nx).EQ.SOLV_PCG_BMG) THEN !BOXMG
             ERROR='>>Must be sparse for calling PCG-BOXMG'
             GOTO 9999
          ENDIF

          SPARSEGKK(nx)=0
        ENDIF
      ENDIF !iotype

C         Options for the sparse LU solvers
      IF(SPARSELU) THEN
        SPARSEGKK(nx)=1
        IDEFLT(1)=2
        FORMAT='('' Specify the LU solver [2]: '''//
     '    '/''   (1) SuperLU'''//
     '    '/''   (2) Umfpack'''//
     '    '/$,''    '',I1)'
C DAH 7-01-2003 Fixing the ;w command option
        IF(IOTYPE.EQ.3) THEN
          SOLVEROPTION(nx)=OLD_SOLVEROPTION
          IDATA(1)=SOLVEROPTION(nx)-SOLV_SPARSELU
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SOLVEROPTION(nx)=IDATA(1)+SOLV_SPARSELU

        IF(SOLVEROPTION(nx).EQ.SOLV_SUPERLU) THEN
          CALL IPSOLU_SUPERLU(SOLV_SLU_PARAM(1,nx),NOQUES,
     '      FILEIP,ERROR,*9999)
        ELSE IF(SOLVEROPTION(nx).EQ.SOLV_UMFPACK4) THEN
          CALL IPSOLU_UMFPACK(SOLV_UMF_PARAM(1,nx),NOQUES,
     '      FILEIP,ERROR,*9999)
        ENDIF

C         Options for the iterative solvers
      ELSE IF(ITERATIVE) THEN
        CALL IPSOLU_ITERATVE(nx,SOLVEROPTION(nx),NOQUES,
     '    FILEIP,ERROR,*9999)
      ENDIF

C         Options to fix broken singular systems (thanks Carey).
C         Switched off to make sure that dodgy systems don't get ignored
      SOLV_FIXIT(nx)=.FALSE.
C         FORMAT='($,'' Do you want to fix singular systems [N]? '',A)'
C         CDEFLT(1)='N'
C         IF(IOTYPE.EQ.3) THEN
C           IF(SOLV_FIXIT(NX)) THEN
C             ADATA(1)='Y'
C           ELSE
C             ADATA(1)='N'
C           ENDIF
C         ENDIF
C         CALL GINOUT(IOTYPE,1,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
C    '      1,ADATA,AYES,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,IMAX,
C    '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C         IF(IOTYPE.NE.3) THEN
C           IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
C             SOLV_FIXIT(NX)=.TRUE.
C           ELSE
C             SOLV_FIXIT(NX)=.FALSE.
C           ENDIF
C         ENDIF

C SGM 26 Oct 2000 grid-based Finite element also
C MLT 29Nov02 grid finite volume also
      IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6
     '  .OR.ITYP4(nr,nx).EQ.7) THEN
        IF(SOLVEROPTION(nx).EQ.1.AND.BEMOVERDETERMINED(nr,nx)) THEN
          CALL ASSERT(SPARSEGKK(nx).EQ.0,
     '      '>>Cannot use sparse LU with an overdetermined system',
     '      ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


