      SUBROUTINE BMG2_SymStd_SOLVE_pcg(
     &                Nx, Ny, 
     &                BMG_PCG_iPARMS, BMG_PCG_rPARMS, BMG_PCG_OUTPUT,
     &                BMG_PCG_RES, MAX_PCG_ITERS, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, NFm, SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_IJK_MAP
     &                )


C ==========================================================================
C
C***BEGIN PROLOGUE  BMG2_SymStd_SOLVE_pcg
C
C***PURPOSE  
C 
C     This subroutine is a preconditioned conjugate gradient code
C     which exclusively uses the black box multigrid code, BOXMG,
C     for solving the precoditioned system.  It is used for discretizations
C     of second order elliptic partial differential equations that 
C     generate, at most, a 9-point stencil on a logically rectangular
C     grid.  It may be applied to other, similarly structured, problems.
C
C***AUTHOR  AUSTIN, TRAVIS
C           GROUP T-7, MAIL STOP B284
C           LOS ALAMOS NATIONAL LABORATORY
C           LOS ALAMOS, NEW MEXICO 87545
C           E-MAIL: taustin@lanl.gov
C
C***LAWYER STUFF (This should suffice for now, but will be updated).
C
C     Copyright, 2000. The Regents of the University of California.
C     This software was produced under a U.S. Government contract
C     (W-7405-ENG-36) by Los Alamos National Laboratory, which is
C     operated by the University of California for the U.S.
C     Department of Energy. The U.S. Government is licensed to use,
C     reproduce, and distribute this software. Permission is
C     granted to the public to copy and use this software without
C     charge, provided that this Notice and any statement of
C     authorship are reproduced on all copies. Neither the
C     Government nor the University makes any warranty, expressed
C     or implied, or assumes any liability or responsibility for
C     the use of this software.
C
C***DESCRIPTION
C
C     BMG2_SymStd_SOLVE_pcg is a preconditioned conjugate gradient 
C     code using BOXMG as a black box preconditioner.  See the 
C     BMG2_SymStd_SOLVE_boxmg code or references cited below for 
C     details on the BOXMG algorithm. 
C
C***PARAMETERS
C
C***INPUT
C
C   Nx      x-dimension of the grid, excluding fictitious points
C   Ny      y-dimension of the grid, excluding fictitious points
C
C   BMG_PCG_iPARMS      Pcg integer parameters
C
C        BMG_PCG_iPARMS(2) --  Maximum number of symmetric V-cycles 
C                              to perform as preconditioner.
C
C        BMG_PCG_iPARMS(3) --  Type of stop test to use for the pcg 
C                              algorithm.  (See PCG_parameters.h)
C
C   BMG_PCG_rPARMS      Pcg real parameters
C
C        BMG_PCG_rPARMS(1) --  Tolerance used to stop PCG iterations.
C
C   BMG_PCG_NORMits     Array of size PCG_iPARMS(1) containing l2
C                       norm of residual after each pcg iteration.
C
C   BMG_iPARMS      See BMG2_boxmg code.
C
C   BMG_rPARMS      See BMG2_boxmg code.
C
C   BMG_IOFLAG      See BMG2_boxmg code.
C
C   NSOm     See INPUT/OUTPUT
C
C   NSORm    See WORK ARRAYS.
C
C   NCIm     See WORK ARRAYS.
C
C   NOGm     Maximum number of grids that can be supported by the pointer
C            array IGRD. It should be set in calling program.
C
C   NXYc     Limit for number of points on coarsest grid.  Coarsification
C            occurs until the number of x or y unknowns minus one is 
C            less than or equal to NXYc.  NXYc should be 3 or greater or 
C            the code will abort.  This parameter is included because 
C            efficiency is a function of machine.   On vector machines, 
C            it frequently pays to take NXYc larger than on scalar machines.
C
C***INPUT/OUTPUT
C
C   SO       User defined real array that contains the coefficient
C            matrix. See above description for format.
C            If core is not set to zero, SO(I)=0., I=1,NSOm should be
C            performed before PUTF is called to set SO.
C
C   NSOm     Dimension of SO, the calling program should have computed
C            this with a call to BMG2_SymStd_SETUP_space.  Clearly,
C            
C            IFD .EQ. 1  =>    NSOm = 3*NFm + 2*NCm 
C            IFD .NE. 1  =>    NSOm = 5*NFm 
C  
C            BOXMG verifies that NSOm is large enough.
C
C   QF       The user defined array that contains the right hand side.
C            It is set to dimension NFm.
C
C   Q        The user defined array that contains The solution vector.
C            It is set to dimension NFm.
C
C***WORK ARRAYS
C
C   ABD      User declared two dimensional real array, which
C            is used to store the coefficient matrix for the coarsest
C            grid. It is then used by the LINPACK routine DPBSL. It
C            should be dimensioned to (NCBWm,NCUm).
C
C   BBD      User declared real array of dimension NCUm for use in the
C            the LINPACK routine DPBSL.
C
C   NCBWm    Maximum first subscript of ABD, which needs to be > or =
C            number of grid points plus 2 in x direction on coarsest
C            grid. BOXMG checks to see if it is large enough.
C
C   NCUm     Maximum second subscript of ABD, which needs to be > or =
C            (x-dimension)*(y-dimension) on coarsest grid. BOXMG
C            checks to see if it is large enough.
C
C   IGRD     A work array. IGRD(.,1),...,IGRD(.,6)
C            are used to store starting locations for the arrays on each
C            grid. IGRD(.,7),...,IGRD(.,9) are used as scratch integer
C            arrays needed by BOXMG. IGRD should be dimensioned to
C            IGRD(NOGm,9) in the calling program.
C
C   SOR      SOR is a real array, which is set to dimension NSORm.
C            It is used as workspace for the tridiagonal factorization 
C            and solves used in the case of line relaxation.  
C
C   NSORm    Dimension of SOR, set in calling program.  Actual requirements
C            need to be cleaned up!
C
C            IRELAX .EQ. 1   =>   NSORm = 0
C            IRELAX .EQ. 2   =>   NSORm = 2*NFm 
C            IRELAX .EQ. 3   =>   NSORm = 4*NFm  (this should be fixed)
C            IRELAX .EQ. 4   =>   NSORm = 4*NFm
C
C            BOXMG checks to see if NSORm is large enough.
C
C   CI       CI is a real array which should be dimensioned
C            to (NCIm) in calling program. It is used to contain
C            the interpolation coefficients.
C
C   NCIm     Dimension of CI, set in calling program.  NCIm = 8*NCm is 
C            enough.  BOXMG checks to see if NCIm is large enough.
C
C***REFERENCES  Dendy, J. E. Jr., "Black Box Multigrid", Journal of 
C                 Computational Physics, Vol. 48, pp. 366-386, 1982
C
C               Dendy, J. E. Jr., "Black Box Multigrid for Nonsymmetric
C                 Problems", Applied Mathematics and Computation, Vol. 13,
C                 pp. 261-283, 1983
C
C               Dendy, J. E. Jr., "Black Box Multigrid for Systems",
C                 Applied Mathematics and Computation, Vol. 19,
C                 pp. 57-74, 1986
C
C               Dendy, J. E. Jr., "Two Multigrid Methods for Three
C                 Dimensional Problems with Discountinuos and Anisotropic
C                 Coefficients", SIAM Journal of Scientific and Satatistical
C                 Computing, Vol. 8, No. 2, September 1987
C
C               Dendy, J. E. Jr., "Black Box Multigrid for Periodic and
C                 Singular Problems", Applied Mathematics and Computation,
C                 Vol. 25, pp. 1-10, 1988
C
C***ROUTINES CALLED  (Someday I'll redo this list)
C
C***REVISION HISTORY  (YYYY/MM/DD)
C
C   2003/07/01  (TMA)
C   - written
C
C***END PROLOGUE  BMG2_PCG
C
C    ************************warning************************************
C   this is an experimental code, and no guarantees are made as to the *
C   validity of the computed solution.                                 *
C    *******************************************************************
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_PCG_parameters.h'
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C -------------------------------
C     Generic Argument Declarations
C     
      INTEGER  MAX_PCG_ITERS, Nx, Ny, NFm
      REAL*8   Q(NFm), QF(NFm)

C ------------------------------------
C     PCG Specific Argument Declarations
C
      INTEGER  BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*8   BMG_PCG_rPARMS(NBMG_PCG_rPARMS),
     &         BMG_PCG_Res(MAX_PCG_ITERS) 

C ------------------------------------
C     BOXMG Specific Argument Declarations
C
      INTEGER  NBMG_iWORKm, NBMG_rWORKm, NSOm, NOG, NOGm
      
      INTEGER BMG_iPARMS(NBMG_iPARMS), BMG_iWORK(NBMG_iWORKm),
     &        BMG_pWORK(NBMG_pWORK), BMG_IJK_MAP(*)
      REAL*8  BMG_rPARMS(NBMG_rPARMS), BMG_rWORK(NBMG_rWORKm),
     &        SO(NSOm)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG), BMG_PCG_OUTPUT

C ----------------------------
C     Local Declarations
C
      INTEGER  J, II, JJ
      INTEGER  NMGCYCLES, STOP_TEST, IFD, IBC
      INTEGER  IRELAX, IRELAX_SYM, NStncl, UPDOWN
      INTEGER  PRECON
      INTEGER  NF, NC, NCI, NSO, NSOR, NCBW, NCU, p_P, p_R, p_Z
      
      REAL*8   ANORM, tol, derr,  RES_L2_0, RES_L2, RHS_L2, SCALED_TOL
      REAL*8   malpha, alpha, delta0, delta1, beta, T, T1, T2

      REAL*8   rTINY
      PARAMETER( rTINY=1E-10 )

C =========================================================================


      T  = rZERO
      T1 = rZERO
      T2 = rZERO

C!!
C!!   We need to add some sanity checks here  !!!
C!!   Maybe a subroutine to check all parameter values??
C!!   

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      !  PCG parameters
      !
      NMGCYCLES    = BMG_PCG_iPARMS(id_BMG_PCG_NMG_CYCLES)
      STOP_TEST    = BMG_PCG_iPARMS(id_BMG_PCG_STOP_TEST)
      !
      PRECON       = BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
      !
      tol          = BMG_PCG_rPARMS(id_BMG_PCG_STOP_TOL)
      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_DIM_NCU)

      !
      !  BoxMG cycle parameters
      !
      IFD          = BMG_iPARMS(id_BMG2_STENCIL)
      IBC          = BMG_iPARMS(id_BMG2_BC) 
      IRELAX       = BMG_iPARMS(id_BMG2_RELAX)
      IRELAX_SYM   = BMG_iPARMS(id_BMG2_RELAX_SYM )

      !
      !  Local pointers
      !
      p_P   = BMG_pWORK(ip_BMG_PCG_P)
      p_R   = BMG_pWORK(ip_BMG_PCG_R)
      p_Z   = BMG_pWORK(ip_BMG_PCG_Z)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  PCG ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      II = Nx+2
      JJ = Ny+2

C -------------- peform zeroth iteration as special case -------------

      IF (IFD.NE.1) THEN
         NStncl=5
      ELSE
         NStncl=3
      ENDIF

C -------------------------
C     Perform Setup:
C -------------------------

      CALL BMG_timer(T1)

      IF ( BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
     &     .EQ. BMG_PCG_BMG_SETUP_all .OR.
     &     BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
     &     .EQ. BMG_PCG_BMG_SETUP_only ) THEN

         CALL BMG_timer(T1)

         IF ( PRECON.eq.BMG_PCG_PRECON_BMG ) THEN
            
            BMG_iPARMS(id_BMG2_SETUP) = BMG_SETUP_only
            
            CALL BMG2_SymStd_PRECON_boxmg( 
     &           Nx, Ny, 
     &           BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &           BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &           SO, NSOm, NOGm,
     &           BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &           BMG_rWORK, NBMG_rWORKm
     &           )
            IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN 
               RETURN  
            ENDIF 
            BMG_iPARMS(id_BMG2_SETUP) = BMG_SETUP_none  

         ENDIF 

         CALL BMG_timer(T2)

         IF( BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP) 
     &        .EQ. BMG_PCG_BMG_SETUP_only ) THEN 
            
            BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)  
     &           = BMG_PCG_BMG_SETUP_none 
            RETURN
            
         ELSE 
            
            BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)  
     &           = BMG_PCG_BMG_SETUP_none 
         ENDIF 
         
      ENDIF 

C --------------------------------------- 
C     Calculate the norm of rhs 
C ---------------------------------------

      CALL BMG2_SymStd_UTILS_norm_l2(
     &          QF, II, JJ, RHS_L2
     &          )

      IF( RHS_L2 .le. rTINY ) THEN
         IF ( BMG_IOFLAG(iBMG2_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) '** NB: Zero RHS implies solution equals zero.'
         ENDIF     
         
         CALL BMG2_SymStd_UTILS_rV_zero( Q, Nx, Ny )
         RETURN
      ENDIF


      IF ( STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2 .OR.
     &     STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2 ) THEN
      
         IF( RHS_L2.LE.1E-4 ) THEN
            SCALED_TOL = rTINY
         ELSE
            SCALED_TOL = tol*RHS_L2
         ENDIF

      ELSE
         
         SCALED_TOL = tol

      ENDIF

C --------------------------------------- 
C     Calculate the residual R = F - A*x
C ---------------------------------------

      J=0

      CALL BMG2_SymStd_residual(
     &          NOG, SO, QF, Q, BMG_rWORK(p_R), II, JJ, 
     &          NOG, IFD, NStncl
     &          )

      CALL BMG2_SymStd_UTILS_norm_l2(
     &          BMG_rWORK(p_R), II, JJ, RES_L2_0 
     &          )

      IF ( RES_L2_0 .LE. SCALED_TOL ) THEN
         IF ( BMG_IOFLAG(iBMG2_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) 'NOTE: Zero initial residual vector.'
            WRITE(*,505) '      Initial guess is the solution.'
         ENDIF
         RETURN
      ENDIF

      IF( BMG_PCG_OUTPUT ) THEN
         WRITE(*,*) 
         WRITE(*,*) 
         WRITE(*,270) ' **** Initial Residual (L2 norm)  = ', RES_L2_0
         WRITE(*,*) 
         
         WRITE(*,*) '   ================================= '
         WRITE(*,*) '     Iteration       Stopping Test   '
         WRITE(*,*) '   ================================= '
      ENDIF

      J=1

C -----------------------------------------------------------
C     Compute Z = M^{-1} R where 
C     M^{-1} is BOXMG symmetric V-cycle
C      BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C                 = BMG_PCG_PRECON_NONE  -> M = I (IDENTITY)
C      BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C                 = BMG_PCG_PRECON_DIAG  -> M = diag(A)
C      BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C                 = BMG_PCG_PRECON_BMG   -> M = N BMG NCYCLES 
C -----------------------------------------------------------

      IF ( PRECON.eq.BMG_PCG_PRECON_NONE ) THEN
         !
         CALL BMG2_SymStd_UTILS_dcopy(
     &             BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ 
     &             )
         !
      ELSEIF ( PRECON.eq.BMG_PCG_PRECON_DIAG ) THEN
         !
         CALL BMG2_SymStd_PRECON_diag(
     &             SO, BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ, NStncl
     &              )
         !
      ELSEIF ( PRECON.eq.BMG_PCG_PRECON_BMG ) THEN
         !
         CALL BMG2_SymStd_UTILS_rV_zero(
     &        BMG_rWORK(p_Z), II, JJ 
     &        )
         !
         CALL BMG2_SymStd_PRECON_boxmg( 
     &             Nx, Ny, 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &             SO, NSOm, NOGm,
     &             BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &             BMG_rWORK, NBMG_rWORKm
     &             )
         !
      ENDIF

C ----------------------------------
C     Copy Z into P.
C ----------------------------------

      CALL BMG2_SymStd_UTILS_dcopy(
     &          BMG_rWORK(p_Z), BMG_rWORK(p_P), II, JJ
     &          )

      CALL BMG2_SymStd_UTILS_dot_l2(
     &          BMG_rWORK(p_Z), BMG_rWORK(p_R), II, JJ, delta0
     &          )

C ===============================
C     Start main loop.
C ================================

      DO J=1, MAX_PCG_ITERS

C ---------------------------------------
C
C        Calculate delta0/(Apj,pj) by
C           (2) Zj     = Apj
C           (3) delta1 = (Zj,Pj)
C           (4) alpha  = delta0/delta1
C
C ---------------------------------------

         CALL BMG2_SymStd_UTILS_matvec(
     &             NOG, SO, BMG_rWORK(p_Z), BMG_rWORK(p_P),
     &             II, JJ, NOG, IFD, NStncl
     &             )
         CALL BMG2_SymStd_UTILS_dot_l2(
     &             BMG_rWORK(p_Z), BMG_rWORK(p_P), II, JJ, delta1
     &             )

         alpha = delta0/delta1

C ----------------------------------
C        Calculate q <- q + alpha*p
C ----------------------------------

         CALL BMG2_SymStd_UTILS_daxpy( 
     &             alpha, BMG_rWORK(p_P), Q, II, JJ
     &             )
         
C ----------------------------------
C        Calculate r <- r - alpha*z
C ----------------------------------
         
         malpha = -1.0*alpha

         CALL BMG2_SymStd_UTILS_daxpy(
     &             malpha, BMG_rWORK(p_Z), BMG_rWORK(p_R), II, JJ
     &             )

         IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2 .OR.
     &        STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2) THEN
            !
            CALL BMG2_SymStd_UTILS_norm_l2(
     &                BMG_rWORK(p_R), II, JJ, RES_L2 
     &                )
            !
            BMG_PCG_Res(J) = RES_L2
            !
            derr = BMG_PCG_Res(J)
            !
            IF( BMG_PCG_OUTPUT ) THEN
               WRITE(*,400) J, derr
            ENDIF

            IF ( derr.LT.SCALED_TOL ) GO TO 200 ! Jump to post processing and return

         ENDIF



C -----------------------------------------------------------
C        Compute Z = M^{-1} R where 
C        M^{-1} is BOXMG symmetric V-cycle
C
C        BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C               = BMG_PCG_PRECON_NONE  -> M = I (IDENTITY)
C        BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C               = BMG_PCG_PRECON_DIAG  -> M = diag(A)
C        BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C               = BMG_PCG_PRECON_BMG   -> M = N BMG NCYCLES
C 
C -----------------------------------------------------------

         IF ( PRECON.eq.BMG_PCG_PRECON_NONE ) then
            !
            CALL BMG2_SymStd_UTILS_dcopy(
     &                BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ
     &                )
            !
         ELSEIF ( PRECON.eq.BMG_PCG_PRECON_DIAG ) THEN
            !
            CALL BMG2_SymStd_PRECON_diag(
     &                SO, BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ, NStncl
     &                )
            !
         ELSEIF ( PRECON.eq.BMG_PCG_PRECON_BMG ) THEN
            !
            CALL BMG2_SymStd_UTILS_rV_zero(
     &                BMG_rWORK(p_Z), II, JJ 
     &                )
            !
            CALL BMG2_SymStd_PRECON_boxmg( 
     &                Nx, Ny, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm
     &                ) 
            !
         ENDIF

C ------------------------------------------------------
C        Calculate delta1 = <R,Z> = <R,inv(M)R>
C ------------------------------------------------------

         CALL BMG2_SymStd_UTILS_dot_l2(
     &             BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ, delta1
     &             )
      
    
C ------------------------------------------------------
C        Check convergence.  If small, jump to the end!
C ------------------------------------------------------

         IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2 .OR.
     &        STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2) THEN
            !
            BMG_PCG_Res(J) = sqrt(delta1)
            !
            derr = BMG_PCG_Res(J)
            !
            IF( BMG_PCG_OUTPUT ) THEN
               WRITE(*,400) J, derr
            ENDIF

            IF ( derr.LT.SCALED_TOL ) GO TO 200 ! Jump to post processing and return

         ENDIF
     
C ----------------------------------------------
C        Calculate beta. 
C ---------------------------------------
     
         beta   = delta1/delta0
         delta0 = delta1

C ----------------------------------
C        Calculate p <- z + beta*p
C ----------------------------------

         CALL BMG2_SymStd_UTILS_dxpby(
     &             beta, BMG_rWORK(p_Z), BMG_rWORK(p_P), II, JJ
     &             )

         !
      ENDDO

 200  CONTINUE      ! <<<<<<<<<<<< OUTSIDE PCG LOOP

      BMG_PCG_rPARMS(id_BMG_PCG_FINAL_TOL) = derr
      BMG_PCG_iPARMS(id_BMG_PCG_NUM_ITERS) = J

      !
      !  Post Processing could be added here
      !

      IF( BMG_PCG_OUTPUT ) THEN
         WRITE(*,*) '   ================================= '
      ENDIF

C =======================================================================

 270  FORMAT (2X,A,1X,E14.7,/)
 400  FORMAT ( 8X,I4,9X,1P,E14.7 )
 505  FORMAT (/,5X,A,/)

C ============================

      RETURN
      END
