      SUBROUTINE BMG2_SymStd_SETUP_parts( 
     &                IFD, IBC, IRELAX, SO, NSO, SOR, NSOR, CI, NCI,
     &                QF, RES, NF,  !! RES is workspace to be zeroed !!
     &                ISKIP, ISTRT,  !! Temporary for periodic
     &                ABD, BBD, NCBW, NCU, IGRD, NOGm, NOG, BMG_IOFLAG,
     &                BMG_iPARMS
     &                )

C ==========================================================================
C
C   BMG2_SymStd_SETUP_parts.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG2_SymStd_SETUP_parts.f loops over all grids calling the
C   routines to construct the operator induced interpolation
C   and the coarse grid operators.  It also factors the operator
C   on the coarsest grid.
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/08/04 (JDM)
C
C ==================================================================
C   INPUT:
C ========================
C
C
C
C ==================================================================
C   OUTPUT:
C ===========================
C
C
C
C ==================================================================
C   LOCAL:
C ========================
C
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 

      INTEGER  NCI, NF, NSO, NSOR, NCBW, NCU, NOGm

      INTEGER  IBC, IFD, IGRD(NOGm,9), IRELAX, ISKIP, ISTRT, NOG
      REAL*8   ABD(NCBW*NCU), BBD(NCU), 
     &         CI(NCI), QF(NF), RES(NF), SO(NSO), SOR(NSOR)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG2_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, K, NSORv, NStncl, NStncl_CG, Nxc, Nxf, Nyf, Nyc, 
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &         INTERP

C ==========================================================================

      !
      ! Sanity check
      !
      IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,1)
         RETURN

      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,2)
         RETURN

      ENDIF

      INTERP = BMG_iPARMS(id_BMG2_INTERP)

      !
      ! Fine Grid stencil size
      !
      IF (IFD.NE.1) THEN
         NStncl=5
      ELSE
         NStncl=3
      ENDIF

      !
      ! Coarse-grid operators are always 9 points.
      !
      NStncl_CG=5

      !
      !  Number of temporary vectors in SOR
      !
      IF ( iRELAX.EQ.BMG_GS_RB_point 
     &    .OR. iRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         NSORv = 2
      ELSE 
         NSORv = 4
      ENDIF

      !
      ! Fine Grid dimensions (Nx,Ny)
      ! 
      CALL BMG2_SymStd_GET_pointers(
     &          NOG, IGRD, NOGm, Nxf, Nyf,
     &          p_U, p_SO, p_SOR, p_CI
     &          )

      !
      ! Interactively examine the fine grid stencil
      !
      IF ( BMG_IOFLAG(iBMG2_BUG_STENCIL_FG) ) THEN
         CALL BMG2_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SO), Nxf, Nyf,
     &             IFD, NStncl, NOG, NOG
     &             )
      ENDIF

      !
      ! Zero workspace
      ! 
C$OMP PARALLEL DO 
      DO i=NStncl*NF+1,NSO
         SO(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1,NSOR
         SOR(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1,NCI
         CI(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1, NCBW*NCU
         ABD(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1,NCU
         BBD(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1,NF
         RES(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C ------------------------------------------
C     Setup operators on all coarse grids.
C ------------------------------------------

      IF ( NOG.GT.1 ) THEN  ! NOG=1 => there isn't a coarse grid
         !
         ! Loop over grids
         !
         DO K = NOG, 2, -1
            !
            ! Determine the number of points in the stencil
            !
            IF (K.NE.NOG .OR. IFD.NE.1) THEN
               NStncl=5
            ELSE
               NStncl=3
            ENDIF
            !
            ! (fake) memory pointers
            !
            CALL BMG2_SymStd_GET_pointers(
     &                K, IGRD, NOGm, Nxf, Nyf,
     &                p_U, p_SO, p_SOR, p_CI
     &                )
            CALL BMG2_SymStd_GET_pointers( 
     &                K-1, IGRD, NOGm, Nxc, Nyc,
     &                p_UC, p_SOC, p_SORC, p_CIC 
     &                )
         
            If ( IBC.EQ.BMG_BCs_definite .OR. 
     &           IBC.EQ.BMG_BCs_indef_nonper ) THEN

               ! ---------------------------------------------------------
               !  Setup for the definite case
               ! ---------------------------------------------------------

               IF ( INTERP.EQ.BMG_INTERP_OI ) THEN
                  !
                  ! Operator Induced Interpolation
                  !
                  CALL BMG2_SymStd_SETUP_interp_OI(
     &                 K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                 Nxf, Nyf, Nxc, Nyc, NOG, IFD, NStncl,
     &                 IBC, IRELAX, BMG_IOFLAG, BMG_iPARMS
     &                 )
                  !
               ELSEIF ( INTERP.EQ.BMG_INTERP_GMG ) THEN
                  !
                  ! Bilinear Interpolation
                  !
                  CALL BMG2_SymStd_SETUP_interp_BI(
     &                 K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                 Nxf, Nyf, Nxc, Nyc, NOG, IFD, NStncl,
     &                 IBC, IRELAX, BMG_IOFLAG, BMG_iPARMS
     &                 )
                  !
               ELSE

                  IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
                     WRITE(*,500) 'This interp operator is invalid!'
                     WRITE(*,520) 
     &                    'HAVE: BMG_iPARMS(id_BMG2_INTERP) =', 
     &                    BMG_iPARMS(id_BMG2_INTERP) 
                  END IF  
         
                  CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,27) 
                  RETURN 

               ENDIF


               IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
                  RETURN
               ENDIF

               IF ( BMG_iPARMS(id_BMG2_CG_CONSTRUCT)
     &            .EQ.BMG_CG_CONS_explicit
     &            ) THEN
                  !
                  ! Coarse-grid operator I^{T}*L*I: explicit construction
                  !
                  CALL BMG2_SymStd_SETUP_ITLI_ex( 
     &                      K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                      Nxf, Nyf, Nxc, Nyc, NOG, IFD, NStncl,
     &                      BMG_IOFLAG, BMG_iPARMS
     &                      )

               ELSE IF ( BMG_iPARMS(id_BMG2_CG_CONSTRUCT)
     &                 .EQ.BMG_CG_CONS_block
     &                 ) THEN
                  !
                  ! Coarse-grid operator I^{T}*L*I:  block construction
                  !
                  CALL BMG2_SymStd_SETUP_ITLI_bl( 
     &                      K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                      Nxf, Nyf, Nxc, Nyc, NOG, IFD, NStncl,
     &                      BMG_IOFLAG, BMG_iPARMS
     &                      )

               ENDIF

               IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
                  RETURN
               ENDIF

               !
            ELSE
               ! ---------------------------------------------------------
               !  Setup for the indefinite and/or periodic case
               ! ---------------------------------------------------------

               !
               !  Operator Induced Interpolation and I^{T}LI
               !
               CALL BMG2_PerSymStd_SETUP_cofp(
     &                   K, K-1, SO(p_SO), SO(p_SOC),
     &                   QF(p_U), QF(p_UC),
     &                   SOR(p_SOR), SOR(p_SORC),
     &                   CI(p_CIC), Nxf, Nyf, Nxc, Nyc,
     &                   NOG, IFD, NStncl, NSORv, IBC,
     &                   IRELAX, ISTRT, ISKIP
     &                   )
               !
            ENDIF

            !
            ! Interactively examine the coarse grid stencil
            !
            IF ( BMG_IOFLAG(iBMG2_BUG_STENCIL_CG) ) THEN
               CALL BMG2_SymStd_DUMP_stencil( 
     &                   BMG_IOFLAG, SO(p_SOC), Nxc, Nyc,
     &                   IFD, NStncl_CG, K-1, NOG
     &                   )
            ENDIF

            !
            ! Interactively examine the restriction operator
            !
            IF ( BMG_IOFLAG(iBMG2_BUG_RESTRICT) ) THEN
               CALL BMG2_SymStd_DUMP_restrict( 
     &                   BMG_IOFLAG, CI(p_CIC), K, K-1, NOG, 
     &                   Nxf, Nyf, Nxc, Nyc
     &                   )
            ENDIF
            
         END DO


         ! 
         ! Relaxation 
         !
         IF ( IBC.EQ.BMG_BCs_definite .OR. 
     &        IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            CALL BMG2_SymStd_SETUP_relax(
     &                IFD, IBC, IRELAX, SO, NSO, SOR, NSOR, 
     &                IGRD, NOGm, NOG, BMG_IOFLAG, BMG_iPARMS
     &                )
            IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF
            !
         ENDIF
         !
      ENDIF


C ------------------------------------------
C     Setup the coarse grid solve
C ------------------------------------------

      IF ( NOG.EQ. 1 ) THEN 
         !
         ! Determine the number of points in the stencil
         !
         IF ( IFD.NE.1 ) THEN
            NStncl=5
         ELSE
            NStncl=3
         ENDIF

      ENDIF

      !
      !  Coarsest grid index is always 1
      !
      K = 1
      CALL BMG2_SymStd_GET_pointers(
     &          K, IGRD, NOGm, Nxc, Nyc,
     &          p_UC, p_SOC, p_SORC, p_CIC
     &          )

      !
      ! Interactively examine the stencil on the coarsest grid
      !
      IF ( BMG_IOFLAG(iBMG2_BUG_STENCIL_CG1) ) THEN
         CALL BMG2_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SOC), Nxc, Nyc,
     &             IFD, NStncl, K, NOG
     &             )
      ENDIF

      !
      ! Setup and factor the matrix for direct solve on grid 1
      !        
      IF ( IBC.EQ.BMG_BCs_definite .OR. 
     &     IBC.EQ.BMG_BCs_indef_nonper ) THEN
         !
         CALL BMG2_SymStd_SETUP_cg_LU(
     &             SO(p_SOC), Nxc, Nyc, NStncl, ABD, NCBW, NCU,
     &             BMG_IOFLAG, BMG_iPARMS
     &             )
         IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
            RETURN
         ENDIF
         !
      ELSE
         !
         CALL BMG2_PerSymStd_SETUP_cg_LU(
     &             SO(p_SOC), Nxc, Nyc, NStncl, ABD, NCBW, NCU, IBC
     &             )
         !
      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_parts.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
