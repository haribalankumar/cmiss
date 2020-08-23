      SUBROUTINE BMG3_SymStd_SETUP_parts( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                ABD, BBD, NABD1, NABD2, IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_parts.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_parts.f loops over all grids calling the
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

C ---------------------------
C    Argument Declarations:
C
      INTEGER  NABD1, NABD2, NBMG_iWORK_PL, NBMG_rWORK_PL, 
     &         NCI, NOGm, NSO, NSOR

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IGRD(NOGm*NBMG_pIGRD), NOG
      REAL*8   abd(nabd1*nabd2), bbd(nabd2),
     &         BMG_rPARMS(NBMG_rPARMS), 
     &         BMG_rWORK_PL(NBMG_rWORK_PL), CI(NCI), SO(NSO),
     &         SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C    Local Declarations:
C
      INTEGER  i, IFD, IBC, IRELAX, INTERP, K, NStncl, NStncl_CG,
     &         Nx, Nxc, Ny, Nyc, Nz, Nzc, 
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &         p_CGTEMP_yo, p_CGTEMP_zo

C ==========================================================================

      !
      ! Unpack
      !
      IFD    = BMG_iPARMS(id_BMG3_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_BC)
      IRELAX = BMG_iPARMS(id_BMG3_RELAX)
      INTERP = BMG_iPARMS(id_BMG3_INTERP)

      !
      ! Coarse-grid operators are always 27 points.
      !
      NStncl_CG=14

      !
      ! Sanity check: NOG
      !
      IF ( NOG.EQ.0 ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,2)
         RETURN
         !
      ELSE IF ( NOG.LT.0 ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,3)
         RETURN
         !
      ENDIF

      !
      ! Sanity check: CG type
      !
      IF ( 
     &   BMG_iPARMS(id_BMG3_CG_TYPE).NE.BMG_CG_ITLI_IzIyIx
     &   .AND.
     &   BMG_iPARMS(id_BMG3_CG_TYPE).NE.BMG_CG_ITLI 
     &   ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The Coarse-Grid Operator Type is invalid!'
            WRITE(*,520) 
     &           'HAVE: BMG_iPARMS(id_BMG3_CG_TYPE) =',
     &           BMG_iPARMS(id_BMG3_CG_TYPE)
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,24)
         RETURN
         !
      ENDIF

      !
      ! Sanity check: CG construction method
      !
      IF ( 
     &     BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &        .NE.BMG_CG_CONS_explicit 
     & 
     &   .AND.
     & 
     &     BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &        .NE.BMG_CG_CONS_block
     &
     &   ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The CG-operator construction is invalid!'
            WRITE(*,520)
     &           'HAVE: BMG_iPARMS(id_BMG3_CG_CONSTRUCT) = ',
     &           BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,25)
         RETURN
         !
      ENDIF

      !
      ! Fine Grid stencil size
      !
      IF ( IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      ! Fine Grid dimensions (Nx,Ny,Nz)
      ! 
      CALL BMG3_SymStd_GET_pointers(
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )

      !
      ! Interactively examine the fine grid stencil
      !
      IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_FG) ) THEN
         CALL BMG3_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SO), Nx, Ny, Nz, 
     &             IFD, NStncl, NOG, NOG
     &             )
      ENDIF

      !
      ! Zero workspace
      ! 
C$OMP PARALLEL DO
      DO i=NStncl*Nx*Ny*Nz+1, NSO
         SO(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1, NSOR
         SOR(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO i=1, NCI
         CI(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO
      
C$OMP PARALLEL DO
      DO i=1, NABD1*NABD2
         ABD(i)=rZERO
      ENDDO
C$OMP END PARALLEL DO
      
C$OMP PARALLEL DO
      DO i=1, NABD2
         BBD(i) = rZERO
      ENDDO
C$OMP END PARALLEL DO

      !
      !   Zero out temporary workspace for operator construction
      !
      DO i=NBMG_iWORK_PL_ptrs+1, NBMG_iWORK_PL
         BMG_iWORK_PL(i) = iZERO
      ENDDO
      
      DO i=1, NBMG_rWORK_PL
         BMG_rWORK_PL(i) = rZERO
      ENDDO

C ------------------------------------------
C     Setup operators on all coarse grids.
C ------------------------------------------

      IF ( NOG.GT.1 ) THEN  ! NOG=1 => there isn't a coarse grid


C         IF ( 
C     &      BMG_iPARMS(id_BMG3_CG_TYPE).EQ.BMG_CG_ITLI_IzIyIx
C     &      ) THEN

            !
            ! This approximation is only formed explicitly, so 
            ! we are ignoring the BMG_iPARMS(id_BMG3_CG_CONSTRUCT).
            ! A warning or flag should be set earlier.
            !

            !
            ! Loop over grids
            !
C            DO K = NOG, 2, -1

               !
               ! Determine the number of points in the stencil
               !
C               IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
C                  NStncl=14
C               ELSE
C                  NStncl=4
C               ENDIF

C               CALL BMG3_SymStd_GET_pointers(
C     &                   K, IGRD, NOGm,
C     &                   p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
C     &                   )
C               CALL BMG3_SymStd_GET_pointers( 
C     &                   K-1, IGRD, NOGm,
C     &                   p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
C     &                   )

               !
               ! Temporary interpolation based on tensor product 
               ! approach is used to compute cg.
               !
C               p_CGTEMP_zo=BMG_iWORK_PL(ip_BMG_CGTEMP_zo)
C               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_CGTEMP_yo)
C               CALL BMG3_SymStd_SETUP_ITLI_Izyx( 
C     &                   K, K-1,
C     &                   SO(p_SO), SO(p_SOC), SOR(p_SOR), CI(p_CIC),
C     &                   Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG,
C     &                   ifd, NStncl, irelax,
C     &                   BMG_rWORK_PL(p_CGTEMP_zo),
C     &                   BMG_rWORK_PL(p_CGTEMP_yo)
C     &                   ) 

C               !
C               ! Interactively examine the coarse grid stencil
C               !
C               IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) ) THEN
C                  CALL BMG3_SymStd_DUMP_stencil( 
C     &                      BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
C     &                      IFD, NStncl_CG, K-1, NOG
C     &                      )
C               ENDIF


C               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_CGTEMP_yo)
C               IF ( INTERP.EQ.BMG_INTERP_OI ) THEN
C                  !
C                  ! Operator Induced Interpolation
C                  !
C                  CALL BMG3_SymStd_SETUP_interp_OI(
C     &                 K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
C     &                 Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
C     &                 ifd, NStncl, irelax, 
C     &                 BMG_rWORK_PL(p_CGTEMP_yo)
C     &                 )
C               ELSEIF ( INTERP.EQ.BMG_INTERP_GMG ) THEN
C                  !
C                  ! Trilinear Interpolation
C                  !   
C                  CALL BMG3_SymStd_SETUP_interp_TI(
C     &                 K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
C     &    BMG3_SymStd_SETUP_ITLI_Izyx             Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
C     &                 ifd, NStncl, irelax, BMG_rWORK_PL(p_CGTEMP_yo),
C     &                 BMG_IOFLAG, BMG_iPARMS
C     &                 )
C               ELSE  
C
C                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
C                     WRITE(*,500) 'This interp operator is invalid!'
C                     WRITE(*,520) 
C     &                    'HAVE: BMG_iPARMS(id_BMG3_INTERP) =', 
C     &                    BMG_iPARMS(id_BMG3_INTERP) 
C                  END IF  
         
C                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,27) 
C                  RETURN 

C               ENDIF 

C               ! 
C               ! Interactively examine the restriction operator 
C               ! 
C               IF ( BMG_IOFLAG(iBMG3_BUG_RESTRICT) ) THEN 
C                  CALL BMG3_SymStd_DUMP_restrict( 
C     &                      BMG_IOFLAG, CI(p_CIC), K, K-1,  
C     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc  
C     &                      ) 
C               ENDIF

C            ENDDO 
 
C         ELSE IF (  
C     &             BMG_iPARMS(id_BMG3_CG_TYPE)
C     &                .EQ.BMG_CG_ITLI 
C     & 
C     &             .AND. 
C     &             
C     &             BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
C     &                .EQ.BMG_CG_CONS_explicit 
C     &
C     &           ) THEN

            ! ---------------------------------------------------------
            !  Explicit construction of the "honest" I^{T}LI operator.
            ! ---------------------------------------------------------

            !
            ! Loop over grids
            !
C            DO K = NOG, 2, -1

               !
               ! Determine the number of points in the stencil
               !
C               IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
C                  NStncl=14
C               ELSE
C                  NStncl=4
C               ENDIF

C               CALL BMG3_SymStd_GET_pointers(
C     &                   K, IGRD, NOGm,
C     &                   p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
C     &                   )
C               CALL BMG3_SymStd_GET_pointers( 
C     &                   K-1, IGRD, NOGm,
C     &                   p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
C     &                   )


C               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_CGTEMP_yo)
C               IF ( INTERP.EQ.BMG_INTERP_OI ) THEN
                  !
                  ! Operator Induced Interpolation
                  !
C                  CALL BMG3_SymStd_SETUP_interp_OI(
C     &                 K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
C     &                 Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
C     &                 ifd, NStncl, irelax, 
C     &                 BMG_rWORK_PL(p_CGTEMP_yo)
C     &                 )
C               ELSEIF ( INTERP.EQ.BMG_INTERP_GMG ) THEN
C                  !
C                  ! Trilinear Interpolation
C                  !                  
C                  CALL BMG3_SymStd_SETUP_interp_TI(
C     &                 K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
C     &                 Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
C     &                 ifd, NStncl, irelax, BMG_rWORK_PL(p_CGTEMP_yo),
C     &                 BMG_IOFLAG, BMG_iPARMS
C     &                 )
C               ELSE  

C                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
C                     WRITE(*,500) 'This interp operator is invalid!'
C                     WRITE(*,520) 
C     &                    'HAVE: BMG_iPARMS(id_BMG3_INTERP) =', 
C     &                    BMG_iPARMS(id_BMG3_INTERP) 
C                  END IF  
         
C                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,27) 
C                  RETURN 

C               ENDIF 

               !
               ! Interactively examine the restriction operator
               !
C               IF ( BMG_IOFLAG(iBMG3_BUG_RESTRICT) ) THEN
C                  CALL BMG3_SymStd_DUMP_restrict( 
C     &                      BMG_IOFLAG, CI(p_CIC), K, K-1, 
C     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
C     &                      )
C               ENDIF


C               IF ( IFD.NE.1 .OR. k.LT.NOG ) then
C                  CALL BMG3_SymStd_SETUP_ITLI27_ex(
C     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
C     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
C     &                      )
C               ELSE
C                  CALL BMG3_SymStd_SETUP_ITLI07_ex(
C     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
C     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
C     &                      )
C               ENDIF

C               !
C               ! Interactively examine the coarse grid stencil
C               !
C               IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) ) THEN
C                  CALL BMG3_SymStd_DUMP_stencil( 
C     &                      BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
C     &                      IFD, NStncl_CG, K-1, NOG
C     &                      )
C               ENDIF

C            ENDDO

         IF ( 
     &             BMG_iPARMS(id_BMG3_CG_TYPE)
     &                .EQ.BMG_CG_ITLI 
     & 
     &             .AND. 
     &             
     &             BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &                .EQ.BMG_CG_CONS_block
     &
     &           ) THEN

            ! ---------------------------------------------------------
            !  Block construction of the "honest" I^{T}LI operator.
            ! ---------------------------------------------------------

            !
            ! Loop over grids
            !
            DO K = NOG, 2, -1

               !
               ! Determine the number of points in the stencil
               !
               IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
                  NStncl=14
               ELSE
                  NStncl=4
               ENDIF

               CALL BMG3_SymStd_GET_pointers(
     &                   K, IGRD, NOGm,
     &                   p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                   )
               CALL BMG3_SymStd_GET_pointers( 
     &                   K-1, IGRD, NOGm,
     &                   p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
     &                   )

               !
               ! Operator Induced Interpolation
               !
               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_CGTEMP_yo)
               IF ( INTERP.EQ.BMG_INTERP_OI ) THEN
                  !
                  ! Operator Induced Interpolation
                  !
                  CALL BMG3_SymStd_SETUP_interp_OI(
     &                 K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                 Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
     &                 ifd, NStncl, irelax
     &                 )
               ELSEIF ( INTERP.EQ.BMG_INTERP_GMG ) THEN
                  !
                  ! Trilinear Interpolation
                  !                  
                  CALL BMG3_SymStd_SETUP_interp_TI(
     &                 K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                 Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
     &                 ifd, NStncl, irelax, BMG_rWORK_PL(p_CGTEMP_yo),
     &                 BMG_IOFLAG, BMG_iPARMS
     &                 )
               ELSE  

                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
                     WRITE(*,500) 'This interp operator is invalid!'
                     WRITE(*,520) 
     &                    'HAVE: BMG_iPARMS(id_BMG3_INTERP) =', 
     &                    BMG_iPARMS(id_BMG3_INTERP) 
                  END IF  
         
                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,27) 
                  RETURN 

               ENDIF 
               !
               ! Interactively examine the restriction operator
               !
               IF ( BMG_IOFLAG(iBMG3_BUG_RESTRICT) ) THEN
                  CALL BMG3_SymStd_DUMP_restrict( 
     &                      BMG_IOFLAG, CI(p_CIC), K, K-1, 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF

               IF ( IFD.NE.1 .OR. k.LT.NOG ) then
                  CALL BMG3_SymStd_SETUP_ITLI27_bl(
     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ELSE
                  CALL BMG3_SymStd_SETUP_ITLI07_bl(
     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF

               !
               ! Interactively examine the coarse grid stencil
               !
C               IF( K.EQ.NOG-1 ) THEN
C                  CALL BMG3_SymStd_PRINT_matrix( 
C     &                 BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
C     &                 IFD, NStncl_CG, K-1, NOG
C     &                 )
C               ENDIF

               IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) ) THEN
                  CALL BMG3_SymStd_DUMP_stencil( 
     &                      BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
     &                      IFD, NStncl_CG, K-1, NOG
     &                      )
               ENDIF

            ENDDO

         ELSE IF ( 
     &      BMG_iPARMS(id_BMG3_CG_TYPE).EQ.BMG_CG_USER
     &      ) THEN

            
            !
            !  Set error trap: USER option not supported yet.
            !

         ELSE
            
            !
            !  Set error trap: Unsupported CG (TYPE,CONSTRUCTION)
            !

         ENDIF


         !
         !  Re-zero workspace
         !
         DO i=1, NBMG_rWORK_PL
            BMG_rWORK_PL(i) = rZERO
         ENDDO

         !
         !  Relaxation
         ! 
         CALL BMG3_SymStd_SETUP_relax(  
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             IFD, IBC, IRELAX, NOGm, NOG, IGRD,  
     &             SO, NSO, SOR, NSOR, CI, NCI, 
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL 
     &             ) 
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF
         

      ENDIF

C ------------------------------------------
C     Setup the coarse grid solve
C ------------------------------------------

      IF ( NOG.EQ. 1 ) THEN 
         !
         ! Determine the number of points in the stencil
         !
         IF ( IFD.NE.BMG_STENCIL_7pt ) THEN
            NStncl=14
         ELSE
            NStncl=4
         ENDIF
      ELSE
         !
         ! There has been at least one coarsening.
         !
         NStncl=14
         !
      ENDIF
      !
      ! Obtain pointers for the stencil on the coarsest grid
      ! 
      K=1          ! Coarsest-grid index is always 1
      CALL BMG3_SymStd_GET_pointers(
     &          K, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )

      !
      ! Interactively examine the stencil on the coarsest grid
      !
      IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG1) ) THEN
         CALL BMG3_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SO), Nx, Ny, Nz, 
     &             IFD, NStncl, K, NOG
     &             )
      ENDIF

      !
      ! Setup and factor the matrix for direct solve on the coarsest grid.
      !
      CALL BMG3_SymStd_SETUP_cg_LU( 
     &          SO(p_SO), Nx, Ny, Nz, NStncl, abd, nabd1, nabd2,
     &          BMG_IOFLAG, BMG_iPARMS 
     &          )
      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF
      

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_parts.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
