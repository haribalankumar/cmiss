      SUBROUTINE BMG3_SymStd_SETUP_interp_TI( 
     &                kgf, kgc, so, soc, ci,
     &                iif, jjf, kkf, iic, jjc, kkc, 
     &                NOG, ifd, NStncl, irelax,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_interp_TI.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_interp_TI.f constructs the trilinear interpolation 
C   operator CI.  It is only being used at the time for problems for which 
C   the opertors are 27-pt operators.  The operator CI interpolates a vector 
C   from the coarse grid KC, to the fine grid, KF. 
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2005/02/21  (TMA)
C   - the core computation was cut from mgcf.f 
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
      INCLUDE 'BMG_parameters.h' 
      INCLUDE 'BMG_stencils.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER NOG, NStncl, ifd, iic, iif, irelax,
     &        jjc, jjf, kgc, kgf, kkc, kkf
      REAL*8  ci(iic,jjc,kkc,26),
     &        so(iif,jjf,kkf,NStncl), soc(iic,jjc,kkc,14)

      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER BMG_iPARMS(NBMG_iPARMS)

C --------------------------
C     Local Declarations:
C
      INTEGER ic, i, iicf, iicf1, iic1, iif1, jc, j, jjcf, jjcf1, 
     &        jjc1, jjf1, kc, k, kk, kkcf, kkc1, kkcf1, kkcp1,  kkf1, 
     &	  kpz 
      REAL*8  eMACH

      REAL*8    rHALF, rFOURTH, rEIGHTH
      PARAMETER (rHALF=0.5D0, rFOURTH=0.25D0, rEIGHTH=0.125D0)

C ==========================================================================

      eMACH = 1.d-13

      iic1 = iic-1
      jjc1 = jjc-1
      kkc1 = kkc-1

      iif1 = iif-1
      jjf1 = jjf-1
      kkf1 = kkf-1

      iicf = (iif-2)/2+3
      jjcf = (jjf-2)/2+3
      kkcf = (kkf-2)/2+3

      iicf1 = iicf-1
      jjcf1 = jjcf-1
      kkcf1 = kkcf-1

      kkcp1 = kkc+1

      IF( kgf.LT.NOG .OR. ifd.NE.1 ) THEN

c ***********************************************************************
c   begin computation of i when k difference operator is 27 point
c ***********************************************************************

         ! 
         ! Set-up interpolant for fine grid points on  
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iicf1,jjc1,kkc1)
	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkc1-1)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+2
            k = 2*(kc-1)
            ci(ic,jc,kc,lxyl)=rHALF
            ci(ic,jc,kc,lxyr)=rHALF
         ENDDO
C$OMP END PARALLEL DO
         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF k-planes that sit on FINE-ONLY x-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iic1,jjcf1,kkc1)
	 DO kk=1,(iic1-1)*(jjcf1-2)*(kkc1-1)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjcf1-2))+2
            k = 2*(kc-1)
            ci(ic,jc,kc,lxya)=rHALF
            ci(ic,jc,kc,lxyb)=rHALF
         ENDDO
C$OMP END PARALLEL DO
         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iic1,jjc1,kkcf1)
	 DO kk=1,(iic1-1)*(jjc1-1)*(kkcf1-2)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjc1-1))+3
            k = 2*(kc-1)
            ci(ic,jc,kc,lxza)=rHALF
            ci(ic,jc,kc,lxzb)=rHALF
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iicf1,jjcf1,kkc1)
	 DO kk=1,(iicf1-2)*(jjcf1-2)*(kkc1-1)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjcf1-2))+2
            k = 2*(kc-1)
            ci(ic,jc,kc,lxynw)=rFOURTH
            ci(ic,jc,kc,lxyne)=rFOURTH
            ci(ic,jc,kc,lxyse)=rFOURTH
            ci(ic,jc,kc,lxysw)=rFOURTH
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on 
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iicf1,jjc1,kkcf1)
 	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkcf1-2)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+3
            k = 2*(kc-1)
            ci(ic,jc,kc,lxznw)=rFOURTH
            ci(ic,jc,kc,lxzne)=rFOURTH
            ci(ic,jc,kc,lxzse)=rFOURTH
            ci(ic,jc,kc,lxzsw)=rFOURTH
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on 
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iic1,jjcf1,kkcf1)
 	 DO kk=1,(iic1-1)*(jjcf1-2)*(kkcf1-2)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjcf1-2))+3
            k = 2*(kc-1)
            ci(ic,jc,kc,lyznw)=rFOURTH
            ci(ic,jc,kc,lyzne)=rFOURTH
            ci(ic,jc,kc,lyzse)=rFOURTH
            ci(ic,jc,kc,lyzsw)=rFOURTH
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, y-lines, and z-lines.
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(ci,iicf1,jjcf1,kkcf1)
 	 DO kk=1,(iicf1-2)*(jjcf1-2)*(kkcf1-2)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjcf1-2))+3
            k = 2*(kc-1)
            ci(ic,jc,kc,ltnw)=rEIGHTH
            ci(ic,jc,kc,ltne)=rEIGHTH
            ci(ic,jc,kc,lbnw)=rEIGHTH
            ci(ic,jc,kc,lbne)=rEIGHTH
            ci(ic,jc,kc,lbsw)=rEIGHTH
            ci(ic,jc,kc,ltsw)=rEIGHTH
            ci(ic,jc,kc,ltse)=rEIGHTH
            ci(ic,jc,kc,lbse)=rEIGHTH
         ENDDO

      ELSE

         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,*) 'ERROR: BMG3_SymStd_SETUP_interp_TI   .... '
            WRITE(*,*) '***** Bilinear interp not designed for 7-pt ops'
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,26)

      ENDIF

C ==========================================================================

      RETURN
      END
