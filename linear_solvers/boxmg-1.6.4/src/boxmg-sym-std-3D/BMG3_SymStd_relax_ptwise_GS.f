      SUBROUTINE BMG3_SymStd_relax_ptwise_GS( 
     &                kg, so, qf, q, RES, sor, ii, jj, kk,
     &                RES_L2, BMG_IOFLAG, NOG, ifd, NStncl, NSORv, 
     &                iRELAX, iRELAX_SYM, UPDOWN,
     &                BMG_IJK_MAP
     &                )

C ==========================================================================
C
C   BMG3_SymStd_relax_GS.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel
C   (with the correct ordering depending on whether we are
C   on the way down, or up in the symmetric cycling case)
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   2000/02/23   - written (M.Berndt)
C                - core computation cut from mgrlx3.f
C   2000/03/06   - updated residual calculation (JDM) 
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
C ======================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER ii, jj, kk, NSORv, NStncl

      INTEGER ifd, iRELAX, iRELAX_SYM, kg, NOG, UPDOWN, 
     &        BMG_IJK_MAP(*)
      REAL*8  q(ii,jj,kk), qf(ii,jj,kk), RES(ii,jj,kk), RES_L2,
     &        so(ii,jj,kk,NStncl), sor(ii,jj,kk,NSORv)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, i1, ibeg, iend, jbeg, jend, kbeg, kend,
     &        j, j1, k, k1, kl, ip, jp, kpp, ibmg_start,
     &        istride, jstride, kstride,
     &        ptstart, ptend, ptstride,
     &        pts

C ==========================================================================


      j1 = jj-1
      i1 = ii-1
      k1 = kk-1

      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN 

         IF ( UPDOWN.EQ.BMG_UP 
     &       .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
            ibeg = i1
            iend = 2
            istride = -1
            jbeg = j1
            jend = 2
            jstride = -1
            kbeg = k1
            kend = 2
            kstride = -1
         ELSE
            ibeg = 2
            iend = i1
            istride = 1
            jbeg = 2
            jend = j1
            jstride = 1
            kbeg = 2
            kend = k1
            kstride = 1
         ENDIF 

C$OMP    PARALLEL DO 
C$OMP&    DEFAULT(NONE)
C$OMP&    PRIVATE(i,j,k)
C$OMP&    SHARED(ibeg,iend,jbeg,jend,kbeg,kend,istride,
C$OMP&           jstride,kstride,q,qf,so,sor)
         DO k=kbeg,kend,kstride
            DO j=jbeg,jend,jstride
               DO i=ibeg,iend,istride

                  q(i,j,k) = ( qf(i,j,k) 
     &                 + so(i,j,k,kpw)*q(i-1,j,k)
     &                 + so(i,j+1,k,kpnw)*q(i-1,j+1,k)
     &                 + so(i,j+1,k,kps)*q(i,j+1,k)
     &                 + so(i+1,j+1,k,kpsw)*q(i+1,j+1,k)
     &                 + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                 + so(i+1,j,k,kpnw)*q(i+1,j-1,k)
     &                 + so(i,j,k,kps)*q(i,j-1,k)
     &                 + so(i,j,k,kpsw)*q(i-1,j-1,k)
     &                 + so(i,j,k,kb)*q(i,j,k-1)
     &                 + so(i,j,k,kbw)*q(i-1,j,k-1)
     &                 + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)
     &                 + so(i,j+1,k,kbn)*q(i,j+1,k-1)
     &                 + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)
     &                 + so(i+1,j,k,kbe)*q(i+1,j,k-1)
     &                 + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)
     &                 + so(i,j,k,kbs)*q(i,j-1,k-1)
     &                 + so(i,j,k,kbsw)*q(i-1,j-1,k-1)
     &                 + so(i,j,k+1,kb)*q(i,j,k+1)
     &                 + so(i,j,k+1,kbe)*q(i-1,j,k+1)
     &                 + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)
     &                 + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)
     &                 + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)
     &                 + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)
     &                 + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)
     &                 + so(i,j,k+1,kbn)*q(i,j-1,k+1)
     &                 + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)
     &                 )*sor(i,j,k,msor)
               ENDDO
            ENDDO
         ENDDO
C$OMP END PARALLEL DO            
         !
      ELSE
         !
         IF ( UPDOWN.EQ.BMG_UP 
     &       .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
            ibeg = i1
            iend = 2
            istride = -1
            jbeg = j1
            jend = 2
            jstride = -1
            kbeg = k1
            kend = 2
            kstride = -1
         ELSE
            ibeg = 2
            iend = i1
            istride = 1
            jbeg = 2
            jend = j1
            jstride = 1
            kbeg = 2
            kend = k1
            kstride = 1
         ENDIF  


C$OMP    PARALLEL DO 
C$OMP&    DEFAULT(NONE)
C$OMP&    PRIVATE(i,j,k)
C$OMP&    SHARED(ibeg,iend,jbeg,jend,kbeg,kend,istride,
C$OMP&           jstride,kstride,q,qf,so,sor)
         DO k=kbeg,kend,kstride
            DO j=jbeg,jend,jstride
               DO i=ibeg,iend,istride

                  q(i,j,k) = ( qf(i,j,k)
     &                 + so(i,j,k,kpw)*q(i-1,j,k)
     &                 + so(i,j+1,k,kps)*q(i,j+1,k)
     &                 + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                 + so(i,j,k,kps)*q(i,j-1,k)
     &                 + so(i,j,k,kb)*q(i,j,k-1)
     &                 + so(i,j,k+1,kb)*q(i,j,k+1)
     &                 )*sor(i,j,k,msor)
               ENDDO
            ENDDO
         ENDDO
C$OMP END PARALLEL DO
         !
      ENDIF

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         CALL BMG3_SymStd_residual( 
     &             kg, NOG, IFD, Q, QF, SO, RES, II, JJ, KK, NStncl,
     &             BMG_IJK_MAP
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES, II, JJ, KK, RES_L2
     &             )
      ENDIF

C =============================

      RETURN
      END

