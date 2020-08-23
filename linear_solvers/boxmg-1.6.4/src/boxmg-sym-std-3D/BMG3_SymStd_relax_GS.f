      SUBROUTINE BMG3_SymStd_relax_GS( 
     &                kg, so, qf, q, RES, sor, ii, jj, kk,
     &                RES_L2, BMG_IOFLAG, NOG, ifd, NStncl, NSORv, 
     &                iRELAX, iRELAX_SYM, UPDOWN, BMG_IJK_MAP
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
     &        BMG_IJK_MAP, p_IJKMAP
      REAL*8  q(ii,jj,kk), qf(ii,jj,kk), RES(ii,jj,kk), RES_L2,
     &        so(ii,jj,kk,NStncl), sor(ii,jj,kk,NSORv)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, i1, ibeg, iend, jbeg, jend, kbeg, kend,
     &        j, j1, k, k1, kl, ip, jp, kpp,
     &        istride, jstride, kstride, maxij, maxijk,
     &        ptstart, ptend, ptstride,
     &        pts


C ==========================================================================

      j1 = jj-1
      i1 = ii-1
      k1 = kk-1

      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN 
         !
         !   27-point relaxation (8-color)
         !

         IF ( UPDOWN.EQ.BMG_UP 
     &       .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
            ptstart = 1
            ptend   = 8
            ptstride = 1
         ELSE
            ptstart = 8
            ptend = 1
            ptstride = -1
         ENDIF         
        
         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

            ibeg = 2+mod(pts-1,2)
            if(mod(ibeg,2).eq.0) then  ! begin on even number
               if(mod(I1,2).eq.0) then
                  iend = I1
               else
                  iend = I1-1
               endif
            else
               if(mod(I1,2).eq.0) then
                  iend = I1-1
               else
                  iend = I1
               endif
            endif

            jbeg = 2+mod(mod((pts-1)/2,2),2)
            if(mod(jbeg,2).eq.0) then  ! begin on even number
               if(mod(J1,2).eq.0) then
                  jend = J1
               else
                  jend = J1-1
               endif
            else
               if(mod(J1,2).eq.0) then
                  jend = J1-1
               else
                  jend = J1
               endif
            endif

            kbeg = 2+mod((pts-1)/4,2)
            if(mod(kbeg,2).eq.0) then  ! begin on even number
               if(mod(K1,2).eq.0) then
                  kend = K1
               else
                  kend = K1-1
               endif
            else
               if(mod(K1,2).eq.0) then
                  kend = K1-1
               else
                  kend = K1
               endif
            endif

            istride = (iend-ibeg)/2+1
            jstride = (jend-jbeg)/2+1
            kstride = (kend-kbeg)/2+1

            maxijk = istride*jstride*kstride
            maxij = istride*jstride

C$OMP PARALLEL DO 
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ip,j,jp,k,kl,kpp)
C$OMP& SHARED(ibeg,istride,jbeg,jstride,
C$OMP&        kbeg,maxij,maxijk,q,qf,so,sor)
            DO kl=0,maxijk-1

               ip  = mod(kl,istride)
               jp  = mod(kl/istride,jstride)
               kpp = kl/maxij

               i=ibeg+2*ip
               j=jbeg+2*jp
               k=kbeg+2*kpp

               q(i,j,k) = ( qf(i,j,k) 
     &              + so(i,j,k,kpw)*q(i-1,j,k)
     &              + so(i,j+1,k,kpnw)*q(i-1,j+1,k)
     &              + so(i,j+1,k,kps)*q(i,j+1,k)
     &              + so(i+1,j+1,k,kpsw)*q(i+1,j+1,k)
     &              + so(i+1,j,k,kpw)*q(i+1,j,k)
     &              + so(i+1,j,k,kpnw)*q(i+1,j-1,k)
     &              + so(i,j,k,kps)*q(i,j-1,k)
     &              + so(i,j,k,kpsw)*q(i-1,j-1,k)
     &              + so(i,j,k,kb)*q(i,j,k-1)
     &              + so(i,j,k,kbw)*q(i-1,j,k-1)
     &              + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)
     &              + so(i,j+1,k,kbn)*q(i,j+1,k-1)
     &              + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)
     &              + so(i+1,j,k,kbe)*q(i+1,j,k-1)
     &              + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)
     &              + so(i,j,k,kbs)*q(i,j-1,k-1)
     &              + so(i,j,k,kbsw)*q(i-1,j-1,k-1)
     &              + so(i,j,k+1,kb)*q(i,j,k+1)
     &              + so(i,j,k+1,kbe)*q(i-1,j,k+1)
     &              + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)
     &              + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)
     &              + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)
     &              + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)
     &              + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)
     &              + so(i,j,k+1,kbn)*q(i,j-1,k+1)
     &              + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)
     &              )*sor(i,j,k,msor)
            END DO
C$OMP END PARALLEL DO            

         ENDDO ! loop over colors
      
      ELSE
         !  
         !  7 point relaxation ( four colors )
         !
         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

            ibeg=mod(j+k+pts,2) + 2
            iend=2*((i1-ibeg)/2)+ibeg

            istride = (iend-ibeg)/2+1
            jstride =  j1-2
            kstride =  k1-2

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k)
C$OMP& SHARED(i1,ibeg,iend,istride,jstride,kstride,pts,q,qf,so,sor)
            DO kl=0,jstride*kstride-1
               !
               j = mod(k,kstride) + 2
               k = (kl/kstride) + 2
               !
               ibeg=mod(j+k+pts,2)+2
               iend=2*((i1-ibeg)/2)+ibeg
               !
               DO i=ibeg,iend,2
 
                  q(i,j,k) = ( qf(i,j,k)
     &                 + so(i,j,k,kpw)*q(i-1,j,k)
     &                 + so(i,j+1,k,kps)*q(i,j+1,k)
     &                 + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                 + so(i,j,k,kps)*q(i,j-1,k)
     &                 + so(i,j,k,kb)*q(i,j,k-1)
     &                 + so(i,j,k+1,kb)*q(i,j,k+1)
     &                 )*sor(i,j,k,msor)
                  
               END DO
            END DO
C$OMP END PARALLEL DO
            !
         END DO ! loop over colors 
         !
      ENDIF

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         p_IJKMAP = BMG_IJK_MAP(kg)
         CALL BMG3_SymStd_residual( 
     &             kg, NOG, IFD, Q, QF, SO, RES, II, JJ, KK, NStncl,
     &             BMG_IJK_MAP(p_IJKMAP)
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES, II, JJ, KK, RES_L2
     &             )
      ENDIF

C =============================

      RETURN
      END

