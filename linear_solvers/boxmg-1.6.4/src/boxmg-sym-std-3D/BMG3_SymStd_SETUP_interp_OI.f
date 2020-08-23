      SUBROUTINE BMG3_SymStd_SETUP_interp_OI( 
     &                kgf, kgc, so, soc, ci,
     &                iif, jjf, kkf, iic, jjc, kkc, 
     &                NOG, ifd, NStncl, irelax
     &                )

C ==========================================================================
C
C   BMG3_SymStd_SETUP_interp_OI.f
C
C   --------------------
C   DESCRIPTION:
C   --------------------
C
C   BMG3_SymStd_SETUP_interp_OI.f constructs the operator-induced 
C   interpolation operator CI from the fine-grid stencil, SO.
C   The operator CI interpolates a vector from the coarse grid
C   KC, to the fine grid, KF. 
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C   
C   Written:    2000/03/06  (JDM)
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

      INCLUDE 'BMG_stencils.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER NOG, NStncl, ifd, iic, iif, irelax,
     &        jjc, jjf, kgc, kgf, kkc, kkf
      REAL*8  ci(iic,jjc,kkc,26),
     &        so(iif,jjf,kkf,NStncl), soc(iic,jjc,kkc,14)
C --------------------------
C     Local Declarations:
C
      INTEGER ic, i, iicf, iicf1, iic1, iif1, jc, j, jjcf, jjcf1, 
     &        jjc1, jjf1, kc, k, kk, kkcf, kkc1, kkcf1, kkcp1,  kkf1, 
     &	      kpz 
      REAL*8  a, b, c, de, dn, dne, dnw, dp, ds, dse, dsw, dw,
     &        d1, d2, d3, eMACH, ep, sum

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

      if(kgf.lt.NOG.or.ifd.ne.1) then

c ***********************************************************************
c   begin computation of i when k difference operator is 27 point
c ***********************************************************************

         ! 
         ! Set-up interpolant for fine grid points on  
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iicf1,jjc1,kkc1,so)
	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkc1-1)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+2
            k = 2*(kc-1)
            a=so(i-1,j+1,k,kpnw)+so(i-1,j,k,kpw)+so(i-1,j,k,kpsw)
     &           +so(i-1,j+1,k,kbnw)+so(i-1,j,k,kbw)
     &           +so(i-1,j,k,kbsw)+so(i-1,j+1,k+1,kbse)
     &           +so(i-1,j,k+1,kbe)+so(i-1,j,k+1,kbne)
            b=so(i,j+1,k,kpsw)+so(i,j,k,kpw)+so(i,j,k,kpnw)
     &           +so(i,j+1,k,kbne)+so(i,j,k,kbe)+so(i,j,k,kbse)
     &           +so(i,j+1,k+1,kbsw)+so(i,j,k+1,kbw)
     &           +so(i,j,k+1,kbnw)
            c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)
     &           +so(i-1,j+1,k,kbn)+so(i-1,j,k,kb)+so(i-1,j,k,kbs)
     &           +so(i-1,j+1,k+1,kbs)+so(i-1,j,k+1,kb)
     &           +so(i-1,j,k+1,kbn)
            ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b/so(i-1,j,k,kp)))
            c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)
     &           -(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc,lxyl)=a/c
            ci(ic,jc,kc,lxyr)=b/c
         ENDDO
C$OMP END PARALLEL DO		
         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF k-planes that sit on FINE-ONLY x-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iic1,jjcf1,kkc1,so)
	 DO kk=1,(iic1-1)*(jjcf1-2)*(kkc1-1)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjcf1-2))+2
            k = 2*(kc-1)
            a=so(i,j,k,kpnw)+so(i,j,k,kps)+so(i+1,j,k,kpsw)
     &           +so(i,j,k,kbnw)+so(i,j,k,kbn)+so(i+1,j,k,kbne)
     &           +so(i,j,k+1,kbse)+so(i,j,k+1,kbs)
     &           +so(i+1,j,k+1,kbsw)
            b=so(i,j-1,k,kpsw)+so(i,j-1,k,kps)+so(i+1,j-1,k,kpnw)
     &           +so(i,j-1,k,kbsw)+so(i,j-1,k,kbs)
     &           +so(i+1,j-1,k,kbse)+so(i,j-1,k+1,kbne)
     &           +so(i,j-1,k+1,kbn)+so(i+1,j-1,k+1,kbnw)
            ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
            c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)
     &           +so(i,j-1,k,kbw)+so(i,j-1,k,kb)+so(i+1,j-1,k,kbe)
     &           +so(i,j-1,k+1,kbe)+so(i,j-1,k+1,kb)
     &           +so(i+1,j-1,k+1,kbw)
            c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)
     &           -(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc,lxya)=a/c
            ci(ic,jc,kc,lxyb)=b/c
         ENDDO
C$OMP END PARALLEL DO
         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iic1,jjc1,kkcf1,so)
	 DO kk=1,(iic1-1)*(jjc1-1)*(kkcf1-2)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjc1-1))+3
            k = 2*(kc-1)
            a=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)+so(i+1,j+1,k,kbsw)
     &           +so(i,j,k,kbe)+so(i,j,k,kb)+so(i+1,j,k,kbw)
     &           +so(i,j,k,kbne)+so(i,j,k,kbn)+so(i+1,j,k,kbnw)
            b=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)
     &           +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)
     &           +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)
     &           +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)
     &           +so(i+1,j,k-1,kbse)
            c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)
     &           +so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)
     &           +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpsw)
     &           +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
            ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
            c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &           -(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc,lxza)=a/c
            ci(ic,jc,kc,lxzb)=b/c
            !PRINT *, ' 3 : ',kk
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
C         PRINT *,' 4 '
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(dnw,dn,dne,dw,de,dsw,ds,dse,ep,dp,
C$OMP&         i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iicf1,jjcf1,kkc1,so)
	 DO kk=1,(iicf1-2)*(jjcf1-2)*(kkc1-1)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjcf1-2))+2
            k = 2*(kc-1)
            dnw=so(i-1,j,k,kpnw)+so(i-1,j,k,kbnw)
     &           +so(i-1,j,k+1,kbse)
            dn=so(i-1,j,k,kps)+so(i-1,j,k,kbn)+so(i-1,j,k+1,kbs)
            dne=so(i,j,k,kpsw)+so(i,j,k,kbne)+so(i,j,k+1,kbsw)
            dw=so(i-1,j-1,k,kpw)+so(i-1,j-1,k,kbw)
     &           +so(i-1,j-1,k+1,kbe)
            de=so(i,j-1,k,kpw)+so(i,j-1,k,kbe)+so(i,j-1,k+1,kbw)
            dsw=so(i-1,j-1,k,kpsw)+so(i-1,j-1,k,kbsw)
     &           +so(i-1,j-1,k+1,kbne)
            ds=so(i-1,j-1,k,kps)+so(i-1,j-1,k,kbs)
     &           +so(i-1,j-1,k+1,kbn)
            dse=so(i,j-1,k,kpnw)+so(i,j-1,k,kbse)
     &           +so(i,j-1,k+1,kbnw)
            ep=MIN(abs((dsw+dw+dnw)/so(i-1,j-1,k,kp)),
     &           abs((dnw+dn+dne)/so(i-1,j-1,k,kp)),
     &           abs((dne+de+dse)/so(i-1,j-1,k,kp)),
     &           abs((dse+ds+dsw)/so(i-1,j-1,k,kp))
     &           )
            dp=dw+dnw+dn+dne+de+dse+ds+dsw
            sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)
     &           -so(i-1,j-1,k+1,kb)
            dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)
     &           /(abs(sum-(rONE+ep)*dp)+eMACH)
            dp=rONE/dp
            ci(ic,jc,kc,lxynw)=dp*(dnw+ci(ic-1,jc,kc,lxya)*dw
     &           +ci(ic,jc,kc,lxyl)*dn)
            ci(ic,jc,kc,lxyne)=dp*(dne+ci(ic,jc,kc,lxyr)*dn
     &           +ci(ic,jc,kc,lxya)*de)
            ci(ic,jc,kc,lxyse)=dp*(dse+ci(ic,jc,kc,lxyb)*de
     &           +ci(ic,jc-1,kc,lxyr)*ds)
            ci(ic,jc,kc,lxysw)=dp*(dsw+ci(ic,jc-1,kc,lxyl)*ds
     &           +ci(ic-1,jc,kc,lxyb)*dw)
            !PRINT *, ' 4 : ',kk
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on 
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         ! 
C         PRINT *,' 5 '
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(dnw,dn,dne,dw,de,dsw,ds,dse,ep,dp,
C$OMP&         i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iicf1,jjc1,kkcf1,so)
 	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkcf1-2)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+3
            k = 2*(kc-1)
            dnw=so(i-1,j+1,k,kbse)+so(i-1,j,k,kbe)
     &           +so(i-1,j,k,kbne)
            dn=so(i-1,j+1,k,kbs)+so(i-1,j,k,kb)+so(i-1,j,k,kbn)
            dne=so(i,j+1,k,kbsw)+so(i,j,k,kbw)+so(i,j,k,kbnw)
            dw=so(i-1,j+1,k-1,kpnw)+so(i-1,j,k-1,kpw)
     &           +so(i-1,j,k-1,kpsw)
            de=so(i,j+1,k-1,kpsw)+so(i,j,k-1,kpw)
     &           +so(i,j,k-1,kpnw)
            dsw=so(i-1,j+1,k-1,kbnw)+so(i-1,j,k-1,kbw)
     &           +so(i-1,j,k-1,kbsw)
            ds=so(i-1,j+1,k-1,kbn)+so(i-1,j,k-1,kb)
     &           +so(i-1,j,k-1,kbs)
            dse=so(i,j+1,k-1,kbne)+so(i,j,k-1,kbe)
     &           +so(i,j,k-1,kbse)
            ep=MIN(abs((dsw+dw+dnw)/so(i-1,j,k-1,kp)),
     &           abs((dnw+dn+dne)/so(i-1,j,k-1,kp)),
     &           abs((dne+de+dse)/so(i-1,j,k-1,kp)),
     &           abs((dse+ds+dsw)/so(i-1,j,k-1,kp)))
            dp=dw+dnw+dn+dne+de+dse+ds+dsw
            sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)
     &           -so(i-1,j,k-1,kps)
            dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)
     &           /(abs(sum-(rONe+ep)*dp)+eMACH)
            dp=rONE/dp
            ci(ic,jc,kc,lxznw)=dp*(dnw+ci(ic-1,jc,kc,lxza)*dw
     &           +ci(ic,jc,kc,lxyl)*dn)
            ci(ic,jc,kc,lxzne)=dp*(dne+ci(ic,jc,kc,lxyr)*dn
     &           +ci(ic,jc,kc,lxza)*de)
            ci(ic,jc,kc,lxzse)=dp*(dse+ci(ic,jc,kc,lxzb)*de
     &           +ci(ic,jc,kc-1,lxyr)*ds)
            ci(ic,jc,kc,lxzsw)=dp*(dsw+ci(ic,jc,kc-1,lxyl)*ds
     &           +ci(ic-1,jc,kc,lxzb)*dw)
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on 
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         ! 
C         PRINT *,' 6 '
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(dnw,dn,dne,dw,de,dsw,ds,dse,ep,dp,
C$OMP&         i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iic1,jjcf1,kkcf1,so)
 	 DO kk=1,(iic1-1)*(jjcf1-2)*(kkcf1-2)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjcf1-2))+3
            k = 2*(kc-1)
            dnw=so(i,j,k,kbse)+so(i,j,k,kbs)+so(i+1,j,k,kbsw)
            dn=so(i,j-1,k,kbe)+so(i,j-1,k,kb)+so(i+1,j-1,k,kbw)
            dne=so(i,j-1,k,kbne)+so(i,j-1,k,kbn)
     &           +so(i+1,j-1,k,kbnw)
            dw=so(i,j,k-1,kpnw)+so(i,j,k-1,kps)
     &           +so(i+1,j,k-1,kpsw)
            de=so(i,j-1,k-1,kpsw)+so(i,j-1,k-1,kps)
     &           +so(i+1,j-1,k-1,kpnw)
            dsw=so(i,j,k-1,kbnw)+so(i,j,k-1,kbn)
     &           +so(i+1,j,k-1,kbne)
            ds=so(i,j-1,k-1,kbw)+so(i,j-1,k-1,kb)
     &           +so(i+1,j-1,k-1,kbe)
            dse=so(i,j-1,k-1,kbsw)+so(i,j-1,k-1,kbs)
     &           +so(i+1,j-1,k-1,kbse)
            ep=MIN(abs((dsw+dw+dnw)/so(i,j-1,k-1,kp)),
     &           abs((dnw+dn+dne)/so(i,j-1,k-1,kp)),
     &           abs((dne+de+dse)/so(i,j-1,k-1,kp)),
     &           abs((dse+ds+dsw)/so(i,j-1,k-1,kp))
     &           )
            dp=dw+dnw+dn+dne+de+dse+ds+dsw
            sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)
     &           -so(i+1,j-1,k-1,kpw)
            dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)
     &           /(abs(sum-(rONE+ep)*dp)+eMACH)
            dp=rONE/dp
            ci(ic,jc,kc,lyznw)=dp*(dnw+ci(ic,jc,kc,lxza)*dw
     &           +ci(ic,jc,kc,lxya)*dn)
            ci(ic,jc,kc,lyzne)=dp*(dne+ci(ic,jc,kc,lxyb)*dn
     &           +ci(ic,jc-1,kc,lxza)*de)
            ci(ic,jc,kc,lyzse)=dp*(dse+ci(ic,jc-1,kc,lxzb)*de
     &           +ci(ic,jc,kc-1,lxyb)*ds)
            ci(ic,jc,kc,lyzsw)=dp*(dsw+ci(ic,jc,kc-1,lxya)*ds
     &           +ci(ic,jc,kc,lxzb)*dw)
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, y-lines, and z-lines.
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(d1,d2,d3,i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iicf1,jjcf1,jjf,kkcf1,so)
 	 DO kk=1,(iicf1-2)*(jjcf1-2)*(kkcf1-2)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjcf1-2))+3
            k = 2*(kc-1)
            d1=so(i-1,j-1,k-1,kpw)+so(i-1,j,k-1,kpnw)
     &           +so(i-1,j,k-1,kps)+so(i,j,k-1,kpsw)
     &           +so(i,j-1,k-1,kpw)
     &           +so(i,j-1,k-1,kpnw)+so(i-1,j-1,k-1,kps)
     &           +so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k-1,kb)
     &           +so(i-1,j-1,k-1,kbw)+so(i-1,j,k-1,kbnw)
     &           +so(i-1,j,k-1,kbn)+so(i,j,k-1,kbne)
     &           +so(i,j-1,k-1,kbe)
     &           +so(i,j-1,k-1,kbse)+so(i-1,j-1,k-1,kbs)
     &           +so(i-1,j-1,k-1,kbsw)+so(i-1,j-1,k,kb)
     &           +so(i-1,j-1,k,kbe)+so(i-1,j,k,kbse)
     &           +so(i-1,j,k,kbs)
     &           +so(i,j,k,kbsw)+so(i,j-1,k,kbw)+so(i,j-1,k,kbnw)
     &           +so(i-1,j-1,k,kbn)+so(i-1,j-1,k,kbne)
            d2=MIN(abs(so(i-1,j-1,k-1,kpw)+so(i-1,j,k-1,kpnw)
     &           + so(i-1,j,k,kbse)+so(i-1,j-1,k,kbe)
     &           + so(i-1,j-1,k,kbne)
     &           + so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k-1,kbsw)
     &           + so(i-1,j-1,k-1,kbw)+so(i-1,j,k-1,kbnw))/
     &             so(i-1,j-1,k-1,kp),
     &            abs(so(i,j-1,k-1,kpw)+so(i,j,k-1,kpsw)
     &           + so(i,j,k,kbsw)+so(i,j-1,k,kbw)
     &           + so(i,j-1,k,kbnw)+so(i,j-1,k-1,kpnw)
     &           + so(i,j-1,k-1,kbse)
     &           + so(i,j-1,k-1,kbe)+so(i,j,k-1,kbne))/
     &             so(i-1,j-1,k-1,kp),
     &            abs(so(i-1,j,k-1,kps)+so(i-1,j,k-1,kpnw)
     &           + so(i-1,j,k,kbse)+so(i-1,j,k,kbs)+so(i,j,k,kbsw)
     &           + so(i,j,k-1,kpsw)+so(i,j,k-1,kbne)
     &           + so(i-1,j,k-1,kbn)
     &           + so(i-1,j,k-1,kbnw)),abs(so(i-1,j-1,k-1,kps)
     &           + so(i-1,j-1,k-1,kpsw)+so(i-1,j-1,k,kbne)
     &           + so(i-1,j-1,k,kbn)+so(i,j-1,k,kbnw)
     &           + so(i,j-1,k-1,kpnw)+so(i,j-1,k-1,kbse)
     &           + so(i-1,j-1,k-1,kbs)+so(i,j-1,k-1,kbse))/
     &             so(i-1,j-1,k-1,kp))
            d2=MIN(d2,abs(so(i-1,j-1,k-1,kb)+so(i-1,j-1,k-1,kbw)
     &            + so(i-1,j,k-1,kbnw)+so(i-1,j,k-1,kbn)
     &            + so(i,j,k-1,kbne)
     &            + so(i,j-1,k-1,kbe)+so(i,j-1,k-1,kbse)
     &            + so(i-1,j-1,k-1,kbs)
     &            + so(i-1,j-1,k-1,kbsw)),
     &           abs(so(i-1,j-1,k,kb)
     &            + so(i-1,j-1,k,kbe)+so(i-1,j,k,kbse)
     &            + so(i-1,j,k,kbs)
     &            + so(i,j,k,kbsw)+so(i,j-1,k,kbw)+so(i,j-1,k,kbnw)
     &            + so(i-1,j-1,k,kbn)+so(i-1,j-1,k,kbne))/
     &           so(i-1,j-1,k-1,kp)
     &            )
            d1=d1+(so(i-1,j-1,k-1,kp)-d1)*MAX(so(i-1,j-1,k-1,kp)
     &           -(rONE+d2)*d1,rZERO)
     &           /(abs(so(i-1,j-1,k-1,kp)-(rONE+d2)
     &           *d1)+eMACH)
            d1=rONE/d1
            ci(ic,jc,kc,ltnw)
     &           =d1*(so(i-1,j,k,kbse)
     &           +ci(ic-1,jc,kc,lyznw)
     &           *so(i-1,j-1,k-1,kpw)+ci(ic-1,jc,kc,lxza)
     &           *so(i-1,j,k-1,kpnw)
     &           +ci(ic,jc,kc,lxznw)*so(i-1,j,k-1,kps)
     &           +ci(ic-1,jc,kc,lxya)
     &           *so(i-1,j-1,k,kbe)+ci(ic,jc,kc,lxyl)
     &           *so(i-1,j,k,kbs)
     &           +ci(ic,jc,kc,lxynw)*so(i-1,j-1,k,kb))
            ci(ic,jc,kc,ltne)
     &           =d1*(so(i,j,k,kbsw)
     &           +ci(ic,jc,kc,lxzne)
     &           *so(i-1,j,k-1,kps)+ci(ic,jc,kc,lxza)
     &           *so(i,j,k-1,kpsw)
     &           +ci(ic,jc,kc,lyznw)*so(i,j-1,k-1,kpw)
     &           +ci(ic,jc,kc,lxyr)
     &           *so(i-1,j,k,kbs)+ci(ic,jc,kc,lxya)
     &           *so(i,j-1,k,kbw)
     &           +ci(ic,jc,kc,lxyne)*so(i-1,j-1,k,kb))
            ci(ic,jc,kc,lbnw)
     &           =d1*(so(i-1,j,k-1,kbnw)
     &           +ci(ic-1,jc,kc-1,lxya)*so(i-1,j-1,k-1,kbw)
     &           +ci(ic,jc,kc-1,lxyl)*so(i-1,j,k-1,kbn)
     &           +ci(ic,jc,kc-1,lxynw)*so(i-1,j-1,k-1,kb)
     &           +ci(ic-1,jc,kc,lyzsw)*so(i-1,j-1,k-1,kpw)
     &           +ci(ic-1,jc,kc,lxzb)*so(i-1,j,k-1,kpnw)
     &           +ci(ic,jc,kc,lxzsw)*so(i-1,j,k-1,kps))
            ci(ic,jc,kc,lbne)
     &           =d1*(so(i,j,k-1,kbne)
     &           +ci(ic,jc,kc-1,lxyne)
     &           *so(i-1,j-1,k-1,kb)+ci(ic,jc,kc-1,lxyr)
     &           *so(i-1,j,k-1,kbn)
     &           +ci(ic,jc,kc-1,lxya)*so(i,j-1,k-1,kbe)
     &           +ci(ic,jc,kc,lxzse)
     &           *so(i-1,j,k-1,kps)+ci(ic,jc,kc,lxzb)
     &           *so(i,j,k-1,kpsw)
     &           +ci(ic,jc,kc,lyzsw)*so(i,j-1,k-1,kpw))
            ci(ic,jc,kc,lbsw)
     &           =d1*(so(i-1,j-1,k-1,kbsw)
     &           +ci(ic-1,jc,kc-1,lxyb)*so(i-1,j-1,k-1,kbw)
     &           +ci(ic,jc,kc-1,lxysw)
     &           *so(i-1,j-1,k-1,kb)+ci(ic,jc-1,kc-1,lxyl)
     &           *so(i-1,j-1,k-1,kbs)+ci(ic-1,jc,kc,lyzse)
     &           *so(i-1,j-1,k-1,kpw)
     &           +ci(ic,jc-1,kc,lxzsw)*so(i-1,j-1,k-1,kps)
     &           +ci(ic-1,jc-1,kc,lxzb)*so(i-1,j-1,k-1,kpsw))
            ci(ic,jc,kc,ltsw)
     &           =d1*(so(i-1,j-1,k,kbne)
     &           +ci(ic-1,jc,kc,lxyb)*so(i-1,j-1,k,kbe)
     &           +ci(ic,jc,kc,lxysw)*so(i-1,j-1,k,kb)
     &           +ci(ic,jc-1,kc,lxyl)*so(i-1,j-1,k,kbn)
     &           +ci(ic-1,jc,kc,lyzne)
     &           *so(i-1,j-1,k-1,kpw)+ci(ic,jc-1,kc,lxznw)
     &           *so(i-1,j-1,k-1,kps)
     &           +ci(ic-1,jc-1,kc,lxza)*so(i-1,j-1,k-1,kpsw))
            ci(ic,jc,kc,ltse)
     &           =d1*(so(i,j-1,k,kbnw)
     &           +ci(ic,jc-1,kc,lxyr)
     &           *so(i-1,j-1,k,kbn)+ci(ic,jc,kc,lxyse)
     &           *so(i-1,j-1,k,kb)
     &           +ci(ic,jc,kc,lxyb)*so(i,j-1,k,kbw)
     &           +ci(ic,jc-1,kc,lxzne)
     &           *so(i-1,j-1,k-1,kps)+ci(ic,jc,kc,lyzne)
     &           *so(i,j-1,k-1,kpw)
     &           +ci(ic,jc-1,kc,lxza)*so(i,j-1,k-1,kpnw))
            ci(ic,jc,kc,lbse)
     &           =d1*(so(i,j-1,k-1,kbse)
     &           +ci(ic,jc-1,kc-1,lxyr)*so(i-1,j-1,k-1,kbs)
     &           +ci(ic,jc,kc-1,lxyse)*so(i-1,j-1,k-1,kb)
     &           +ci(ic,jc,kc-1,lxyb)*so(i,j-1,k-1,kbe)
     &           +ci(ic,jc-1,kc,lxzse)*so(i-1,j-1,k-1,kps)
     &           +ci(ic,jc,kc,lyzse)*so(i,j-1,k-1,kpw)
     &           +ci(ic,jc-1,kc,lxzb)*so(i,j-1,k-1,kpnw))
         ENDDO
C$OMP END PARALLEL DO

c     end of computation of i when kf diference operator is 27 point
c******************************

      else ! if kgf.ge.NOG.and.ifd.eq.1

c ***********************************************************************
c     begin computation of i when kf difference operator is seven point.
c *********************************************************************** 

         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF k-planes that sit on FINE-ONLY y-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iicf1,jjc1,kkc1,so)
	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkc1-1)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+2
            k = 2*(kc-1)
            a=so(i-1,j,k,kpw)
            b=so(i,j,k,kpw)
            ep=MIN(abs(a/so(i-1,j,k,kp)),abs(b))
            c=a+b+so(i-1,j,k,kps)+so(i-1,j+1,k,kps)
     &           +so(i-1,j,k,kb)+so(i-1,j,k+1,kb)
            c=a+b+(so(i-1,j,k,kp)-c)*MAX(so(i-1,j,k,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i-1,j,k,kp)
     &           -(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc,lxyl)=a/c
            ci(ic,jc,kc,lxyr)=b/c
         ENDDO
C$OMP END PARALLEL DO

         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF k-planes that sit on FINE-ONLY x-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iic1,jjcf1,kkc1,so)
	 DO kk=1,(iic1-1)*(jjcf1-2)*(kkc1-1)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjcf1-2))+2
            k = 2*(kc-1)
            a=so(i,j,k,kps)
            b=so(i,j-1,k,kps)
            c=a+b+so(i,j-1,k,kpw)+so(i+1,j-1,k,kpw)
     &           +so(i,j-1,k,kb)+so(i,j-1,k+1,kb)
            ep=MIN(abs(a/so(i,j-1,k,kp)),abs(b/so(i,j-1,k,kp)))
            c=a+b+(so(i,j-1,k,kp)-c)*MAX(so(i,j-1,k,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i,j-1,k,kp)
     &           -(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc,lxya)=a/c
            ci(ic,jc,kc,lxyb)=b/c
         ENDDO
C$OMP END PARALLEL DO
         ! 
         ! Set-up interpolant for fine grid points on 
         ! CF j-planes that sit on FINE-ONLY x-lines
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iic1,jjc1,kkcf1,so)
	 DO kk=1,(iic1-1)*(jjc1-1)*(kkcf1-2)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjc1-1))+3
            k = 2*(kc-1)
            a=so(i,j,k,kb)
            b=so(i,j,k-1,kb)
            c=a+b+so(i,j,k-1,kpw)+so(i+1,j,k-1,kpw)
     &           +so(i,j+1,k-1,kps)+so(i,j,k-1,kps)
            ep=MIN(abs(a/so(i,j,k-1,kp)),abs(b/so(i,j,k-1,kp)))
            c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &           -(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc,lxza)=a/c
            ci(ic,jc,kc,lxzb)=b/c
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on
         ! CF k-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY y-lines.
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(dnw,dn,dne,dw,de,dsw,ds,dse,ep,dp,
C$OMP&         i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iicf1,jjcf1,kkc1,so)
	 DO kk=1,(iicf1-2)*(jjcf1-2)*(kkc1-1)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjcf1-2))+2
            k = 2*(kc-1)
            dn=so(i-1,j,k,kps)
            dw=so(i-1,j-1,k,kpw)
            de=so(i,j-1,k,kpw)
            ds=so(i-1,j-1,k,kps)
            dp=dw+dn+de+ds
            sum=so(i-1,j-1,k,kp)-so(i-1,j-1,k,kb)
     &           -so(i-1,j-1,k+1,kb)
            ep=MIN(abs(dw/so(i-1,j-1,k,kp)),
     &           abs(dn/so(i-1,j-1,k,kp)),
     &           abs(de/so(i-1,j-1,k,kp)),
     &           abs(ds/so(i-1,j-1,k,kp)) 
     &           )
            dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)
     &           /(abs(sum-(rONE+ep)*dp)+eMACH)
            dp=rONE/dp
            ci(ic,jc,kc,lxynw)=dp*(ci(ic-1,jc,kc,lxya)*dw
     &           +ci(ic,jc,kc,lxyl)*dn)
            ci(ic,jc,kc,lxyne)=dp*(ci(ic,jc,kc,lxyr)*dn
     &           +ci(ic,jc,kc,lxya)*de)
            ci(ic,jc,kc,lxyse)=dp*(ci(ic,jc,kc,lxyb)*de
     &           +ci(ic,jc-1,kc,lxyr)*ds)
            ci(ic,jc,kc,lxysw)=dp*(ci(ic,jc-1,kc,lxyl)*ds
     &           +ci(ic-1,jc,kc,lxyb)*dw)
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on 
         ! CF j-planes that sit on FINE-ONLY x-lines
         ! and FINE-ONLY z-lines.
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(dnw,dn,dne,dw,de,dsw,ds,dse,ep,dp,
C$OMP&         i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iicf1,jjc1,kkcf1,so)
 	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkcf1-2)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+3
            k = 2*(kc-1)
            dn=so(i-1,j,k,kb)
            dw=so(i-1,j,k-1,kpw)
            de=so(i,j,k-1,kpw)
            ds=so(i-1,j,k-1,kb)
            dp=dw+dn+de+ds
            sum=so(i-1,j,k-1,kp)-so(i-1,j+1,k-1,kps)
     &           -so(i-1,j,k-1,kps)
            ep=MIN(abs(dw/so(i-1,j,k-1,kp)),
     &           abs(dn/so(i-1,j,k-1,kp)),
     &           abs(de/so(i-1,j,k-1,kp)),
     &           abs(ds/so(i-1,j,k-1,kp)) 
     &           )
            dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)
     &           /(abs(sum-(rONE+ep)*dp)+eMACH)
            dp=rONE/dp
            ci(ic,jc,kc,lxznw)=dp*(ci(ic-1,jc,kc,lxza)*dw
     &           +ci(ic,jc,kc,lxyl)*dn)
            ci(ic,jc,kc,lxzne)=dp*(ci(ic,jc,kc,lxyr)*dn
     &           +ci(ic,jc,kc,lxza)*de)
            ci(ic,jc,kc,lxzse)=dp*(ci(ic,jc,kc,lxzb)*de
     &           +ci(ic,jc,kc-1,lxyr)*ds)
            ci(ic,jc,kc,lxzsw)=dp*(ci(ic,jc,kc-1,lxyl)*ds
     &           +ci(ic-1,jc,kc,lxzb)*dw)
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for fine grid points on 
         ! CF i-planes that sit on FINE-ONLY y-lines
         ! and FINE-ONLY z-lines.
         ! 
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(dn,dw,de,ds,dp,ep,
C$OMP&         i,ic,j,jc,k,kc,kk,sum)
C$OMP& SHARED(eMACH,ci,iic1,jjcf1,kkcf1,so)
 	 DO kk=1,(iic1-1)*(jjcf1-2)*(kkcf1-2)
            ic = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iic1-1)*(jjcf1-2))+3
            k = 2*(kc-1)
            dn=so(i,j-1,k,kb)
            dw=so(i,j,k-1,kps)
            de=so(i,j-1,k-1,kps)
            ds=so(i,j-1,k-1,kb)
            dp=dw+dn+de+ds
            sum=so(i,j-1,k-1,kp)-so(i,j-1,k-1,kpw)
     &           -so(i+1,j-1,k-1,kpw)
            ep=MIN(abs(dw/so(i,j-1,k-1,kp)),
     &           abs(dn/so(i,j-1,k-1,kp)),
     &           abs(de/so(i,j-1,k-1,kp)),
     &           abs(ds/so(i,j-1,k-1,kp)) 
     &           )
            dp=dp+(sum-dp)*MAX(sum-(rONE+ep)*dp,rZERO)
     &           /(abs(sum-(rONE+ep)*dp)+eMACH)
            dp=rONE/dp
            ci(ic,jc,kc,lyznw)=dp*(ci(ic,jc,kc,lxza)*dw
     &           +ci(ic,jc,kc,lxya)*dn)
            ci(ic,jc,kc,lyzne)=dp*(ci(ic,jc,kc,lxyb)*dn
     &           +ci(ic,jc-1,kc,lxza)*de)
            ci(ic,jc,kc,lyzse)=dp*(ci(ic,jc-1,kc,lxzb)*de
     &           +ci(ic,jc,kc-1,lxyb)*ds)
            ci(ic,jc,kc,lyzsw)=dp*(ci(ic,jc,kc-1,lxya)*ds
     &           +ci(ic,jc,kc,lxzb)*dw)
         ENDDO
C$OMP END PARALLEL DO
         !
         ! Set-up interpolant for those fine grid points that
         ! sit on FINE-ONLY x-lines, FINE-ONLY y-lines, and 
         ! FINE-ONLY z-lines.
         !
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(ep,dp,i,ic,j,jc,k,kc,kk)
C$OMP& SHARED(eMACH,ci,iicf1,jjcf1,kkcf1,so)
 	 DO kk=1,(iicf1-2)*(jjcf1-2)*(kkcf1-2)
            ic = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iicf1-2)*(jjcf1-2))+3
            k = 2*(kc-1)
            dp=so(i-1,j-1,k-1,kpw)+so(i-1,j,k-1,kps)
     &           +so(i,j-1,k-1,kpw)+so(i-1,j-1,k-1,kps)
     &           +so(i-1,j-1,k-1,kb)+so(i-1,j-1,k,kb)
            ep=MIN(abs(so(i-1,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),
     &           abs(so(i-1,j,k-1,kps)/so(i-1,j-1,k-1,kp)),
     &           abs(so(i,j-1,k-1,kpw)/so(i-1,j-1,k-1,kp)),
     &           abs(so(i-1,j-1,k-1,kps)/so(i-1,j-1,k-1,kp)),
     &           abs(so(i-1,j-1,k-1,kb)/so(i-1,j-1,k-1,kp)),
     &           abs(so(i-1,j-1,k,kb)/so(i-1,j-1,k-1,kp)) 
     &           )
            dp=(so(i-1,j-1,k-1,kp)-dp)*MAX(so(i-1,j-1,k-1,kp)
     &           -(rONE+ep)*dp,rZERO)/(abs(so(i-1,j-1,k-1,kp)
     &           -(rONE+ep)*dp)+eMACH)+dp
            dp=rONE/dp
            ci(ic,jc,kc,ltnw)=dp*(ci(ic-1,jc,kc,lyznw)
     &           *so(i-1,j-1,k-1,kpw)
     &           +ci(ic,jc,kc,lxznw)*so(i-1,j,k-1,kps)
     &           +ci(ic,jc,kc,lxynw)*so(i-1,j-1,k,kb))
            ci(ic,jc,kc,ltne)=dp*(ci(ic,jc,kc,lxzne)
     &           *so(i-1,j,k-1,kps)
     &           +ci(ic,jc,kc,lyznw)*so(i,j-1,k-1,kpw)
     &           +ci(ic,jc,kc,lxyne)*so(i-1,j-1,k,kb))
            ci(ic,jc,kc,lbnw)=dp*(ci(ic,jc,kc-1,lxynw)
     &           *so(i-1,j-1,k-1,kb)
     &           +ci(ic-1,jc,kc,lyzsw)*so(i-1,j-1,k-1,kpw)
     &           +ci(ic,jc,kc,lxzsw)*so(i-1,j,k-1,kps))
            ci(ic,jc,kc,lbne)=dp*(ci(ic,jc,kc-1,lxyne)
     &           *so(i-1,j-1,k-1,kb)
     &           +ci(ic,jc,kc,lxzse)*so(i-1,j,k-1,kps)
     &           +ci(ic,jc,kc,lyzsw)*so(i,j-1,k-1,kpw))
            ci(ic,jc,kc,lbsw)=dp*(ci(ic,jc,kc-1,lxysw)
     &           *so(i-1,j-1,k-1,kb)
     &           +ci(ic-1,jc,kc,lyzse)*so(i-1,j-1,k-1,kpw)
     &           +ci(ic,jc-1,kc,lxzsw)*so(i-1,j-1,k-1,kps))
            ci(ic,jc,kc,ltsw)=dp*(ci(ic,jc,kc,lxysw)
     &           *so(i-1,j-1,k,kb)
     &           +ci(ic-1,jc,kc,lyzne)*so(i-1,j-1,k-1,kpw)
     &           +ci(ic,jc-1,kc,lxznw)*so(i-1,j-1,k-1,kps))
            ci(ic,jc,kc,ltse)=dp*(ci(ic,jc,kc,lxyse)
     &           *so(i-1,j-1,k,kb)
     &           +ci(ic,jc-1,kc,lxzne)*so(i-1,j-1,k-1,kps)
     &           +ci(ic,jc,kc,lyzne)*so(i,j-1,k-1,kpw))
            ci(ic,jc,kc,lbse)=dp*(ci(ic,jc,kc-1,lxyse)
     &           *so(i-1,j-1,k-1,kb)
     &           +ci(ic,jc-1,kc,lxzse)*so(i-1,j-1,k-1,kps)
     &           +ci(ic,jc,kc,lyzse)*so(i,j-1,k-1,kpw))
         enddo
C$OMP END PARALLEL DO

      endif ! of if(kgf.lt.NOG.or.ifd.ne.1)

C ==========================================================================

      RETURN
      END
