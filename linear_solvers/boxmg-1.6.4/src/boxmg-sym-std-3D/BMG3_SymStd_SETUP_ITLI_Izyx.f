      SUBROUTINE BMG3_SymStd_SETUP_ITLI_Izyx( 
     &                kgf, kgc, so, soc, sor, ci,
     &                iif, jjf, kkf, iic, jjc, kkc, 
     &                NOG, ifd, NStncl, irelax, zo, yo 
     &                )

C ========================================================================

C***BEGIN PROLOGUE  MGCF
C***REFER TO  BMG3D
C
C***DESCRIPTION
C             MGCF computes the interpolation operator
C             from the coarse grid to the fine grid. It also calculates
C             the difference operator on the coarse grid and computes
C             the right hand side for the coarse grid if required. If
C             NOG.gt.4, plane relaxation is needed, and the
C             plane relaxation needs line relaxation, for which lu
C             decompositions are computed and stored.
C
C***PARAMETERS
C***INPUT
C   KGF       Fine grid number.
C   KGC       Coarse grid number.
C   IIF       Number. of grid points in x direction on fine grid,
C             including two fictitious points.
C   JJF       Number. of grid points in y direction on fine grid,
C             including two fictitious points.
C   KKF       Number. of grid points in z direction on fine grid,
C             including two fictitious points.
C   IIC       Number. of grid points in x direction on coarse grid,
C             including two fictitious points.
C   JJC       Number. of grid points in y direction on coarse grid,
C             including two fictitious points.
C   KKC       Number. of grid points in z direction on coarse grid,
C             including two fictitious points.
C   NOG     Refer to BMG3D.
C   IFD       Refer to BMG3D.
C   IRELAX    Refer to BMG3D.
C   SO        Refer to BMG3D.
C   SOR       Refer to BMG3D.
C   CI        Refer to BMG3D.
C***OUTPUT
C   SOC       SO for coarse grid
C   SORC      SOR for coarse grid
C***ROUTINES CALLED (NONE)
C***CHANGES 
C  991207     code indented, almost all gotos replaced by if-then-else
C             replaced variables (if,jf,kf) by (i,j,k), by M.Berndt
C***END PROLOGUE  MGCOEF
      IMPLICIT NONE

      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C     ARGUMENTS
      INTEGER ifd, iic, iif, irelax, jjc, jjf, kgc, kgf, kkc, kkf,
     &        NOG, NStncl

      REAL*8  ci(iic,jjc,kkc,26),
     &        so(iif,jjf,kkf,NStncl), soc(iic,jjc,kkc,14), 
     &        sor(iif,jjf,kkf,2), 
     &        yo(iif,jjc,kkc+1,14), zo(iif,jjf,kkc,14)

C     LOCAL VARIABLES
      INTEGER  ic, i, iicf, iicf1, iic1, iif1, jc, j, jjcf, jjcf1, 
     &         jjc1, jjf1, kc, k, kk, kkcf, kkc1, kkcf1, kkcp1, kkf1, 
     &         kpz
      REAL*8   a, b, c, eMACH, ep 

C ========================================================================

C***  FIRST EXECUTABLE STATEMENT  MGCF

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

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk,kpz)
C$OMP& SHARED(iif,jjc,kkcp1,yo)
      DO kk=1,iif*jjc*kkcp1*14
         i=mod(kk-1,iif)+1
         j=mod((kk-1)/iif,jjc)+1
         k=mod((kk-1)/(iif*jjc),2)+1
         kpz=(kk-1)/(iif*jjc*kkcp1)+1
         yo(i,j,k,kpz) = rZERO
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kk,kpz)
C$OMP& SHARED(iif,jjf,kkc,zo)
      DO kk=1,iif*jjf*kkc*14
         i=mod(kk-1,iif)+1
         j=mod((kk-1)/iif,jjf)+1
         k=mod((kk-1)/(iif*jjf),2)+1
         kpz=(kk-1)/(iif*jjf*kkc)+1
         zo(i,j,k,kpz) = rZERO
      ENDDO

C      do kpz=1,14
C         do k=1,kkcp1
C            do j=1,jjc
C               do i=1,iif
C                  yo(i,j,k,kpz)= rZERO
C               enddo
C            enddo
C         enddo
C      enddo

C      do kpz=1,14
C         do k=1,kkc
C            do j=1,jjf
C               do i=1,iif
C                  zo(i,j,k,kpz)= rZERO
C               enddo
C            enddo
C         enddo
C      enddo


c
c   compute interpolation operator iz, from grid kgfz to grid kgf, where
c   grid kgfz is obtained form grid kgf by coarsening in the z-direction
c   only.
c
      if(kgf.lt.NOG.or.ifd.ne.1) then
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,j,k,kc,kk)
C$OMP& SHARED(eMACH,iif1,jjf1,kkcf1,so,sor)
	 DO kk=1,(iif1-1)*(jjf1-1)*(kkcf1-2)
            i = mod((kk-1),(iif1-1))+2
            j = mod((kk-1)/(iif1-1),(jjf1-1))+2
            kc = (kk-1)/((iif1-1)*(jjf1-1))+3
            k = 2*(kc-1)
            a=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)
     &           +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)
     &           +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)
     &           +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)
     &           +so(i+1,j,k-1,kbse)
            b=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)
     &           +so(i+1,j+1,k,kbsw)+so(i,j,k,kbe)+so(i,j,k,kb)
     &           +so(i+1,j,k,kbw)+so(i,j,k,kbne)+so(i,j,k,kbn)
     &           +so(i+1,j,k,kbnw)
            c=a+b+so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)
     &           +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpw)
     &           +so(i+1,j,k-1,kpw)+so(i,j,k-1,kpsw)
     &           +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
            ep=MIN(abs(a),abs(b),rONE)
            c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &           -(rONE+ep)*c)+eMACH)
            sor(i,j,k-1,mtot)=a/c
            sor(i,j,k,mtot)=b/c
         ENDDO            
C$OMP END PARALLEL DO
      else ! if kgf.ge.NOG.and.ifd.eq.1   
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,j,k,kc,kk)
C$OMP& SHARED(eMACH,iif1,jjf1,kkcf1,so,sor)
	 DO kk=1,(iif1-1)*(jjf1-1)*(kkcf1-2)
            i = mod((kk-1),(iif1-1))+2
            j = mod((kk-1)/(iif1-1),(jjf1-1))+2
            kc = (kk-1)/((iif1-1)*(jjf1-1))+3
            k = 2*(kc-1)
            a=so(i,j,k-1,kb)
            b=so(i,j,k,kb)
            c=a+b+so(i,j+1,k-1,kps)+so(i,j,k-1,kpw)
     &           +so(i+1,j,k-1,kpw)+so(i,j,k-1,kps)
            ep=MIN(abs(a),abs(b),rONE)
            c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &           -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &           -(rONE+ep)*c)+eMACH)
            sor(i,j,k-1,mtot)=a/c
            sor(i,j,k,mtot)=b/c
         enddo
C$OMP END PARALLEL DO
      endif
c
c   compute lz = iz(transpose) l iz, where l is the difference operator
c   on grid kgf.
c
      
      if(kgf.lt.NOG.or.ifd.ne.1) then    
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kc,kk)
C$OMP& SHARED(eMACH,iif1,jjf1,kkc1,so,zo,sor)
	 DO kk=1,(iif1-1)*(jjf1-1)*(kkc1-1)
            i = mod((kk-1),(iif1-1))+2
            j = mod((kk-1)/(iif1-1),(jjf1-1))+2
            kc = (kk-1)/((iif1-1)*(jjf1-1))+2
            k = 2*(kc-1)
            zo(i,j,kc,kpw)=
     &           so(i,j,k,kpw)
     &           +so(i,j,k+1,kbe)*sor(i-1,j,k+1,mtot)
     &           +so(i,j,k,kbw)*sor(i-1,j,k,mtot)
     &           +sor(i,j,k+1,mtot)
     &           *( so(i,j,k+1,kpw)*sor(i-1,j,k+1,mtot)
     &           +so(i,j,k+1,kbw))
     &           +sor(i,j,k,mtot) 
     &           *( so(i,j,k-1,kpw)*sor(i-1,j,k,mtot)
     &           +so(i,j,k,kbe))

            zo(i,j,kc,kpsw)=so(i,j,k,kpsw)+so(i,j,k+1,kbne)
     &           *sor(i-1,j-1,k+1,mtot)+so(i,j,k,kbsw)
     &           *sor(i-1,j-1,k,mtot)+sor(i,j,k,mtot)
     &           *(so(i,j,k-1,kpsw)*sor(i-1,j-1,k,mtot)
     &           +so(i,j,k,kbne))+sor(i,j,k+1,mtot)
     &           *(so(i,j,k+1,kpsw)*sor(i-1,j-1,k+1,mtot)
     &           +so(i,j,k+1,kbsw))

            zo(i,j+1,kc,kpnw)=so(i,j+1,k,kpnw)+so(i,j+1,k,kbnw)
     &           *sor(i-1,j+1,k,mtot)+so(i,j+1,k+1,kbse)
     &           *sor(i-1,j+1,k+1,mtot)+sor(i,j,k,mtot)
     &           *(so(i,j+1,k-1,kpnw)*sor(i-1,j+1,k,mtot)
     &           +so(i,j+1,k,kbse))+sor(i,j,k+1,mtot)
     &           *(so(i,j+1,k+1,kpnw)*sor(i-1,j+1,k+1,mtot)
     &           +so(i,j+1,k+1,kbnw))

            zo(i,j,kc,kps)=so(i,j,k,kps)+so(i,j,k+1,kbn)
     &           *sor(i,j-1,k+1,mtot)+so(i,j,k,kbs)
     &           *sor(i,j-1,k,mtot)+sor(i,j,k,mtot)
     &           *(so(i,j,k-1,kps)*sor(i,j-1,k,mtot)
     &           +so(i,j,k,kbn))+sor(i,j,k+1,mtot)
     &           *(so(i,j,k+1,kps)*sor(i,j-1,k+1,mtot)
     &           +so(i,j,k+1,kbs))

            zo(i,j,kc,kp)=so(i,j,k,kp)-so(i,j,k,kb)
     &           *sor(i,j,k,mtot)-so(i,j,k+1,kb)
     &           *sor(i,j,k+1,mtot)-sor(i,j,k,mtot)
     &           *(-so(i,j,k-1,kp)*sor(i,j,k,mtot)
     &           +so(i,j,k,kb))-sor(i,j,k+1,mtot)
     &           *(-so(i,j,k+1,kp)*sor(i,j,k+1,mtot)
     &           +so(i,j,k+1,kb))

            zo(i,j,kc,kb)=so(i,j,k,kb)*sor(i,j,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i,j,k-1,kb)
     &           -so(i,j,k-1,kp)*sor(i,j,k-1,mtot))

            zo(i,j,kc,kbw)=so(i,j,k,kbw)*sor(i-1,j,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i,j,k-1,kpw)
     &           *sor(i-1,j,k-1,mtot)+so(i,j,k-1,kbw))

            zo(i,j+1,kc,kbnw)=so(i,j+1,k,kbnw)
     &           *sor(i-1,j+1,k-1,mtot)+sor(i,j,k,mtot)
     &           *(so(i,j+1,k-1,kpnw)*sor(i-1,j+1,k-1,mtot)
     &           +so(i,j+1,k-1,kbnw))

            zo(i,j+1,kc,kbn)=so(i,j+1,k,kbn)
     &           *sor(i,j+1,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i,j+1,k-1,kps)
     &           *sor(i,j+1,k-1,mtot)+so(i,j+1,k-1,kbn))

            zo(i+1,j+1,kc,kbne)=so(i+1,j+1,k,kbne)
     &                 *sor(i+1,j+1,k-1,mtot)+sor(i,j,k,mtot)
     &                 *(so(i+1,j+1,k-1,kpsw)*sor(i+1,j+1,k-1,mtot)
     &                 +so(i+1,j+1,k-1,kbne))

            zo(i+1,j,kc,kbe)=so(i+1,j,k,kbe)
     &           *sor(i+1,j,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i+1,j,k-1,kpw)
     &           *sor(i+1,j,k-1,mtot)+so(i+1,j,k-1,kbe))

            zo(i+1,j,kc,kbse)=so(i+1,j,k,kbse)
     &           *sor(i+1,j-1,k-1,mtot)+sor(i,j,k,mtot)
     &           *(so(i+1,j,k-1,kpnw)*sor(i+1,j-1,k-1,mtot)
     &           +so(i+1,j,k-1,kbse))

            zo(i,j,kc,kbs)=so(i,j,k,kbs)*sor(i,j-1,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i,j,k-1,kps)
     &           *sor(i,j-1,k-1,mtot)+so(i,j,k-1,kbs))

            zo(i,j,kc,kbsw)=so(i,j,k,kbsw)
     &           *sor(i-1,j-1,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i,j,k-1,kpsw)
     &           *sor(i-1,j-1,k-1,mtot)+so(i,j,k-1,kbsw))

         ENDDO
C$OMP END PARALLEL DO

      else                      ! if kgf.ge.NOG.and.ifd.eq.1                
         
C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,j,k,kc,kk)
C$OMP& SHARED(eMACH,iif1,jjf1,kkc1,so,zo,sor)
	 DO kk=1,(iif1-1)*(jjf1-1)*(kkc1-1)
            i = mod((kk-1),(iif1-1))+2
            j = mod((kk-1)/(iif1-1),(jjf1-1))+2
            kc = (kk-1)/((iif1-1)*(jjf1-1))+2
            k = 2*(kc-1)

            zo(i,j,kc,kpw)=so(i,j,k,kpw)+sor(i,j,k+1,mtot)
     &           *so(i,j,k+1,kpw)*sor(i-1,j,k+1,mtot)
     &           +sor(i,j,k,mtot)*so(i,j,k-1,kpw)
     &           *sor(i-1,j,k,mtot)
            zo(i,j,kc,kpsw)= rZERO
            zo(i,j+1,kc,kpnw)= rZERO
            zo(i,j,kc,kps)=so(i,j,k,kps)+sor(i,j,k,mtot)
     &           *so(i,j,k-1,kps)*sor(i,j-1,k,mtot)
     &           +sor(i,j,k+1,mtot)*so(i,j,k+1,kps)
     &           *sor(i,j-1,k+1,mtot)
            zo(i,j+1,kc,kpnw)= rZERO
            zo(i,j,kc,kp)=so(i,j,k,kp)-sor(i,j,k,mtot)
     &           *(-so(i,j,k-1,kp)*sor(i,j,k,mtot)
     &           +so(i,j,k,kb))-sor(i,j,k+1,mtot)
     &           *(-so(i,j,k+1,kp)*sor(i,j,k+1,mtot)
     &           +so(i,j,k+1,kb))-so(i,j,k,kb)*sor(i,j,k,mtot)
     &           -so(i,j,k+1,kb)*sor(i,j,k+1,mtot)

            zo(i,j,kc,kb)=so(i,j,k,kb)*sor(i,j,k-1,mtot)
     &           +sor(i,j,k,mtot)*(so(i,j,k-1,kb)
     &           -so(i,j,k-1,kp)*sor(i,j,k-1,mtot))

            zo(i,j,kc,kbw)=sor(i,j,k,mtot)*so(i,j,k-1,kpw)
     &           *sor(i-1,j,k-1,mtot)
            zo(i,j+1,kc,kbnw)= rZERO

            zo(i,j+1,kc,kbn)=sor(i,j,k,mtot)*so(i,j+1,k-1,kps)
     &           *sor(i,j+1,k-1,mtot)
            zo(i+1,j+1,kc,kbne)= rZERO


            zo(i+1,j,kc,kbe)=sor(i,j,k,mtot)*so(i+1,j,k-1,kpw)
     &           *sor(i+1,j,k-1,mtot)
            zo(i+1,j,kc,kbse)= rZERO

            zo(i,j,kc,kbs)=sor(i,j,k,mtot)*so(i,j,k-1,kps)
     &           *sor(i,j-1,k-1,mtot)
            zo(i,j,kc,kbsw)= rZERO

         ENDDO
C$OMP END PARALLEL DO

      endif                     ! of if (kgf.lt.NOG.or.ifd.ne.1)


c     compute iy, the interpolation operator form grid kgzy to grid
c     kgz, where grid kgz is obtained from grid kgz by coarsening
c     in the y-direction only.
c     

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,j,jc,kc,kk)
C$OMP& SHARED(eMACH,iif1,jjcf1,kkc,so,zo,sor)
	 DO kk=1,(iif1-1)*(jjcf1-2)*(kkc-2)
            i  = mod((kk-1),(iif1-1))+2
            jc = mod((kk-1)/(iif1-1),(jjcf1-2))+3
            j = 2*(jc-1)
            kc = (kk-1)/((iif1-1)*(jjcf1-2))+3

            a=zo(i,j-1,kc-1,kbsw)+zo(i,j-1,kc-1,kbs)
     &           +zo(i+1,j-1,kc-1,kbse)+zo(i,j-1,kc-1,kpsw)
     &           +zo(i,j-1,kc-1,kps)+zo(i+1,j-1,kc-1,kpnw)
     &           +zo(i,j-1,kc,kbne)+zo(i,j-1,kc,kbn)
     &           +zo(i+1,j-1,kc,kbnw)
            b=zo(i,j,kc-1,kbnw)+zo(i,j,kc-1,kbn)
     &           +zo(i+1,j,kc-1,kbne)+zo(i,j,kc-1,kpnw)
     &           +zo(i,j,kc-1,kps)+zo(i+1,j,kc-1,kpsw)
     &           +zo(i,j,kc,kbse)+zo(i,j,kc,kbs)
     &           +zo(i+1,j,kc,kbsw)
            c=zo(i,j-1,kc-1,kbw)+zo(i,j-1,kc-1,kb)
     &           +zo(i+1,j-1,kc-1,kbe)+zo(i,j-1,kc-1,kpw)
     &           +a+b+zo(i+1,j-1,kc-1,kpw)
     &           +zo(i,j-1,kc,kbe)+zo(i,j-1,kc,kb)
     &           +zo(i+1,j-1,kc,kbw)
            ep=MIN(abs(a),abs(b),rONE)
            c=a+b+(zo(i,j-1,kc-1,kp)-c)*MAX(zo(i,j-1,kc-1,kp)
     &           -(rONE+ep)*c,rZERO)
     &           /(abs(zo(i,j-1,kc-1,kp)-(rONE+ep)*c)+eMACH)
            sor(i,j-1,kc,msor)=a/c
            sor(i,j,kc,msor)=b/c
            
         ENDDO
C$OMP END PARALLEL DO

c     
c     compute ly = iy(transpose) lz iy.
c     

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,j,jc,kc,kk)
C$OMP& SHARED(eMACH,iif1,jjc1,kkc,sor,yo,zo)
	 DO kk=1,(iif1-1)*(jjc1-1)*(kkc-2)
            i  = mod((kk-1),(iif1-1))+2
            jc = mod((kk-1)/(iif1-1),(jjc1-1))+2
            j = 2*(jc-1)
            kc = (kk-1)/((iif1-1)*(jjc1-1))+3
            yo(i,jc,kc,kp)=zo(i,j,kc-1,kp)-zo(i,j,kc-1,kps)
     &           *sor(i,j,kc,msor)-zo(i,j+1,kc-1,kps)
     &           *sor(i,j+1,kc,msor)-sor(i,j,kc,msor)
     &           *(-zo(i,j-1,kc-1,kp)*sor(i,j,kc,msor)
     &           +zo(i,j,kc-1,kps))-sor(i,j+1,kc,msor)
     &           *(-zo(i,j+1,kc-1,kp)*sor(i,j+1,kc,msor)
     &           +zo(i,j+1,kc-1,kps))

            yo(i,jc,kc,kpw)=zo(i,j,kc-1,kpw)+zo(i,j,kc-1,kpsw)
     &           *sor(i-1,j,kc,msor)+zo(i,j+1,kc-1,kpnw)
     &           *sor(i-1,j+1,kc,msor)+sor(i,j,kc,msor)
     &           *(zo(i,j-1,kc-1,kpw)*sor(i-1,j,kc,msor)
     &           +zo(i,j,kc-1,kpnw))+sor(i,j+1,kc,msor)
     &           *(zo(i,j+1,kc-1,kpw)*sor(i-1,j+1,kc,msor)
     &           +zo(i,j+1,kc-1,kpsw))

            yo(i,jc,kc,kpnw)=zo(i,j-1,kc-1,kpnw)
     &           *sor(i-1,j,kc,msor)
     &           +sor(i,j-1,kc,msor)*(zo(i,j,kc-1,kpnw)
     &           +zo(i,j-1,kc-1,kpw)*sor(i-1,j,kc,msor))
            
            yo(i,jc,kc,kps)=zo(i,j,kc-1,kps)*sor(i,j-1,kc,msor)
     &              +sor(i,j,kc,msor)*(zo(i,j-1,kc-1,kps)
     &              -zo(i,j-1,kc-1,kp)*sor(i,j-1,kc,msor))

            yo(i,jc,kc,kpsw)=zo(i,j,kc-1,kpsw)
     &           *sor(i-1,j-1,kc,msor)
     &           +sor(i,j,kc,msor)*(zo(i,j-1,kc-1,kpsw)
     &           +zo(i,j-1,kc-1,kpw)*sor(i-1,j-1,kc,msor))

            yo(i,jc,kc,kb)=zo(i,j,kc-1,kb)+sor(i,j,kc,msor)
     &           *(zo(i,j,kc-1,kbn)+zo(i,j-1,kc-1,kb)
     &           *sor(i,j,kc-1,msor))+sor(i,j+1,kc,msor)
     &           *(zo(i,j+1,kc-1,kbs)+zo(i,j+1,kc-1,kb)
     &           *sor(i,j+1,kc-1,msor))+zo(i,j,kc-1,kbs)
     &           *sor(i,j,kc-1,msor)+zo(i,j+1,kc-1,kbn)
     &           *sor(i,j+1,kc-1,msor)

            yo(i,jc,kc,kbw)=zo(i,j,kc-1,kbw)+zo(i,j,kc-1,kbsw)
     &           *sor(i-1,j,kc-1,msor)+zo(i,j+1,kc-1,kbnw)
     &           *sor(i-1,j+1,kc-1,msor)+sor(i,j,kc,msor)
     &           *(zo(i,j,kc-1,kbnw)+zo(i,j-1,kc-1,kbw)
     &           *sor(i-1,j,kc-1,msor))+sor(i,j+1,kc,msor)
     &           *(zo(i,j+1,kc-1,kbsw)+zo(i,j+1,kc-1,kbw)
     &           *sor(i-1,j+1,kc-1,msor))
            
            yo(i,jc,kc,kbnw)=zo(i,j-1,kc-1,kbnw)
     &           *sor(i-1,j,kc-1,msor)+sor(i,j-1,kc,msor)
     &           *(zo(i,j-1,kc-1,kbw)*sor(i-1,j,kc-1,msor)
     &           +zo(i,j,kc-1,kbnw))

            yo(i,jc,kc,kbn)=zo(i,j-1,kc-1,kbn)
     &           *sor(i,j,kc-1,msor)
     &           +sor(i,j-1,kc,msor)*(zo(i,j,kc-1,kbn)
     &           +zo(i,j-1,kc-1,kb)*sor(i,j,kc-1,msor))

            yo(i,jc,kc,kbn)=zo(i,j-1,kc-1,kbn)
     &           *sor(i,j,kc-1,msor)
     &           +sor(i,j-1,kc,msor)*(zo(i,j,kc-1,kbn)
     &           +zo(i,j-1,kc-1,kb)*sor(i,j,kc-1,msor))
            
            yo(i+1,jc,kc,kbne)=zo(i+1,j-1,kc-1,kbne)
     &           *sor(i+1,j,kc-1,msor)+sor(i,j-1,kc,msor)
     &           *(zo(i+1,j,kc-1,kbne)+zo(i+1,j-1,kc-1,kbe)
     &           *sor(i+1,j,kc-1,msor))

            yo(i+1,jc,kc,kbe)=zo(i+1,j,kc-1,kbe)
     &           +zo(i+1,j,kc-1,kbse)
     &           *sor(i+1,j,kc-1,msor)+zo(i+1,j+1,kc-1,kbne)
     &           *sor(i+1,j+1,kc-1,msor)+sor(i,j,kc,msor)
     &           *(zo(i+1,j,kc-1,kbne)+zo(i+1,j-1,kc-1,kbe)
     &           *sor(i+1,j,kc-1,msor))+sor(i,j+1,kc,msor)
     &           *(zo(i+1,j+1,kc-1,kbse)+zo(i+1,j+1,kc-1,kbe)
     &           *sor(i+1,j+1,kc-1,msor))

            yo(i+1,jc,kc,kbse)=zo(i+1,j,kc-1,kbse)
     &           *sor(i+1,j-1,kc-1,msor)+sor(i,j,kc,msor)
     &           *(zo(i+1,j-1,kc-1,kbse)+zo(i+1,j-1,kc-1,kbe)
     &           *sor(i+1,j-1,kc-1,msor))

            yo(i,jc,kc,kbs)=zo(i,j,kc-1,kbs)
     &           *sor(i,j-1,kc-1,msor)
     &           +sor(i,j,kc,msor)*(zo(i,j-1,kc-1,kbs)
     &           +zo(i,j-1,kc-1,kb)*sor(i,j-1,kc-1,msor))

            yo(i,jc,kc,kbsw)=zo(i,j,kc-1,kbsw)
     &           *sor(i-1,j-1,kc-1,msor)+sor(i,j,kc,msor)
     &           *(zo(i,j-1,kc-1,kbsw)+zo(i,j-1,kc-1,kbw)
     &           *sor(i-1,j-1,kc-1,msor))

         ENDDO
C$OMP END PARALLEL DO

c     
c     compute ix, the interpolation operator form grid kgzy to grid kgc
c     

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(a,b,c,ep,i,ic,jc,kc,kk)
C$OMP& SHARED(ci,eMACH,iicf1,jjc1,kkcp1,yo)
	 DO kk=1,(iicf1-2)*(jjc1-1)*(kkcp1-3)
            ic  = mod((kk-1),(iicf1-2))+3
            i = 2*(ic-1)
            jc = mod((kk-1)/(iicf1-2),(jjc1-1))+2
            kc = (kk-1)/((iicf1-2)*(jjc1-1))+4
            a=yo(i-1,jc+1,kc-1,kbnw)+yo(i-1,jc,kc-1,kbw)
     &           +yo(i-1,jc,kc-1,kbsw)+yo(i-1,jc+1,kc-1,kpnw)
     &           +yo(i-1,jc,kc-1,kpw)+yo(i-1,jc,kc-1,kpsw)
     &           +yo(i-1,jc+1,kc,kbse)+yo(i-1,jc,kc,kbe)
     &           +yo(i-1,jc,kc,kbne)
            b=yo(i,jc+1,kc-1,kbne)+yo(i,jc,kc-1,kbe)
     &           +yo(i,jc,kc-1,kbse)
     &           +yo(i,jc+1,kc-1,kpsw)+yo(i,jc,kc-1,kpw)
     &           +yo(i,jc,kc-1,kpnw)+yo(i,jc+1,kc,kbsw)
     &           +yo(i,jc,kc,kbw)+yo(i,jc,kc,kbnw)
            c=yo(i-1,jc+1,kc-1,kbn)+yo(i-1,jc,kc-1,kb)
     &           +yo(i-1,jc,kc-1,kbs)+yo(i-1,jc+1,kc-1,kps)+a+b
     &           +yo(i-1,jc,kc-1,kps)+yo(i-1,jc+1,kc,kbs)
     &           +yo(i-1,jc,kc,kb)+yo(i-1,jc,kc,kbn)
            ep=MIN(abs(a),abs(b),rONE)
            c=a+b+(yo(i-1,jc,kc-1,kp)-c)*MAX(yo(i-1,jc,kc-1,kp)
     &           -(rONE+ep)*c,rZERO)
     &           /(abs(yo(i-1,jc,kc-1,kp)-(rONE+ep)*c)+eMACH)
            ci(ic,jc,kc-2,lxyl)=a/c
            ci(ic,jc,kc-2,lxyr)=b/c
         ENDDO
C$OMP END PARALLEL DO

c               
c     compute lc = ix (transpose) ly ix.
c     

C$OMP PARALLEL DO
C$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(i,ic,jc,kc,kk)
C$OMP& SHARED(ci,eMACH,iic1,jjc1,kkcp1,soc,yo)
	 DO kk=1,(iic1-1)*(jjc1-1)*(kkcp1-3)
            ic  = mod((kk-1),(iic1-1))+2
            i = 2*(ic-1)
            jc = mod((kk-1)/(iic1-1),(jjc1-1))+2
            kc = (kk-1)/((iic1-1)*(jjc1-1))+4
            soc(ic,jc,kc-2,kp)=yo(i,jc,kc-1,kp)-yo(i,jc,kc-1,kpw)
     &           *ci(ic,jc,kc-2,lxyr)-yo(i+1,jc,kc-1,kpw)
     &           *ci(ic+1,jc,kc-2,lxyl)-ci(ic,jc,kc-2,lxyr)
     &           *(-yo(i-1,jc,kc-1,kp)*ci(ic,jc,kc-2,lxyr)
     &           +yo(i,jc,kc-1,kpw))-ci(ic+1,jc,kc-2,lxyl)
     &           *(-yo(i+1,jc,kc-1,kp)*ci(ic+1,jc,kc-2,lxyl)
     &           +yo(i+1,jc,kc-1,kpw))

            soc(ic,jc,kc-2,kpw)=yo(i,jc,kc-1,kpw)
     &           *ci(ic,jc,kc-2,lxyl)
     &           +ci(ic,jc,kc-2,lxyr)*(yo(i-1,jc,kc-1,kpw)
     &           -yo(i-1,jc,kc-1,kp)*ci(ic,jc,kc-2,lxyl))

            soc(ic,jc+1,kc-2,kpnw)=yo(i,jc+1,kc-1,kpnw)
     &           *ci(ic,jc+1,kc-2,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i-1,jc+1,kc-1,kpnw)+yo(i-1,jc+1,kc-1,kps)
     &           *ci(ic,jc+1,kc-2,lxyl))

            soc(ic,jc,kc-2,kpsw)=yo(i,jc,kc-1,kpsw)
     &           *ci(ic,jc-1,kc-2,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i-1,jc,kc-1,kpsw)+yo(i-1,jc,kc-1,kps)
     &           *ci(ic,jc-1,kc-2,lxyl)) 

            soc(ic,jc,kc-2,kps)=yo(i,jc,kc-1,kps)+yo(i,jc,kc-1,kpsw)
     &           *ci(ic,jc-1,kc-2,lxyr)+yo(i+1,jc,kc-1,kpnw)
     &           *ci(ic+1,jc-1,kc-2,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i,jc,kc-1,kpnw)+yo(i-1,jc,kc-1,kps)
     &           *ci(ic,jc-1,kc-2,lxyr))+ci(ic+1,jc,kc-2,lxyl)
     &           *(yo(i+1,jc,kc-1,kpsw)+yo(i+1,jc,kc-1,kps)
     &           *ci(ic+1,jc-1,kc-2,lxyl))

            soc(ic,jc,kc-2,kb)=yo(i,jc,kc-1,kb)
     &           +ci(ic,jc,kc-2,lxyr)
     &           *(yo(i,jc,kc-1,kbe)+yo(i-1,jc,kc-1,kb)
     &           *ci(ic,jc,kc-3,lxyr))+ci(ic+1,jc,kc-2,lxyl)
     &           *(yo(i+1,jc,kc-1,kbw)+yo(i+1,jc,kc-1,kb)
     &           *ci(ic+1,jc,kc-3,lxyl))+yo(i,jc,kc-1,kbw)
     &           *ci(ic,jc,kc-3,lxyr)+yo(i+1,jc,kc-1,kbe)
     &           *ci(ic+1,jc,kc-3,lxyl)

            soc(ic,jc,kc-2,kbw)=yo(i,jc,kc-1,kbw)
     &           *ci(ic,jc,kc-3,lxyl)
     &           +ci(ic,jc,kc-2,lxyr)*(yo(i-1,jc,kc-1,kbw)
     &           +yo(i-1,jc,kc-1,kb)*ci(ic,jc,kc-3,lxyl))

            soc(ic,jc+1,kc-2,kbnw)=yo(i,jc+1,kc-1,kbnw)
     &           *ci(ic,jc+1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i-1,jc+1,kc-1,kbnw)+yo(i-1,jc+1,kc-1,kbn)
     &           *ci(ic,jc+1,kc-3,lxyl))

            soc(ic,jc+1,kc-2,kbnw)=yo(i,jc+1,kc-1,kbnw)
     &           *ci(ic,jc+1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i-1,jc+1,kc-1,kbnw)+yo(i-1,jc+1,kc-1,kbn)
     &           *ci(ic,jc+1,kc-3,lxyl))

            soc(ic,jc+1,kc-2,kbn)=yo(i,jc+1,kc-1,kbn)
     &              +yo(i,jc+1,kc-1,kbnw)*ci(ic,jc+1,kc-3,lxyr)
     &              +yo(i+1,jc+1,kc-1,kbne)*ci(ic+1,jc+1,kc-3,lxyl)
     &              +ci(ic,jc,kc-2,lxyr)*(yo(i,jc+1,kc-1,kbne)
     &              +yo(i-1,jc+1,kc-1,kbn)*ci(ic,jc+1,kc-3,lxyr))
     &              +ci(ic+1,jc,kc-2,lxyl)*(yo(i+1,jc+1,kc-1,kbnw)
     &              +yo(i+1,jc+1,kc-1,kbn)*ci(ic+1,jc+1,kc-3,lxyl))
           
            soc(ic,jc+1,kc-2,kbne)=yo(i-1,jc+1,kc-1,kbne)
     &           *ci(ic,jc+1,kc-3,lxyr)+ci(ic,jc,kc-2,lxyl)
     &           *(yo(i,jc+1,kc-1,kbne)+yo(i-1,jc+1,kc-1,kbn)
     &           *ci(ic,jc+1,kc-3,lxyr))

            soc(ic,jc,kc-2,kbe)=yo(i-1,jc,kc-1,kbe)
     &           *ci(ic,jc,kc-3,lxyr)+ci(ic,jc,kc-2,lxyl)
     &           *(yo(i,jc,kc-1,kbe)+yo(i-1,jc,kc-1,kb)
     &           *ci(ic,jc,kc-3,lxyr)) 

            soc(ic,jc,kc-2,kbse)=yo(i-1,jc,kc-1,kbse)
     &              *ci(ic,jc-1,kc-3,lxyr)+ci(ic,jc,kc-2,lxyl)
     &              *(yo(i,jc,kc-1,kbse)+yo(i-1,jc,kc-1,kbs)
     &              *ci(ic,jc-1,kc-3,lxyr))

            soc(ic,jc,kc-2,kbs)=yo(i,jc,kc-1,kbs)+yo(i,jc,kc-1,kbsw)
     &           *ci(ic,jc-1,kc-3,lxyr)+yo(i+1,jc,kc-1,kbse)
     &           *ci(ic+1,jc-1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i,jc,kc-1,kbse)+yo(i-1,jc,kc-1,kbs)
     &           *ci(ic,jc-1,kc-3,lxyr))+ci(ic+1,jc,kc-2,lxyl)
     &           *(yo(i+1,jc,kc-1,kbsw)+yo(i+1,jc,kc-1,kbs)
     &           *ci(ic+1,jc-1,kc-3,lxyl))  

            soc(ic,jc,kc-2,kbsw)=yo(i,jc,kc-1,kbsw)
     &           *ci(ic,jc-1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &           *(yo(i-1,jc,kc-1,kbsw)+yo(i-1,jc,kc-1,kbs)
     &           *ci(ic,jc-1,kc-3,lxyl))

         ENDDO
C$OMP END PARALLEL DO

      return
      end






