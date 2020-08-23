      SUBROUTINE SHEARSTRESS(NBJ,NEELEM,NPNE,NVJE,XP,ERROR,*)

      IMPLICIT NONE
      INCLUDE 'b00.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'lung00.cmn'  
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      INTEGER nb,ne,noelem,np1,np2,nv1,nv2
      REAL*8 radius,velocity

      CALL ENTERS('SHEARSTRESS',*9999)

C  Calculate the shear stress
      nj_Pe=9
      nj_mu=10
      
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np1=NPNE(1,nb,ne) !end node
        np2=NPNE(2,nb,ne) !end node
        nv1=NVJE(1,nb,nj_radius,ne) !should be using nj_flow
        nv2=NVJE(2,nb,nj_radius,ne)

C.......Shear stress = dynamic_viscosity*4*flow/(pi*r^3)
        XP(1,nv1,nj_mu,np1)=GAS_VISCOSITY*4.d0*XP(1,nv1,nj_flow,np1)/
     &    (PI*XP(1,nv1,nj_radius,np1)**3) !Pa
        XP(1,nv2,nj_mu,np2)=GAS_VISCOSITY*4.d0*XP(1,nv2,nj_flow,np2)/
     &    (PI*XP(1,nv2,nj_radius,np2)**3) !Pa

C.......Pe = Re * Sc = velocity * diameter / diffusivity
c        XP(1,nv1,nj_Pe,np1)=2.d0*XP(1,nv1,nj_flow,np1)/
c     &    (PI*XP(1,nv2,nj_radius,np2))
c        XP(1,nv2,nj_Pe,np2)=2.d0*XP(1,nv2,nj_flow,np2)/
c     &    (PI*XP(1,nv2,nj_radius,np2))

      ENDDO

      CALL EXITS('SHEARSTRESS')
      RETURN
 9999 CALL ERRORS('SHEARSTRESS',ERROR)
      CALL EXITS('SHEARSTRESS')
      RETURN 1
      END

