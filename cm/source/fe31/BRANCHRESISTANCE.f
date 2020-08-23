      SUBROUTINE BRANCHRESISTANCE(NBJ,NEELEM,NPNE,NVJE,CE,
     &  XP,ERROR,*)
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn' 
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'pulm00.cmn' 
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 alpha,CE(NMM,NEM),
     &  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

      INTEGER nb,ne,noelem,np,nv1,nv2
      REAL*8 C_ZETA,gamma,length,Re,R1,R2,radius,Vdot,zeta
      LOGICAL END1
      
      REAL*8 LENGTH_1D,RADIUS_1D

      CALL ENTERS('BRANCHRESISTANCE',*9999)

      alpha = 4.d0/3.d0
      C_ZETA=1.85d0 !constant for turbulence correction factor
      IF(XP(1,1,nj_flow,1).GE.0.d0)THEN
        gamma = 0.357d0 !inspiration
      ELSE
        gamma = 0.460d0 !expiration
      ENDIF

      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        np=NPNE(2,nb,ne)
        length=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP) !mm
        nv1=NVJE(1,nb,nj_radius,ne)
        nv2=NVJE(2,nb,nj_radius,ne)
        R1=XP(1,nv1,nj_radius,NPNE(1,nb,ne))
        R2=XP(1,nv2,nj_radius,NPNE(2,nb,ne))
        radius=0.5d0*(R1+R2) !mm
C.......Resistance in units of Pa.s.mm-3        
        CE(nm_R,ne)=8.d0*GAS_VISCOSITY*length/(PI*radius**4) !laminar resistance
        Vdot = XP(1,1,nj_flow,np) !flow in mm^3
        Re=DABS(Vdot*2.d0*GAS_DENSITY/(PI*radius*GAS_VISCOSITY))
        zeta = DSQRT(2.d0*radius*Re/length)*gamma
        IF(zeta.LT.1.d0)zeta=1.d0
        IF(GAS_VISCOSITY.LT.1.d-3)THEN
          CE(nm_R,ne)=CE(nm_R,ne)*zeta !Turbulent Resistance
        ENDIF

        CE(nm_Rt,ne)=CE(nm_R,ne)

      ENDDO

      CALL EXITS('BRANCHRESISTANCE')
      RETURN
 9999 CALL ERRORS('BRANCHRESISTANCE',ERROR)
      CALL EXITS('BRANCHRESISTANCE')
      RETURN 1
      END
