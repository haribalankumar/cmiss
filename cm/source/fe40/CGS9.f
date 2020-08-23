      SUBROUTINE CGS9(nr,nw,nx,CG,EG,RT,TG,ZG,ERROR,*)

C#### Subroutine: CGS9
C###  Description:
C###    CGS9 evaluates strain (EG) and stress (TG) for 3D elasticity.

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp40.cmn'
!     Parameter List
      INTEGER nr,nw,nx
      REAL*8 CG(NMM),EG(3,3),RT(3,3),TG(3,3),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k,NU1(0:3)
      REAL*8 CC(6,6),C1,C2,C3,
     '  DENS,POIS1,POIS2,S_thermal(3),
     '  YMOD1,YMOD2

      DATA NU1/1,2,4,7/

      CALL ENTERS('CGS9',*9999)

! Compute strains referred to r.c. coords. Note: ZG has derivs wrt
C GMH 10/1/96 When material coords, use extended version to get
C             strains (include dX_k/dnu_j).
C      DO ni1=1,3
C        DO ni2=1,3
C          EG(ni1,ni2)=0.5d0*(ZG(ni1,NU1(ni2))+ZG(ni2,NU1(ni1)))
C        ENDDO !ni2
C      ENDDO !ni1
      DO i=1,3
        DO j=1,3
          EG(i,j)=0.0d0
          DO k=1,3
            EG(i,j)=EG(i,j)+0.5d0*
     '        (RT(k,i)*ZG(k,NU1(j))+RT(k,j)*ZG(k,NU1(i)))
          ENDDO !k
        ENDDO !j
      ENDDO !i

!     EG(1,1)=dX_dNu(1,1)*ZG(1,NU1(1))+dX_dNu(2,1)*ZG(2,NU1(1))
!     EG(2,2)=dX_dNu(1,2)*ZG(1,NU1(2))+dX_dNu(2,2)*ZG(2,NU1(2))
!     EG(1,2)=0.5d0*(dX_dNu(1,1)*ZG(1,NU1(2))+dX_dNu(2,1)*ZG(2,NU1(2))
!    '             + dX_dNu(1,2)*ZG(1,NU1(1))+dX_dNu(2,2)*ZG(2,NU1(1)))
!     EG(2,1)=EG(1,2)

! Isotropic
      IF(IMT(nw).EQ.1.OR.IMT(nw).EQ.5) THEN
        CALL CGC9(nr,nw,nx,CC,C1,C2,C3,CG,DENS,S_thermal,ERROR,*9999)
        TG(1,1)=C1*EG(1,1)+C2*EG(2,2)+C2*EG(3,3)
        TG(2,2)=C2*EG(1,1)+C1*EG(2,2)+C2*EG(3,3)
        TG(3,3)=C2*EG(1,1)+C2*EG(2,2)+C1*EG(3,3)
        TG(1,2)=C3*EG(1,2)
        TG(2,3)=C3*EG(2,3)
        TG(3,1)=C3*EG(3,1)
        TG(2,1)=TG(1,2)
        TG(3,2)=TG(2,3)
        TG(1,3)=TG(3,1)
C GMH 5/1/96 Old code - use CGC9 to get coefficients
C        YMOD=CG(1)  !is Young's modulus
C        POIS=CG(2)  !is Poisson's ratio nu
C        RLAMDA=YMOD*POIS/((1.d0+POIS)*(1.d0-2.d0*POIS))
C        RMU   =YMOD/(2.d0*(1.d0+POIS))
C        DIL=EG(1,1)+EG(2,2)+EG(3,3)
C        TG(1,1)=RLAMDA*DIL+2.d0*RMU*EG(1,1)
C        TG(2,2)=RLAMDA*DIL+2.d0*RMU*EG(2,2)
C        TG(3,3)=RLAMDA*DIL+2.d0*RMU*EG(3,3)
C        TG(1,2)=2.d0*RMU*EG(1,2)
C        TG(2,3)=2.d0*RMU*EG(2,3)
C        TG(3,1)=2.d0*RMU*EG(3,1)
C        TG(2,1)=TG(1,2)
C        TG(3,2)=TG(2,3)
C        TG(1,3)=TG(3,1)

! Transversly isotropic (fibre in 1 dir.n)
! Note: this is copied from CGS11 and needs modifying for 3D
      ELSE IF(IMT(nw).EQ.2) THEN
        CALL CGC11(nr,nw,nx,CC,CG,DENS,ERROR,*9999)
        TG(1,1)=CC(1,1)*EG(1,1)+CC(1,2)*EG(2,2)+CC(1,3)*EG(1,2)
        TG(2,2)=CC(2,1)*EG(1,1)+CC(2,2)*EG(2,2)+CC(2,3)*EG(1,2)
        TG(1,2)=CC(3,1)*EG(1,1)+CC(3,2)*EG(2,2)+CC(3,3)*EG(1,2)
        YMOD1=CG(1) !is Young's modulus E1
        YMOD2=CG(2) !is Young's modulus E2
C        SMOD1=CG(3) !is shear modulus G1
        POIS1=CG(4) !is Poisson's ratio nu1
        POIS2=CG(5) !is Poisson's ratio nu2
        EG(3,3)=-POIS1/YMOD1*TG(1,1)-POIS2/YMOD2*TG(2,2)

! Transversly isotropic (fibre in 2 dir.n)
! Note: this is copied from CGS11 and needs modifying for 3D
      ELSE IF(IMT(nw).EQ.3) THEN
        CALL CGC11(nr,nw,nx,CC,CG,DENS,ERROR,*9999)
        TG(1,1)=CC(1,1)*EG(1,1)+CC(1,2)*EG(2,2)+CC(1,3)*EG(1,2)
        TG(2,2)=CC(2,1)*EG(1,1)+CC(2,2)*EG(2,2)+CC(2,3)*EG(1,2)
        TG(1,2)=CC(3,1)*EG(1,1)+CC(3,2)*EG(2,2)+CC(3,3)*EG(1,2)
        YMOD1=CG(1) !is Young's modulus E1
        YMOD2=CG(2) !is Young's modulus E2
C        SMOD1=CG(3) !is shear modulus G1
        POIS1=CG(4) !is Poisson's ratio nu1
        POIS2=CG(5) !is Poisson's ratio nu2
        EG(3,3)=-POIS1/YMOD1*TG(1,1)-POIS2/YMOD2*TG(2,2)

! Orthotropic (fibre in 1 dir.n)
      ELSE IF(IMT(nw).EQ.4) THEN

      ENDIF
      IF(KTYP43.GT.0) THEN !include thermal stress effects
        DO i=1,3
          TG(i,i)=TG(i,i)+S_thermal(i)
        ENDDO !i
      ENDIF !ktyp43

      CALL EXITS('CGS9')
      RETURN
 9999 CALL ERRORS('CGS9',ERROR)
      CALL EXITS('CGS9')
      RETURN 1
      END


