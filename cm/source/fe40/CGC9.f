      SUBROUTINE CGC9(nr,nw,nx,CC,C1,C2,C3,CG,DENS,S_thermal,ERROR,*)

C#### Subroutine: CGC9
C###  Description:
C###    <HTML> <PRE>
C###    CGC9 evaluates coefficents for 3D elasticity.
C###    CC(i,j) is the 6*6 matrix of elastic material coefficients:
C###       1 isotropic,
C###       2 transversely isotropic wrt Xi(1),
C###       3 transversely isotropic wrt Xi(2),
C###       4 orthotropic.
C###       5 anisotropic (special case).
C###    Note: 2d_stress_vector = CC * 2d_strain_vector.
C###    If KTYP43.gt.1 thermal stress is calculated in S_thermal.
C###    </PRE> </HTML>

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'ktyp40.cmn'
!     Parameter List
      INTEGER nr,nw,nx
      REAL*8 CC(6,*),C1,C2,C3,CG(NMM),DENS,S_thermal(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j
      REAL*8 ALFA,POIS,POIS1,POIS2,SMOD1,TEMPERATURE,
     '  YMOD,YMOD1,YMOD2,YY

      CALL ENTERS('CGC9',*9999)

! Zero material matrix
      DO i=1,6
        DO j=1,6
          CC(i,j)=0.d0
        ENDDO
      ENDDO

! Isotropic
      IF(IMT(nw).EQ.1) THEN
        YMOD=CG(1) !is Young's modulus (kPa)
        POIS=CG(2) !is Poisson's ratio
        DENS=CG(ILT(nw,nr,nx))*1.D-3 !is density in Mg/m^3
        YY=YMOD/((1.d0+POIS)*(1.d0-2.d0*POIS))
        C1=YY*(1.d0-POIS)
        C2=YY*POIS
C GMH 3/1/96 mu=YY*(0.5d0-POIS) c3=2*mu
        C3=YMOD/(1.0D0+POIS)

        IF(KTYP43.GT.0) THEN !include thermal stress effects
          ALFA=CG(3)
          TEMPERATURE=CG(4)
          DO i=1,3
            S_thermal(i)=-ALFA*TEMPERATURE*YMOD/(3.d0*(1.d0-2.d0*POIS))
          ENDDO
        ENDIF !ktyp43


        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,
     '      '('' YMOD='',E11.3,'' POIS='',E11.3,'' YY='',E11.3)')
     '      YMOD,POIS,YY
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' C1='',E11.3,'' C2='',E11.3,'' C3='',E11.3)')
     '      C1,C2,C3
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

! Transversly isotropic (fibre in 1 dir.n)
! Note: this is copied from CGC11 and needs modifying for 3D
      ELSE IF(IMT(nw).EQ.2) THEN
        YMOD1=CG(1) !is Young's modulus E1 (GPa)
        YMOD2=CG(2) !is Young's modulus E2 (GPa)
        SMOD1=CG(3) !is shear modulus G1   (GPa)
        POIS1=CG(4) !is Poisson's ratio nu1
        POIS2=CG(5) !is Poisson's ratio nu2
        DENS=CG(ILT(nw,nr,nx))*1.0D-3 !is density (Mg/m^3)
        IF(nw.EQ.9) THEN
          YY=1.d0/(1.d0-POIS1*POIS2)
          CC(1,1)=YMOD1*YY
          CC(1,2)=YMOD1*POIS2*YY
          CC(1,3)=0.d0
c         CC(2,1)=YMOD2*POIS1*YY
          CC(2,1)=CC(1,2)
          CC(2,2)=YMOD2*YY
          CC(2,3)=0.d0
          CC(3,1)=0.d0
          CC(3,2)=0.d0
          CC(3,3)=SMOD1
        ENDIF

! Transversly isotropic (fibre in 2 dir.n)
! Note: this is copied from CGC11 and needs modifying for 3D
      ELSE IF(IMT(nw).EQ.3) THEN
        YMOD1=CG(1) !is Young's modulus E1 (GPa)
        YMOD2=CG(2) !is Young's modulus E2 (GPa)
        SMOD1=CG(3) !is shear modulus G1   (GPa)
        POIS1=CG(4) !is Poisson's ratio nu1
        POIS2=CG(5) !is Poisson's ratio nu2
        DENS=CG(ILT(nw,nr,nx))*1.0D-3 !is density (Mg/m^3)
        IF(nw.EQ.9) THEN
          YY=1.d0/(1.d0-POIS1*POIS2)
          CC(1,1)=YMOD1*YY
          CC(1,2)=YMOD1*POIS2*YY
          CC(1,3)=0.d0
          CC(2,1)=YMOD2*POIS1*YY
          CC(2,2)=YMOD2*YY
          CC(2,3)=0.d0
          CC(3,1)=0.d0
          CC(3,2)=0.d0
          CC(3,3)=SMOD1
        ENDIF

! Orthotropic (fibre in 1 dir.n)
      ELSE IF(IMT(nw).EQ.4) THEN

! Anisotropic (special case)
      ELSE IF(IMT(nw).EQ.5) THEN
        C1=CG(1)
        C2=CG(2)
        C3=CG(3)
        IF(KTYP43.GT.0) THEN !include thermal stress effects
          !calculate YMOD and POIS based on c1 and c2
          !calculate as YMOD/(1-2*POIS) as this is the term needed
          ALFA=CG(4)
          TEMPERATURE=CG(5)
          DO i=1,3
            S_thermal(i)=-ALFA*TEMPERATURE*(C1+2.0D0*C2)/3.0D0
          ENDDO
        ENDIF !ktyp43
      ENDIF !anisotropy

      CALL EXITS('CGC9')
      RETURN
 9999 CALL ERRORS('CGC9',ERROR)
      CALL EXITS('CGC9')
      RETURN 1
      END


