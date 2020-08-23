      SUBROUTINE BASE_XI(IBT,IDO,INP,NBJ,G,XE,XI,ERROR,*)

C#### Subroutine: BASE_XI
C###  Description:
C###    BASE_XI computes normalized Xi coord base vectors
C###    g1=G(1,nj),g2=G(2,nj),g3=G(3,nj).

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM)
      REAL*8 G(3,3),XE(NSM,NJM),XI(3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ni,nj,NU1(3)
      REAL*8 COSH_LAMDA,COS_MU,COS_THETA,DX(3,3),PXI,RLAMDA,RMU,
     '  SINH_LAMDA,SIN_MU,SIN_THETA,SUM,THETA

      DATA NU1/2,4,7/

      CALL ENTERS('BASE_XI',*9999)
      DO nj=1,NJT !calc. derivs of Xj (curvilinear) wrt Xi
        nb=NBJ(nj)
        DO ni=1,NIT(nb)
          DX(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      NU1(ni),XI,XE(1,nj))
        ENDDO !ni
      ENDDO !nj

      IF(ITYP10(1).EQ.1) THEN      !rectangular cartesian coords

      ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar coords

      ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar coords

      ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal coords
        nb=NBJ(1)
        RLAMDA     = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '                   XE(1,1))
        sinh_LAMDA = DSINH(RLAMDA)
        cosh_LAMDA = DCOSH(RLAMDA)
        nb=NBJ(2)
        RMU        = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '                   XE(1,2))
        sin_MU     = DSIN(RMU)
        cos_MU     = DCOS(RMU)
        nb=NBJ(3)
        THETA      = PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '                   XE(1,3))
        sin_THETA  = DSIN(THETA)
        cos_THETA  = DCOS(THETA)
        DO ni=1,NIT(NBJ(1))
          G(ni,1)=FOCUS*(sinh_LAMDA*DX(1,ni)*cos_MU-cosh_LAMDA
     '                  *sin_MU*DX(2,ni))
          G(ni,2)=FOCUS*(cosh_LAMDA*DX(1,ni)*sin_MU*cos_THETA
     '                  +sinh_LAMDA*cos_MU*DX(2,ni)*cos_THETA
     '                  -sinh_LAMDA*sin_MU*sin_THETA*DX(3,ni))
          G(ni,3)=FOCUS*(cosh_LAMDA*DX(1,ni)*sin_MU*sin_THETA
     '                  +sinh_LAMDA*cos_MU*DX(2,ni)*sin_THETA
     '                  +sinh_LAMDA*sin_MU*cos_THETA*DX(3,ni))
        ENDDO

      ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal coords
      ENDIF

      DO ni=1,NIT(NBJ(1)) !normalize base vectors
        SUM=DSQRT(G(ni,1)*G(ni,1)+G(ni,2)*G(ni,2)+G(ni,3)*G(ni,3))
        DO nj=1,NJT
          G(ni,nj)=G(ni,nj)/SUM
        ENDDO
      ENDDO

      CALL EXITS('BASE_XI')
      RETURN
 9999 CALL ERRORS('BASE_XI',ERROR)
      CALL EXITS('BASE_XI')
      RETURN 1
      END


