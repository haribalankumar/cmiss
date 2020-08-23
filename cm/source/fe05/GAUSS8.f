      SUBROUTINE GAUSS8(NAN,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS8
C###  Description:
C###    GAUSS8 defines the Gaussian quadrature coords XIGG and weights
C###    WG and evaluates the basis function Gauss point array PG for
C###    auxiliary type basis function nb.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NAN(NIM,NAM),nb,NGAP(NIM)
      REAL*8  PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,na,ng,ng1,ng2,ng3,ni,nu
      REAL*8 D(7,7),PSI8,W(7,7),XI(3),XIGG(7,7,7,3)

      DATA D/7*0.0D0,
     '      -0.2886751345948130D0,         0.2886751345948130D0,5*0.0D0,
     '      -0.3872983346207410D0, 0.0D0,  0.3872983346207410D0,4*0.0D0,
     '      -0.4305681557970260D0,        -0.1699905217924280D0,
     '       0.1699905217924280D0,         0.4305681557970260D0,3*0.0D0,
     '      -0.4530899229693320D0,        -0.2692346550528410D0,  0.0D0,
     '       0.2692346550528410D0,         0.4530899229693320D0,2*0.0D0,
     '      -0.4662347571015760D0,        -0.3306046932331330D0,
     '      -0.1193095930415990D0,         0.1193095930415990D0,
     '       0.3306046932331330D0,         0.4662347571015760D0,  0.0D0,
     '      -0.4745539561713800D0,        -0.3707655927996970D0,
     '      -0.2029225756886990D0, 0.0D0,  0.2029225756886990D0,
     '       0.3707655927996970D0,         0.4745539561713800D0/
      DATA W/1.0D0,6*0.0D0,2*0.50D0,5*0.0D0,0.2777777777777780D0,
     '       0.4444444444444440D0,         0.2777777777777780D0,4*0.0D0,
     '       0.1739274225687270D0,         0.3260725774312730D0,
     '       0.3260725774312730D0,         0.1739274225687270D0,3*0.0D0,
     '       0.1184634425280940D0,         0.2393143352496830D0,
     '       0.2844444444444440D0,         0.2393143352496830D0,
     '       0.1184634425280940D0, 2*0.0D0,
     '       0.0856622461895850D0,         0.1803807865240700D0,
     '       0.2339569672863460D0,         0.2339569672863460D0,
     '       0.1803807865240700D0,         0.0856622461895850D0,  0.0D0,
     '       0.0647424830844350D0,         0.1398526957446390D0,
     '       0.1909150252525600D0,         0.2089795918367350D0,
     '       0.1909150252525600D0,         0.1398526957446390D0,
     '       0.064742483084435D0/

      CALL ENTERS('GAUSS8',*9999)
      ng1=NGAP(1)
      ng2=1
      ng3=1
      IF(NIT(nb).GT.1) ng2=NGAP(2)
      IF(NIT(nb).GT.2) ng3=NGAP(3)
      DO 20 K=1,ng3
      DO 20 J=1,ng2
      DO 20 I=1,ng1
        XIGG(I,J,K,1)=0.50D0+D(I,ng1)
        XIGG(I,J,K,2)=0.50D0+D(J,ng2)
        XIGG(I,J,K,3)=0.50D0+D(K,ng3)
        ng=I+(J-1+(K-1)*ng2)*ng1
        WG(ng)=W(I,ng1)*W(J,ng2)*W(K,ng3)
        DO ni=1,NIT(nb)
          XI(ni)=XIGG(I,J,K,ni)
          XIG(ni,ng)=XI(ni)
        ENDDO
        DO na=1,NAT(nb)
          DO nu=1,NUT(nb)
C MPN 7/4/93 - Separating out NST(nb) and NAT(nb)
c            PG(na+NST(nb)-NAT(nb),nu,ng)=PSI8(NAN(1,na),nb,nu,XI)
            PG(na+NST(nb),nu,ng)=PSI8(NAN(1,na),nb,nu,XI)
          ENDDO
        ENDDO
 20   CONTINUE

      NBASEF(nb,0)=1  !number of children in family
      NBASEF(nb,1)=nb !parent number
      NFBASE(1,nb)=nb !family number of global basis number nbbem
      NFBASE(2,nb)=1  !local child number in family

      CALL EXITS('GAUSS8')
      RETURN
 9999 CALL ERRORS('GAUSS8',ERROR)
      CALL EXITS('GAUSS8')
      RETURN 1
      END


