      SUBROUTINE GAUSS7(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS7
C###  Description:
C###    GAUSS7 - for Extended Lagrange basis.
C###    Defines the Gaussian quadrature coords XIG and weights WG
C###    and evaluates the basis function Gauss point array PG for
C###    Lagrange / Hermite  tensor product type basis function nb.
C###    Also calculates the position of the gauss points on the
C###    boundary of the element.
C**** 22/1/92:   Changed to have even spacing of gauss/grid pts. GBS

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),NGAP(NIM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,nb,ng,ng1,ng2,ng3,ni,nk,nn,nu
      REAL*8 D(9,9),PSI1,W(9,9),XI(3),XIGG(9,9,9,3)

C      DATA D/                          0.0,                          8*0.0,
C     ' -0.5,                                               0.5, 7*0.0,
C     ' -0.5,                      0.0,                     0.5, 6*0.0,
C     ' -0.5, -0.2886751345948130,      0.2886751345948130, 0.5, 5*0.0,
C     ' -0.5, -0.3872983346207410, 0.0, 0.3872983346207410, 0.5, 4*0.0,
C     ' -0.5, -0.4305681557970260,     -0.1699905217924280,
C     '        0.1699905217924280,      0.4305681557970260, 0.5, 3*0.0,
C     ' -0.5, -0.4530899229693320,     -0.2692346550528410, 0.0,
C     '             0.2692346550528410,      0.4530899229693320, 0.5, 2*0.0,
C     '      -0.5, -0.4662347571015760,     -0.3306046932331330,
C     '            -0.1193095930415990,      0.1193095930415990,
C     '             0.3306046932331330,      0.4662347571015760, 0.5,   0.0,
C     '      -0.5, -0.4745539561713800,     -0.3707655927996970,
C     '            -0.2029225756886990, 0.0, 0.2029225756886990,
C     '             0.3707655927996970,      0.4745539561713800, 0.5/
C
C      DATA W/ 1.0D0,                                                   8*0.0,
C     '        0.0,                                              0.0, 7*0.0,
C     '        0.0, 1.0D0,                                         0.0, 6*0.0,
C     '        0.0, 2*0.50,                                      0.0, 5*0.0,
C     '        0.0, 0.2777777777777780,
C     '             0.4444444444444440,      0.2777777777777780, 0.0, 4*0.0,
C     '        0.0, 0.1739274225687270,      0.3260725774312730,
C     '             0.3260725774312730,      0.1739274225687270, 0.0, 3*0.0,
C     '        0.0, 0.1184634425280940,      0.2393143352496830,
C     '             0.2844444444444440,      0.2393143352496830,
C     '                                      0.1184634425280940, 0.0, 2*0.0,
C     '        0.0, 0.0856622461895850,      0.1803807865240700,
C     '             0.2339569672863460,      0.2339569672863460,
C     '             0.1803807865240700,      0.0856622461895850, 0.0,   0.0,
C     '        0.0, 0.0647424830844350,      0.1398526957446390,
C     '             0.1909150252525600,      0.2089795918367350,
C     '             0.1909150252525600,      0.1398526957446390,
C     '                                      0.0647424830844350, 0.0/

      DATA D/                         0.0D0,
     '8*0.0D0,
     '   -0.5D0,                                               0.5D0,
     '7*0.0D0,
     '   -0.5D0,                      0.0D0,                   0.5D0,
     '6*0.0D0,
     '   -0.5D0,-0.1666666666666667D0,      0.1666666666666667D0,0.5D0,
     '5*0.0D0,
     '   -0.5D0,-0.2500000000000000D0,0.0D0,0.2500000000000000D0,0.5D0,
     '4*0.0D0,
     '   -0.5D0,-0.3000000000000000D0,     -0.1000000000000000D0,
     '           0.1000000000000000D0,      0.3000000000000000D0,0.5D0,
     '3*0.0D0,
     '   -0.5D0,-0.3333333333333333D0,     -0.1666666666666667D0,0.0D0,
     '           0.1666666666666667D0,      0.3333333333333333D0,0.5D0,
     '2*0.0D0,
     '   -0.5D0,-0.3571428571428571D0,     -0.2142857142857142D0,
     '          -0.0714285714285714D0,      0.0714285714285714D0,
     '           0.2142857142857142D0,      0.3571428571428571D0,0.5D0,
     '    0.0D0,
     '   -0.5D0,-0.3750000000000000D0,     -0.2500000000000000D0,
     '          -0.1250000000000000D0,0.0D0,0.1250000000000000D0,
     '           0.2500000000000000D0,      0.3750000000000000D0,0.5D0/

      DATA W/ 1.0D0,                          8*0.0D0,
     '        0.0D0,                            0.0D0, 7*0.0D0,
     '        0.0D0, 1.0D0,                     0.0D0, 6*0.0D0,
     '        0.0D0, 2*0.50D0,                  0.0D0, 5*0.0D0,
     '        0.0D0, 3*0.33333333333333D0,      0.0D0, 4*0.0D0,
     '        0.0D0, 4*0.25D0,                  0.0D0, 3*0.0D0,
     '        0.0D0, 5*0.2D0,                   0.0D0, 2*0.0D0,
     '        0.0D0, 6*0.16666666666667D0,      0.0D0,   0.0D0,
     '        0.0D0, 7*0.14285714285714D0,      0.0D0 /

      CALL ENTERS('GAUSS7',*9999)
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
        DO 10 nk=1,NKT(0,nb)
        DO 10 nn=1,NNT(nb)
        DO 10 nu=1,NUT(nb)
          PG(nk+(nn-1)*NKT(0,nb),nu,ng)=PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)
 10     CONTINUE
 20   CONTINUE

      NBASEF(nb,0)=1 !number of children in family
      NBASEF(nb,1)=nb !parent number
      NFBASE(1,nb)=nb !family number of global basis number nbbem
      NFBASE(2,nb)=1  !local child number in family

      CALL EXITS('GAUSS7')
      RETURN
 9999 CALL ERRORS('GAUSS7',ERROR)
      CALL EXITS('GAUSS7')
      RETURN 1
      END


