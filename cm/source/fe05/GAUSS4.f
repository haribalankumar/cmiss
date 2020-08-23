      SUBROUTINE GAUSS4(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS4
C###  Description:
C###    GAUSS4 defines the Gaussian quadrature coords XIG and weights WG
C###    and evaluates the basis function Gauss point array PG for
C###    Lagrange/Hermite/Fourier tensor product type basis function nb.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),NGAP(NIM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k,l,nb,ng,ng1,ng2,ng3,ng4,ni,nk,nn,nu
      REAL*8 D(7,7),PSI1,W(7,7),XI(4),XIGG(7,7,7,7,4)

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

      CALL ENTERS('GAUSS4',*9999)
      ng1=NGAP(1)
      ng2=1
      ng3=1
      ng4=1
      IF(NIT(nb).GT.1) ng2=NGAP(2)
      IF(NIT(nb).GT.2) ng3=NGAP(3)
      IF(NIT(nb).GT.3) ng4=NGAP(4)
      DO l=1,ng4
        DO k=1,ng3
          DO j=1,ng2
            DO i=1,ng1
              XIGG(i,j,k,l,1)=0.50D0+D(i,ng1)
              XIGG(i,j,k,l,2)=0.50D0+D(j,ng2)
              XIGG(i,j,k,l,3)=0.50D0+D(k,ng3)
              XIGG(i,j,k,l,4)=0.50D0+D(l,ng4)
              ng=i+(j-1+(k-1+(l-1)*ng3)*ng2)*ng1
              WG(ng)=W(i,ng1)*W(j,ng2)*W(k,ng3)*W(l,ng4)
              DO ni=1,NIT(nb)
                XI(ni)=XIGG(i,j,k,l,ni)
                XIG(ni,ng)=XI(ni)
              ENDDO
              DO nk=1,NKT(0,nb)
                DO nn=1,NNT(nb)
                  DO nu=1,NUT(nb)
                    PG(nk+(nn-1)*NKT(0,nb),nu,ng)=
     '                PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      NBASEF(nb,0)=1 !number of children in family
      NBASEF(nb,1)=nb !parent number
      NFBASE(1,nb)=nb !family number of global basis number nbbem
      NFBASE(2,nb)=1  !local child number in family

      CALL EXITS('GAUSS4')
      RETURN
 9999 CALL ERRORS('GAUSS4',ERROR)
      CALL EXITS('GAUSS4')
      RETURN 1
      END


