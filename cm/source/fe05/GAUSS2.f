      SUBROUTINE GAUSS2(IBT,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS2
C###  Description:
C###    GAUSS2 defines the Gaussian quadrature coords XIGG and weights
C###    WG and evaluates the basis function Gauss point array PG for
C###    Simplex/Serendipity/Sector  type  basis function nb.

C**** XL(1..3) are area coordinates.
C**** XI(1..3) are Xi coordinates.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),INP(NNM,NIM),nb,NGAP(NIM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I(3,3),li,li2,MU,ng,ng1,NGPP(7),ni,nn,nu
      REAL*8 PSI2,PSI2_XI,WLG(7,4),XI(4),XL(4),XLG(3,7,4),
     '  XLG3D(4)

      DATA NGPP/1,0,2,3,0,0,4/
      DATA XLG/0.3333333333D0,0.3333333333D0,0.3333333333D0,
     '         18*0.0D0,
     '         0.5D0         ,0.5D0         ,0.0D0         ,
     '         0.0D0         ,0.5D0         ,0.5D0         ,
     '         0.5D0         ,0.0D0         ,0.5D0         ,
     '         12*0.0D0,
     '         0.3333333333D0,0.3333333333D0,0.3333333333D0,
     '         0.6D0         ,0.2D0         ,0.2D0         ,
     '         0.2D0         ,0.6D0         ,0.2D0         ,
     '         0.2D0         ,0.2D0         ,0.6D0         ,
     '         9*0.0D0,
     '         0.3333333333D0,0.3333333333D0,0.3333333333D0,
     '         0.0597158717D0,0.4701420641D0,0.4701420641D0,
     '         0.4701420641D0,0.0597158717D0,0.4701420641D0,
     '         0.4701420641D0,0.4701420641D0,0.0597158717D0,
     '         0.7974269853D0,0.1012865073D0,0.1012865073D0,
     '         0.1012865073D0,0.7974269853D0,0.1012865073D0,
     '         0.1012865073D0,0.1012865073D0,0.7974269853D0/
      DATA XLG3D/0.25D0,0.25D0,0.25D0,0.25D0/ !Gauss points in terms of
                                              !Xi coordinates for 3D
                                              !only one gauss point so
      DATA WLG/1.0D0,6*0.0D0,                 !far
     '         0.3333333333D0,0.3333333333D0,0.3333333333D0,
     '         4*0.0D0,
     '         -0.5625D0     ,0.5208333333D0,0.5208333333D0,
     '         0.5208333333D0,3*0.0D0,
     '         0.225D0       ,0.1323941527D0,0.1323941527D0,
     '         0.1323941527D0,0.1259391805D0,0.1259391805D0,
     '         0.1259391805D0/
      DATA I/10,13,6,6,10,13,9,15,14/

      CALL ENTERS('GAUSS2',*9999)

      IF(IBT(1,1).EQ.3) THEN !Simplex
        IF(IBT(3,1).EQ.1) THEN !Area coordinates
          ng1=NGPP(NGAP(1))
          DO ng=1,NGAP(1)
            DO li=1,3
              XL(li)=XLG(li,ng,ng1)
              XIG(li,ng)=XL(li)
              WG(ng)=WLG(ng,ng1)*0.5D0
            ENDDO
            DO nn=1,NNT(nb)
              PG(nn,1,ng)=PSI2(IBT,INP,nb,1,nn,XL)
              DO ni=1,NIT(nb)
                mu=1+ni*(ni+1)/2
                nu=2+ni*(ni+3)/2
                PG(nn,MU  ,ng)=-PSI2(IBT,INP,nb,2   ,nn,XL)
     '                         +PSI2(IBT,INP,nb,nu  ,nn,XL)
                PG(nn,MU+1,ng)= PSI2(IBT,INP,nb,3,   nn,XL)
     '                         +PSI2(IBT,INP,nb,nu+1,nn,XL)
     '                   -2.0D0*PSI2(IBT,INP,nb,1+ni*(11-ni)/2,nn,XL)
                mu=1+ni*(6-ni)
                PG(nn,MU  ,ng)= PSI2(IBT,INP,nb,3      ,nn,XL)
     '                         -PSI2(IBT,INP,nb,I(ni,1),nn,XL)
     '                         -PSI2(IBT,INP,nb,I(ni,2),nn,XL)
     '                         +PSI2(IBT,INP,nb,I(ni,3),nn,XL)
              ENDDO
            ENDDO
          ENDDO
C CS 12/8/97 new Xi coordinate basis functions for simplex elements

        ELSE IF(NIT(nb).EQ.3) THEN !3D simplex
! At the moment only considering a linear basis function - will
! do the quadratic etc stuff later
          ng1=NGPP(NGAP(1))
          DO ng=1,NGAP(1)
            DO li = 1,3
              XI(li)  = XLG3D(li)
              XIG(li,ng) = XI(li)
            ENDDO
            WG(1) = 1.0D0
            IF(NBSC(2,nb).EQ.1) THEN ! Linear basis function only so far
              DO nn = 1,NNT(nb)
                IF(nn.EQ.1) THEN
                  PG(nn,1,ng) = 1.0d0 - XI(1) - XI(2) - XI(3)
                ELSE IF(nn.EQ.2) THEN
                  PG(nn,1,ng) = XI(1)
                ELSE IF(nn.EQ.3) THEN
                  PG(nn,1,ng) = XI(2)
                ELSE
                  PG(nn,1,ng) = XI(3)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ELSE !Xi coordinates 2D
          ng1=NGPP(NGAP(1))
          DO ng=1,NGAP(1)
            li2=0
            DO li=1,2
              XI(li)=XLG(li+li2,ng,ng1)
              XIG(li,ng)=XI(li)
              WG(ng)=WLG(ng,ng1)*0.5D0
              li2=li2+1
            ENDDO
            DO nn=1,NNT(nb)
              DO nu=1,NUT(nb)
                PG(nn,nu,ng)=PSI2_XI(nb,nu,nn,XI)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ELSE IF(IBT(1,1).EQ.4) THEN !Serendipity

      ELSE IF(IBT(1,1).EQ.5) THEN !Sector

      ENDIF

      NBASEF(nb,0)=1 !number of children in family
      NBASEF(nb,1)=nb !parent number
      NFBASE(1,nb)=nb  !family number of global basis number nbbem
      NFBASE(2,nb)=1  !local child number in family

      CALL EXITS('GAUSS2')
      RETURN
 9999 CALL ERRORS('GAUSS2',ERROR)
      CALL EXITS('GAUSS2')
      RETURN 1
      END


