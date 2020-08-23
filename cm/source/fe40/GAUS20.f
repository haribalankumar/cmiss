      SUBROUTINE GAUS20(IDO,INP,nb,NGP,D3PG,ERROR,*)

C#### Subroutine: GAUS20
C###  Description:
C###    GAUS20 defines the Gaussian quadrature coords XIG_LOCAL and weights WG
C###    and evaluates the basis function Gauss point array D3PG for
C###    Lagrange / Hermite  tensor product type basis function nb.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,NGP(*)
      REAL*8 D3PG(16,4,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,k,nd,ng,NG1,NG2,NG3,ni,nk,nn
      REAL*8 D(5,5),PSI20,XI(3),XIG_LOCAL(5,5,5,3)

      DATA D/5*0.d0,-0.2886751345948130d0,0.2886751345948130d0,3*0.d0,
     '       -0.3872983346207410d0,0.d0,0.3872983346207410d0,2*0.d0,
     '       -0.4305681557970260d0,    -0.1699905217924280d0,
     '        0.1699905217924280d0,     0.4305681557970260d0,0.d0,
     '       -0.4530899229693320d0,-0.2692346550528410d0,0.d0,
     '        0.2692346550528410d0,0.453089922969332d0/


      CALL ENTERS('GAUS20',*9999)
      NG1=NGP(1)
      NG2=NGP(2)
C      NG3=NGP(3)  !by MPN 29/1/92 - doesn't make sense assigning NG3 twice!
      NG3=1
      DO 20 k=1,NG3
      DO 20 j=1,NG2
      DO 20 i=1,NG1
        XIG_LOCAL(i,j,k,1)=0.50d0+D(i,NG1)
        XIG_LOCAL(i,j,k,2)=0.50d0+D(j,NG2)
        XIG_LOCAL(i,j,k,3)=0.50d0+D(k,NG3)
        ng=i+(j-1+(k-1)*NG2)*NG1
        DO ni=1,NIT(nb)
          XI(ni)=XIG_LOCAL(i,j,k,ni)
        ENDDO
        DO 10 nk=1,NKT(0,nb)
        DO 10 nn=1,NNT(nb)
        DO 10 nd=1,4
          D3PG(nk+(nn-1)*NKT(0,nb),nd,ng)
     '      =PSI20(IDO,INP,nb,nd,nk,nn,XI)
 10     CONTINUE
 20   CONTINUE

      CALL EXITS('GAUS20')
      RETURN
 9999 CALL ERRORS('GAUS20',ERROR)
      CALL EXITS('GAUS20')
      RETURN 1
      END


