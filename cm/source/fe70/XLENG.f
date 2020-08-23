      REAL*8 FUNCTION XLENG(X1,X2)

C#### Function: XLENG
C###  Type: REAL*8
C###  Description:
C###    XLENG returns integrated length of cubic in x,y or x,z plane.
C###    X1(nk,j),X2(nk,j);nk=1,2;j=1,2 are coordinates of two nodes.

      IMPLICIT NONE
!     Parameter List
      REAL*8 X1(2,2),X2(2,2)
!     Local Variables
      INTEGER J,ng
      REAL*8 DXDXI(2),PH3,WG_LOCAL(4),XI,XIG_LOCAL(4)
      DATA XIG_LOCAL/0.0694318442029D0,0.3300094782075D0,
     '         0.6699905217924D0,0.9305681557970D0/
      DATA  WG_LOCAL/0.1739274225687D0,0.3260725774313D0,
     '         0.3260725774313D0,0.1739274225687D0/

      XLENG=0.0d0
      DO ng=1,4
        XI=XIG_LOCAL(ng)
        DO j=1,2
          DXDXI(j)=PH3(1,1,2,XI)*X1(1,J)+PH3(1,2,2,XI)*X1(2,J)
     '            +PH3(2,1,2,XI)*X2(1,J)+PH3(2,2,2,XI)*X2(2,J)
        ENDDO
        XLENG=XLENG+WG_LOCAL(ng)*SQRT(DXDXI(1)**2+DXDXI(2)**2)
      ENDDO

      RETURN
      END


