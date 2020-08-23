      SUBROUTINE CALC_LATT_GDOTN(NITB,AQ,GDOTN,PROPQ,ERROR,*)

C#### Subroutine: CALC_LATT_GDOTN
C###  Description:
C###    CALC_LATT_GDOTN calculates the product of the conductivity
C###    tensor G with the normal tensor XNLOCAL and returns
C###    this value in the vector GDOTN.       

      IMPLICIT NONE
      
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'tol00.cmn'
      
!     Parameter List
      INTEGER NITB
      REAL*8 AQ(NMAQM),GDOTN(3),PROPQ(3,3)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER maq,ni,nj
      REAL*8 XNLOCAL(3)
      
      CALL ENTERS('CALC_LATT_GDOTN',*9999)

C     Calculate a normal vector at nq
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_X,ERROR,*9999)
      XNLOCAL(1)=AQ(maq)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Y,ERROR,*9999)
      XNLOCAL(2)=AQ(maq)
      CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Z,ERROR,*9999)
      XNLOCAL(3)=AQ(maq)

C     Find the coefficient
      DO ni=1,NITB
        GDOTN(ni)=0.0d0
        DO nj=1,NITB
          GDOTN(ni)=GDOTN(ni)+(PROPQ(nj,ni)*XNLOCAL(nj))
        ENDDO
      ENDDO
      
      CALL EXITS('CALC_LATT_GDOTN')
      RETURN
 9999 CALL ERRORS('CALC_LATT_GDOTN',ERROR)
      CALL EXITS('CALC_LATT_GDOTN')
      RETURN 1
      END        


