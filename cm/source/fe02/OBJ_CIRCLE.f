      SUBROUTINE OBJ_CIRCLE(ISEGM,ITOT,ZD,WD,ERROR,*)

C#### Subroutine: OBJ_CIRCLE
C###  Description:
C###    OBJ_CIRCLE defines object from GKS circle created by label
C###    command.

      IMPLICIT NONE
      INCLUDE 'draw00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEGM,ITOT
      REAL*8 WD(NJM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I
      REAL*8 CIRCUM,RADIUS,THETA,THETA1,THETA2,XCENTR,YCENTR

      CALL ENTERS('OBJ_CIRCLE',*9999)

      XCENTR=DBLE(XSEGMENT_DATA(1,ISEGM))
      YCENTR=DBLE(YSEGMENT_DATA(1,ISEGM))
      RADIUS=DBLE(XSEGMENT_DATA(2,ISEGM))
      THETA1=DBLE(XSEGMENT_DATA(3,ISEGM))
      THETA2=DBLE(YSEGMENT_DATA(3,ISEGM))
      CIRCUM=(THETA2-THETA1)*RADIUS
      ITOT=INT(100.0D0*CIRCUM/DBLE(DIAG))
      CALL ASSERT(NDT+ITOT.LE.NDM,'NDM too small',ERROR,*9999)
      NTOBJE=NTOBJE+1
      OBJECT_NAME(NTOBJE)='CIRCLE'
      NSOBJE(1,NTOBJE)=ISEGM    !is segment number of object
      NSOBJE(2,NTOBJE)=1        !is number of parts to object
      NSOBJE(3,NTOBJE)=NDT+1    !is 1st nd in object
      NSOBJE(4,NTOBJE)=NDT+ITOT !is 2nd nd in object
      DO I=1,ITOT
        THETA=THETA1+DBLE(I-1)/DBLE(ITOT-1)*(THETA2-THETA1)
        ZD(1,NDT+I)=XCENTR+RADIUS*DCOS(THETA)
        ZD(2,NDT+I)=YCENTR+RADIUS*DSIN(THETA)
        WD(1,NDT+I)=1.0D0
        WD(2,NDT+I)=1.0D0
      ENDDO

      CALL EXITS('OBJ_CIRCLE')
      RETURN
 9999 CALL ERRORS('OBJ_CIRCLE',ERROR)
      CALL EXITS('OBJ_CIRCLE')
      RETURN 1
      END


