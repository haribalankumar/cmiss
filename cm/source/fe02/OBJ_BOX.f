      SUBROUTINE OBJ_BOX(ISEGM,ITOT,ZD,WD,ERROR,*)

C#### Subroutine: OBJ_BOX
C###  Description:
C###    OBJ_BOX defines object from GKS box created by label command.

      IMPLICIT NONE
      INCLUDE 'draw00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEGM,ITOT
      REAL*8 WD(NJM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I
      REAL*8 XPTS(4),XSIDE,XWC,YPTS(4),YSIDE,YWC

      CALL ENTERS('OBJ_BOX',*9999)

      CALL ASSERT(NDT+40.LE.NDM,'>>NDM too small',ERROR,*9999)
      NTOBJE=NTOBJE+1
      OBJECT_NAME(NTOBJE)='BOX'
      ITOT=40
      NSOBJE(1,NTOBJE)=ISEGM  !is segment number of object
      NSOBJE(2,NTOBJE)=4      !is number of parts to object
      NSOBJE(3,NTOBJE)=NDT+1  !is 1st nd in object
      NSOBJE(4,NTOBJE)=NDT+10 !is 2nd nd in object
      NSOBJE(5,NTOBJE)=NDT+20 !is 3rd nd in object
      NSOBJE(6,NTOBJE)=NDT+30 !is 4th nd in object
      NSOBJE(7,NTOBJE)=NDT+40 !is 5th nd in object
      XWC  =DBLE(XSEGMENT_DATA(1,ISEGM))
      YWC  =DBLE(YSEGMENT_DATA(1,ISEGM))
      XSIDE=DBLE(XSEGMENT_DATA(2,ISEGM))
      YSIDE=DBLE(YSEGMENT_DATA(2,ISEGM))
      XPTS(1)=XWC-XSIDE
      YPTS(1)=YWC-YSIDE
      XPTS(2)=XWC+XSIDE
      YPTS(2)=YPTS(1)
      XPTS(3)=XPTS(2)
      YPTS(3)=YWC+YSIDE
      XPTS(4)=XPTS(1)
      YPTS(4)=YPTS(3)
      DO I=1,10
        ZD(1,NDT+I)=XPTS(1)+DBLE(I-1)/9.0D0*(XPTS(2)-XPTS(1))
        ZD(2,NDT+I)=YPTS(1)
      ENDDO
      DO I=11,20
        ZD(1,NDT+I)=XPTS(2)
        ZD(2,NDT+I)=YPTS(2)+DBLE(I-10)/10.0D0*(YPTS(3)-YPTS(2))
      ENDDO
      DO I=21,30
        ZD(1,NDT+I)=XPTS(3)-DBLE(I-20)/10.0D0*(XPTS(3)-XPTS(4))
        ZD(2,NDT+I)=YPTS(3)
      ENDDO
      DO I=31,40
        ZD(1,NDT+I)=XPTS(4)
        ZD(2,NDT+I)=YPTS(4)-DBLE(I-30)/10.0D0*(YPTS(4)-YPTS(1))
      ENDDO
      DO I=1,40
        WD(1,NDT+I)=1.0D0
        WD(2,NDT+I)=1.0D0
      ENDDO

      CALL EXITS('OBJ_BOX')
      RETURN
 9999 CALL ERRORS('OBJ_BOX',ERROR)
      CALL EXITS('OBJ_BOX')
      RETURN 1
      END


