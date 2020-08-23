      SUBROUTINE OBJ_LINE(ISEGM,ITOT,ZD,WD,ERROR,*)

C#### Subroutine: OBJ_LINE
C###  Description:
C###    OBJ_LINE defines object from GKS line created by label command.

      IMPLICIT NONE
      INCLUDE 'draw00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEGM,ITOT
      REAL*8 WD(NJM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,N1PTS,N2PTS,NDPTS,nopts,NTPTS
      REAL*8 XPTS(4),YPTS(4)

      CALL ENTERS('OBJ_LINE',*9999)

      CALL ASSERT(NDT+40.LE.NDM,'NDM too small',ERROR,*9999)
      NTOBJE=NTOBJE+1
      OBJECT_NAME(NTOBJE)='LINE'
      NTPTS=ISEGMENT_DATA(1,ISEGM)
      NSOBJE(1,NTOBJE)=ISEGM !is segment number of object
      NSOBJE(2,NTOBJE)=NTPTS !is number of parts to object
      NSOBJE(3,NTOBJE)=NDT+1 !is 1st nd in object
      NDPTS=INT(40.0D0/DBLE(NTPTS))
      XPTS(1)=DBLE(XSEGMENT_DATA(1,ISEGM))
      YPTS(1)=DBLE(YSEGMENT_DATA(1,ISEGM))
      ZD(1,NDT+1)=XPTS(1)
      ZD(2,NDT+1)=YPTS(1)
      N1PTS=1
      DO nopts=1,NTPTS
        XPTS(2)=DBLE(XSEGMENT_DATA(nopts+1,ISEGM))
        YPTS(2)=DBLE(YSEGMENT_DATA(nopts+1,ISEGM))
        N2PTS=N1PTS+NDPTS
        DO I=N1PTS+1,N2PTS
          ZD(1,NDT+I)=XPTS(1)+DBLE(I-N1PTS)/DBLE(NDPTS)
     '                       *(XPTS(2)-XPTS(1))
          ZD(2,NDT+I)=YPTS(1)+DBLE(I-N1PTS)/DBLE(NDPTS)
     '                       *(YPTS(2)-YPTS(1))
        ENDDO
        NSOBJE(3+nopts,NTOBJE)=NDT+N2PTS !is next nd in object
        N1PTS=N2PTS
        XPTS(1)=XPTS(2)
        YPTS(1)=YPTS(2)
      ENDDO
      ITOT=N2PTS
      DO I=1,ITOT
        WD(1,NDT+I)=1.0D0
        WD(2,NDT+I)=1.0D0
      ENDDO

      CALL EXITS('OBJ_LINE')
      RETURN
 9999 CALL ERRORS('OBJ_LINE',ERROR)
      CALL EXITS('OBJ_LINE')
      RETURN 1
      END


