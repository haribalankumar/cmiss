      SUBROUTINE OBJ_BEZIER(ISEGM,ITOT,ZD,WD,ERROR,*)

C#### Subroutine: OBJ_BEZIER
C###  Description:
C###    OBJ_BEZIER defines object from GKS bezier curve created by
C###    label command.

      IMPLICIT NONE
      INCLUDE 'draw00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      REAL*8 ZD(NJM,NDM),WD(NJM,NDM),XPTS(4),YPTS(4)
      CHARACTER ERROR*(*)
!     Lo cal Variables
      INTEGER I,ISEGM,ITOT,N1PTS,N2PTS,NDPTS,NOPT1,nopts,NTPTS
      REAL*8 XI,XI2,XI3

      CALL ENTERS('OBJ_BEZIER',*9999)

      NTPTS=ISEGMENT_DATA(1,ISEGM)
      CALL ASSERT(NDT+10*NTPTS.LE.NDM,'NDM too small',ERROR,*9999)
      NTOBJE=NTOBJE+1
      OBJECT_NAME(NTOBJE)='BEZIER'
      NSOBJE(1,NTOBJE)=ISEGM !is segment number of object
      NSOBJE(2,NTOBJE)=NTPTS !is number of parts to object
      NSOBJE(3,NTOBJE)=NDT+1 !is 1st nd in object
      NDPTS=10
      XPTS(1)=DBLE(XSEGMENT_DATA(1,ISEGM))
      YPTS(1)=DBLE(YSEGMENT_DATA(1,ISEGM))
      ZD(1,NDT+1)=XPTS(1)
      ZD(2,NDT+1)=YPTS(1)
      N1PTS=1
      DO nopts=1,NTPTS
        NOPT1=3*(nopts-1)+1
        XPTS(2)=DBLE(XSEGMENT_DATA(NOPT1+1,ISEGM))
        YPTS(2)=DBLE(YSEGMENT_DATA(NOPT1+1,ISEGM))
        N2PTS=N1PTS+NDPTS
        DO I=N1PTS+1,N2PTS
          XI=DBLE(I-N1PTS)/DBLE(NDPTS)
          XI2=XI*XI
          XI3=XI2*XI
          ZD(1,NDT+I)=(1.0D0-3.0D0*XI+3.0D0*XI2-XI3)
     '               * DBLE(XSEGMENT_DATA(NOPT1  ,ISEGM)) +
     '                (3.0D0*XI-6.0D0*XI2+3.0D0*XI3)
     '               * DBLE(XSEGMENT_DATA(NOPT1+1,ISEGM)) +
     '                (3.0D0*XI2-3.0D0*XI3      )
     '               * DBLE(XSEGMENT_DATA(NOPT1+2,ISEGM)) +
     '                 XI3
     '               * DBLE(XSEGMENT_DATA(NOPT1+3,ISEGM))
          ZD(2,NDT+I)=(1.0D0-3.0D0*XI+3.0D0*XI2-XI3)
     '               * DBLE(YSEGMENT_DATA(NOPT1  ,ISEGM)) +
     '                (3.0D0*XI-6.0D0*XI2+3.0D0*XI3)
     '               * DBLE(YSEGMENT_DATA(NOPT1+1,ISEGM)) +
     '                (3.0D0*XI2-3.0D0*XI3      )
     '               * DBLE(YSEGMENT_DATA(NOPT1+2,ISEGM)) +
     '                 XI3
     '               * DBLE(YSEGMENT_DATA(NOPT1+3,ISEGM))
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

      CALL EXITS('OBJ_BEZIER')
      RETURN
 9999 CALL ERRORS('OBJ_BEZIER',ERROR)
      CALL EXITS('OBJ_BEZIER')
      RETURN 1
      END


