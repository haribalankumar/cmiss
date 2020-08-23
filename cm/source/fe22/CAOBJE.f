      SUBROUTINE CAOBJE(ISEG,ISOBJE,STRING,ERROR,*)

C#### Subroutine: CAOBJE
C###  Description:
C###    CAOBJE cancels object segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'obje00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISOBJE(NWM,NGRSEGM,NGRSEGM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INSTAT,IPICK,ISEGM,iw,IWK(6),noiw,noobje,
     '  nopart,NTIW,NTPART
      LOGICAL ABBREV,MOUSE,SEGME

      CALL ENTERS('CAOBJE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel object;s
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to cancel the
C###    objects on.
C###  Description:
C###    Cancel object segment on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM cancel object;m
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation (GX window) to cancel the
C###    objects on.
C###  Description:
C###    Cancel object segment with mouse on specified workstation.

        OP_STRING(1)=BLANK(1:15) //';m'
        OP_STRING(2)=BLANK(1:15) //'<on WS#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAOBJE',ERROR,*9999)
      ELSE
        CALL CHECKQ('SM',noco,1,CO,COQU,STRING,*1)

        SEGME=.FALSE.
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
        ENDIF
        IF(SEGME.OR.MOUSE) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ENDIF

        IF(SEGME) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO noobje=1,NTOBJE
c                NOSEGM=NSOBJE(1,noobje) !is segment number of complete object
                NTPART=NSOBJE(2,noobje) !is number of parts to object
                DO nopart=1,NTPART
                  IF(ISOBJE(iw,noobje,nopart).GT.0) THEN
                    CALL DELETE_SEGMENT(ISOBJE(iw,noobje,nopart),ISEG,
     '                IW,ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO

        ELSE IF(MOUSE) THEN
          iw=IWK(1)
          CALL ACWK(iw,0,ERROR,*9999)
          DO noobje=1,NTOBJE
            NTPART=NSOBJE(2,noobje) !is number of parts to object
            DO nopart=1,NTPART
              IF(ISEG(ISOBJE(iw,noobje,nopart)).GT.0) THEN
                CALL DETECT(iw,ISEG,ISOBJE(iw,noobje,nopart),
     '            'DETECTABLE',ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
          WRITE(OP_STRING,'('' >>Pick object parts on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL PICK(iw,'EVENT',INSTAT,ISEGM,IPICK,ERROR,*9999)
C          CONTINUE=.TRUE.
C          DO WHILE(CONTINUE)
c           CALL EVENT(ID_WS,ID_DEVICE,INPUT_STATUS,CLASS,IDATA,
c    '        R4DATA,SDATA,ERROR,*9999)
c            IF(DOP) THEN
c              WRITE(OP_STRING,*) ' INPUT_CLASS=',CLASS
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            ENDIF
c            IF(CLASS(1:4).EQ.'PICK') THEN
c              IF(INPUT_STATUS.EQ.1) THEN
c                ISEGM=IDATA(1)
c                WRITE(OP_STRING,'('' Picked object segment '',I3)')
c     '           ISEGM
c                 CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                NAOBJE=IFROMC(CSEG(ISEGM)(8:9))   !is current object number
c                NAPART=IFROMC(CSEG(ISEGM)(56:57)) !is current part number
c                nd1=NSOBJE(2+NAPART,NAOBJE) !is beginning of part
c                ND2=NSOBJE(3+NAPART,NAOBJE) !is end of part
c                NDIFF=ND2-nd1
c                IF(NAPART.EQ.1) THEN
c                  ND3=nd1
c                  ND4=NDT-NDIFF
c                ELSE IF(NAPART.GT.1) THEN
c                  ND3=nd1+1
c                  ND4=NDT-NDIFF
c                ENDIF
c                DO nd=ND3,ND4 !to restack data arrays
c                  DO nj=1,NJT
c                    ZD(nj,nd)=ZD(nj,nd+NDIFF)
c                    WD(nj,nd)=WD(nj,nd+NDIFF)
c                  ENDDO
c                ENDDO
c                CALL DELETE_SEGMENT(ISEGM,ISEG,iw,ERROR,*9999)
c                IF(NAPART.LT.NSOBJE(2,NAOBJE)) THEN
c                  DO nopart=NAPART,NSOBJE(2,NAOBJE) !to restack object array
c                    NSOBJE(2+nopart,NAOBJE)=NSOBJE(3+nopart,NAOBJE)
c                    ISOBJE(iw,NAOBJE,nopart)=ISOBJE(iw,NAOBJE,1+nopart)
c                  ENDDO
c                ENDIF
c                NSOBJE(2,NAOBJE)=NSOBJE(2,NAOBJE)-1 !to reduce no of parts
c                DO noobje=1,NTOBJE !to reallocate nd's in object arrays
c                  NTPART=NSOBJE(2,noobje)
c                  DO nopart=1,NTPART
c                    nd=NSOBJE(2+nopart,noobje)
c                    IF(nd.GE.nd1) THEN
c                      NSOBJE(2+nopart,noobje)=
c     '                  NSOBJE(2+nopart,noobje)-NDIFF
c                      IF(nopart.EQ.NTPART) THEN
c                        NSOBJE(3+nopart,noobje)=
c     '                    NSOBJE(3+nopart,noobje)-NDIFF
c                      ENDIF
c                    ENDIF
c                  ENDDO
c                ENDDO
c                NDT=NDT-NDIFF
c
c              ELSE
c               CALL INPUT_MODE(iw,LD1,'PICK','REQUEST',ERROR,*9999)
c                CONTINUE=.FALSE.
c              ENDIF
c            ENDIF
C          ENDDO
          CALL DAWK(iw,0,ERROR,*9999)
          CALL ACWK(iw,1,ERROR,*9999)
          DO noobje=1,NTOBJE
            NTPART=NSOBJE(2,noobje) !is number of parts to object
            DO nopart=1,NTPART
              IF(ISEG(ISOBJE(iw,noobje,nopart)).GT.0) THEN
                CALL DETECT(iw,ISEG,ISOBJE(iw,noobje,nopart),
     '            'UNDETECTABLE',ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('CAOBJE')
      RETURN
 9999 CALL ERRORS('CAOBJE',ERROR)
      CALL EXITS('CAOBJE')
      RETURN 1
      END


