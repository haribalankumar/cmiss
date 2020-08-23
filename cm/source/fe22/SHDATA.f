      SUBROUTINE SHDATA(ISDANO,ISDAPR,ISDATA,ISDATR,ISEG,NEELEM,
     '  STRING,ERROR,*)

C#### Subroutine: SHDATA
C###  Description:
C###    SHDATA shows data segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISDANO(NWM,NEM),ISDAPR(NWM,NEM),ISDATA(NWM,NGRSEGM),
     '  ISDATR(NWM,NEM),ISEG(*),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,ISD,iw,IWK(6),N3CO,ne,nodata,noelem,
     '  noiw,nr,NTIW
      CHARACTER TYPE*11
      LOGICAL CBBREV

      CALL ENTERS('SHDATA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show data
C###  Description:
C###    Make the specified data segments, or their numbers, projections
C###    or trace visible on the specified workstation.
C###  Parameter:      <(positions/numbers/projections/trace)[position]>
C###    Specify what type of data information to make visible.
C###  Parameter:      <at (DATA_SET#/last)[last]>
C###    Specify which data sets are to be made visible. The options are
C###    all data sets or just the last data set.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify which windows to show the data sets on. The
C###    default is to show data in all windows.

        OP_STRING(1)=STRING(1:IEND)
     '    //' <(positions/numbers/projections/trace)[position]>'
        OP_STRING(2)=BLANK(1:15) //'<at (DATA_SET#/last)[last]>'
        OP_STRING(3)=BLANK(1:15)  //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHDATA',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'POSITIONS',2,noco+1,NTCO,N3CO)) THEN
          TYPE(1:8)='POSITION'
        ELSE IF(CBBREV(CO,'NUMBERS',2,noco+1,NTCO,N3CO)) THEN
          TYPE(1:7)='NUMBERS'
        ELSE IF(CBBREV(CO,'PROJECTIONS',2,noco+1,NTCO,N3CO)) THEN
          TYPE(1:11)='PROJECTIONS'
        ELSE IF(CBBREV(CO,'TRACE',2,noco+1,NTCO,N3CO)) THEN
          TYPE(1:5)='TRACE'
        ELSE
          TYPE(1:8)='POSITION'
        ENDIF
        IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
          nodata=IFROMC(CO(N3CO+1))
        ELSE
          nodata=NTDATA
        ENDIF

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(TYPE(1:8).EQ.'POSITION') THEN
              ISD=ISDATA(iw,nodata)
              IF(ISD.GT.0.AND.ISEG(ISD).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISD,'VISIBLE',ERROR,*9999)
              ELSE IF(ISD.EQ.0.OR.ISEG(ISD).EQ.0) THEN
                WRITE(OP_STRING,
     '            '('' >>Data at '',I3,'' is not defined on '',I1)')
     '            nodata,iw
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE IF(TYPE(1:7).EQ.'NUMBERS') THEN
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  ISD=ISDANO(iw,ne)
                  IF(ISD.GT.0.AND.ISEG(ISD).EQ.1) THEN
                    CALL VISIB(iw,ISEG,ISD,'VISIBLE',ERROR,*9999)
                  ELSE IF(ISD.EQ.0.OR.ISEG(ISD).EQ.0) THEN
                    WRITE(OP_STRING,'('' >>Data number in '
     '                //'element '',I4,'' is not defined on '',I1)')
     '                ne,iw
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
            ELSE IF(TYPE(1:11).EQ.'PROJECTIONS') THEN
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  ISD=ISDAPR(iw,ne)
                  IF(ISD.GT.0.AND.ISEG(ISD).EQ.1) THEN
                    CALL VISIB(iw,ISEG,ISD,'VISIBLE',ERROR,*9999)
                  ELSE IF(ISD.EQ.0.OR.ISEG(ISD).EQ.0) THEN
                    WRITE(OP_STRING,
     '                '('' >>Data projection in element '',I4,'
     '                //''' is not defined on '',I1)') ne,iw
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
            ELSE IF(TYPE(1:5).EQ.'TRACE') THEN
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  ISD=ISDATR(iw,ne)
                  IF(ISD.GT.0.AND.ISEG(ISD).EQ.1) THEN
                    CALL VISIB(iw,ISEG,ISD,'VISIBLE',ERROR,*9999)
                  ELSE IF(ISD.EQ.0.OR.ISEG(ISD).EQ.0) THEN
                    WRITE(OP_STRING,'('' >>Data trace in '
     '                //'element '',I4,'' is not defined on '',I1)')
     '                ne,iw
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHDATA')
      RETURN
 9999 CALL ERRORS('SHDATA',ERROR)
      CALL EXITS('SHDATA')
      RETURN 1
      END


