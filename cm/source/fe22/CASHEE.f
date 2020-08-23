      SUBROUTINE CASHEE(ISEG,ISSHEE,NEELEM,STRING,ERROR,*)

C#### Subroutine: CASHEE
C###  Description:
C###    CASHEE cancels sheet segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSHEE(NWM,NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noelem,noiw,noshee,
     '  nr,NTIW
      LOGICAL ABBREV,SEGME

      CALL ENTERS('CASHEE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel sheets;s
C###  Parameter:      <(all/last)[all]>
C###    Specify either 'all' lines or just the 'last' line drawn.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Cancel sheet data on specified workstation.

        OP_STRING(1)=STRING(1:IEND) //';s'
        OP_STRING(2)=BLANK(1:15) //'<(all/last)[all]>'
        OP_STRING(3)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CASHEE',ERROR,*9999)
      ELSE
C MPN 12-Oct-94 Adding assert
        CALL ASSERT(NJ_LOC(NJL_FIBR,0,0).EQ.3,
     '    '>>Sheets not defined',ERROR,*9999)
        CALL CHECKQ(' S',noco,1,CO,COQU,STRING,*1)
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ELSE
          SEGME=.FALSE.
        ENDIF

        IF(SEGME) THEN
          IF(ABBREV(CO(noco+1),'LAST',1)) THEN
            DO noiw=1,NTIW
              iw=IWK(noiw)
              IF(IWKS(iw).GT.0) THEN
                CALL ACWK(iw,1,ERROR,*9999)
                DO nr=1,NRT
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    IF(ISSHEE(iw,ne,NTSHEE).GT.0) THEN
                      CALL DELETE_SEGMENT(ISSHEE(iw,ne,NTSHEE),ISEG,iw,
     '                  ERROR,*9999)
                    ENDIF
                  ENDDO
                ENDDO
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
            NTSHEE=NTSHEE-1
          ELSE
            DO noiw=1,NTIW
              iw=IWK(noiw)
              IF(IWKS(iw).GT.0) THEN
                CALL ACWK(iw,1,ERROR,*9999)
                DO noshee=1,NTSHEE
                  DO nr=1,NRT
                    DO noelem=1,NEELEM(0,nr)
                      ne=NEELEM(noelem,nr)
                      IF(ISSHEE(iw,ne,noshee).GT.0) THEN
                        CALL DELETE_SEGMENT(ISSHEE(iw,ne,noshee),ISEG,
     '                    iw,ERROR,*9999)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
                CALL DAWK(iw,1,ERROR,*9999)
              ENDIF
            ENDDO
            NTSHEE=0
          ENDIF

        ELSE IF(.NOT.SEGME) THEN
          CALL_DATA_SHEET=.FALSE.
        ENDIF
      ENDIF

      CALL EXITS('CASHEE')
      RETURN
 9999 CALL ERRORS('CASHEE',ERROR)
      CALL EXITS('CASHEE')
      RETURN 1
      END


