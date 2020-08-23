      SUBROUTINE CASURF(ISEG,ISSURF,NEELEM,STRING,ERROR,*)

C#### Subroutine: CASURF
C###  Description:
C###    CASURF cancels surface segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSURF(NWM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noiw,noelem,nr,NTIW
      LOGICAL ABBREV,SEGME

      CALL ENTERS('CASURF',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel surface;s
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###    Cancel rendered surface structure on specified workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CASURF',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ELSE
          SEGME=.FALSE.
        ENDIF

        IF(SEGME) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(ISSURF(iw,ne).GT.0) THEN
                    CALL DELETE_SEGMENT(ISSURF(iw,ne),ISEG,iw,ERROR,
     '                *9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('CASURF')
      RETURN
 9999 CALL ERRORS('CASURF',ERROR)
      CALL EXITS('CASURF')
      RETURN 1
      END


