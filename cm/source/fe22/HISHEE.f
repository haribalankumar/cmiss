      SUBROUTINE HISHEE(ISEG,ISSHEE,NEELEM,STRING,ERROR,*)

C#### Subroutine: HISHEE
C###  Description:
C###    HISHEE hides sheet segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSHEE(NWM,NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),ne,noelem,noiw,noshee,nr,NTIW

      CALL ENTERS('HISHEE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide sheets
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Description:
C###     Hide sheet segments on specified workstations

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HISHEE',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO noshee=1,NTSHEE
              DO nr=1,NRT
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(ISEG(ISSHEE(iw,ne,noshee)).EQ.2) THEN
                    CALL VISIB(iw,ISEG,ISSHEE(iw,ne,noshee),'INVISIBLE',
     '                ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HISHEE')
      RETURN
 9999 CALL ERRORS('HISHEE',ERROR)
      CALL EXITS('HISHEE')
      RETURN 1
      END


