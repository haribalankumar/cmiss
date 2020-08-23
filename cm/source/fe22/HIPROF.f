      SUBROUTINE HIPROF(ISEG,ISPROF,STRING,ERROR,*)

C#### Subroutine: HIPROF
C###  Description:
C###    HIPROF hides profile segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPROF(2)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,noprof

      CALL ENTERS('HIPROF',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide profile
C###  Description:
C###    Hide profile.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIPROF',ERROR,*9999)
      ELSE
        DO noprof=1,3
          iw=30+noprof
          CALL ACWK(iw,1,ERROR,*9999)
          IF(ISEG(ISPROF(noprof)).EQ.2) THEN
            CALL VISIB(iw,ISEG,ISPROF(noprof),'INVISIBLE',ERROR,*9999)
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('HIPROF')
      RETURN
 9999 CALL ERRORS('HIPROF',ERROR)
      CALL EXITS('HIPROF')
      RETURN 1
      END


