      SUBROUTINE SHPROF(ISEG,ISPROF,STRING,ERROR,*)

C#### Subroutine: SHPROF
C###  Description:
C###    SHPROF shows profile segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPROF(2)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,noprof

      CALL ENTERS('SHPROF',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show profile
C###  Description:
C###    Make the profile visible.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHPROF',ERROR,*9999)
      ELSE
        DO noprof=1,3
          iw=30+noprof
          CALL ACWK(iw,1,ERROR,*9999)
          IF(ISEG(ISPROF(noprof)).EQ.1) THEN
            CALL VISIB(iw,ISEG,ISPROF(noprof),'VISIBLE',ERROR,*9999)
          ELSE IF(ISEG(ISPROF(noprof)).EQ.0) THEN
            WRITE(OP_STRING,
     '        '('' >>Profile is not defined at '',I2)') noprof
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('SHPROF')
      RETURN
 9999 CALL ERRORS('SHPROF',ERROR)
      CALL EXITS('SHPROF')
      RETURN 1
      END


