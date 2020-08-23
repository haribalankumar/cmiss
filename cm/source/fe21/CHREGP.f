      SUBROUTINE CHREGP(REG_PARAMETER,STRING,ERROR,*)

C#### Subroutine: CHREGP
C###  Description:
C###    CHREGP changes specific epicardial inverse regularisation
C###    parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
!     Parameter List
      REAL*8 REG_PARAMETER(0:NTSM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N3CO,nts
      REAL*8 NEW_PARAMETER
      LOGICAL CBBREV

!     Functions
      INTEGER IFROMC
      REAL*8 RFROMC

      CALL ENTERS('CHREGP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change reg-parameters
C###  Parameter <index #[1]>
C###    Specify the time/equation index of the regularisation parameter.
C###  Parameter <all>
C###    Change all the regularisation parameters to be the same value
C###  Parameter <value #[0.0]>
C###    Specify the new regularisation parameter.
C###  Description:
C###    Changes a specific epicardial inverse regularisation parameter.
C###    The current regularisation parameters for all time/equations can
C###    be seen from the call to DRREGP. Once all the necessary
C###    regularisation parameters have been modified the subsequent
C###    call should be made to EVINVE to re-evaluate the
C###     epicardial inverse.

        OP_STRING(1)=STRING(1:IEND)//' <index #[1]>'
        OP_STRING(2)=BLANK(1:15)//'<all>'
        OP_STRING(3)=BLANK(1:15)//'<value #[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHREGP',ERROR,*9999)
      ELSE
        CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '    ERROR,*9999)


        IF(CBBREV(CO,'VALUE',1,noco+1,NTCO,N3CO)) THEN
          NEW_PARAMETER=RFROMC(CO(N3CO+1))
        ELSE
          NEW_PARAMETER=0.0d0
        ENDIF


        IF(CBBREV(CO,'ALL',1,noco+1,NTCO,N3CO)) THEN
          DO nts=1,INT(REG_PARAMETER(0))
            REG_PARAMETER(nts)=NEW_PARAMETER
          ENDDO
        ELSE
          IF(CBBREV(CO,'INDEX',1,noco+1,NTCO,N3CO)) THEN
            nts=IFROMC(CO(N3CO+1))
            CALL ASSERT((nts.GE.1).AND.(nts.LE.
     '        INT(REG_PARAMETER(0))),'>>Non valid index',ERROR,*9999)
          ELSE
            nts=1
          ENDIF

          REG_PARAMETER(nts)=NEW_PARAMETER

        ENDIF

      ENDIF

      CALL EXITS('CHREGP')
      RETURN
 9999 CALL ERRORS('CHREGP',ERROR)
      CALL EXITS('CHREGP')
      RETURN 1
      END


