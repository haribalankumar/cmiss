      SUBROUTINE FUNCT5(MODE,n,X,F,G,NSTATE,NPROB,Z,NWCORE,
     '  IUSER,USER)

C#### Subroutine: FUNCT5
C###  Description:
C###    FUNCT5 evaluates the objective function, F, for MINOS
C###    optimisation routines.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IUSER(*),MODE,n,NPROB,NSTATE,NWCORE
      REAL*8 F,X(N),G(N),Z(NWCORE),USER(*)
!     Local Variables
      INTEGER IBEG,IEND,noopti
      CHARACTER ERROR*(MXCH),STRING*(MXCH)
      LOGICAL END

      CALL ENTERS('FUNCT5',*9999)

      IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN !Data Fitting
C **    Objective function is sum of squares of all the orthogonal
C **    projection distances of the data points
        DO noco=1,5
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,N
          USER(OS_PAOPTI+noopti-1)=X(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='OBJECTIVE'
        CO(4)='WITH'
        CO(5)='OPTIMISATION_PARAMETERS'
        NTCO=5
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        F=FUNC
        IF(MODE.EQ.2) THEN
          DO noopti=1,n
            G(noopti)=USER(OS_FGRADM+noopti-1)
          ENDDO
        ENDIF

      ELSE
        ERROR=' Invalid optimisation KTYP''s'
        GOTO 9999
      ENDIF

      CALL EXITS('FUNCT5')
      RETURN
 9999 CALL ERRORS('FUNCT5',ERROR)
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      ERROR(1:)=' '
      MODE=-9999
 9998 CALL EXITS('FUNCT5')
      RETURN
      END



