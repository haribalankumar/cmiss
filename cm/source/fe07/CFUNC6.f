      SUBROUTINE CFUNC6(MODE,m,n,NJAC,X,F,NSTATE,NPROB,
     '  IUSER,USER)

C#### Subroutine: CFUNC6
C###  Description:
C###    CFUNC6 sets up the nonlinear constraints, C, for MINOS
C###    optimisation routines for a sparse constraint jacobian.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IUSER(*),m,MODE,n,NJAC,NPROB,NSTATE
      REAL*8 F(M),X(N),USER(*)
!     Local Variables
      INTEGER IBEG,IEND,nocont,noopti
      CHARACTER ERROR*(MXCH),STRING*(MXCH)
      LOGICAL END

      CALL ENTERS('CFUNC6',*9999)
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(CFUNC6_1)
        WRITE(OP_STRING,'('' MODE='',I3,'' M='',I4,'' N='',I4,'
     '    //''' NJAC='',I4,'' NSTATE='',I3,'' NPROB='',I3)')
     '    MODE,m,n,NJAC,NSTATE,NPROB
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(CFUNC6_1)
      ENDIF

      IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN !Data Fitting
        DO noco=1,5
          NTCOQU(noco)=0
        ENDDO
        DO noopti=1,NTOPTI
          USER(OS_PAOPTI+noopti-1)=X(noopti)
        ENDDO

        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='CONSTRAINTS'
        CO(4)='WRT'
        CO(5)='DATA_FITTING'
        NTCO=5
        noco=1
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

        DO nocont=1,m
          F(nocont)=USER(OS_CM+nocont-1)
        ENDDO
        IF(MODE.EQ.2) THEN
C          DO noopti=1,n
C            DO nocont=1,m
C              G(nocont,noopti)=USER(OS_CJACM+nocont-1+
C     '          (noopti-1)*NCOM)
C            ENDDO
C          ENDDO
        ENDIF

      ELSE
        ERROR=' Invalid optimisation KTYP''s'
        GOTO 9999
      ENDIF

      CALL EXITS('CFUNC6')
      RETURN
 9999 CALL ERRORS('CFUNC6',ERROR)
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      ERROR(1:)=' '
      MODE=-9999
 9998 CALL EXITS('CFUNC6')
      RETURN
      END


