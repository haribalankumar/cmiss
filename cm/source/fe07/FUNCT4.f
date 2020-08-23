      SUBROUTINE FUNCT4(MODE,NTRES,NT_OPTI,LDFJ,XC,FC,FJAC,NSTATE,
     '  IUSER,USER)

C#### Subroutine: FUNCT4
C###  Description:
C###    FUNCT4 evaluates a set of residuals, FC, for optimisation
C###    routine MINLSSQP, during a SOLVE.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IUSER(*),LDFJ,MODE,NT_OPTI,NSTATE,NTRES
      REAL*8 FC(*),FJAC(LDFJ,*),USER(*),XC(*)
!     Local Variables
      INTEGER IBEG,IEND,noopti,nores
      CHARACTER ERROR*(MXCH),STRING*(MXCH)
      LOGICAL END

      CALL ENTERS('FUNCT4',*9999)
      DO noco=1,16
        NTCOQU(noco)=0
      ENDDO
      DO noopti=1,NT_OPTI
        USER(OS_PBOPTI+noopti-1)=XC(noopti)
      ENDDO

      CO(1)='FEM'
      CO(2)='UPDATE'
      CO(3)='OPTIMISATION'
      CO(4)='SOLVE' !new MPN 20-Apr-95: for NONLIN solns with MINLSSQP
      NTCO=4
      noco=1
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP   CRITICAL(FUNCT4_1)
        WRITE(OP_STRING,'('' fem update optimisation..'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP   END CRITICAL(FUNCT4_1)
      ENDIF
      CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)

      IF(MODE.EQ.0.OR.MODE.EQ.2) THEN !evaluate residuals
        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='WRT'
        CO(5)='GEOM_PARAMS'
        CO(6)='NOVIEW'
        NTCO=6
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT4_2)
          WRITE(OP_STRING,'('' fem eval resid wrt geom noview..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT4_2)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)
        DO nores=1,NTRES
          FC(nores)=USER(OS_RESID+nores-1)
        ENDDO
      ENDIF

      IF((MODE.EQ.1.OR.MODE.EQ.2)   !evaluate first derivatives
     '   .AND.KTYP1D.GT.0) THEN     !KTYP1D is the derivative level set
        CO(1)='FEM'
        CO(2)='EVALUATE'
        CO(3)='RESIDUALS'
        CO(4)='ANALYT_DERIVS'
        CO(5)='WRT'
        CO(6)='GEOM_PARAMS'
        CO(7)='NOVIEW'
        NTCO=7
        noco=1
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP     CRITICAL(FUNCT4_3)
          WRITE(OP_STRING,'('' fem eval resid analyt wrt geom '
     '      //'noview..'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP     END CRITICAL(FUNCT4_3)
        ENDIF
        CALL FEM(ISEG_TEMP,CSEG_TEMP,END,STRING,IUSER,USER,ERROR,*9999)
        DO nores=1,NTRES
          DO noopti=1,NT_OPTI
            FJAC(nores,noopti)=USER(OS_RESJAC+nores-1+
     '       (noopti-1)*NREM)
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('FUNCT4')
      RETURN
 9999 CALL ERRORS('FUNCT4',ERROR)
      CALL STRING_TRIM(ERROR,IBEG,IEND)
      WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)
      CALL WRITES(IOER,OP_STRING,ERROR,*9998)
      ERROR(1:)=' '
      MODE=-9999
 9998 CALL EXITS('FUNCT4')
      RETURN
      END


