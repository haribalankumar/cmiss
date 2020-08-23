      SUBROUTINE NAGMIN(INTWORK,nr,SQUID_CONFIG,BASELINE,
     '  REALWORK,TIME,WARMSTART,ERROR,*)

C#### Subroutine: NAGMIN
C###  Description:
C###    NAGMIN is a buffer routine to NAGMINA the minimisation
C###    routine using MINLSSQP.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER INTWORK(NIWM),nr,SQUID_CONFIG
      REAL*8 BASELINE,REALWORK(NRWM),TIME
      CHARACTER ERROR*(*)
      LOGICAL WARMSTART
!     Local Variables
      INTEGER i,j

      CALL ENTERS('NAGMIN',*9999)

C CPB 2/9/93 ZEROING JACOBIAN ARRAYS etc
      IF(.NOT.WARMSTART) THEN
C KAT 15Feb01: corrected initialization. LDFJ = NREM,
C       although perhaps could be NT_RES.
        DO j=0,NTOPTI-1
          DO i=0,NT_RES-1
            REALWORK(OS_RESJAC+i+NREM*j)=0.0d0
          ENDDO
        ENDDO
C KAT 15Feb01: corrected initialization. LDCJ = NCOM,
C       although perhaps could be NTCNTR.
        DO j=0,NTOPTI-1
          DO i=0,NTCNTR-1
            REALWORK(OS_CONJAC+i+NCOM*j)=0.0d0
          ENDDO
        ENDDO
        DO I=0,NT_RES-1
          REALWORK(OS_RESID+I)=0.0d0
        ENDDO
        DO I=0,NTCNTR-1
          REALWORK(OS_CONTR+I)=0.0d0
        ENDDO
      ELSE
        DO I=1,NOPM*NREM
          REALWORK(OS_RESJAC+I-1)=0.0d0
        ENDDO
        DO I=1,NOPM*NCOM
          REALWORK(OS_CONJAC+I-1)=0.0d0
        ENDDO
        DO I=1,NREM
          REALWORK(OS_RESID+I-1)=0.0d0
        ENDDO
        DO I=1,NCOM
          REALWORK(OS_CONTR+I-1)=0.0d0
        ENDDO
      ENDIF

C*** LKC 8-MAY-2003 storing the TIME and REGION into integer and real
C     work arrays for magnetic dipole inverses
C    LKC 9-JUL-2003 baseline for magnetic calculations
      REALWORK(OS_WORK)=TIME
      REALWORK(OS_WORK+1)=BASELINE
      INTWORK(OS_IWORK)=nr
      INTWORK(OS_IWORK+1)=SQUID_CONFIG

C*** SEN 19/2/03 changed so that calling is consistant between fitting and
C     other calls to nagmina.
      CALL NAGMINA(INTWORK(OS_ISTATE),REALWORK(OS_PAOPTI),
     '  REALWORK(OS_PMIN),REALWORK(OS_PMAX),REALWORK(OS_R),
     '  REALWORK(OS_F),REALWORK(OS_FJAC),REALWORK(OS_XC),
     '  INTWORK,REALWORK,.TRUE.,ERROR,*9999)

      CALL EXITS('NAGMIN')
      RETURN
 9999 CALL ERRORS('NAGMIN',ERROR)
      CALL EXITS('NAGMIN')
      RETURN 1
      END


