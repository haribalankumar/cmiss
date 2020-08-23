      SUBROUTINE MINOSOPTIA(HA,HS,KA,A,PAOPTI,PI,PMIN,PMAX,RC,
     '  X_N,Z,IUSER,USER,ERROR,*)

C#### Subroutine: MINOSOPTIA
C###  Description:
C###    MINOSOPTIA is a minimisation routine using MINOS routine.

C**** KTYP26=2 KTYP27=5: (Geometric data fitting)

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'time02.cmn'
!     Parameter List
      INTEGER KA(*),IUSER(*)
      INTEGER*4 HA(*),HS(*)
      REAL*8 A(*),PAOPTI(*),PI(*),PMIN(*),PMAX(*),RC(*),
     '  USER(*),X_N(*),Z(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NUMNAMES
      PARAMETER (NUMNAMES=1)
      INTEGER i,IBEG,IEND,
     '  INFORM,IOBJ,ISPEC,IPRINT,
     '  ISUMM,m,MINCOR,n,NAME1(NUMNAMES),NAME2(NUMNAMES),
     '  nb,ne,NINF,NNAME,NNCON,NNJAC,NNOBJ,
     '  noopti,ns,NWCORE
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 OBJADD,SINF,OBJ
      CHARACTER NAMES(NUMNAMES)*8
      EXTERNAL CFUNC5,CFUNC6,FUNCT5,MATMOD

      CALL ENTERS('MINOSOPTIA',*9999)

      IF(KTYP26.EQ.2.AND.KTYP27.EQ.5) THEN !Data Fitting

C
C Set the initial parameter values and bounds on the variables
C
        DO noopti=1,NTOPTI
          X_N(noopti)=PAOPTI(noopti)
        ENDDO

C
C Set the output/input file unit numbers
C
        IPRINT=IOOP
        ISPEC=0
        IF(SUMMFILE) THEN
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          CALL OPENF(IOFILE6,'DISK',FILE00(IBEG:IEND)//'.misumm',
     '      'UNKNOWN','SEQUEN','FORMATTED',132,ERROR,*9999)
          ISUMM=IOFILE6
        ELSE
          ISUMM=0
        ENDIF
C
C Set the problem specification options
C
        CALL MIOPT(  'Defaults            ',IPRINT,ISUMM,INFORM)
        CALL MIOPTI( 'Print Level         ',IPPLEV,IPRINT,ISUMM,INFORM)
        CALL MIOPTI( 'Derivative Level    ',IPDLEV,IPRINT,ISUMM,INFORM)
        CALL MIOPTI( 'Verify Level        ',IPVLEV,IPRINT,ISUMM,INFORM)
        CALL MIOPTR( 'Optimality Tolerance',OPTTOL,IPRINT,ISUMM,INFORM)
        CALL MIOPTR( 'Row Tolerance       ',NLFTOL,IPRINT,ISUMM,INFORM)
        CALL MIOPTR( 'Linesearch Tolerance',LNSTOL,IPRINT,ISUMM,INFORM)
        CALL MIOPTR( 'Function Precision  ',1.0D-12,IPRINT,ISUMM,INFORM)
        IF(SPARSEJAC) THEN
          CALL MIOPT('Jacobian Sparse     ',IPRINT,ISUMM,INFORM)
        ELSE
          CALL MIOPT('Jacobian Dense      ',IPRINT,ISUMM,INFORM)
        ENDIF
        CALL MIOPTI( 'Debug Level         ',DBGLEV,IPRINT,ISUMM,INFORM)
        CALL MIOPTI( 'Iterations Limit    ',ITERLM,IPRINT,ISUMM,INFORM)
C
C Set the problem dimensions
C
        m=NTCNTR
        n=NTOPTI
        nb=m+n
        NNAME=NUMNAMES
        NNCON=NTCNTR
        NNOBJ=NTOPTI
        NNJAC=NTOPTI
        NWCORE=L_ZM
        IOBJ=0
        OBJADD=0.0D0
C
C Set up the jacobian pointer arrays to define the matrix
C
        IF(SPARSEJAC) THEN
C         Do nothing at the moment
          ne=1
        ELSE
          ne=NTOPTI*NTCNTR
          KA(1)=1
          KA(NTOPTI+1)=ne+1
          DO noopti=2,NTOPTI
            KA(noopti)=KA(noopti-1)+NTCNTR
            DO i=1,NTCNTR
              HA(i+KA(noopti-1)-1)=i
            ENDDO
          ENDDO
          DO i=1,NTCNTR
             HA(i+KA(NTOPTI)-1)=i
          ENDDO
        ENDIF
C
C Solve the problem
C
        CALL CPU_TIMER(CPU_USER,TIME_START)

        IF(SPARSEJAC) THEN
          CALL MINOS('Cold',M,N,nb,ne,NNAME,NNCON,NNOBJ,NNJAC,IOBJ,
     '      OBJADD,NAMES,A,HA,KA,PMIN,PMAX,NAME1,NAME2,HS,X_N,PI,RC,
     '      INFORM,MINCOR,ns,NINF,SINF,OBJ,Z,NWCORE,ISPEC,IPRINT,
     '      ISUMM,CFUNC6,FUNCT5,MATMOD,IUSER,USER)
        ELSE
          CALL MINOS('Cold',M,N,nb,ne,NNAME,NNCON,NNOBJ,NNJAC,IOBJ,
     '      OBJADD,NAMES,A,HA,KA,PMIN,PMAX,NAME1,NAME2,HS,X_N,PI,RC,
     '      INFORM,MINCOR,ns,NINF,SINF,OBJ,Z,NWCORE,ISPEC,IPRINT,
     '      ISUMM,CFUNC5,FUNCT5,MATMOD,IUSER,USER)
        ENDIF

        CALL CPU_TIMER(CPU_USER,TIME_STOP)
        ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

        IF(SUMMFILE) THEN
          CALL CLOSEF(IOFILE6,ERROR,*9999)
        ENDIF

        IF(inform.ne.0) THEN
          IF(INFORM.EQ.42) THEN
            WRITE(OP_STRING,'('' NWCORE='',I7,'' MINCOR='',I7)')
     '        NWCORE,MINCOR
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE
          WRITE(OP_STRING,'(/'' Objective function = '',D12.4)') OBJ
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Solution time:'',D11.4,'' seconds'')')
     '      ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

      ENDIF

      CALL EXITS('MINOSOPTIA')
      RETURN
 9999 CALL ERRORS('MINOSOPTIA',ERROR)
      CALL EXITS('MINOSOPTIA')
      RETURN 1
      END


