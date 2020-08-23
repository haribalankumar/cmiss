      SUBROUTINE OPTRSF(NPLIST3,NPLIST4,NPLIST5,ERROR,*)

C#### Subroutine: OPTRSF
C###  Description:
C###    OPTRSF outputs transfer matrix (from first surface/epicardium
C###    to a second surface/body surface) parameters.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER NPLIST3(0:NPM),NPLIST4(0:NPM),NPLIST5(0:NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nolist
      REAL*8 DUMMY(1)

      CALL ENTERS('OPTRSF',*9999)

      DUMMY(1)=0.0d0

      IF(KTYP95.EQ.1) THEN
        WRITE(OP_STRING,'('' Single layer Transfer Matrix'')')
      ELSEIF(KTYP95.EQ.2) THEN
        WRITE(OP_STRING,'('' Double layer Transfer Matrix'')')
      ENDIF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(KTYP94.EQ.1) THEN
        WRITE(OP_STRING,'('' Transfer matrix computed directly'')')
      ELSEIF(KTYP94.EQ.2) THEN
        WRITE(OP_STRING,'('' Transfer matrix computed algebraicly'')')
      ENDIF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
CC AJPs 11-11-97
      WRITE(OP_STRING,'('' The number of first surface nodes = '',I5)')
     '  NPLIST3(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' These nodes belong to region '',I1)')
     '  TRSF_NR_FIRST
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' The number of second surface nodes = '',I5)')
     '  NPLIST4(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' These nodes belong to region '',I1)')
     '  TRSF_NR_SECOND
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(KTYP95.EQ.3) THEN !double layer from transmembrane to epi
        WRITE(OP_STRING,'('' The number of outer surface nodes '
     '    //'= '',I5)')NPLIST5(0)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' These nodes belong to region '',I1)')
     '    TRSF_NR_OUTER
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'('' Total number of regions involved in the '
     '  //'transfer matrix calculation :'',I1)') TRSF_NRLIST(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' These regions are '',7(I1,2X))')
     '  (TRSF_NRLIST(nolist),nolist=1,TRSF_NRLIST(0))
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(/'' First/inner surface nodes: '')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C cpb 7/6/01 Use write_long to handle record overflow. Also nodes should be
C written out with I5 format.
C      WRITE(OP_STRING,'('' '',10I5,/:(1X,10I5))')
C     '  (NPLIST3(nolist),nolist=1,NPLIST3(0))
C      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      CALL WRITE_LONG(INTTYPE,1,1,IOFI,NPLIST3(0),10,10,NPLIST3(1),
     '  DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      WRITE(OP_STRING,'(/'' Second surface nodes: '')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C cpb 7/6/01 Use write_long to handle record overflow. Also nodes should be
C written out with I5 format.
C      WRITE(OP_STRING,'('' '',10I5,/:(1X,10I5))')
C     '  (NPLIST4(nolist),nolist=1,NPLIST4(0))
C      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      CALL WRITE_LONG(INTTYPE,1,1,IOFI,NPLIST4(0),10,10,NPLIST4(1),
     '  DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      IF(KTYP95.EQ.3) THEN !double layer from transmembrane to epi
        WRITE(OP_STRING,'(/'' Outer surface nodes: '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C cpb 7/6/01 Use write_long to handle record overflow. Also nodes should be
C written out with I5 format.
C        WRITE(OP_STRING,'('' '',10I5,/:(1X,10I5))')
C     '    (NPLIST5(nolist),nolist=1,NPLIST4(0))
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        CALL WRITE_LONG(INTTYPE,1,1,IOFI,NPLIST5(0),10,10,NPLIST5(1),
     '    DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      ENDIF
CC AJPe

      WRITE(OP_STRING,'('' The sampling frequency is :'',F11.4)')
     '  TRSF_FREQUENCY
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(KTYP95.EQ.2.OR.KTYP95.EQ.3) THEN !Double layer

        IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function

          WRITE(OP_STRING,'('' The type of activation wavefront is a '
     '      //'Heavyside step function'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The transmembrane resting potential '
     '      //'is :'',F11.4)') TRSF_ACTN_WAVE_REST
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The transmembrane jump is :'',F11.4)')
     '      TRSF_ACTN_WAVE_JUMP
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function

          WRITE(OP_STRING,'('' The type of activation wavefront is a '
     '      //'sigmoidal function'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The transmembrane resting potential '
     '      //'is :'',F11.4)') TRSF_ACTN_WAVE_REST
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The transmembrane jump is :'',F11.4)')
     '      TRSF_ACTN_WAVE_JUMP
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
           WRITE(OP_STRING,'('' The wavefront width is :'',F11.4)')
     '      TRSF_ACTN_WAVE_WIDTH
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !Arctan function

          WRITE(OP_STRING,'('' The type of activation wavefront is '
     '      //'an arctangent function'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The transmembrane resting potential '
     '      //'is :'',F11.4)') TRSF_ACTN_WAVE_REST
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' The transmembrane jump is :'',F11.4)')
     '      TRSF_ACTN_WAVE_JUMP
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
           WRITE(OP_STRING,'('' The wavefront width is :'',F11.4)')
     '      TRSF_ACTN_WAVE_WIDTH
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ENDIF

      ENDIF

      CALL EXITS('OPTRSF')
      RETURN
 9999 CALL ERRORS('OPTRSF',ERROR)
      CALL EXITS('OPTRSF')
      RETURN 1
      END


