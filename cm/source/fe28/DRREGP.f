      SUBROUTINE DRREGP(ISEG,ISPLOTXY,REG_PARAMETER,WK5_INV,CSEG,STRING,
     '  ERROR,*)

C#### Subroutine: DRREGP
C###  Description:
C###    DRREGP draws the regularisation parameters over time.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'plot02.cmn'
      INCLUDE 'mach00.inc'
!     Parameter list
      INTEGER ISEG(*),ISPLOTXY(2)
      REAL*8 REG_PARAMETER(0:NTSM),WK5_INV(*)
      CHARACTER CSEG(*)*(*),STRING*(MXCH),ERROR*(*)
!     Local variables
      INTEGER i,IBEG,IEND,INDEX,iw,N3CO
      INTEGER*4 D_PTS_PTR
      LOGICAL CBBREV
      EXTERNAL POLYMARKER_DYNAM
      PARAMETER (iw=10)
!     Functions
      INTEGER INDEX_POLYMARKER

      CALL ENTERS('DRREGP',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C-----------------------------------------------------------------------

C#### Command: FEM draw reg-parameters
C###  Parameter:    <rgb=RGB[blue]>
C###    Specify the colour. The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Description:
C###    Draw the regularisation parameters on worksation 10.

        OP_STRING(1)=STRING(1:IEND)//' <rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C-----------------------------------------------------------------------

      ELSEIF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRREGP',ERROR,*9999)
      ELSE
        CALL ASSERT(EVALUATE_INVERSE,'>>Evaluate inverse first',
     '    ERROR,*9999)
        CALL ASSERT(USE_GRAPHICS.EQ.1,'>>Set USE_GRAPHICS=1',ERROR,
     '    *9999)

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYMARKER(0,'CIRCLE','SIZE1','BLUE')
        ENDIF
        DO i=1,INT(REG_PARAMETER(0))
          WK5_INV(i)=DBLE(i)
        ENDDO

        ! Initialisation
        NEXTPLOT='REPLACE'
        ISHOLD=.FALSE.
        XSCALE='LINEAR'
        XGRID=.TRUE.
        YGRID=.TRUE.
        TITLE='Regularisation Parameters'

        IF(ISTABILISE.EQ.2) THEN
          ! TGSVD
          YSCALE='LINEAR'
          XLABEL='time'
          YLABEL='k'
        ELSEIF(ISTABILISE.EQ.3) THEN
          ! Tikhonov
          YSCALE='LOG'
          XLABEL='time'
          YLABEL='lambda'
        ENDIF

        IF(ICOUPLING.EQ.2) THEN
          ! Greensite
          XLABEL='equation'
        ENDIF

        ! Plot the regularisation parameters
        XMAX=REG_PARAMETER(0)
        XMIN=1
        YMAX=REG_PARAMETER(1)
        YMIN=YMAX
        DO i=2,INT(REG_PARAMETER(0))
          YMAX=MAX(YMAX,REG_PARAMETER(i))
          YMIN=MIN(YMIN,REG_PARAMETER(i))
        ENDDO
        D_PTS_PTR=0
        CALL ALLOCATE_MEMORY(3*INT(REG_PARAMETER(0)),6,DPTYPE,D_PTS_PTR,
     '    MEM_INIT,ERROR,*9999)
        CALL ACWK(iw,0,ERROR,*9999)
        CALL SGPLOTXY(POLYMARKER_DYNAM,INDEX,ISEG,1,ISPLOTXY(1),iw,
     '    INT(REG_PARAMETER(0)),%VAL(D_PTS_PTR),WK5_INV,
     '    REG_PARAMETER(1),CSEG,ERROR,*9999)
        CALL SGPLOTXY(%VAL(0),INDEX,ISEG,2,ISPLOTXY(2),iw,%VAL(0),
     '    %VAL(D_PTS_PTR),%VAL(0),%VAL(0),CSEG,ERROR,*9999)
        CALL DAWK(iw,0,ERROR,*9999)
        CALL REFRESH_GRAPHICS(0.0,ERROR,*9999)
        CALL FREE_MEMORY(D_PTS_PTR,ERROR,*9999)
      ENDIF

      CALL EXITS('DRREGP')
      RETURN
 9999 CALL ERRORS('DRREGP',ERROR)
      CALL EXITS('DRREGP')
      RETURN 1
      END

