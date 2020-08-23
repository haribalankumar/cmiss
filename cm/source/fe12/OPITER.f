      SUBROUTINE OPITER(ERROR,*)

C#### Subroutine: OPITER
C###  Description:
C###    OPITER outputs iteration parameters.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'iter00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG1,IBEG2,IBEG3,IEND1,IEND2,IEND3
      CHARACTER TITLE1(2)*20,TITLE2(3)*40,TITLE3(3)*40,
     '  TITLE4(3)*20,TITLE5(2)*40

      DATA TITLE1/'solution',                         ! ktyp20=1
     '            'optimisation'/                     ! ktyp20=2

      DATA TITLE2/'forward problem calculations',     ! ktyp21=1,ktyp20=1
     '            'surface potential fitting ',       ! ktyp21=2    "
     '            ' '/                                ! ktyp21=3

      DATA TITLE3/'data fitting',                     ! ktyp21=1,ktyp20=2
     '            ' ',                                ! ktyp21=2    "
     '            ' '/                                ! ktyp21=3    "

      DATA TITLE4/'ipfiles',                          ! informat_code = 1
     '            'map3d',                            !               = 2
     '            'emap'/                             !               = 3

      DATA TITLE5/'files',                            ! init_code = 1
     '            'fitting electrode potentials'/     !           = 2

      CALL ENTERS('OPITER',*9999)

      IF(KTYP20.EQ.1) THEN
        CALL STRING_TRIM(TITLE1(KTYP20),IBEG1,IEND1)
        CALL STRING_TRIM(TITLE2(KTYP21),IBEG2,IEND2)
        WRITE(OP_STRING,'('' Iteration on '',A,'' wrt '',A)')
     '    TITLE1(KTYP20)(IBEG1:IEND1),TITLE2(KTYP21)(IBEG2:IEND2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP21.EQ.1) THEN !Forward problem calculations
          WRITE(OP_STRING,'('' Number of iterations : '',I4)')
     '      NUM_ITS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Base filename : '',A)') IT_FNAME
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL STRING_TRIM(TITLE4(INFORMAT_CODE),IBEG1,IEND1)
          OP_STRING(1)=' Input file format is '//
     '     TITLE4(INFORMAT_CODE)(IBEG1:IEND1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL STRING_TRIM(TITLE4(OUTFORMAT_CODE),IBEG1,IEND1)
          OP_STRING(1)=' Output file format is '//
     '     TITLE4(OUTFORMAT_CODE)(IBEG1:IEND1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL STRING_TRIM(TITLE5(INIT_CODE),IBEG1,IEND1)
          OP_STRING(1)=' Initial conditions obtained from '//
     '      TITLE5(INIT_CODE)(IBEG1:IEND1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Initial condition output '','
     '      //'''frequency : '',I4)') OUTPUT_FREQ(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Nodal solution output '','
     '      //'''frequency : '',I4)') OUTPUT_FREQ(2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP21.EQ.2) THEN !surface potential fitting
          WRITE(OP_STRING,'('' Number of iterations : '',I4)') NUM_ITS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Base filename : '',A)') IT_FNAME
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          CALL STRING_TRIM(TITLE4(INFORMAT_CODE),IBEG1,IEND1)
          OP_STRING(1)=' Input file format is '//
     '     TITLE4(INFORMAT_CODE)(IBEG1:IEND1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Fitted field values output '','
     '      //'''frequency : '',I4)') OUTPUT_FREQ(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP21.EQ.3) THEN
        ENDIF
      ELSE IF(KTYP20.EQ.2) THEN
        CALL STRING_TRIM(TITLE1(KTYP20),IBEG1,IEND1)
        CALL STRING_TRIM(TITLE3(KTYP21),IBEG3,IEND3)
        WRITE(OP_STRING,'('' Iteration on '',A,'' wrt '',A)')
     '    TITLE1(KTYP20)(IBEG1:IEND1),TITLE3(KTYP21)(IBEG3:IEND3)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP21.EQ.1) THEN
          WRITE(OP_STRING,'('' Maximum number of iterations : '',I4)')
     '      NUM_ITS
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Data RMS error tolerance : '',D11.4)')
     '      ITER_TOL(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Scale factor tolerance : '',D11.4)')
     '      ITER_TOL(2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP12.GT.0) THEN
            WRITE(OP_STRING,'('' Smooting update factor : '',D11.4)')
     '        ITER_STORE(1)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' Maximum data error : '',D11.4)')
     '      ITER_STORE(2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(WARMOPTI) THEN
            WRITE(OP_STRING,'('' Warm optimisations are used'')')
          ELSE
            WRITE(OP_STRING,'('' Warm optimisations are not used'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Base filename : '',A)') IT_FNAME
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Fitted nodal values output '','
     '      //'''frequency : '',I4)') OUTPUT_FREQ(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Data projection point output '','
     '      //'''frequency : '',I4)') OUTPUT_FREQ(2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP12.GT.0) THEN
            WRITE(OP_STRING,'('' Smooting weight update '','
     '        //'''frequency : '',I4)') OUTPUT_FREQ(3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(KTYP21.EQ.2) THEN
        ENDIF
      ELSEIF(KTYP20.EQ.3) THEN !extracellular solution iteration

        WRITE(OP_STRING,'('' Iteration type is              : '
     '    //'Extracellular Solution'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'('' Maximum potential tolerance is : '',F8.6)')
     '    ITER_TOL(1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Maximum flux tolerance is      : '',F8.6)')
     '    ITER_TOL(2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'('' Maximum number of iterations   : '',I4)')
     '    NUM_ITS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ENDIF

      CALL EXITS('OPITER')
      RETURN
 9999 CALL ERRORS('OPITER',ERROR)
      CALL EXITS('OPITER')
      RETURN 1
      END


