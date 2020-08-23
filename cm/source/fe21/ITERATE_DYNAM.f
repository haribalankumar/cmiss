      SUBROUTINE ITERATE_DYNAM(ISEG,RESID,RESIDM,
     '  CSEG,ECHO,END,STRING,INTWORK,REALWORK,ERROR,*)

C#### Subroutine: ITERATE_DYNAM
C###  Description:
C###    ITERATE_DYNAM handles all the iteration calls to FEM
C###    except coupled bidomain iterations.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'iter00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER ISEG(*),INTWORK(*)
      REAL*8 RESID(*),RESIDM(*),REALWORK(*)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL ECHO,END
!     Local Variables
      INTEGER I,IBEG1,IBEG2,IBEG3,IBEG4,IEND1,IEND2,IEND3,IEND4,
     '  it_count,nores
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 SUMSQUARE
      CHARACTER CHAR1*4,FILENAME*25

      CALL ENTERS('ITERATE_DYNAM',*9999)

      IF(KTYP20.EQ.1) THEN !Iterate on solution
        IF(KTYP21.EQ.1) THEN !Forward problem calaculation
C******** ALTERED FOR FITTED FIELDS. AJP 11-5-94
C          CALL ASSERT(KTYP8.GT.0,'>>Define fit first',ERROR,*9999)
C******** END FIRST ALTERATION
C ajp 5-3-94   IF(INIT_CODE.EQ.2) THEN
C            CALL ASSERT(CALL_INIT,'>>Define b.c.s first',ERROR,*9999)
C          ENDIF
          CALL STRING_TRIM(IT_FNAME,IBEG1,IEND1)

C KAT 2001-04-09: UPDATE_MATRIX never used
C          UPDATE_MATRIX(1)=.TRUE.
          DO it_count=1,NUM_ITS
            IF(ECHO) THEN
              WRITE(OP_STRING,'('' ### Begining iteration '','
     '          //'''number : '',I4)') it_count
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(CHAR1(1:4),'(I4)') it_count
            CALL STRING_TRIM(CHAR1,IBEG2,IEND2)

C*** obtain the initial conditions

            IF(INIT_CODE.EQ.1) THEN !obtain initial conditions from files

C***    read in the initial conditions from a file

              FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '          CHAR1(IBEG2:IEND2)
              CO(1)='FEM'
              CO(2)='DEFINE'
              CO(3)='INITIAL'
              COQU(3,1)='r'
              COQU(3,2)=FILENAME
              NTCO=3
              NTCOQU(3)=2
              noco=1
              IF(ECHO) THEN
                OP_STRING(1)=' ### FEM DEFINE INITIAL;r;'//FILENAME
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
              NTCOQU(3)=0
              COQU(3,1)=' '
              COQU(3,2)=' '

            ELSE IF(INIT_CODE.EQ.2) THEN !obtain init. cond. from fitting

C 24/2/97 LC removed section :
C                   AJP 11-5-94.  Altered to read in fitted field files.

C*** NEW TO READ IN FITTED FIELD INFO
              IF(INFORMAT_CODE.EQ.1) THEN !ipfiles
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
              ELSE IF(INFORMAT_CODE.EQ.2) THEN !map3d
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
              ENDIF
              CO(1)='FEM'
              CO(2)='DEFINE'
              CO(3)='FIELD'
              COQU(3,1)='r'
              COQU(3,2)=FILENAME
              NTCO=3
              NTCOQU(3)=2
              noco=1
              IF(ECHO) THEN
                CALL STRING_TRIM(FILENAME,IBEG3,IEND3)
                OP_STRING(1)=' ### FEM DEFINE FIELD;r;'//
     '            FILENAME(IBEG3:IEND3)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
              NTCOQU(3)=0
              COQU(3,1)=' '
              COQU(3,2)=' '
C END NEW


C***    set up fix array (changed in fit)

              CO(1)='FEM'
              CO(2)='DEFINE'
              CO(3)='INITIAL'
              COQU(3,1)='r'
c              COQU(3,2)='tank' !Needs generalising
              COQU(3,2)='tank_ref' !Needs generalising
              NTCO=3
              NTCOQU(3)=2
              noco=1
              IF(ECHO) THEN
c                OP_STRING(1)=' ### FEM DEFINE INITIAL;r;tank'
                OP_STRING(1)=' ### FEM DEFINE INITIAL;r;tank_ref'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)
              NTCOQU(3)=0
              COQU(3,1)=' '
              COQU(3,2)=' '

C AJP 11-5-94.  Altered to read in fitted field files.
cC***    update the initial conditions from the field
c
c              CO(1)='FEM'
c              CO(2)='UPDATE'
c              CO(3)='INITIAL'
c              CO(4)='FROM'
c              CO(5)='ITERATION'
c              CO(6)='REGION'
c              CO(7)='1'
c              NTCO=7
c              noco=1
c              IF(ECHO) THEN
c                OP_STRING(1)=' ### FEM UPDATE INITIAL FROM '//
c     '            'ITERATION FOR 1'
c                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c              ENDIF
c              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
c     '          ERROR,*9999)
c
c

C*** NEW TO READ IN FITTED FIELD INFO
C***    update the initial conditions from the field

              CO(1)='FEM'
              CO(2)='UPDATE'
              CO(3)='INITIAL'
              CO(4)='FROM'
              CO(5)='FIELD'
              NTCO=5
              noco=1
              IF(ECHO) THEN
                OP_STRING(1)=' ### FEM UPDATE INITIAL FROM '//
     '            'FIELD'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)

C*** END NEW

C***    write out the initial conditions if required

              IF(OUTPUT_FREQ(1).ne.0) THEN
                IF(MOD(it_count,OUTPUT_FREQ(1)).EQ.0) THEN
                  FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '              CHAR1(IBEG2:IEND2)
                  CO(1)='FEM'
                  CO(2)='DEFINE'
                  CO(3)='INITIAL'
                  COQU(3,1)='w'
                  COQU(3,2)=FILENAME
                  NTCO=3
                  NTCOQU(3)=2
                  noco=1
                  IF(ECHO) THEN
                    OP_STRING(1)=' ### FEM DEFINE INITIAL;w;'//
     '                FILENAME
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '              ERROR,*9999)
                  NTCOQU(3)=0
                  COQU(3,1)=' '
                  COQU(3,2)=' '

                 ENDIF
              ENDIF
            ENDIF

CC***    list initial for debug
C
C              CO(1)='FEM'
C              CO(2)='LIST'
C              CO(3)='INITIAL'
C              CO(4)='FULL'
C              NTCO=4
C              noco=1
C              IF(ECHO) THEN
C                OP_STRING(1)=' ### FEM LIST INITIAL FULL '
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              ENDIF
C              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
C     '          ERROR,*9999)
C
C*** obtain solve parameters

            FILENAME=IT_FNAME(IBEG1:IEND1)
            CO(1)='FEM'
            CO(2)='DEFINE'
            CO(3)='SOLVE'
            COQU(3,1)='r'
            COQU(3,2)=FILENAME
            NTCO=3
            NTCOQU(3)=2
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM DEFINE SOLVE;r;'//FILENAME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '        ERROR,*9999)
            NTCOQU(3)=0
            COQU(3,1)=' '
            COQU(3,2)=' '

CC***    list solve for debug
C
C              CO(1)='FEM'
C              CO(2)='LIST'
C              CO(3)='SOLVE'
C              NTCO=3
C              noco=1
C              IF(ECHO) THEN
C                OP_STRING(1)=' ### FEM LIST SOLVE '
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              ENDIF
C              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
C     '          ERROR,*9999)

C*** solve the problem

            CO(1)='FEM'
            CO(2)='SOLVE'
            NTCO=2
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM SOLVE'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '        ERROR,*9999)

C*** output solution values if required

            IF(OUTPUT_FREQ(2).ne.0) THEN
              IF(MOD(it_count,OUTPUT_FREQ(2)).EQ.0) THEN
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
                CO(1)='FEM'
                CO(2)='DEFINE'
                CO(3)='NODE'
                COQU(3,1)='w'
                COQU(3,2)=FILENAME
                NTCO=3
                NTCOQU(3)=2
                noco=1
                IF(ECHO) THEN
                  OP_STRING(1)=' ### FEM DEFINE NODE;w;'//FILENAME
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '            ERROR,*9999)
                NTCOQU(3)=0
                COQU(3,1)=' '
                COQU(3,2)=' '
C cpb 11/2/94 putting in this next command just for debugging
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
                CO(1)='FEM'
                CO(2)='LIST'
                CO(3)='NODE'
                COQU(3,1)=FILENAME
                CO(4)='SOLUTION'
                NTCO=4
                NTCOQU(3)=1
                noco=1
                IF(ECHO) THEN
                  CALL STRING_TRIM(FILENAME,IBEG3,IEND3)
                  OP_STRING(1)=' ### FEM LIST NODE;'//
     '              FILENAME(IBEG3:IEND3)//' SOLUTION'
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '            ERROR,*9999)
                NTCOQU(3)=0
                COQU(3,1)=' '
              ENDIF
            ENDIF

C***    update potential from solution

              CO(1)='FEM'
              CO(2)='UPDATE'
              CO(3)='POTENTIAL'
              CO(4)='FROM'
              CO(5)='INITIAL'
              CO(6)='DATASET'
              CO(7)=CHAR1
              NTCO=7
              noco=1
              IF(ECHO) THEN
                OP_STRING(1)=' ### FEM UPDATE POTENTIAL FROM '//
     '            'INITIAL DATASET '//CHAR1(IBEG2:IEND2)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '          ERROR,*9999)

C KAT 2001-04-09: UPDATE_MATRIX never used
C            UPDATE_MATRIX(1)=.FALSE.
          ENDDO

        ELSE IF(KTYP21.EQ.2) THEN !Surface potential fitting
          CALL ASSERT(KTYP8.EQ.4,'>>Define potential fit first',
     '      ERROR,*9999)
          CALL STRING_TRIM(IT_FNAME,IBEG1,IEND1)

C KAT 2001-04-09: UPDATE_MATRIX never used
C          UPDATE_MATRIX(1)=.TRUE.
          DO it_count=1,NUM_ITS
            IF(ECHO) THEN
              WRITE(OP_STRING,'('' ### Begining iteration '','
     '          //'''number : '',I4)') it_count
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(CHAR1(1:4),'(I4)') it_count
            CALL STRING_TRIM(CHAR1,IBEG2,IEND2)

C***    read the potentials values from the electrode dataset

            IF(INFORMAT_CODE.EQ.1) THEN !ipfiles
              FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '          CHAR1(IBEG2:IEND2)
              CO(5)='IPFILE'
            ELSE IF(INFORMAT_CODE.EQ.2) THEN !map3d
              FILENAME=IT_FNAME(IBEG1:IEND1)
              CO(5)='MAP3D'
            ELSE IF(INFORMAT_CODE.EQ.3) THEN !emap
              CO(5)='EMAP'
            ENDIF
            CO(1)='FEM'
            CO(2)='DEFINE'
            CO(3)='POTENTIAL'
            COQU(3,1)='r'
            COQU(3,2)=FILENAME
            CO(4)='DATA_FORMAT'
            NTCO=5
            NTCOQU(3)=2
            noco=1
            IF(ECHO) THEN
              CALL STRING_TRIM(FILENAME,IBEG3,IEND3)
              CALL STRING_TRIM(CO(5),IBEG4,IEND4)
              OP_STRING(1)=' ### FEM DEFINE POTENTIAL;r;'//
     '          FILENAME(IBEG3:IEND3)//' '//
     '          ' DATA_FORMAT '//CO(5)(IBEG4:IEND4)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '        ERROR,*9999)
            NTCOQU(3)=0
            COQU(3,1)=' '
            COQU(3,2)=' '

C***    fit a field to the electrode data

            CO(1)='FEM'
            CO(2)='FIT'
            CO(3)='POTENTIAL'
            NTCO=3
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM FIT POTENTIAL'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '        ERROR,*9999)

            IF(OUTPUT_FREQ(1).ne.0) THEN
              IF(MOD(it_count,OUTPUT_FREQ(1)).EQ.0) THEN
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
                CO(1)='FEM'
                CO(2)='DEFINE'
                CO(3)='FIELD'
                COQU(3,1)='w'
                COQU(3,2)=FILENAME
                NTCO=3
                NTCOQU(3)=2
                noco=1
                IF(ECHO) THEN
                  OP_STRING(1)=' ### FEM DEFINE FIELD;w;'//FILENAME
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '            ERROR,*9999)
                NTCOQU(3)=0
                COQU(3,1)=' '
                COQU(3,2)=' '
              ENDIF
            ENDIF

C KAT 2001-04-09: UPDATE_MATRIX never used
C            UPDATE_MATRIX(1)=.FALSE.
          ENDDO

        ENDIF

C 25/2/97 LC removed section : OLD??? AJP 1-6-94

      ELSE IF(KTYP20.EQ.2) THEN !Iterate on optimisation
        IF(KTYP21.EQ.1) THEN !Data fitting by optimisation iteration
          CALL ASSERT(KTYP26.EQ.2.AND.KTYP27.EQ.5,
     '      '>>Define data fitting by optimisation first',ERROR,*9999)

          it_count=0
          ITER_CONVERGED(1)=.FALSE.
          ITER_CONVERGED(2)=.FALSE.
          CALL STRING_TRIM(IT_FNAME,IBEG1,IEND1)

          CALL CPU_TIMER(CPU_USER,TIME_START)

C Calculate begining residuals

          CO(1)='FEM'
          CO(2)='EVALUATE'
          CO(3)='RESIDUALS'
          CO(4)='WRT'
          CO(5)='DATA'
          NTCO=5
          noco=1
          IF(ECHO) THEN
            OP_STRING(1)=' ### FEM EVALUATE RESIDUALS WRT DATA'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

          SUMSQUARE=0.0d0
          IF(KTYP29.EQ.1) THEN
            DO nores=1,NT_RES-KTYP12
              SUMSQUARE=SUMSQUARE+RESID(nores)
            ENDDO
          ELSE IF(KTYP29.EQ.2) THEN
            DO nores=1,NT_RES-KTYP12
              SUMSQUARE=SUMSQUARE+RESIDM(nores)
            ENDDO
          ENDIF
          ITER_ERR(1)=DSQRT(SUMSQUARE/DBLE(NT_RES-KTYP12))

          IF(ECHO) THEN
            WRITE(OP_STRING,'('' ### Begining data RMS error = '','
     '        //'D11.4)') ITER_ERR(1)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          ITER_ERR(2)=0.0d0

          DO WHILE(it_count.LT.NUM_ITS.AND.(.NOT.ITER_CONVERGED(1).OR.
     '      .NOT.ITER_CONVERGED(2)))

            it_count=it_count+1

            IF(ECHO) THEN
              WRITE(OP_STRING,'('' ### Begining iteration '','
     '          //'''number : '',I4)') it_count
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(CHAR1(1:4),'(I4)') it_count
            CALL STRING_TRIM(CHAR1,IBEG2,IEND2)

C Perform the optimisation

            CO(1)='OPTIMISE'
            IF(it_count.GT.1.AND.WARMOPTI) THEN
              CO(2)='WARM'
              NTCO=2
              IF(ECHO) THEN
                OP_STRING(1)=' ### OPTIMISE WARM'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ELSE
              NTCO=1
              IF(ECHO) THEN
                OP_STRING(1)=' ### OPTIMISE'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            CALL SYNTAX(ISEG,CSEG,END,ERROR,*9999)

C Update the scale factors

            CO(1)='FEM'
            CO(2)='UPDATE'
            CO(3)='SCALE_FACTORS'
            CO(4)='ITERATE'
            NTCO=4
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM UPDATE SCALE_FACTORS ITERATE'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

C If required write out the nodes

            IF(OUTPUT_FREQ(1).ne.0) THEN
              IF(MOD(it_count,OUTPUT_FREQ(1)).EQ.0) THEN
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
                CO(1)='FEM'
                CO(2)='DEFINE'
                CO(3)='NODE'
                COQU(3,1)='w'
                COQU(3,2)=FILENAME
                NTCO=3
                NTCOQU(3)=2
                noco=1
                IF(ECHO) THEN
                  OP_STRING(1)=' ### FEM DEFINE NODE;w;'//FILENAME
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '            ERROR,*9999)
                NTCOQU(3)=0
                COQU(3,1)=' '
                COQU(3,2)=' '
              ENDIF
            ENDIF

C Update the data point projections

            CO(1)='FEM'
            CO(2)='DEFINE'
            CO(3)='DATA'
            COQU(3,1)='c'
            CO(4)='XI'
            CO(5)='OR'
            CO(6)='OLD'
            NTCO=6
            NTCOQU(3)=1
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM DEFINE DATA;c XI ORTHOGONAL OLD'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)
            NTCOQU(3)=0
            COQU(3,1)=' '

C Cancel any data point projections over the maximum error limit

            CO(1)='FEM'
            CO(2)='CANCEL'
            CO(3)='DATA'
            CO(4)='ERROR'
            CO(5)='GREATER'
            WRITE(CO(6),'(D12.4)') ITER_STORE(2)
            NTCO=6
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM CANCEL DATA ERROR '
     '          //'GREATER '//CO(6)(1:12)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

C If required write out the data point projections

            IF(OUTPUT_FREQ(2).ne.0) THEN
              IF(MOD(it_count,OUTPUT_FREQ(2)).EQ.0) THEN
                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
     '            CHAR1(IBEG2:IEND2)
                CO(1)='FEM'
                CO(2)='DEFINE'
                CO(3)='DATA'
                COQU(3,1)='w'
                COQU(3,2)=FILENAME
                CO(4)='XI'
                NTCO=4
                NTCOQU(3)=2
                noco=1
                IF(ECHO) THEN
                  CALL STRING_TRIM(FILENAME,IBEG3,IEND3)
                  OP_STRING(1)=' ### FEM DEFINE DATA;w;'//
     '              FILENAME(IBEG3:IEND3)//' XI'
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '            ERROR,*9999)
                NTCOQU(3)=0
                COQU(3,1)=' '
                COQU(3,2)=' '
              ENDIF
            ENDIF

C Define optimisation parameters

            FILENAME=IT_FNAME(IBEG1:IEND1)
            CO(1)='FEM'
            CO(2)='DEFINE'
            CO(3)='OPTIMISE'
            COQU(3,1)='r'
            COQU(3,2)=FILENAME
            NTCO=3
            NTCOQU(3)=2
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM DEFINE OPTIMISE;r;'//FILENAME
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)
            NTCOQU(3)=0
            COQU(3,1)=' '
            COQU(3,2)=' '

C Calculate current residuals

            CO(1)='FEM'
            CO(2)='EVALUATE'
            CO(3)='RESIDUALS'
            CO(4)='WRT'
            CO(5)='DATA'
            NTCO=5
            noco=1
            IF(ECHO) THEN
              OP_STRING(1)=' ### FEM EVALUATE RESIDUALS WRT DATA'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,ERROR,*9999)

C Check convergence

            SUMSQUARE=0.0d0
            IF(KTYP29.EQ.1) THEN
              DO nores=1,NT_RES-KTYP12
                SUMSQUARE=SUMSQUARE+RESID(nores)
              ENDDO
            ELSE IF(KTYP29.EQ.2) THEN
              DO nores=1,NT_RES-KTYP12
                SUMSQUARE=SUMSQUARE+RESIDM(nores)
              ENDDO
            ENDIF
            ITER_ERR_OLD(1)=ITER_ERR(1)
            ITER_ERR(1)=DSQRT(SUMSQUARE/DBLE(NT_RES-KTYP12))

            IF(ECHO) THEN
              WRITE(OP_STRING,'('' ### Iteration data RMS error = '','
     '          //'D11.4)') ITER_ERR(1)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

            DO I=1,2
              ITER_RELERR(I)=DABS(ITER_ERR(I)-ITER_ERR_OLD(I))/
     '          (1.d0+ITER_ERR_OLD(I))
              ITER_CONVERGED(I)=ITER_RELERR(I).LE.ITER_TOL(I)
            ENDDO

            IF(ECHO) THEN
              WRITE(OP_STRING,'('' ### Error tolerances:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ### Data RMS relative error    : '','
     '          //'''actual ='',D11.4,'' required = '',D11.4)')
     '          ITER_RELERR(1),ITER_TOL(1)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ### Scale factor relative error: '','
     '          //'''actual ='',D11.4,'' required = '',D11.4)')
     '          ITER_RELERR(2),ITER_TOL(2)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

C If required update the smooting weights

            IF(KTYP12.GT.0) THEN
              IF(OUTPUT_FREQ(3).NE.0) THEN
                IF(MOD(it_count,OUTPUT_FREQ(3)).EQ.0) THEN
                  CO(1)='FEM'
                  CO(2)='UPDATE'
                  CO(3)='SOBOLEV'
                  CO(4)='FROM'
                  CO(5)='FACTOR'
                  CO(6)='VALUE'
                  WRITE(CO(7),'(D12.4)') ITER_STORE(1)
                  NTCO=7
                  noco=1
                  IF(ECHO) THEN
                    OP_STRING(1)=' ### FEM UPDATE SOBOLEV '
     '                //'FROM FACTOR VALUE '//CO(7)(1:12)
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
     '              ERROR,*9999)
                ENDIF
              ENDIF
            ENDIF

          ENDDO !it_count

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

          IF(.NOT.ITER_CONVERGED(1).OR..NOT.ITER_CONVERGED(2)) THEN
            IF(.NOT.ECHO) THEN
              WRITE(OP_STRING,'('' ### Error tolerances:'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ### Data RMS relative error:     '','
     '          //'''actual ='',D11.4,'' required = '',D11.4)')
     '          ITER_RELERR(1),ITER_TOL(1)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' ### Scale factor relative error: '','
     '          //'''actual ='',D11.4,'' required = '',D11.4)')
     '          ITER_RELERR(2),ITER_TOL(2)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(OP_STRING,'('' *** Maximum number of iterations '','
     '        //'''reached. Convergence not achieved'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          WRITE(OP_STRING,'(/,'' ### Total number of iterations: '','
     '      //'I4)') it_count
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ### Solution time         : '',D11.4,'
     '      //''' seconds'')') ELAPSED_TIME
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' ### Average iteration time: '',D11.4,'
     '      //''' seconds'')') DBLE(ELAPSED_TIME)/DBLE(it_count)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('ITERATE_DYNAM')
      RETURN
 9999 CALL ERRORS('ITERATE_DYNAM',ERROR)
      CALL EXITS('ITERATE_DYNAM')
      RETURN 1
      END


