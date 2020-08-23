      SUBROUTINE IPITER(ERROR,*)

C#### Subroutine: IPITER
C###  Description:
C###    IPITER inputs iteration parameters.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iter00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter list
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER IBEG,ICHAR,IEND,INFO,NOQUES
      CHARACTER BASE_FNAME*15
      LOGICAL FILEIP

      CALL ENTERS('IPITER',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='('' Specify whether iterating on [1]:'''//
     '  '/''   (1) Solution               '''//
     '  '/''   (2) Optimisation           '''//
     '  '/''   (3) Circular Bidomain      '''//
     '  '/''   (4) Simple Slice Bidomain  '''//
     '  '/''  *(5) Full Slice Bidomain    '''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP20
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP20=IDATA(1)

      IF(KTYP20.LT.3) THEN
        IF(KTYP20.EQ.1) THEN      !iterating on solution
          FORMAT='('' Specify type of solution iteration [1]: '''//
     '      '/''   (1) Forward problem calculation'''//
     '      '/''   (2) Surface potential fitting'''//
     '      '/''   (3) Unused'''//
     '        '/$,''    '',I1)'
        ELSE IF(KTYP20.EQ.2) THEN !iterating on optimisation
          FORMAT='('' Specify type of optimisation iteration [1]: '''//
     '      '/''   (1) Data fitting'''//
     '      '/''   (2) Unused'''//
     '      '/''   (3) Unused'''//
     '      '/$,''    '',I1)'
        ENDIF
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP21
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP21=IDATA(1)
      ENDIF

      IF(KTYP20.EQ.1) THEN
        IF(KTYP21.EQ.1) THEN !forward problem calculations

          IDEFLT(1)=1
          FORMAT='($,'' Enter the number of iterations [1]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_ITS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_ITS=IDATA(1)

          CALL STRING_TRIM(FILE00,IBEG,IEND)
          CDEFLT(1)=FILE00(IBEG:IEND)
          BASE_FNAME(1:15)=' '
          BASE_FNAME(15-(IEND-IBEG):15)=FILE00(IBEG:IEND)
          FORMAT='($,'' Enter the base filename ['//
     '      BASE_FNAME(1:15)//']: '',A)'
          IF(IOTYPE.EQ.3) CDATA(1)=IT_FNAME
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,15,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c         CALL CINOUT(IOTYPE,IOIP,IFILE,FORMAT,1,CDATA,CDEFLT,15,
c    '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IT_FNAME=CDATA(1)(1:15)

          FORMAT='('' Specify the input format [1]:'''//
     '      '/''   (1) IP files'''//
     '      '/''   (2) Map3d'''//
     '      '/''   (3) Emap'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=INFORMAT_CODE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) INFORMAT_CODE=IDATA(1)

          FORMAT='('' Specify the output format [1]:'''//
     '      '/''   (1) IP files'''//
     '      '/''   (2) Map3d'''//
     '      '/''   (3) Emap'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=OUTFORMAT_CODE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) OUTFORMAT_CODE=IDATA(1)

          FORMAT='('' Specify whether initial conditions are '','
     '      //'''obtained from [1]:'''
     '      //'/''   (1) Files'''
     '      //'/''   (2) Fitting electrode data'''
     '      //'/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=INIT_CODE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) INIT_CODE=IDATA(1)

          IF(INIT_CODE.EQ.2) THEN
            IDEFLT(1)=1
            FORMAT='($,'' Enter the initial condition output '','
     '        //'''frequency [1]: '',I4)'
            IF(IOTYPE.EQ.3) IDATA(1)=OUTPUT_FREQ(1)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IDEFLT,0,9999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '        INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) OUTPUT_FREQ(1)=IDATA(1)
          ENDIF

          IDEFLT(1)=1
          FORMAT='($,'' Enter the nodal solution output '','
     '      //'''frequency [1]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=OUTPUT_FREQ(2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) OUTPUT_FREQ(2)=IDATA(1)

        ELSE IF(KTYP21.EQ.2) THEN !Surface potential fitting

          IDEFLT(1)=1
          FORMAT='($,'' Enter the number of iterations [1]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_ITS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_ITS=IDATA(1)

          CALL STRING_TRIM(FILE00,IBEG,IEND)
          CDEFLT(1)=FILE00(IBEG:IEND)
          BASE_FNAME(1:15)=' '
          BASE_FNAME(15-(IEND-IBEG):15)=FILE00(IBEG:IEND)
          FORMAT='($,'' Enter the base filename ['//
     '      BASE_FNAME(1:15)//']: '',A)'
          IF(IOTYPE.EQ.3) CDATA(1)=IT_FNAME
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,15,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c         CALL CINOUT(IOTYPE,IOIP,IFILE,FORMAT,1,CDATA,CDEFLT,15,
c    '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IT_FNAME=CDATA(1)(1:15)

          FORMAT='('' Specify the input format [1]:'''//
     '      '/''   (1) IP files'''//
     '      '/''   (2) Map3d'''//
     '      '/''   (3) Emap'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=INFORMAT_CODE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) INFORMAT_CODE=IDATA(1)

          IDEFLT(1)=1
          FORMAT='($,'' Enter the fitted field values output '','
     '      //'''frequency [1]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=OUTPUT_FREQ(1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,0,9999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) OUTPUT_FREQ(1)=IDATA(1)

        ELSE IF(KTYP21.EQ.3) THEN
        ENDIF
      ELSE IF(KTYP20.EQ.2) THEN
        IF(KTYP21.EQ.1) THEN !Data fitting by optimisation
          IDEFLT(1)=2
          FORMAT='($,'' Enter the max. num. of iterations [2]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_ITS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_ITS=IDATA(1)

          RDEFLT(1)=1.0d-2
          FORMAT='($,'' Specify the data RMS error tolerance '','
     '      //'''[0.01]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=ITER_TOL(1)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,1.0d0,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITER_TOL(1)=RDATA(1)

          RDEFLT(1)=1.0d-2
          FORMAT='($,'' Specify the scale factor tolerance '','
     '      //'''[0.01]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=ITER_TOL(2)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,1.0d0,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITER_TOL(2)=RDATA(1)

          IF(KTYP12.GT.0) THEN
            RDEFLT(1)=1.0d0
            FORMAT='($,'' Specify the smoothing update factor '','
     '        //'''[1.0]: '',D11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=ITER_STORE(1)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,1000.0d0,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) ITER_STORE(1)=RDATA(1)
          ENDIF

          RDEFLT(1)=25.0d0
          FORMAT='($,'' Specify the maximum data error '','
     '      //'''[25.0]: '',D11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=ITER_STORE(2)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,1000.0d0,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITER_STORE(2)=RDATA(1)

          FORMAT='($,'' Use warm optimisations [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(WARMOPTI) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(2)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            WARMOPTI=(ADATA(1).EQ.'Y').OR.(ADATA(1).EQ.'y')
          ENDIF

          CALL STRING_TRIM(FILE00,IBEG,IEND)
          CDEFLT(1)=FILE00(IBEG:IEND)
          BASE_FNAME(1:15)=' '
          BASE_FNAME(15-(IEND-IBEG):15)=FILE00(IBEG:IEND)
          FORMAT='($,'' Enter the base filename ['//
     '      BASE_FNAME(1:15)//']: '',A)'
          IF(IOTYPE.EQ.3) CDATA(1)=IT_FNAME
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,15,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c         CALL CINOUT(IOTYPE,IOIP,IFILE,FORMAT,1,CDATA,CDEFLT,15,
c    '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IT_FNAME=CDATA(1)(1:15)

          IDEFLT(1)=0
          FORMAT='($,'' Enter the fitted nodal values output '','
     '      //'''frequency [0]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=OUTPUT_FREQ(1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,0,9999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) OUTPUT_FREQ(1)=IDATA(1)

          IDEFLT(1)=0
          FORMAT='($,'' Enter the data projection point output '','
     '      //'''frequency [0]: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=OUTPUT_FREQ(2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,0,9999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) OUTPUT_FREQ(2)=IDATA(1)

          IF(KTYP12.GT.0) THEN
            IDEFLT(1)=0
            FORMAT='($,'' Enter the smooting weight update '','
     '        //'''frequency [0]: '',I4)'
            IF(IOTYPE.EQ.3) IDATA(1)=OUTPUT_FREQ(3)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '        IDEFLT,0,9999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '        INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) OUTPUT_FREQ(3)=IDATA(1)
          ENDIF

        ELSE IF(KTYP21.EQ.2) THEN
        ENDIF
      ELSEIF(KTYP20.GE.3) THEN !extracellular solution iteration
        !MLB 17-Feb-1998

        !Set the convergence tolerance for potential
        RDEFLT(1)=0.1d0
        FORMAT='($,'' Specify the maximum potential '
     '    //'tolerance [0.1]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=ITER_TOL(1)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITER_TOL(1)=RDATA(1)

        !Set the convergence tolerance for flux
        RDEFLT(1)=0.1d0
        FORMAT='($,'' Specify the maximum flux '
     '    //'tolerance [0.1]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=ITER_TOL(2)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITER_TOL(2)=RDATA(1)

        !Specify a maximum number of iterations
        IDEFLT(1)=10
        FORMAT='($,'' Enter the maximum number of '
     '    //'iterations [10]: '',I4)'
        IF(IOTYPE.EQ.3) IDATA(1)=NUM_ITS
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_ITS=IDATA(1)

        !Specify upinit alpha value
        RDEFLT(1)=1.0d0
        FORMAT='($,'' Specify the update initial alpha '
     '   //'value [1.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=ITER_ALPHA1
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-0.1d0,1.1d0,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITER_ALPHA1=RDATA(1)

        !Specify solve iterate alpha value
        RDEFLT(1)=0.0d0
        FORMAT='($,'' Specify the solve iterate alpha '
     '    //'value [0.0]: '',D11.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=ITER_ALPHA2
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-0.1d0,1.1d0,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITER_ALPHA2=RDATA(1)
      ENDIF

      CALL EXITS('IPITER')
      RETURN
 9999 CALL ERRORS('IPITER',ERROR)
      CALL EXITS('IPITER')
      RETURN 1
      END


