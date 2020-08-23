      SUBROUTINE IPINVE(ERROR,*)

C#### Subroutine: IPINVE
C###  Description:
C###    IPINVE inputs regularisation options for inverting the transfer
C###    matrix.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPINVE',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C LKC 1-DEC-1999 Adding new question
      IDEFLT(1)=1
      FORMAT='('' Specify inverse approach [1]:'''//
     '    '/''   (1) Potential imaging'''//
     '    '/''   (2) Activation imaging'''//
     '    '/''   (3) Unused'''//
     '    '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) THEN
        IDATA(1)=INV_APPROACH
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) INV_APPROACH=IDATA(1)

C*** Potential Imaging approach
      IF(INV_APPROACH.EQ.1) THEN
        EVALUATE_INVERSE=.FALSE.

        IDEFLT(1)=1
        FORMAT='(/'' Specify whether [1]:'''//
     '    '/''   (1) Inverse transfer matrix calculated explicitly'''//
     '    '/''   (2) Inverse solution calculated from matrix eqtns'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=ICALC_TRANSFER
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ICALC_TRANSFER=IDATA(1)

CC JMB 07-FEB-2000
        IDEFLT(1)=4
        IF(ICALC_TRANSFER.EQ.1) THEN
          FORMAT='(/'' Enter the type of inversion stabilisation '//
     '      'scheme [4]:'''//
     '      '/''   (1) None'''//
     '      '/''   (2) Truncated SVD'''//
     '      '/''  *(3) Twomey Regularisation'''//
     '      '/''   (4) Zero-order Tikhonov'''//
     '      '/''  *(5) First-order Tikhonov'''//
     '      '/''  *(6) Second-order Tikhonov'''//
     '      '/''  *(7) Local regularisation'''//
     '      '/''   (8) Unused'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=IREGULARISE
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,8,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IREGULARISE=IDATA(1)
          IF(IREGULARISE.EQ.2) THEN
            RDEFLT(1)=1.0d-3
            FORMAT='($,'' Enter the cutoff-level for the'
     '        //' singular values: '',D11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=RMIN_SVD
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) RMIN_SVD=RDATA(1)
          ELSE IF(IREGULARISE.GE.4.AND.IREGULARISE.LE.6) THEN !Tikhonov
            RDEFLT(1)=1.0d-3
            FORMAT='($,'' Enter the Tikhonov regularisation parameter'
     '        //': '',D11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=TIKH_VALUE
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) TIKH_VALUE=RDATA(1)
          ENDIF
          IF(IREGULARISE.EQ.3.OR.IREGULARISE.GE.5) THEN
            ERROR='>>That option not implemented yet'
            GOTO 9999
          ENDIF

CC JMB 07-FEB-200 Adding new section
        ELSEIF(ICALC_TRANSFER.EQ.2) THEN
          IDEFLT(1)=3
          FORMAT='(/'' Enter the type of inversion stabilisation '
     '      //'scheme [3]:'''//
     '      '/''   (1) None'''//
     '      '/''   (2) Truncated SVD'''//
     '      '/''   (3) Tikhonov'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=ISTABILISE
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISTABILISE=IDATA(1)
          IF(ISTABILISE.GT.3) THEN
            ERROR='Unknown stabilisation scheme'
          ENDIF
          IF(ISTABILISE.GT.1) THEN
            IDEFLT(1)=1
            FORMAT='(/'' Specify any additional contraints [1]:'''//
     '      '/''   (1) None'''//
     '      '/''  *(2) Surface Gradient Operator'''//
     '      '/''  *(3) Surface Laplacian Operator'''//
     '      '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=ICONSTRAINT
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ICONSTRAINT=IDATA(1)
            IF(ICONSTRAINT.GT.3) THEN
              ERROR='Unknown additional constraints'
              GOTO 9999
            ENDIF
            IDEFLT(1)=1
            FORMAT='(/'' Specify any additional coupling [1]:'''//
     '      '/''   (1) None'''//
     '      '/''   (2) Greensite'''//
     '      '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=ICOUPLING
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) ICOUPLING=IDATA(1)
            IF(ICOUPLING.GT.2) THEN
              ERROR='Unknown additional coupling'
              GOTO 9999
            ENDIF
          ENDIF
          IF(ISTABILISE.EQ.2) THEN
            IDEFLT(1)=1
            FORMAT='(/'' Enter the regularisation scheme [1]:'''//
     '      '/''   (1) Generalised Cross Validation (GCV) Criterion'''//
     '      '/''   (2) L-curve Criterion'''//
     '      '/''   (3) Picard Criterion'''//
     '      '/''   (4) Quasi-optimality Criterion'''//
     '      '/''   (5) Optimal Criterion'''//
     '      '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=IREGULARISE
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) IREGULARISE=IDATA(1)
            IF(IREGULARISE.GT.5) THEN
              ERROR='Unknown regularisation scheme'
              GOTO 9999
            ENDIF
          ELSEIF(ISTABILISE.EQ.3) THEN
            IDEFLT(1)=3
            FORMAT='(/'' Enter the regularisation scheme [1]:'''//
     '      '/''   (1) Generalised Cross Validation (GCV) Criterion'''//
     '      '/''   (2) L-curve Criterion'''//
     '      '/''   (3) Picard Criterion'''//
     '      '/''   (4) Quasi-optimality Criterion'''//
     '      '/''   (5) Optimal Criterion'''//
     '      '/''   (6) CRESO Criterion'''//
     '      '/''   (7) Zero-crossing Criterion'''//
     '      '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=IREGULARISE
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) IREGULARISE=IDATA(1)
            IF(IREGULARISE.GT.7) THEN
             ERROR='Unknown regularisation scheme'
              GOTO 9999
            ENDIF
          ENDIF
        ENDIF

C*** Activation Imaging approach
      ELSEIF(INV_APPROACH.EQ.2) THEN
        IDEFLT(1)=1
        FORMAT='(/'' Specify the imaging approach [1]:'''//
     '    '/''   (1) Zero-crossing'''//
     '    '/''   (2) Unused'''//
     '    '/''   (3) Unused'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=ICALC_TRANSFER
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ICALC_TRANSFER=IDATA(1)
        IF(ICALC_TRANSFER.GE.2) THEN
          ERROR='Unknown imaging approach'
          GOTO 9999
        ENDIF

        IDEFLT(1)=1
        FORMAT='(/'' Enter the regularisation scheme [1]:'''//
     '    '/''   (1) None'''//
     '    '/''   (2) Surface Laplacian'''//
     '    '/''   (3) Unused'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=IREGULARISE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) IREGULARISE=IDATA(1)

        IF(IREGULARISE.EQ.2) THEN
          IDEFLT(1)=1
          FORMAT='(/'' Specify regularisation parameter choice [1]:'''//
     '      '/''   (1) Constant'''//
     '      '/''   (2) Ratio'''//
     '      '/''   (3) L-curve zero'''//
     '      '/''   (4) Maximum curvature'''//
     '      '/''   (5) Minimum product of norms'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=ICONSTRAINT
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ICONSTRAINT=IDATA(1)

          RDEFLT(1)=5.0d-1
          IF(ICONSTRAINT.EQ.1) THEN
            FORMAT='($,'' Enter the regularisation parameter'
     '        //' [0.5]: '',D11.4)'
          ELSE
            FORMAT='($,'' Enter the initial regularisation parameter'
     '        //' [0.5]: '',D11.4)'
          ENDIF
          IF(IOTYPE.EQ.3) RDATA(1)=REG_PARAM_LAPLACE
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) REG_PARAM_LAPLACE=RDATA(1)
        ENDIF !IREGULARISE = 2

      ELSE
        ERROR=' Unknown inverse approach'
        GOTO 9999
      ENDIF !INV_APPROACH

      CALL_INVE=.TRUE.

      CALL EXITS('IPINVE')
      RETURN
 9999 CALL ERRORS('IPINVE',ERROR)
      CALL EXITS('IPINVE')
      RETURN 1
      END


