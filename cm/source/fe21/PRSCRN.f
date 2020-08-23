      SUBROUTINE PRSCRN(STRING,ERROR,*)

C#### Subroutine: PRSCRN
C###  Description:
C###    PRSCRN prints a GX graphics window to a file (postscript
C###    or ppm format).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER DEPTH,HEIGHT,IBEG,IBEG1,IEND,IEND1,IFROMC,N3CO,
     '  ORIENT,PS_TYPE,WIDTH
      CHARACTER FILE*(MXCH),TYPE*10
      LOGICAL ABBREV,CBBREV

      CALL ENTERS('PRSCRN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: print<;FILENAME>[file] postscript
C###  Parameter:   <colour/monochrome>[colour]
C###  Parameter:   <portrait/landscape>[portrait]
C###  Parameter:   <eps/ps>[eps]
C###  Description:
C###    Prints a GX graphics window to a postscript format file.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//'] postscript'
        OP_STRING(2)=BLANK(1:IEND)//'<colour/monochrome>[colour]'
        OP_STRING(3)=BLANK(1:IEND)//'<portrait/landscape>[portrait]'
        OP_STRING(4)=BLANK(1:IEND)//'<eps/ps>[eps]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: print<;FILENAME>[file] portable_pix_map
C###  Parameter:   height=<HEIGHT>[100]
C###  Parameter:   width=<WIDTH>[100]
C###  Description:
C###    Prints a GX graphics window to a ppm format file.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//'] portable_pix_map'
        OP_STRING(2)=BLANK(1:IEND)//'height=<HEIGHT>[100]'
        OP_STRING(3)=BLANK(1:IEND)//'width=<WIDTH>[100]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','PRSCRN',ERROR,*9999)
      ELSE
        TYPE='NONE'
        IF(NTCOQU(noco).EQ.0) THEN !no filename given so use default
          FILE=FILE00
        ELSE
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        ENDIF
        CALL STRING_TRIM(FILE,IBEG,IEND)
        IF(ABBREV(CO(noco+1),'POSTSCRIPT',3)) THEN
          TYPE='POSTSCRIPT'
          IF(CBBREV(CO,'LANDSCAPE',3,noco+1,NTCO,N3CO)) THEN
            ORIENT=2
          ELSE !Portrait
            ORIENT=1
          ENDIF
          IF(CBBREV(CO,'MONOCHROME',3,noco+1,NTCO,N3CO)) THEN
            DEPTH=2
          ELSE !Colour
            DEPTH=1
          ENDIF
C MPN 26Mar2002: I don't think this option does anything since
C                it doesn't appear to be implemented in GX???
          IF(CBBREV(CO,'PS',2,noco+1,NTCO,N3CO)) THEN
            PS_TYPE=1
C MPN 26Mar2002: use standard postscript extension
            FILE=FILE(IBEG:IEND)//'.ps'
C             FILE=FILE(IBEG:IEND)//'.psc'
          ELSE !EPS
            PS_TYPE=2
            FILE=FILE(IBEG:IEND)//'.eps'
         ENDIF
        ELSE IF(ABBREV(CO(noco+1),'PORTABLE_PIX_MAP',3)) THEN
          TYPE='PORTABLE'
          IF(CBBREV(CO,'HEIGHT',1,noco+1,NTCO,N3CO)) THEN
            HEIGHT=IFROMC(CO(N3CO+1))
          ELSE
            HEIGHT=100
          ENDIF
          IF(CBBREV(CO,'WIDTH',1,noco+1,NTCO,N3CO)) THEN
            WIDTH=IFROMC(CO(N3CO+1))
          ELSE
            WIDTH=100
          ENDIF
          FILE=FILE(IBEG:IEND)//'.ppm'
        ENDIF

        CALL PRINT_SCREEN_TO_FILE(DEPTH,HEIGHT,ORIENT,PS_TYPE,WIDTH,
     '    FILE,TYPE,ERROR,*9999)

      ENDIF

      CALL EXITS('PRSCRN')
      RETURN
 9999 CALL ERRORS('PRSCRN',ERROR)
      CALL EXITS('PRSCRN')
      RETURN 1
      END


