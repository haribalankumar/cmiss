      SUBROUTINE DECUST(NXLIST,STRING,ERROR,*)
C  SMAR009 18/01/99 removed NRLIST,from list
C#### Subroutine: DECUST
C###  Description:
C###    DECUST defines parameters for customising a given mesh
C###    with prompted input or from filename.ipcust.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NXLIST(0:NXM)
C  SMAR009 18/01/99 removed NRLIST(0:NRM),
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,nxc !SMAR009 22/12/98 ,no_nrlist
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DECUST',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define customisation;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define mesh customisation parameters - mainly for use in customising
C###    a generic surface mesh to a given set of meaurements on length and
C###    circumference. Parameters can be read
C###    from or written to the file FILENAME.ipcust, with $current
C###    specifing the current default file.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to read the history
C###    data into.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECUST',ERROR,*9999)
      ELSE

        CALL PARSE_QUALIFIERS('DLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
C        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        nxc=nxc !routine not currently class dependent

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          IPFILE=2
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'cust',STATUS,
     '        ERR,ERROR,*9999)
C            DO no_nrlist=1,NRLIST(0)
C MLB 19 Feb 1997 nr not used
C             nr=NRLIST(no_nrlist)
C              CALL IPCUST(ERROR,*9999)
C            ENDDO !no_nrlist
C
            CALL IPCUST(ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF !filio
      ENDIF

      CALL EXITS('DECUST')
      RETURN
 9999 CALL ERRORS('DECUST',ERROR)
      CALL EXITS('DECUST')
      RETURN 1
      END


