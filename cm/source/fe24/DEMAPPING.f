      SUBROUTINE DEMAPPING(NBJ,NEELEM,NKJ,NPNODE,NPNY,
     '  NRLIST,NVJP,NXLIST,NYNP,NYNY,CYNY,XP,ZP,
     '  STRING,ERROR,*)

C#### Subroutine: DEMAPPING
C###  Description:
C###    DEMAPPING defines ny->ny mapping

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),
     '  NRLIST(0:NRM),NVJP(NJM,NPM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNY(0:NYYM,NYM,NRM,NXM)
      REAL*8 CYNY(0:NYYM,NYM,NRM,NXM),
     '  XP(NKM,NVM,NJM,NPM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,no_nrlist,nr,nx,nxc,ny
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEMAPPING',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define mapping;l/p/r<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define a mapping of one mesh dof to another
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEMAPPING',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number
        CALL PARSE_QUALIFIERS('LPR',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        nx=NXLIST(1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C       CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C       CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(USE_MAPS.NE.0) THEN
          DO ny=1,NYM
            DO nr=1,NRM
              DO nx=1,NXM
                NYNY(0,ny,nr,nx)=0
              ENDDO
            ENDDO
          ENDDO
        ELSE
          CALL ASSERT(.FALSE.,'>>USE_MAPS must be'
     '      //' defined',ERROR,*9999)
        ENDIF

        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'map',
     '        STATUS,ERR,ERROR,*9999)

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)
              CALL IPMAPPING(NBJ,NEELEM,NKJ,NPNODE,NPNY,
     '          nr,NVJP,nxc,NYNP,NYNY,CYNY,XP,ZP,ERROR,*9999)
            ENDDO !no_nrlist
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF
      ENDIF

      CALL_MAPPING=.TRUE.
C*** 22/10/08 JHC Added a warning to re-define solve, 
C                 if the mapping is defined after defining solve.
      IF(CALL_SOLV) THEN
C       Need to 'define solve' again to set mappings.
        WRITE(OP_STRING,'('' >>Warning: Need to '
     '    //'define solve again.'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('DEMAPPING')
      RETURN
 9999 CALL ERRORS('DEMAPPING',ERROR)
      CALL EXITS('DEMAPPING')
      RETURN 1
      END


