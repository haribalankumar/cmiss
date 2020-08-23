      SUBROUTINE DETRSF(NELIST3,NHP,NKH,NPLIST3,NPLIST4,NPLIST5,NPNY,
     '  NVHP,NXLIST,NYNP,NYNR,YP,FIX,STRING,ERROR,*)
CC AJPe

C#### Subroutine: DETRSF
C###  Description:
C###    DETRSF defines parameters for transfer matrix from one
C###    surface to another surface.  The current application is for
C###    trasferring potentials between epicardial and body surface.

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
!     Parameter List
CC AJPs 11-11-97
      INTEGER NELIST3(0:NEM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPLIST3(0:NPM),NPLIST4(0:NPM),NPLIST5(0:NPM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
CC AJPe
      REAL*8 YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,nx,nxc
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE
      CHARACTER FILE*(MXCH),STATUS*3

      CALL ENTERS('DETRSF',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C#### Command: FEM define transfer;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    This command defines the two surfaces between which the
C###    transfer matrix is to be constructed.
C###  Parameter: <class #[1]>
C###    Specify the class number (of solve type) which the
C###    transfer matrix corresponds to.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<class #>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DETRSF',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS(' DLMPRW',NOCO,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        IF(FILIO) CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)
        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

C LKC 28-JUL-1999
        CALL ASSERT(USE_TRANSFER.EQ.1,'Set USE_TRANSFER to 1',ERROR,
     '    *9999)

        IF(FILIO) THEN
          IPFILE=2
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'trsf',
     '        STATUS,ERR,ERROR,*9999)
            CALL IPTRSF(NELIST3,NHP(1,0,nx),NKH(1,1,1,0),NPLIST3,
     '        NPLIST4,NPLIST5,NPNY,NVHP(1,1,1,0),NYNP,
     '        NYNR(0,0,1,0,nx),YP(1,1,nx),FIX(1,1,nx),ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF !filio
      ENDIF

      CALL_TRANSFER=.TRUE.
      CALL EXITS('DETRSF')
      RETURN
 9999 CALL ERRORS('DETRSF',ERROR)
      CALL EXITS('DETRSF')
      RETURN 1
      END


