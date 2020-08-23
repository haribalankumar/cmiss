      SUBROUTINE DEACTI(NBH,NBJ,NEELEM,NELIST,NGLIST,NRLIST,NXLIST,
     '  FEXT,STRING,ERROR,*)

C#### Subroutine: DEACTI
C###  Description:
C###    DEACTI defines parameters for active muscle contraction
C###    properties with prompted input or from filename.ipacti.
C###    Note: the mouse option stuff here was copied from DEOXS1 in FE17

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
      INTEGER NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NGLIST(0:NGM),NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 FEXT(NIFEXTM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,no_nrlist,nr,nxc,nx
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DEACTI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define active;l/m/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define active muscle properties. Active muscle properties are
C###    read from or written to the file FILENAME.ipacti in the directory
C###    specified by PATH with $current specifing the current default file.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/m/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEACTI',ERROR,*9999)
      ELSE

        CALL PARSE_QUALIFIERS('DLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL ASSERT(CALL_EQUA,'>>Equation not defined',ERROR,*9999)
        CALL ASSERT(CALL_MATE,'>>Material parameters not defined',
     '    ERROR,*9999)

C CPB 8/6/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(FILIO) THEN
          IPFILE=2
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'acti',STATUS,
     '        ERR,ERROR,*9999)
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL IPACTI(NBH,NBJ,NEELEM,NELIST,NGLIST,nr,nx,FEXT,
     '          ERROR,*9999)
            ENDDO !no_nrlist
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ENDIF !filio
        CALL_ACTI=.TRUE.
      ENDIF

      CALL EXITS('DEACTI')
      RETURN
 9999 CALL ERRORS('DEACTI',ERROR)
      CALL EXITS('DEACTI')
      RETURN 1
      END


