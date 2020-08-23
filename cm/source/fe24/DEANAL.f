      SUBROUTINE DEANAL(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  NBJ,NDIPOLES,NEELEM,NENP,NKJE,NKH,
     '  NPF,NP_INTERFACE,NPNE,NPNODE,NRE,NRLIST,NVHP,NVJE,NXLIST,
     '  NYNP,CE,DIPOLE_CEN,DIPOLE_DIR,SE,XA,XE,XP,YP,YQ,STRING,ERROR,*)

C#### Subroutine: DEANAL
C###  Description:
C###    DEANAL defines analytic solution for region nr and class nxc.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXLIST(0:NXM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,no_nrlist,nr,
     '  nx,nxc,N3CO
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ACTIVATION,ALL_REGIONS,CALCU,CBBREV,FILEIO,GENER,MOUSE


      CALL ENTERS('DEANAL',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define analytic;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines and calculates the analytic solution for region nr and
C###    class nxc. The analytic solution parameters are read from or
C###    written to the file FILENAME.ipanal, with the current default
C###    filename being $current.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The all value
C###    specifies all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to read the history
C###    data into.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #1>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C#### Command: FEM define analytic;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]> activation
C###  Description:
C###    Defines and calculates an analytic activation sequence for
C###    region nr and class nxc.
C###    The analytic activation parameters are read from or
C###    written to the file FILENAME.ipanal, with the current default
C###    filename being $current.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The all value
C###    specifies all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to read the history
C###    data into.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'activation'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #1>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEANAL',ERROR,*9999)
      ELSE

        CALL ASSERT(CALL_EQUA,'>>Equation not defined',ERROR,*9999)
        CALL ASSERT(CALL_MATE,'>>Material parameters not defined',
     '    ERROR,*9999)

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILEIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'ACTIVATION',1,noco+1,NTCO,N3CO)) THEN
          ACTIVATION=.TRUE.
        ELSE
          ACTIVATION=.FALSE.
        ENDIF


        IF(FILEIO) THEN
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          IPFILE=1 !is input file version number on 24-Mar-1995
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'anal',STATUS,
     '      ERR,ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)
            CALL IPANAL(DIPOLE_CEN_NTIME(1,1,nx),
     '        DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBJ,NDIPOLES(1,nx),
     '        NEELEM,NENP(1,0,nr),NKJE,
     '        NKH,NPF,NP_INTERFACE,NPNE,NPNODE,nr,
     '        NRE,NVHP(1,1,1,nr),NVJE,nx,NYNP,CE(1,1,nx),
     '        DIPOLE_CEN(1,0,1,1,nx),DIPOLE_DIR(1,0,1,1,nx),SE,XA,XE,
     '        XP,YP(1,1,nx),YQ(1,1,1,nx),ACTIVATION,ERROR,*9999)
          ENDDO !no_nrlist
          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF !fileio

        CALL_ANAL=.TRUE.

      ENDIF

      CALL EXITS('DEANAL')
      RETURN
 9999 CALL ERRORS('DEANAL',ERROR)
      CALL EXITS('DEANAL')
      RETURN 1
      END


