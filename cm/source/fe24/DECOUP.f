      SUBROUTINE DECOUP(NEELEM,NP_INTERFACE,NXLIST,
     '  NXQ,STRING,XQ,ERROR,*)

C#### Subroutine: DECOUP
C###  Description:
C###    DECOUP defines coupling between dependent variables of
C###    different regions or problem classes.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 XQ(NJM,NQM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,IFROMC,
     '  IWK(6),no_nx,NTIW,nx,nxc,N3CO,TISSUE_REG,TREE_REG
      LOGICAL ALL_REGIONS,CALCU,FILIO,GENER,MOUSE
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL CBBREV

      CALL ENTERS('DECOUP',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define coupling;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define coupling between regions. Coupling properties are
C###    read from or written to the file FILENAME (with extension
C###    .ipcoup) in the directory specified by PATH.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define coupling;c
C###  Description:
C###    Defines the coupling between 1D trees and 2/3D regions.
C###    The command loops over the grid points in the tree and finds
C###    which points do not have a neighbour in the Xi 1 direction.
C###    These points are coupled to the nearest point on the 2/3D
C###    region.
C###  Parameter:      <tree_region #[2]>
C###    Specify the tree (1D) region to calculate the coupling with.
C###  Parameter:      <tissue_region #[1]>
C###    Specify the tissue (2/3D) region to calculate the coupling with.

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<tree_region #[2]>'
        OP_STRING(3)=BLANK(1:15)//'<tissue_region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DECOUP',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990
        CALL PARSE_QUALIFIERS('DLMPRWC',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO no_nx=1,NXLIST(0)
          nxc=NXLIST(no_nx)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ENDDO

        IF(CALCU) THEN
          IF(CBBREV(CO,'tree_region',2,noco+1,NTCO,N3CO)) THEN
            TREE_REG=IFROMC(CO(N3CO+1))
          ELSE
            TREE_REG=2
          ENDIF
          IF(CBBREV(CO,'tissue_region',2,noco+1,NTCO,N3CO)) THEN
            TISSUE_REG=IFROMC(CO(N3CO+1))
          ELSE
            TISSUE_REG=1
          ENDIF
          CALL IPCOUP(CALCU,NEELEM,NP_INTERFACE,nx,NXQ,
     '      TISSUE_REG,TREE_REG,XQ,ERROR,*9999)
        ELSE
          TISSUE_REG=0
          TREE_REG=0
          IF(FILIO) THEN
            ALL_REGIONS=.FALSE. !when parse_regions not called
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'coup',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPCOUP(CALCU,NEELEM,NP_INTERFACE,nx,NXQ,
     '        TISSUE_REG,TREE_REG,XQ,ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
          ENDIF
        ENDIF

      ENDIF

      CALL EXITS('DECOUP')
      RETURN
 9999 CALL ERRORS('DECOUP',ERROR)
      CALL EXITS('DECOUP')
      RETURN 1
      END


