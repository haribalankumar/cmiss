      SUBROUTINE DESOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  NAN,NBH,NDIPOLES,NEELEM,NENQ,NRLIST,NQET,NQLIST,NQNE,NQS,NQSCNB,
     '  NQXI,NWQ,NXLIST,NXQ,AQ,CG,CQ,DIPOLE_CEN,DIPOLE_DIR,DXDXIQ,GD,PG,
     '  PROPQ,WG,XE,XG,XQ,YQ,STRING,ERROR,*)

C#### Subroutine: DESOUR
C###  Description:
C###    DESOUR defines sources. Currently defines dipole
C###    sources for EEG studies.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'sour00.cmn'

!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NENQ(0:8,NQM),NRLIST(0:NRM),NQET(NQSCM),NQLIST(0:NQM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),
     '  NWQ(8,0:NQM,NAM),NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CG(NMM,NGM),CQ(NMM,NQM,NXM),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),DXDXIQ(3,3,NQM),
     '  GD(NZ_GD_M),PG(NSM,NUM,NGM,NBM),PROPQ(3,3,4,2,NQM,NXM),
     '  WG(NGM,NBM),XE(NSM,NJM),XG(NJM,NUM),XQ(NJM,NQM),
     '  YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,IPFILE,
     '  no_nrlist,nr,nr_grid,nx,nxc,nx_dipole,nx_grid,N3CO
      REAL*8 RFROMC,SETTIME
      CHARACTER FILE*(MXCH),STATUS*3,TYPE*4
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,FIRST_TIME,GENER,MOUSE

      CALL ENTERS('DESOUR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define source;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines dipole source parameters. Parameters can be read from or
C###    written to the file FILENAME.ipsour, with $current specifying the
C###    current default file.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The all value
C###    specifies all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to read the history
C###    data into.
C###  Parameter:      <scale #[1.0]>
C###    Scale the direction of the dipole vector

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<scale      #[1.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define source;d/l/p/r/w/c<;FILENAME[$current]><;(PATH/example)[$current]> grid
C###  Description:
C###    Defines dipole source parameters. Parameters can be read from or
C###    written to the file FILENAME.ipsour, with $current specifing the
C###    current default file. Dipole sources are calculated from
C###    the gradient of the transmembrane potential which is obtained
C###    from a grid solution.
C###  Parameter:      <grregion   #[2]>
C###    Specify the region which contains the grid points.
C###  Parameter:      <grclass    #[2]>
C###    Specify the class of the transmembrane potential problem.
C###  Parameter:      <diregion   #[1]>
C###    Specify the region in which the dipoles are to be made.
C###  Parameter:      <diclass    #[1]>
C###    Specify the class number of the problem the dipoles will
C###    be used as sources for.
C###  Parameter:      <time       #[1.0]>
C###    Specify the time to store the dipole at if the dipole
C###    center or vector is time dependent.
C###  Parameter:      <scale      #[1.0]>
C###    Specify a scale factor to multiply all calculated dipoles by.
C###    This only works with the calculate qualifier.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w/c'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> grid'
        OP_STRING(2)=BLANK(1:15)//'<grregion   #[2]>'
        OP_STRING(3)=BLANK(1:15)//'<grclass    #[2]>'
        OP_STRING(4)=BLANK(1:15)//'<diregion   #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<diclass    #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<time       #[1.0]>'
        OP_STRING(7)=BLANK(1:15)//'<scale      #[1.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DESOUR',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 27-Jul-1994
        CALL PARSE_QUALIFIERS('DLMPRWC',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C LKC 29-APR-2003
C       If no nx then assign a temporary nx class - might want to
C       read in a source with no equations setup for simple maniputation
C       or exporting to CMGUI.
C
C        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
C     '    ERROR,*9999)
        IF(nx.EQ.0) THEN
          nx=1
          WRITE(OP_STRING(1),
     '      '(''>>Warning - No nx defined for this solve class'')') 
          WRITE(OP_STRING(2),
     '      '(''>>          Setting nx=1'')') 
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        

        
C LKC 2-FEB-1999
        CALL ASSERT(USE_DIPOLE.EQ.1,'>>Set USE_DIPOLE to 1',
     '    ERROR,*9999)

        IF(CBBREV(CO,'GRID',3,noco+1,NTCO,N3CO)) THEN
          TYPE='GRID'
        ELSE
          TYPE=' '
        ENDIF

        IF(CALCU) THEN
          IF(TYPE(1:4).NE.'GRID') THEN
            CALL ASSERT(.FALSE.,
     '        '>>Grid must be specified for calculate option',
     '        ERROR,*9999)
          ENDIF
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C LKC 29-APR-2003 propagate the scaling factor parameter
C       to standard dipoles which are not calculated by the grids,
C       ie those which are read in.
        IF(CBBREV(CO,'SCALE',3,noco+1,NTCO,N3CO)) THEN
          DIPOLE_SCALE_FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          DIPOLE_SCALE_FACTOR=1.0d0
        ENDIF
        
        IF(TYPE(1:4).EQ.'GRID') THEN
          IF(CBBREV(CO,'GRREGION',3,noco+1,NTCO,N3CO)) THEN
            nr_grid=IFROMC(CO(N3CO+1))
          ELSE
            nr_grid=2
          ENDIF
          IF(CBBREV(CO,'GRCLASS',3,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_grid,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_grid.NE.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)
          IF(CBBREV(CO,'DIREGION',3,noco+1,NTCO,N3CO)) THEN
            nr=IFROMC(CO(N3CO+1))
          ELSE
            nr=1
          ENDIF
          IF(CBBREV(CO,'DICLASS',3,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_dipole,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_dipole.NE.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)
          nx=nx_dipole
          IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
            SETTIME=RFROMC(CO(N3CO+1))
          ELSE
            SETTIME=1.0d0
          ENDIF

C LKC 29-APR-2003 propagate the scaling factor parameter
C       to standard dipoles which are not calculated by the grids,
C       ie those which are read in.
C         IF(CBBREV(CO,'SCALE',3,noco+1,NTCO,N3CO)) THEN
C           DIPOLE_SCALE_FACTOR=RFROMC(CO(N3CO+1))
C         ELSE
C           DIPOLE_SCALE_FACTOR=1.0d0
C         ENDIF
          
          NRLIST(0)=1
          NRLIST(1)=nr
        ELSE
          nr_grid=0
          nx_grid=1
        ENDIF

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'sour',STATUS,
     '        ERR,ERROR,*9999)

            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              IF(ALL_REGIONS) CALL PROMPT_REGION_ALL(nr,ERROR,*9999)

              CALL IPSOUR(DIPOLE_CEN_NTIME(1,1,nx),
     '          DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NAN,NBH,
     '          NDIPOLES(1,nx),NEELEM,NENQ,nr,nr_grid,NQET,NQLIST,NQNE,
     '          NQS,NQSCNB,NQXI,NWQ(1,0,1),nx_grid,NXQ,AQ,CG,
     '          CQ(1,1,nx_grid),DIPOLE_CEN(1,0,1,1,nx),
     '          DIPOLE_DIR(1,0,1,1,nx),DXDXIQ,PG,
     '          PROPQ(1,1,1,1,1,nx_grid),SETTIME,WG,XE,XG,XQ,
     '          YQ(1,1,1,nx_grid),CALCU,ERROR,*9999)
            ENDDO
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ELSE
          CALL IPSOUR(DIPOLE_CEN_NTIME(1,1,nx),DIPOLE_DIR_NTIME(1,1,nx),
     '      IBT,IDO,INP,NAN,NBH,NDIPOLES(1,nx),NEELEM,NENQ,nr,nr_grid,
     '      NQET,NQLIST,NQNE,NQS,NQSCNB,NQXI,NWQ(1,0,1),nx_grid,NXQ,AQ,
     '      CG,CQ(1,1,nx_grid),DIPOLE_CEN(1,0,1,1,nx),
     '      DIPOLE_DIR(1,0,1,1,nx),DXDXIQ,PG,PROPQ(1,1,1,1,1,nx_grid),
     '      SETTIME,WG,XE,XG,XQ,YQ(1,1,1,nx_grid),CALCU,ERROR,*9999)
        ENDIF
      ENDIF

      !Initialise GD here if dipoles are being used.
      DO i=1,NZ_GD_M
        GD(i)=0.0d0
      ENDDO

      CALL EXITS('DESOUR')
      RETURN
 9999 CALL ERRORS('DESOUR',ERROR)
      CALL EXITS('DESOUR')
      RETURN 1
      END


