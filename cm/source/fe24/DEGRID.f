      SUBROUTINE DEGRID(NAQ,NBJ,NEELEM,NELIST,NENP,NENQ,NLL,
     '  NLQ,NNB,NPNE,NQET,NQLIST,NQSCNB,NQXI,NRLIST,NWQ,NXI,NXQ,
     '  NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,NQNE_PTR,NQNLAT_PTR,NQS_PTR,
     '  STRING,DL,XIQ,ERROR,*)


C#### Subroutine: DEGRID
C###  Description:
C###    DEGRID defines finite difference grid.
C***  Martin Buist, June 1997

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER NAQ(NQM,NAM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NLL(12,NEM),NLQ(NQM),
     '  NNB(4,4,4,NBFM),NPNE(NNM,NBFM,NEM),NQET(NQSCM),NQLIST(0:NQM),
     '  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),
     '  NRLIST(0:NRM),NWQ(8,0:NQM,NAM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      INTEGER*4 NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,
     &  NQNE_PTR,NQNLAT_PTR,NQS_PTR
      REAL*8 DL(3,NLM),XIQ(NIM,NQM)
      CHARACTER ERROR*(*),STRING*(*)

!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,
     '  IEND1,IEND2,IPFILE,N3CO,NITB,noelem,nr
      CHARACTER FILE*(MXCH),STATUS*3,TYPE*8
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,GENER,MOUSE

      CALL ENTERS('DEGRID',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define grid;d/l/p/r/w;<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define grid point information. This includes defining the
C###    type of grid scheme to be used in each element. The material
C###    coordinates and the grid point connectivity are calculated.
C###    Grid point properties are read from or written to the file
C###    FILENAME (with extension .ipgrid) in the directory specified
C###    by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w;'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define grid;d/l/p/r/w;<;FILENAME[$current]><;(PATH/example)[$current]> gauss
C###  Description:
C###    Define grid point information. This includes defining the
C###    type of grid scheme to be used in each element. The material
C###    coordinates and the grid point connectivity are calculated.
C###    Grid point properties are read from or written to the file
C###    FILENAME (with extension .ipgrid) in the directory specified
C###    by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w;'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define grid;d/l/p/r/w;<;FILENAME[$current]><;(PATH/example)[$current]> coronary
C###  Description:
C###    Define grid point information for coronary trees. This works
C###    only for 1D meshes. Two points are defined at each bifurcation.
C###    Grid point properties are read from or written to the file
C###    FILENAME (with extension .ipgrid) in the directory specified
C###    by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w;'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> coronary'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define grid;d/l/p/r/w;<;FILENAME[$current]><;(PATH/example)[$current]> adaptive_level
C###  Description:
C###    Define adaptive level of interpolation for each grid point.
C###    Level na=1 is the highest resolution.
C###    Grid point levels are read from or written to the file
C###    FILENAME (with extension .ipgrlv) in the directory specified by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w;'
     '   //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '   //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> adaptive_level'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C NEWS JHC 22-NOV-2004 Added to define grid at gauss pts
C#### Command: FEM define grid;d/l/p/r/w;<;FILENAME[$current]><;(PATH/example)[$current]> adaptive_level
C###  Description:
C###    Define grid points at gauss points.
C###    FILENAME (with extension .ipgrlv) in the directory specified by PATH.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w;'
     '   //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '   //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> gauss'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C NEWE
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEGRID',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRID.EQ.1,'Must set USEGRID to 1',ERROR,*9999)
        IPFILE=1
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,CALCU,FILIO,
     '    GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL ASSERT(CALL_ELEM,'>>Must define elements first',
     '    ERROR,*9999)
        IF(CBBREV(CO,'ADAPTIVE_LEVEL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='ADAPTIVE'
        ELSE IF(CBBREV(CO,'CORONARY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='CORONARY'
        ELSE
          TYPE='DEFAULT'
        ENDIF

        IF(CBBREV(CO,'GAUSS',2,noco+1,NTCO,N3CO)) THEN
          KTYP3B=2
        ELSE
          KTYP3B=1
        ENDIF

        IF(FILIO) THEN
          CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

          IF(TYPE(1:7).EQ.'DEFAULT'.OR.TYPE(1:8).EQ.'CORONARY') THEN
            NEQM=NEELEM(0,0)
            DO nr=1,NRM !NRLIST(0)
C              nr=NRLIST(nrr)
              DO noelem=1,NEELEM(0,nr)
                IF(NEELEM(noelem,nr).GT.NEQM) NEQM=NEELEM(noelem,nr)
              ENDDO !ne
            ENDDO !nr


            IF(.NOT.CALL_GRID) THEN
              NQNE_PTR=0
              CALL ALLOCATE_MEMORY(NEQM*NQEM,0,INTTYPE,
     &          NQNE_PTR,MEM_INIT,ERROR,*9999)
              NQS_PTR=0
              CALL ALLOCATE_MEMORY(NEQM,0,INTTYPE,
     &          NQS_PTR,MEM_INIT,ERROR,*9999)
            ENDIF

            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'grid',STATUS,
     &        ERR,ERROR,*9999)
            CALL IPGRID(TYPE,NBJ,NEELEM,NELIST,NENP,NENQ,NLL,NNB,NPNE,
     &        NQET,%VAL(NQNE_PTR),%VAL(NQS_PTR),NQSCNB,NQXI,NRLIST,NWQ,
     &        NXI,NXQ,NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,NQNLAT_PTR,DL,
     &        XIQ,ERROR,*9999) 

            CALL_GRID=.TRUE.

            IF(NMGT.GT.1) THEN !more than one grid level

C LKC 2-JUN-1999 this is wrong ! see PJH
              NITB=NIT(NQSCNB(1)) !temporary
              nr=1                !temporary
              CALL ConstructNAQ(NITB,NAQ,NXQ,NWQ(1,1,3),ERROR,*9999)
              CALL ConstructNXQ(NAQ,nr,NWQ,NXQ,ERROR,*9999)
            ENDIF

          ELSE IF(TYPE(1:8).EQ.'ADAPTIVE') THEN
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'grlv',STATUS,
     '        ERR,ERROR,*9999)
            CALL IPGRLV(NLQ,NQLIST,ERROR,*9999)
          ENDIF !TYPE

          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF !FILIO
      ENDIF

      CALL EXITS('DEGRID')
      RETURN
 9999 CALL ERRORS('DEGRID',ERROR)
      CALL EXITS('DEGRID')
      RETURN 1
      END


