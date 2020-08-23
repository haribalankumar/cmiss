      SUBROUTINE IMGRID(NBJ,NEELEM,NENQ,NLQ,NQET,NQNE_PTR,NQS_PTR,
     '  NQSCNB,NQXI,NRLIST,NWQ,NXQ,NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,
     &  NQNLAT_PTR,AQ,DNUDXQ,XQ,STRING,ERROR,*)

C#### Subroutine: IMGRID
C###  Description:
C###    IMGRID imports a finite difference grid from a file containing
C###    grid point location and geometry
C***  Greg Sands, July 2003

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),
     &  NLQ(NQM),NQET(NQSCM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),
     &  NRLIST(0:NRM),NWQ(8,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      INTEGER*4 NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,NQNE_PTR,NQNLAT_PTR,
     &  NQS_PTR
      REAL*8 AQ(NMAQM,NQM),DNUDXQ(3,3,NQM),XQ(NJM,NQM)
      CHARACTER ERROR*(*),STRING*(*)

!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,ILIST(3),N3CO,ni,
     &  nij1,nij2,nq,NQTOTAL,nr
      CHARACTER FILE*(MXCH),FILENAME*(MXCH),LINE*(MXCH),TYPE*8
      LOGICAL ALL_REGIONS,CBBREV,COMPUTE_NXQ

      CALL ENTERS('IMGRID',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM import grid;FILENAME[$current]<;(PATH/example)[$current]>
C###  Description:
C###    Import grid point information, by default into the Lattice array 
C###    structures.
C###  Parameter:      <TrueGrid>
C###    Input file is in TrueGrid format.
C###  Parameter:      <nxq>
C###    Populate NXQ array.  This will only be accurate for a grid scheme
C###    that has consistent xis.
C###  Parameter:      <region (#)[1]>
C###    Specify the region number to use.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)
     &    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     &    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']> <TrueGrid>'
        OP_STRING(2)=BLANK(1:15)//'<nxq>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','IMGRID',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRID.EQ.1,'Must set USEGRID to 1',ERROR,*9999)
        
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)

        IF(CBBREV(CO,'TRUEGRID',2,noco+1,NTCO,N3CO)) THEN
          TYPE='TRUEGRID'
        ELSE
          TYPE='DEFAULT'
        ENDIF

        IF(CBBREV(CO,'NXQ',2,noco+1,NTCO,N3CO)) THEN
          COMPUTE_NXQ=.TRUE.
        ELSE
          COMPUTE_NXQ=.FALSE.
        ENDIF

        CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(TYPE(1:8).EQ.'TRUEGRID') THEN
          CALL STRING_TRIM(FILE,IBEG,IEND)
          FILENAME=FILE(IBEG:IEND)//'.tg'
          CALL OPENF(IFILE,'DISK',FILENAME(IBEG:IEND+3),
     &      'OLD','SEQUEN','FORMATTED',132,ERROR,*9999)

          READ(UNIT=IFILE,FMT='(A)') LINE
          CALL STRING_TRIM(LINE(1:50),IBEG,IEND)
          WRITE(OP_STRING,'('' Importing TrueGrid file: '',A)')
     &      LINE(IBEG:IEND)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

          READ(UNIT=IFILE,FMT=*) ILIST
          NQTOTAL=ILIST(2)
          WRITE(OP_STRING,'('' Number of grid points: '',I8)') NQTOTAL
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          CALL ASSERT(NQTOTAL.LE.NQM,'>>Increase NQM',ERROR,*9999)

          NEQM=ILIST(3)
          WRITE(OP_STRING,'('' Number of elements:    '',I8)') NEQM
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

          CALL ALLOCATE_LATTICE(NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,
     &      NQNLAT_PTR,ERROR,*9999)
          NQNE_PTR=0
          CALL ALLOCATE_MEMORY(NEQM*NQEM,0,INTTYPE,
     &      NQNE_PTR,MEM_INIT,ERROR,*9999)
          NQS_PTR=0
          CALL ALLOCATE_MEMORY(NEQM,0,INTTYPE,
     &      NQS_PTR,MEM_INIT,ERROR,*9999)

C***      Skip 3 lines (all zeros in sample file)
          READ(UNIT=IFILE,FMT='(///)')

          NQT=NQTOTAL
          NQR(1,nr)=1
          NQR(2,nr)=NQTOTAL
          IF(NEELEM(0,0).GE.1) THEN
            NQSCT=1
            CALL ASSERT(NQSCT.LE.NQSCM,'>>Increase NQSCM',ERROR,*9999)
            NQSCNB(1)=NBJ(1,1)
            NQXI(0,1)=3
            NQET(1)=1
            DO ni=1,NIM
              NQXI(ni,1)=2
              NQET(1)=NQET(1)*NQXI(ni,1)
            ENDDO !ni
            NMGT=1
          ENDIF

          CALL IMGRID_TRUEGRID(NENQ,%VAL(NLATNE_PTR),%VAL(NLATNQ_PTR),
     &      %VAL(NLATPNQ_PTR),NLQ,%VAL(NQNE_PTR),%VAL(NQNLAT_PTR),
     &      %VAL(NQS_PTR),NWQ,NXQ,AQ,XQ,COMPUTE_NXQ,ERROR,*9999)

          DO nq=1,NQT
            DO nij1=1,NJT
              DO nij2=1,NJT
                IF(nij1.EQ.nij2) THEN
                  DNUDXQ(nij1,nij2,nq)=1.0d0
                ELSE
                  DNUDXQ(nij1,nij2,nq)=0.0d0
                ENDIF !diag
              ENDDO !nij2
            ENDDO !nij1
          ENDDO

          CALL_GRID=.TRUE.
          USE_LAT=1
        ENDIF !TYPE
        CALL CLOSEF(IFILE,ERROR,*9999)
      ENDIF

      CALL EXITS('IMGRID')
      RETURN
 9999 CALL ERRORS('IMGRID',ERROR)
      CALL EXITS('IMGRID')
      RETURN 1
      END


