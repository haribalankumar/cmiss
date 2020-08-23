C#### Module: CMISS_ARCHIVE3
C###  Description:
C###    Contains archived code from modules FE20 -> FE29

CFE21 Subroutine CHCLOC    changes clock position
C     Subroutine CHCOLO    changes colours
C     Subroutine CHLUT     changes colour lookup table
C     Subroutine NODE3D   
C     Subroutine: RECFID    receives fiducial marker information 
C###  Routine: REFINE_TRIANGLE dynamic memory wrapper for REFINE_TRIANGLE_DYNAM
C###  Routine: REFINE_TRIANGLE_DYNAM refine mesh created using Triangle
C###  Routine: SENDFID   sends fiducial marker information 
C     Subroutine UPGRID    update grid parameters

CFE22
C###  Routine: CAELEC   cancel electrodes
C###  Routine: HIELEC   hide electrodes
C###  Routine: SHELEC   show electrodes
C###  Routine: CAIMAG   cancel image
C###  Routine: HIIMAG   hide image
C###  Routine: SHIMAG   show image
C     Subroutine: SHPMAR   show polymarker
C###  Routine: CATEXT   cancel text
C###  Routine: HITEXT   hide text
C###  Routine: SHTEXT   show text
C###  Routine: CATRAC   cancel trace window
C     Subroutine CAUSER   cancel user defined names

CFE23 Subroutine DRAW       draws polyline on workstation
C     Subroutine PLOT       plots vectors
C     Subroutine SG_ROW_NUMBER  defines row number segments

CFE24
C###  Routine: DEIMAG   define image
C###  Routine: DEMACR   define macro command key
C     Subroutine DETEXT     define text
C     Subroutine DEVECT     define vector
C     Subroutine DEWIND     old VMS window options
C     Subroutine DEXI       define Xi coordinates (rewritten)
C     Subroutine DEXI       so we can remove large comment blocks
C     Subroutine DEGRID     define grd parameters
C###  Routine: DEPLANE  define image or tag plane   
C     Subroutine DEPOTE        define potential   
C     Subroutine DETIME     define time variables

CFE25
C###  Routine: SGELEC   create new electrode segment
C###  Routine: SGIMAG   create segment for image
      
CFE26.f
C     Subroutine DIBASE
C     Subroutine DITRAC
C     Subroutine EXGEOM


CFE27 Subroutine LIFIBR     list fibres
C###  Routine: LIMACR   list macro commands
C     Subroutine LIMATR     list matrix
C###  Routine: LITEXT   list text
C     Subroutine LITIME     list time variable information
C     SUBROUTINE LITRAN     lists Phigs transformation
C###  Routine: LIVSAE   list VSaero parameters

CFE28
C###  Routine: DRELEC   draw electrodes
C###  Routine: DRIMAG   draw image
C     Subroutine DRTEXT     draw text

CFE29 Subroutine UPGRID     upgrid at 30-3-1998 with Gregs stuff

Module FE20
=========== 

C LC 24/2/97 removed from :

C#### Subroutine: SYNTAX
C###  Description:
C###    Determines whether CO is a valid command. If CO is known the 
C###    command qualifiers are checked. When the qualifiers are 
C###    understood the command is executed.

!old MPN   7-Apr-95: unused?
!            IF(CO(noco+1).EQ.'?') THEN
!              CALL TRIM(STRING,IBEG,IEND)
!              IF(MAPOPN) THEN
!                WRITE(OP_STRING,*) STRING(1:IEND)//' off'
!              ELSE IF(.NOT.MAPOPN) THEN
!                WRITE(OP_STRING,*) STRING(1:IEND)//' on'
!     '            //' <imaging/data_fitting>[imaging]'
!              ENDIF
!              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
!            ELSE
!              IF(ABBREV(CO(noco+1),'ON',2)) THEN
!                IF(CBBREV(CO,'IMAGING',1,noco+1,NTCO,N3CO)) THEN
!                  MPLUN=MBOPN(INF,WT,TIMEOUT,LOAD,LFREE,LMSIZ,
!     '              MFREE,MMSIZ,'ESV1$DKB200:[USER.CMISS.FEM]FASTIMP.LO',
!     '              ' ',GETALT)
!                ELSE IF(CBBREV(CO,'DATA_FITTING',1,noco+1,NTCO,N3CO)) THEN
!                  MPLUN=MBOPN(INF,WT,TIMEOUT,LOAD,LFREE,LMSIZ,
!     '              MFREE,MMSIZ,'ESV1$DKB200:[USER.CMISS.FEM]FASTFEM.LO',
!     '              ' ',GETALT)
!                ENDIF
!                WRITE(OP_STRING,*) 'Bytes available = ',LFREE,MFREE
!      	        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!                WRITE(OP_STRING,*) 'Bytes in memory = ',LMSIZ,MMSIZ
!      	        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
!                MAPOPN=.TRUE.
C MLB 20-Nov-1996 these 3 logicals have been removed from cspi00.cmn
!                PUTPICS=.FALSE.
!                PUTMAPIO=.FALSE.
!                PUTMASKMAP=.FALSE.

!              ELSE IF(ABBREV(CO(noco+1),'OFF',2)) THEN
!                MAPOPN=.FALSE.
!              ENDIF
!            ENDIF

C KAT 7Nov00

      ELSE IF(ABBREV(CO(noco),'WRITE',1)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
C#### Command: write
C###  Parameter:  Buffer
C###  Parameter:  RGB 
C###  Description:
C###    
          OPTION_1( 1)='Buffer'
          OPTION_1( 2)='RGB'
          NTCH1=2
          CALL LIST_COMMANDS(2,NTCH1,OPTION_1,ERROR,*9999)
        ELSE
          IF(ABBREV(CO(noco),'BUFFER',1)) THEN
            IUNIT=IOFILE2
            CALL OPENF(IUNIT,'DISK','buffer.com','NEW','SEQUEN',
     '        'FORMATTED',132,ERROR,*9999)
            DO no_buffer=BUFFER_COUNT,2,-1
              WRITE(IUNIT,'(A)') BUFFER(no_buffer)
     '          (4:BUFFER_LENGTH(no_buffer))
            ENDDO
            CALL CLOSEF(IUNIT,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco),'RGB',1)) THEN
            write(*,'('' not yet implemented'')')
          ENDIF
        ENDIF


Module FE21
=========== 

C 22/2/97 LC removed section from : 
C#### Subroutine: FIT
C###  Description:
C###    FIT fits geometry or field parameters to data points defined 
C###    in ZD or YG.

C cpb 20/10/95 Merging all the different fit routines
C        IF(TYPE(1:6).EQ.'FIBRES'.OR.  
C     '    (TYPE(1:8).EQ.'GEOMETRY'.AND.ITYP6(nr,nx).EQ.1).OR.
C     '    TYPE(1:5).EQ.'FIELD'.OR.
C     '    TYPE(1:8).EQ.'POTENTIAL') THEN !fit fibre/sheet,geometry,field
CC                                         or potential
C          CALL FITFLD(IBT,IDO,INP,IPIVOT,ISC_GK,ISC_GKK,ISC2_GKK,
C     '      ISC_GQ,IWK_ITERATIVE(OS_ISC_PRECON_ITERATIVE),ISR_GK,
C     '      ISR_GKK,ISR2_GKK,ISR_GQ,
C     '      IWK_ITERATIVE(OS_ISR_PRECON_ITERATIVE),IWK1,IWK2,
C     '      IWK_ITERATIVE(OS_IWK_PRECON_ITERATIVE),LGE,LN,NBH,NBJ,NDDL,
C     '      NDLT,NJE,NKE,NKH,NONY(0,1,1,0,nx),NPF,NPNE,NPNODE,
C     '      NPNY(0,1,0,nx),NQE,nr,NRE,NVHE,NVJE,NVHP(1,1,1,nr),
C     '      NW,nx,NYNE,
C     '      NYNO(0,1,1,0,nx),NYNP,NYNR(0,0,1,0,nx),CONY(0,1,1,0,nx),
C     '      CYNO(0,1,1,0,nx),EDD,ER,ES,GK,GKK,GR,GRR,GQ,PG,
C     '      WK_ITERATIVE(OS_PRECON_ITERATIVE),RG,SE,VE,WD,WDL,WG,
C     '      WK1,WK_ITERATIVE(OS_WK1_ITERATIVE),
C     '      WK_ITERATIVE(OS_WK2_ITERATIVE),
C     '      WK_ITERATIVE(OS_WK_PRECON_ITERATIVE),WU,XA,XE,XG,XID,XIDL,
C     '      XO(1,nx),XP,ZA,ZD,ZDL,ZE,ZP,FIRSTA,FIX(1,1,nx),ERROR,*9999)
C
C        ELSE IF(TYPE(1:8).EQ.'GEOMETRY'.AND.ITYP6(nr,nx).EQ.2) THEN
C
C          WRITE(OP_STRING,'('' *** FITGEO Moved to archive ***'')')
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC cpb 7/10/94 Moved to archive
CC          CALL FITGEO(IBT,IDO,INP,IWK1,LD,LDR,LGE,LN,NXI,
CC     '      NBJ,NDDL,NDLT,NFF,NJE,NJP,NKE,NKJ,NLL,NNL,
CC     '      NONY(0,1,1,nr,nx),
CC     '      NPF,NPL,NPNE,NPO,NQE,nr,nx,CONY(0,1,1,nr,nx),DL,
CC     '      EDD,ER,ES,SCALE,SE,SQ,VE,
CC     '      WD,WDL,XA,XE,XID,XIDL,XP,YP(1,1,nx),ZD,ZDL,
CC     '      GKK,GRR,WK1,WK2,FIX(1,1,nx),ERROR,*9999)
C
C        ELSE IF(TYPE(1:7).EQ.'FOURIER'.AND.ITYP6(nr,nx).EQ.1) THEN
C
Cc cpb 26/12/94 to archive
CC          CALL FITFOU(IBT,ISC_GKK,ISC2_GKK,IDO,INP,ISR_GKK,ISR2_GKK,
CC     '      IWK1,IWK2,LN,NBJ,NDDL,NDLT,NJE,NKE,NKJ,NONY(0,1,1,nr,nx),
CC     '      NPF,NPNE,NPO,NQE,nr,NRE,NVJE,NVJP,nx,
CC     '      NYNO(0,1,1,nr,nx),CONY(0,1,1,nr,nx),
CC     '      ER,ES,GK,PG,SE,WD,WDL,WG,WU,XA,XE,XID,XIDL,XP,
CC     '      YP(1,1,nx),ZD,ZDL,GKK,GRR,WK1,XO(1,nx),ERROR,*9999)
C
C        ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
C          CALL FITGAU(IBT,IDISP,IDO,INP,ISC_GKK,ISC2_GKK,ISR_GKK,
C     '      ISR2_GKK,IWK1,IWK2,LGE,LN,NBH,NBJ,NKE,NKH,NKJ,
C     '      NONY(0,1,1,0,nx),NPF,NPNE,NPNODE,NPNY(0,1,0,nx),NQE,nr,
C     '      NRE,NVHE,NVJE,NVHP(1,1,1,nr),NVJP,nx,NYNE,
C     '      NYNO(0,1,1,0,nx),NYNP,
C     '      NYNR(0,0,1,0,nx),CONY(0,1,1,0,nx),CYNO(0,1,1,0,nx),ER,ES,
C     '      GK,GKK,GR,GRR,PG,SE,WG,WK1,WU,XA,XE,XIG,XO(1,nx),XP,YG,
C     '      FIX(1,1,nx),ERROR,*9999)
C
C        ENDIF

C 25/2/97 LC removed section from : 
C              AJP 11-5-94.  Altered to read in fitted field files.

C#### Subroutine: ITERATEA
C###  Description:
C###    ITERATEA handles all the iteration calls to FEM.

C AJP 11-5-94.  Altered to read in fitted field files.
cC***    read the field values from the electrode dataset
c
c              IF(INFORMAT_CODE.EQ.1) THEN !ipfiles
c                FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
c     '            CHAR1(IBEG2:IEND2)
c                CO(7)='IPFILE'
c              ELSE IF(INFORMAT_CODE.EQ.2) THEN !map3d
c                FILENAME=IT_FNAME(IBEG1:IEND1)
c                CO(7)='MAP3D'              
c              ENDIF
c              CO(1)='FEM'
c              CO(2)='DEFINE'
c              CO(3)='DATA'
c              COQU(3,1)='r'
c              COQU(3,2)=FILENAME
c              CO(4)='ELECTRODE'
c              CO(5)='POTENTIAL'
c              CO(6)='FORMAT'
c              CO(8)='DATASET'
c              CO(9)=CHAR1
c              NTCO=9
c              NTCOQU(3)=2
c              noco=1
c              IF(ECHO) THEN
c                CALL TRIM(FILENAME,IBEG3,IEND3)
c                CALL TRIM(CO(7),IBEG4,IEND4)
c                OP_STRING(1)=' ### FEM DEFINE DATA;r;'//
c     '            FILENAME(IBEG3:IEND3)//' ELECTRODE POTENTIAL'//
c     '            ' FORMAT '//CO(7)(IBEG4:IEND4)//' DATASET '//
c     '            CHAR1(IBEG2:IEND2)
c                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c              ENDIF
c              CALL FEM(ISEG,END,STRING,INTWORK,REALWORK,
c     '          ERROR,*9999)       
c              NTCOQU(3)=0
C              COQU(3,1)=' '
C              COQU(3,2)=' '
c
cC***    fit a field to the electrode data
c
c              CO(1)='FEM'
c              CO(2)='FIT'
c              CO(3)='POTENTIAL'
c              NTCO=3
c              noco=1
c              IF(ECHO) THEN
c                OP_STRING(1)=' ### FEM FIT POTENTIAL'
c                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c              ENDIF
c              CALL FEM(ISEG,cseg,END,STRING,INTWORK,REALWORK,
c     '          ERROR,*9999)       
C END Part of altereation to read in fitted data

C 25/2/97 LC removed section from :

C OLD??? AJP 1-6-94
C          CALL ASSERT(KTYP8.GT.0,'>>Define fit first',ERROR,*9999)
C          CALL TRIM(IT_FNAME,IBEG1,IEND1)
C
C          UPDATE_MATRIX(1)=.TRUE.
C          DO it_count=1,NUM_ITS
C            IF(ECHO) THEN
C              WRITE(OP_STRING,'('' ### Begining iteration '',
C     '          ''number : '',I4)') it_count
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C            WRITE(CHAR1(1:4),'(I4)') it_count
C            CALL TRIM(CHAR1,IBEG2,IEND2)
C
CC***    read the field values from the electrode dataset
C
C            IF(INFORMAT_CODE.EQ.1) THEN !ipfiles
C              FILENAME=IT_FNAME(IBEG1:IEND1)//'_'//
C     '          CHAR1(IBEG2:IEND2)
C              CO(7)='IPFILE'
C            ELSE IF(INFORMAT_CODE.EQ.2) THEN !map3d
C              FILENAME=IT_FNAME(IBEG1:IEND1)
C              CO(7)='MAP3D'              
C            ELSE IF(INFORMAT_CODE.EQ.3) THEN !emap
C            ENDIF
C            CO(1)='FEM'
C            CO(2)='DEFINE'
C            CO(3)='DATA'
C            COQU(3,1)='r'
C            COQU(3,2)=FILENAME
C            CO(4)='ELECTRODE'
C            CO(5)='POTENTIAL'
C            CO(6)='FORMAT'
C            CO(8)='DATASET'
C            CO(9)=CHAR1
C            NTCO=9
C            NTCOQU(3)=2
C            noco=1
C            IF(ECHO) THEN
C              CALL TRIM(FILENAME,IBEG3,IEND3)
C              CALL TRIM(CO(7),IBEG4,IEND4)
C              OP_STRING(1)=' ### FEM DEFINE DATA;r;'//
C     '          FILENAME(IBEG3:IEND3)//' ELECTRODE POTENTIAL'//
C     '          ' FORMAT '//CO(7)(IBEG4:IEND4)//' DATASET '//
C     '          CHAR1(IBEG2:IEND2)
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
C     '        ERROR,*9999)       
C            NTCOQU(3)=0
C            COQU(3,1)=' '
C            COQU(3,2)=' '
C
CC***    fit a field to the electrode data
C
C            CO(1)='FEM'
C            CO(2)='FIT'
C            CO(3)='POTENTIAL'
C            NTCO=3
C            noco=1
C            IF(ECHO) THEN
C              OP_STRING(1)=' ### FEM FIT POTENTIAL'
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
C     '        ERROR,*9999)       
C
CC***    update potentials from the field            
C
C            CO(1)='FEM'
C            CO(2)='UPDATE'
C            CO(3)='POTENTIAL'
C            CO(4)='FROM'
C            CO(5)='FIELD'
C            CO(6)='DATASET'
C            CO(7)=CHAR1(IBEG2:IEND2)
C            NTCO=7
C            noco=1
C            IF(ECHO) THEN
C              OP_STRING(1)=' ### FEM UPDATE POTENTIAL FROM '//
C     '          'FIELD DATASET '//CHAR1(IBEG2:IEND2)
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C            CALL FEM(ISEG,CSEG,END,STRING,INTWORK,REALWORK,
C     '        ERROR,*9999)       
C
C            UPDATE_MATRIX(1)=.FALSE.
C          ENDDO
C
C        ENDIF
C END OLD ???

C 2001-12-13 KAT removed from
C#### Subroutine: READF
C###  Description:
C###    READF reads .COM, HISTORY, .IOD, MATRIX, .XP, or .SIGNAL files.

C#### Command: FEM read signal;<FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:   <time #[0.0]>
C###  Description:
C###    Not used - should be archived.

        OP_STRING(1)=STRING(1:IEND)
     '    //' signal;<FILENAME['//FILE00(IBEG1:IEND1)//']>' 
        OP_STRING(2)=BLANK(1:15)//'<time #[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

        ELSE IF(ABBREV(CO(noco+1),'SIGNAL',2)) THEN
          TYPE='SIGNAL'


        ELSE IF(TYPE(1:6).EQ.'SIGNAL') THEN

          IF(CBBREV(CO,'TIME',1,noco+1,NTCO,N3CO)) THEN
            ITIME=IFROMC(CO(N3CO+1))
          ELSE
            ITIME=1
          ENDIF

         IUNIT=8
C          FIDUCALC=.FALSE.
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.header','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          READ(IUNIT,'(I4)') NUMSIGNALS
          READ(IUNIT,'(I5)') NUMSAMPLES
          READ(IUNIT,'(I5)') NUMSAMPLESPERSECOND
          DO I=1,NUMSIGNALS
            ACCEPT(I)=.TRUE.
          ENDDO
          DO I=1,NUMSIGNALS
            READ(IUNIT,'(18X,A5)') SIGHEADER(1,I)
            READ(IUNIT,'(9X,A5)') SIGHEADER(2,I)
            READ(IUNIT,'(A)') CHARS
            IF(CHARS(1:9).EQ.' REJECTED') THEN
              ACCEPT(I)=.FALSE.
            ENDIF
            READ(IUNIT,'(9X,E12.5)') TEMPVAR
            READ(IUNIT,'(9X,E12.5)') TEMPVAR
            READ(IUNIT,*)
          ENDDO
          CALL CLOSEF(IUNIT,ERROR,*9999)

          OPEN(UNIT=IUNIT,FILE=FILE(IBEG:IEND)//'.SIGNAL',STATUS='OLD',
     '      FORM='UNFORMATTED',IOSTAT=IOS)
          IF(ios.ne.0)THEN
            WRITE(CHAR1,'(I3)')IOS
            ERROR=' File open error: iostat='//CHAR1(1:3)
            GOTO 9999
          ENDIF
          CALL SIGNAL_READ(I2P(1,1,ITIME),IUNIT,NUMSIGNALS,NUMSAMPLES)
          CALL CLOSEF(IUNIT,ERROR,*9999)
      
      ELSE IF(TYPE(1:6).EQ.'SIGNAL') THEN
        CLOSE(UNIT=8)

      

C 25/2/97 LC removed section from : Not used (i hope). GBS 27-OCT-1994

C#### Subroutine: WRITEF
C###  Description:
C###    WRITEF writes .DAT, .IOD, MATRIX, .OPG, .OPS, .TRACE or 
C###    .VSAERO file.  
C###    CPB 14/2/93 Writes map3d .data files

C!!!! Not used (i hope). GBS 27-OCT-1994
C          IUNIT=IOFILE1
C          IOFI=IOFILE1
C          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.opg','NEW',
C     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
C          WRITE(IUNIT,'(/'' General Parameters''/1X,18(''=''))')
C          CALL OPCOOR(ERROR,*9999)
C          WRITE(IUNIT,'(//'' Basis Functions''/,1X,15(''=''))')
C          CALL OPBASE(IBT,IDO,INP,NAN,0,NGAP,PG,XIG,ERROR,*9999)
C          WRITE(IUNIT,'(//,'' Element Parameters''/,1X,18(''=''))')
C          DO nr=1,NRT
C            CALL OPELEM(NBJ,NBH,NCO,NFF,NGAP,NJE,
C     '        NJP,NKE,NKJ,NLF,NLL,NNF,NNL,NPF,NPL,NPNE,NQE,nr,NW,
C     '        DF,DL,PG,RG,SE,VE,WG,XA,XE,XG,XP,ERROR,*9999)
C          ENDDO
C          WRITE(IUNIT,'(//,'' Global Parameters''/1X,17(''=''))')
C          DO nr=1,NRT
C            CALL OPNODE(NBJ,NHP(1,nr,nx),NJE,NJP,NKH(1,1,1,nr),NKJ,NPNODE,nr,
C     '        NVJE,NVJP,.FALSE.,' ',XA,XP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C          ENDDO
C          CALL OPLINE(NPL,DL,ERROR,*9999)
C          CALL OPFACE(NLF,NPF,NPNE,DF,SE,ERROR,*9999)
C          CALL CLOSEF(IUNIT,ERROR,*9999)
C          IOFI=IOOP


C 25/2/97 LC removed section from : Not used. PJH and AJP  26-4-93

C#### Subroutine: WRITEF
C###  Description:
C###    WRITEF writes .DAT, .IOD, MATRIX, .OPG, .OPS, .TRACE or 
C###    .VSAERO file.  
C###    CPB 14/2/93 Writes map3d .data files

C Not used. PJH and AJP  26-4-93
C          IUNIT=IOFILE1
C          IOFI=IOFILE1
C          nr=1 !modify later
C          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.ops','NEW',
C     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
C          CALL OPEQUA(NHE(1,nx),NPB,nr,NW,ERROR,*9999)
C          WRITE(IUNIT,'(//'' Equation parameters''/1X,19(''=''))')
C          IF(ITYP1(nr,nx).EQ.3) THEN
C            CALL OPMAT3(NJE,NW,CE,CP,ERROR,*9999)
C          ELSE IF(ITYP1(nr,nx).EQ.4) THEN
C            CALL OPMAT4(NJE,NW,CE,CP,ERROR,*9999)
C          ELSE IF(ITYP1(nr,nx).EQ.5) THEN
C            CALL OPMAT5(NJE,NW,CE,CP,ERROR,*9999)
C          ELSE IF(ITYP1(nr,nx).EQ.9) THEN
C            CALL OPMAT3(NJE,NW,CE,CP,ERROR,*9999)
C          ENDIF
C          WRITE(IUNIT,'(//'' Boundary Conditions''/,1X,19(''=''))')
C          IF(ITYP1(nr,nx).EQ.3) THEN
C            CALL OPINI3(ITHRES,NBH,NBJ,NHE(1,nx),NHP(1,nr,nx),NJE,
C     '        NKH(1,1,1,nr),NPB,NVJE,NVJP,NW,
C     '        NYNP,FIX(1,1,nx),PE,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C          ELSE IF(ITYP1(nr,nx).EQ.4) THEN
C            CALL OPINI4(NBH,NHE(1,nx),NHP(1,nr,nx),NJE,NKH(1,1,1,nr),NPB,NW,
C     '        FIX(1,1,nx),PE,
C     '        YP(1,1,nx),ZA,ZP,ERROR,*9999)
C          ELSE IF(ITYP1(nr,nx).EQ.5) THEN
C            CALL OPINI5(NBH,NHE(1,nx),NHP(1,nr,nx),NJE,NKH(1,1,1,nr),NPB,NW,
C     '        FIX(1,1,nx),PE,
C     '        YP(1,1,nx),ZA,ZP,ERROR,*9999)
C          ELSE IF(ITYP1(nr,nx).EQ.9) THEN
C            CALL OPINI9(NBH,NBJ,NHE(1,nx),NHP(1,nr,nx),NJE,NKH(1,1,1,nr),NPB,
C     '        NPNE,nr,NVJE,NVJP,NW,
C     '        NYNP,FIX(1,1,nx),PE,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C          ENDIF
C
C          IF(ITYP2(nr,nx).EQ.1) THEN
C            WRITE(IUNIT,'(//'' Nodal values''/,1X,12(''=''))')
C            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C            CALL OPNODE(NHP(1,nr,nx),NJP,NKH(1,1,1,nr),NKJ,NPNODE,nr,NVJE,NVJP,
C     '        .TRUE.,'Solution',XA,XP,YP(1,1,nx),ZP,ERROR,*9999)
C            WRITE(IUNIT,'(//'' Nodal reactions''/,1X,15(''=''))')
C            CALL YPZP(5,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C            CALL OPNODE(NHP(1,nr,nx),NJP,NKH(1,1,1,nr),NKJ,NPNODE,nr,NVJE,NVJP,
C     '        .TRUE.,'Reaction',XA,XP,YP(1,1,nx),ZP,ERROR,*9999)
C
C          ELSE IF(ITYP2(nr,nx).EQ.14.OR.ITYP2(nr,nx).EQ.15) THEN
C            WRITE(IUNIT,'(//'' Nodal positions''/,1X,15(''=''))')
C            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C            CALL OPNODE(NHP(1,nr,nx),NJP,NKH(1,1,1,nr),NKJ,NPNODE,nr,NVJE,NVJP,
C     '        .TRUE.,'Solution',XA,XP,YP(1,1,nx),ZP,ERROR,*9999)
C            WRITE(IUNIT,'(//'' Nodal displacements''/,1X,
C     '        19(''=''))')
C            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C            DO nonode=1,NPNODE(0,nr)
C              np=NPNODE(nonode,nr)
C              DO nj=1,NJT
C                DO nv=1,NVJP(nj,np)
C                  DO nk=1,NKJ(nj,np)
C                    ZP(nk,nv,nj,np,nc)=ZP(nk,nv,nj,np,nc)
C     '                -XP(nk,nv,nj,np)
C                  ENDDO               
C                ENDDO
C              ENDDO
C            ENDDO
C            CALL OPNODE(NHP(1,nr,nx),NJP,NKH(1,1,1,nr),NKJ,NPNODE,nr,NVJE,NVJP,
C     '        .TRUE.,'Displacement',XA,XP,YP(1,1,nx),ZP,ERROR,*9999)
C            WRITE(IUNIT,'(//'' Nodal reactions''/,1X,15(''=''))')
C            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C            CALL OPNODE(NHP(1,nr,nx),NJP,NKH(1,1,1,nr),NKJ,NPNODE,nr,NVJE,NVJP,
C     '        .TRUE.,'Reaction',XA,XP,YP(1,1,nx),ZP,ERROR,*9999)
C            WRITE(IUNIT,'(//'' Element stresses''/,1X,16(''=''))')
C            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,
C     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
C            ENERGY_TOTAL=0.0
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
C     '            NPNE(1,1,ne),NQE(1,1,ne),nr,NVJE(1,1,1,ne),
C     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
C                CALL ZPZE(NBH(1,1,ne),nc,
C     '            NHE(ne,nx),NKE(1,1,1,ne),NPF(1,1),
C     '            NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),
C     '            SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
C                IF(ITYP1(nr,nx).EQ.4) THEN
C                  CALL CPCG(NW(ne,1),NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
C     '              CE(1,ne),CG,CP,PG,ERROR,*9999)
C                  CALL OPST40(NBH(1,1,ne),NBJ(1,ne),ne,
C     '              NHE(ne,nx),NJE(ne),NW(ne,1),nx,
C     '              CG,ENERGY,PG,TYPE,VE(1,1,ne),XE,XG,
C     '              YG(1,1,ne),ZE,ZG,ERROR,*9999)
C                  ENERGY_TOTAL=ENERGY_TOTAL+ENERGY
C                ELSE IF(ITYP1(nr,nx).EQ.5) THEN
C                  CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
C     '              CE(1,ne),CG,CP,PG,ERROR,*9999)
C                  nb=NBH(NH_LOC(1,nx),1,ne)
C                  CALL OPC50(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),ne,
C     '              NHE(ne,nx),NJE(ne),NKE(1,1,1,ne),NPNE(1,1,ne),
C     '              nr,nx,NW(ne,1),CE(1,ne),CG,CP,FEXT(1,1,ne),
C     '              PG,RG,VE(1,1,ne),
C     '              XA,XE,XF,XG,XIG(1,1,nb),ZA(1,1,1,ne),ZE,ZF,ZG,
C     '              ERROR,*9999)
C                ENDIF
C              ENDDO
C            WRITE(IUNIT,'('' Total energy = '',D12.4)') ENERGY_TOTAL
C          ENDIF
C          CALL CLOSEF(IUNIT,ERROR,*9999)
C          IOFI=IOOP
C

C 3 October 2000:
        ELSE IF(ABBREV(CO(noco+1),'VSAERO',2)) THEN
          TYPE='VSAERO'

        ELSE IF(TYPE(1:6).EQ.'VSAERO') THEN
          IUNIT=IOFILE2
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.vsaero','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          CALL IOVSAE(IBT,IDO,INP,NBJ,NKE,NPF,NPNE,
     '      NRE,NVJE,
     '      SE,XA,XE,XP,ERROR,*9999)
          CALL CLOSEF(IUNIT,ERROR,*9999)

      SUBROUTINE CHCLOC(ISCLOC,ISEG,noco,NTCO,CO,CSEG,STRING,ERROR,*)

C#### Subroutine: CHCLOC
C###  Description:
C###    Change clock position with mouse.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cloc00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISCLOC(*),ISEG(*),noco,NTCO
      CHARACTER CO(*)*(*),CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INSTAT,iw,N3CO
      REAL*8 XWC,YWC
      LOGICAL CBBREV

      CALL ENTERS('CHCLOC',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)//';m <on WS>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHCLOC',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          iw=IFROMC(CO(N3CO+1))
        ELSE
          iw=1
        ENDIF
        CALL ACWK(iw,0,ERROR,*9999)
        INSTAT=1
        DO WHILE(INSTAT.EQ.1)
          WRITE(OP_STRING,'('' >>Relocate clock on '',I1,'':'')') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.d0,XWC,0.d0,YWC,
     '      ERROR,*9999)
          IF(INSTAT.EQ.1) THEN
            X_CLOCK(1)=XWC
            X_CLOCK(2)=YWC
          ENDIF
          CALL SGCLOC(0,ISCLOC(iw),ISEG,iw,CSEG,0.d0,ERROR,*9999)
        ENDDO
        CALL DAWK(iw,0,ERROR,*9999)
      ENDIF

      CALL EXITS('CHCLOC')
      RETURN
 9999 CALL ERRORS('CHCLOC',ERROR)
      CALL EXITS('CHCLOC')
      RETURN 1
      END


      SUBROUTINE CHCOLO(ISEG,noco,NTCO,CO,CSEG,STRING,ERROR,*)

C#### Subroutine: CHCOLO
C###  Description:
C###    Change colours.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),noco,NTCO
      CHARACTER CO(*)*(*),CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,ILUT_BLU,ILUT_RED,INDEX_OLD,INSTAT,
     '  IPICK,ISEG0,ISEG1,ISEG2,ISEGM,iw,LD1,N3CO
      REAL*8 XWC,YWC
      REAL*8 PT(3,2),XCOL_BLU,XCOL_RED,YCOL_BLU,YCOL_RED
      LOGICAL CBBREV
      DATA LD1/1/

      CALL ENTERS('CHCOLO',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)//';m <on WS>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHCOLO',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          iw=IFROMC(CO(N3CO+1))
        ELSE
          iw=1
        ENDIF

        IF(iw.EQ.4) PROJEC='RECTANGULAR'
        CALL ACWK(iw,0,ERROR,*9999)
        ILUT_RED=3
        ILUT_BLU=249
        CALL SET_COLOUR_REP(iw,ILUT_RED,ILUT_BLU,ERROR,*9999)
        IF(iw.EQ.1) THEN
          CALL SET_COLOUR_REP(91,ILUT_RED,ILUT_BLU,ERROR,*9999)
          XCOL_RED=DBLE(XMAX)-0.22d0*(DBLE(XMAX)-DBLE(XMIN))
          XCOL_BLU=XCOL_RED
          YCOL_RED=DBLE(YMIN)+0.9d0*(DBLE(YMAX)-DBLE(YMIN))
          YCOL_BLU=DBLE(YMIN)+0.1d0*(DBLE(YMAX)-DBLE(YMIN))
        ELSE IF(iw.EQ.4) THEN
          CALL SET_COLOUR_REP(94,ILUT_RED,ILUT_BLU,ERROR,*9999)
          XCOL_RED= 0.87d0
          XCOL_BLU=-0.9d0
          YCOL_RED= 0.6d0
          YCOL_BLU=YCOL_RED
        ENDIF
        ISEG0=0
        CALL OPEN_SEGMENT(ISEG0,ISEG,iw,'COLOUR_MARKER',1,INDEX_OLD,
     '    0,1,CSEG,ERROR,*9999)
        PT(1,1)=XCOL_RED
        PT(2,1)=YCOL_RED
        PT(1,2)=XCOL_BLU
        PT(2,2)=YCOL_BLU
        CALL POLYLINE(1,iw,2,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISEG0,iw,ERROR,*9999)
        ISEG1=0
        CALL OPEN_SEGMENT(ISEG1,ISEG,iw,'COLOUR_MARKER',1,INDEX_OLD,
     '    1,1,CSEG,ERROR,*9999)
        PT(1,1)=XCOL_RED
        PT(2,1)=YCOL_RED
        CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISEG1,iw,ERROR,*9999)
        CALL DETECT(iw,ISEG,ISEG1,'DETECTABLE',ERROR,*9999)
        ISEG2=0
        CALL OPEN_SEGMENT(ISEG2,ISEG,iw,'COLOUR_MARKER',1,INDEX_OLD,
     '    2,1,CSEG,ERROR,*9999)
        PT(1,2)=XCOL_BLU
        PT(2,2)=YCOL_BLU
        CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
        CALL CLOSE_SEGMENT(ISEG2,iw,ERROR,*9999)
        CALL DETECT(iw,ISEG,ISEG2,'DETECTABLE',ERROR,*9999)
        INSTAT=1
        DO WHILE(INSTAT.EQ.1)
          WRITE(OP_STRING,'('' >>Pick colour marker on scale'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL PICK(iw,LD1,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,*9999)
          IF(INSTAT.EQ.1) THEN
            WRITE(OP_STRING,'('' >>Relocate colour marker'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL LOCATOR(1,iw,INSTAT,'REQUEST',3,0.d0,XWC,0.d0,YWC,
     '        ERROR,*9999)
            IF(INSTAT.EQ.1) THEN
              IF(iw.EQ.1) THEN
                IF(YWC.GT.YCOL_RED) THEN
                  YWC=YCOL_RED
                ELSE IF(YWC.LT.YCOL_BLU) THEN
                  YWC=YCOL_BLU
                ENDIF
              ELSE IF(iw.EQ.4) THEN
                IF(XWC.GT.XCOL_RED) THEN
                  XWC=XCOL_RED
                ELSE IF(XWC.LT.XCOL_BLU) THEN
                  XWC=XCOL_BLU
                ENDIF
              ENDIF
              CALL OPEN_SEGMENT(ISEGM,ISEG,iw,'COLOUR_MARKER',1,
     '          INDEX_OLD,1,1,CSEG,ERROR,*9999)
              IF(ISEGM.EQ.ISEG1) THEN
                IF(iw.EQ.1) THEN
                  PT(1,1)=XCOL_RED
                  PT(2,1)=YWC
                  ILUT_RED=3+INT((YCOL_RED-YWC)/(YCOL_RED-YCOL_BLU)
     '                                         *246.d0)
                ELSE IF(iw.EQ.4) THEN
                  PT(1,1)=XWC
                  PT(2,1)=YCOL_RED
                  ILUT_RED=3+INT((XCOL_RED-XWC)/(XCOL_RED-XCOL_BLU)
     '                                         *246.d0)
                ENDIF
                CALL POLYMARKER(1,iw,1,PT(1,1),ERROR,*9999)
              ELSE IF(ISEGM.EQ.ISEG2) THEN
                IF(iw.EQ.1) THEN
                  PT(1,2)=XCOL_BLU
                  PT(2,2)=YWC
                  ILUT_BLU=249-INT((YWC-YCOL_BLU)/(YCOL_RED-YCOL_BLU)
     '                                         *246.d0)
                ELSE IF(iw.EQ.4) THEN
                  PT(1,2)=XWC
                  PT(2,2)=YCOL_BLU
                  ILUT_BLU=249-INT((XWC-XCOL_BLU)/(XCOL_RED-XCOL_BLU)
     '                                         *246.d0)
                ENDIF
                CALL POLYMARKER(1,iw,1,PT(1,2),ERROR,*9999)
              ENDIF
              CALL CLOSE_SEGMENT(ISEGM,iw,ERROR,*9999)
              CALL SET_COLOUR_REP(iw,ILUT_RED,ILUT_BLU,ERROR,*9999)
              IF(iw.EQ.1) THEN
                CALL SET_COLOUR_REP(91,ILUT_RED,ILUT_BLU,ERROR,*9999)
              ELSE IF(iw.EQ.4) THEN
                CALL SET_COLOUR_REP(94,ILUT_RED,ILUT_BLU,ERROR,*9999)
              ENDIF
              CALL DETECT(iw,ISEG,ISEGM,'DETECTABLE',ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
        CALL DELETE_SEGMENT(ISEG0,ISEG,iw,ERROR,*9999)
        CALL DELETE_SEGMENT(ISEG1,ISEG,iw,ERROR,*9999)
        CALL DELETE_SEGMENT(ISEG2,ISEG,iw,ERROR,*9999)
        CALL DAWK(iw,0,ERROR,*9999)
      ENDIF

      CALL EXITS('CHCOLO')
      RETURN
 9999 CALL ERRORS('CHCOLO',ERROR)
      CALL EXITS('CHCOLO')
      RETURN 1
      END


      SUBROUTINE CHLUT(noco,NTCO,CO,STRING,ERROR,*)

C#### Subroutine: CHLUT
C###  Description:
C###    Change colour lookup table.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER noco,NTCO
      CHARACTER CO(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX_FINISH,INDEX_START,iw,N3CO
      CHARACTER RGB_FINISH*10,RGB_START*10
      LOGICAL CBBREV,DEFAULT

      CALL ENTERS('CHLUT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)//' <default>'
        OP_STRING(2)=BLANK(1:15)//'<on WS>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<from START_%>[0]'
        OP_STRING(3)=BLANK(1:15)//'<to FINISH_%>[100]'
        OP_STRING(4)=BLANK(1:15)//'<start START_COLOUR>[black]'
        OP_STRING(5)=BLANK(1:15)//'<end FINISH_COLOUR>[white]'
        OP_STRING(6)=BLANK(1:15)//'<on WS>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHLUT',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'DEFAULT',1,noco+1,NTCO,N3CO)) THEN
          DEFAULT=.TRUE.
        ELSE
          DEFAULT=.FALSE.
          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            INDEX_START=IFROMC(CO(N3CO+1))
          ELSE
            INDEX_START=0
          ENDIF
          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            INDEX_FINISH=IFROMC(CO(N3CO+1))
          ELSE
            INDEX_FINISH=100
          ENDIF
          IF(CBBREV(CO,'START',1,noco+1,NTCO,N3CO)) THEN
            RGB_START=CO(N3CO+1)
          ELSE
            RGB_START='BLACK'
          ENDIF
          IF(CBBREV(CO,'END',1,noco+1,NTCO,N3CO)) THEN
            RGB_FINISH=CO(N3CO+1)
          ELSE
            RGB_FINISH='WHITE'
          ENDIF
          IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
            iw=IFROMC(CO(N3CO+1))
          ELSE
            iw=1
          ENDIF
        ENDIF
        IF(DEFAULT) THEN
          CALL SET_COLOUR_LUT(COLOUR_LUT,ERROR,*9999)
        ELSE
          CALL SET_COLOUR_LUT_RANGE(INDEX_START,INDEX_FINISH,COLOUR_LUT,
     '      RGB_START,RGB_FINISH,ERROR,*9999)
        ENDIF
        CALL SET_COLOUR_REP(iw,3,249,ERROR,*9999)
      ENDIF

      CALL EXITS('CHLUT')
      RETURN
 9999 CALL ERRORS('CHLUT',ERROR)
      CALL EXITS('CHLUT')
      RETURN 1
      END


      SUBROUTINE NODE3D(XWC,YWC,ZWC)

C**** Not used anymore? - check image o/p
                    
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:fbgr00.cmn'
      INCLUDE 'cmiss$reference:trans00.cmn'
!     Parameter List
      REAL XWC,YWC,ZWC
!     Local Variables
      REAL Z1(3),Z2(3)

C     CALL ENTERS('NODE3D',*9999)
      IF(RHTRAN) THEN
        CALL FBIMGPT(5,XWC,YWC,ZWC,Z1(1),Z1(2))
        XWC=Z1(1)
        YWC=Z1(2)
      ELSE IF(LFTRAN) THEN
        CALL FBIMGPT(6,XWC,YWC,ZWC,Z1(1),Z1(2))
        XWC=Z1(1)
        YWC=Z1(2)
      ELSE
        Z1(1)=XWC
        Z1(2)=YWC
        Z1(3)=ZWC
        CALL ZZ(Z1,Z2,TRANS)
        XWC=Z2(1)
        YWC=Z2(2)
        ZWC=Z2(3)
      ENDIF

C     CALL EXITS('NODE3D')
      RETURN
      END

      SUBROUTINE RECFID(STRING,ERROR,*)
      
C#### Subroutine: RECFID
C###  Description:
C###    RECFID recieves information about the fiducial markers from 
C###    a Motif front end, through the socket interface.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('RECFID',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','RECFID',ERROR,*9999)
      ELSE
        IF(USE_SOCKET) THEN !Socket interface to front-end is being used
           CALL SOCKET_READ_SIGNAL(ERROR,*9999)
        ELSE     !Socket interface not being used
          !Cannot receive info
        ENDIF
      ENDIF

      CALL EXITS('RECFID')
      RETURN
 9999 CALL ERRORS('RECFID',ERROR)
      CALL EXITS('RECFID')
      RETURN 1
      END


      SUBROUTINE REFINE_TRIANGLE(IBT,IDO,INP,ISC_GK,ISR_GK,LGE,
     '  NAN,NBH,NBHF,NBJ,NBJF,
     '  NEELEM,NEL,NELIST,NENP,NFF,NFFACE,NHE,NHP,
     '  NKEF,NKH,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NNB,NNF,NNL,
     '  NONY,NPF,NP_INTERFACE,NPL,NPLIST,NPNE,NPNF,NPNODE,NPNY,nr,
     '  NRE,NRLIST,NVHE,NVHP,NVJE,NVJF,NVJL,NVJP,
     '  NW,NWP,NXI,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '  PLSG_EDGE_LIST,PLSG_EDGE_MARKER_LIST,
     '  PLSG_POINT_MARKER_LIST,
     '  PLSG_SEGMENT_LIST,PLSG_SEGMENT_MARKER_LIST,CE,CONY,CYNO,
     '  DF,DL,PG,
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,PLSG_PARTITION_LIST,
     '  RG,SE,SF,SP,WG,XA,XE,XG,XP,YP,ZA,FIX,ERROR,*)

C#### Subroutine: REFINE_TRIANGLE
C###  Description:
C###    REFINE_TRIANGLE wraps REFINE_TRIANGLE_DYNAM dynamically 
C###    allocating and freeing the arrays required in
C###    REFINE_TRIANGLE_DYNAM
C**** Created by Carey Stevens 18 September 1997

C#### Variable: PLSG_ELEMENT_AREA(ne) 
C###  Type: REAL*8
C###  Set_up: REFINE_TRIANGLE
C###  Description: 
C###    PLSG_ELEMENT_AREA is an array of triangle area constraints.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),
     '  nr,NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     '  NVJP(NJM,NPM),NW(NEM,3,NXM),NWP(NPM,2),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM),
     '  PLSG_EDGE_LIST(NLM*2),PLSG_EDGE_MARKER_LIST(NLM),
     '  PLSG_POINT_MARKER_LIST(NPM),PLSG_SEGMENT_LIST(NLM*2),
     '  PLSG_SEGMENT_MARKER_LIST(NLM)
      REAL*8 CE(NMM,NEM,NXM),CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),
     '  PLSG_ELEMENT_AREA(NEM),PLSG_HOLE_LIST(NRM*2),
     '  PLSG_PARTITION_LIST(NEM*2),
     '  RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  SP(NKM,NBFM,NPM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER*4 NUM_NODE_VERSIONS_PTR,POINT_ATTRIBUTE_LIST_PTR,
     '  TRIANGLE_ATTRIBUTE_LIST_PTR

      CALL ENTERS('REFINE_TRIANGLE',*9999)

      NUM_NODE_VERSIONS_PTR=0
      POINT_ATTRIBUTE_LIST_PTR=0
      TRIANGLE_ATTRIBUTE_LIST_PTR=0

      CALL ALLOCATE_MEMORY(NPM,1,1,NUM_NODE_VERSIONS_PTR,.TRUE.,
     '  ERROR,*9999)
      CALL ALLOCATE_MEMORY(NPM*NAM,1,3,POINT_ATTRIBUTE_LIST_PTR,.TRUE.,
     '  ERROR,*9999)
      CALL ALLOCATE_MEMORY(NEM*NAM,1,3,TRIANGLE_ATTRIBUTE_LIST_PTR,
     '  .TRUE.,ERROR,*9999)

      CALL REFINE_TRIANGLE_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,LGE,NAN,
     '  NBH,NBHF,NBJ,NBJF,NEELEM,NEL,NELIST,NENP,NFF,NFFACE,
     '  NHE,NHP,NKEF,NKH,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NNB,NNF,NNL,NONY,
     '  NPF,NP_INTERFACE,NPL,NPLIST,NPNE,NPNF,NPNODE,NPNY,nr,NRE,NRLIST,
     '  %VAL(NUM_NODE_VERSIONS_PTR),NVHE,NVHP,NVJE,NVJF,NVJL,
     '  NVJP,NW,NWP,NXI,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '  PLSG_EDGE_LIST,PLSG_EDGE_MARKER_LIST,
     '  PLSG_POINT_MARKER_LIST,PLSG_SEGMENT_LIST,
     '  PLSG_SEGMENT_MARKER_LIST,CE,CONY,CYNO,DF,DL,PG,
     '  %VAL(POINT_ATTRIBUTE_LIST_PTR), 
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,PLSG_PARTITION_LIST,RG,SE,SF,
     '  SP,%VAL(TRIANGLE_ATTRIBUTE_LIST_PTR),WG,XA,XE,XG,XP,YP,ZA,FIX,
     '  ERROR,*9999)

      CALL FREE_MEMORY(NUM_NODE_VERSIONS_PTR,ERROR,*9999)
      CALL FREE_MEMORY(POINT_ATTRIBUTE_LIST_PTR,ERROR,*9999)
      CALL FREE_MEMORY(TRIANGLE_ATTRIBUTE_LIST_PTR,ERROR,*9999)

      CALL EXITS('REFINE_TRIANGLE')
      RETURN
 9999 CALL ERRORS('REFINE_TRIANGLE',ERROR)
      CALL EXITS('REFINE_TRIANGLE')
      RETURN 1
      END


      SUBROUTINE REFINE_TRIANGLE_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,
     '  LGE,NAN,NBH,NBHF,
     '  NBJ,NBJF,NEELEM,NEL,NELIST,NENP,NFF,NFFACE,
     '  NHE,NHP,NKEF,
     '  NKH,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NNB,NNF,NNL,NONY,NPF,
     '  NP_INTERFACE,NPL,NPLIST,NPNE,NPNF,NPNODE,NPNY,nr,NRE,NRLIST,
     '  NUM_NODE_VERSIONS,NVHE,NVHP,NVJE,NVJF,NVJL,
     '  NVJP,NW,NWP,NXI,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '  PLSG_EDGE_LIST,PLSG_EDGE_MARKER_LIST,
     '  PLSG_POINT_MARKER_LIST,PLSG_SEGMENT_LIST,
     '  PLSG_SEGMENT_MARKER_LIST,CE,CONY,CYNO,
     '  DF,DL,PG,POINT_ATTRIBUTE_LIST, 
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,PLSG_PARTITION_LIST,RG,SE,SF,
     '  SP,TRIANGLE_ATTRIBUTE_LIST,WG,XA,XE,XG,XP,YP,ZA,FIX,ERROR,*)

C#### Subroutine: REFINE_TRIANGLE_DYNAM
C###  Description:
C###    Refines a triangulation created using IPMESH9
C###  SEE-ALSO: IPMESH9_DYNAM
C**** Created by Carey Stevens 18 September 1997

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b13.cmn'
      INCLUDE 'cmiss$reference:cbfe01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'cmiss$reference:mesh00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:solv00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NLF(4,NFM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNB(4,4,4,NBFM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPL(5,0:3,NLM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),nr,
     '  NRE(NEM),NRLIST(0:NRM),
     '  NUM_NODE_VERSIONS(NPM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),NVJL(4,NJM,NLM),
     '  NVJP(NJM,NPM),NW(NEM,3,NXM),NWP(NPM,2),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM),
     '  PLSG_EDGE_LIST(NLM*2),
     '  PLSG_EDGE_MARKER_LIST(NLM),PLSG_POINT_MARKER_LIST(NPM),
     '  PLSG_SEGMENT_LIST(NLM*2),PLSG_SEGMENT_MARKER_LIST(NLM)
      REAL*8 CE(NMM,NEM,NXM),CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),
     '  PLSG_ELEMENT_AREA(NEM),PLSG_HOLE_LIST(NRM*2),
     '  PLSG_PARTITION_LIST(NEM*2),POINT_ATTRIBUTE_LIST(NPM*NAM),
     '  RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),SP(NKM,NBFM,NPM),
     '  TRIANGLE_ATTRIBUTE_LIST(NEM*NAM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER basis1,ELEMENT_TYPE,ERROR_FLAG,
     '  ERROR_STRINGC(500),ie,il,LENGTH,nb,nb1,NBH_LOCAL(4),nc,ne,nh,
     '  nhx,nhx_MAX,NLV(12,3),nj,njj,njj1,njj2,nk,
     '  NKJP(3*NJ_LOC_MX),nn,ns,np,noelem,nonode,NUM_POINTS,
     '  NUM_TRIANGLE_ATTRIBUTES,nv,nvar,nx,nxc
      INTEGER*4 WORK_PTR
      CHARACTER ERROR_STRINGF*100
      LOGICAL CALC_SPARSITY

      DATA NLV/2,3,2,1,3,3,3,3,1, 1, 2, 2, !njt=1
     '         2,3,2,1,3,3,3,3,2, 1, 2, 2, !njt=2
     '         2,3,3,1,3,3,3,3,3, 1, 2, 2/ !njt=3

      CALL ENTERS('REFINE_TRIANGLE_DYNAM',*9999)

      CALL ASSERT(USE_TRIANGLE.EQ.1,
     '  '>>USE_TRIANGLE must be set to 1',ERROR,*9999)
      CALL ASSERT(NJT.EQ.2,
     '  '>>The number of global coordinates must be 2',ERROR,*9999)

C     Find the largest NKJ for each field nj
      DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
        nj=NJ_LOC(NJL_FIEL,njj,nr)
        NKJP(nj)=0
        DO nonode=1,NPNODE(0,nr)
          IF(NKJ(nj,NPNODE(nonode,nr)).GT.NKJP(nj)) THEN
            NKJP(nj)=NKJ(nj,NPNODE(nonode,nr))
          ENDIF
        ENDDO
      ENDDO

      nb=NBJ(1,NEELEM(1,nr))
      nxc=1 !?????????????
      nhx=1 !?????????????
      CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
      ELEMENT_TYPE=NW(NEELEM(1,nr),1,nx) 
      IF(nx.ne.0) THEN
        basis1=NBH(NH_LOC(nhx,nx),1,NEELEM(1,nr))
      ELSE
        basis1=1
      ENDIF

      NUM_POINTS=0     
      DO np=1,NPNODE(0,nr)
        NUM_NODE_VERSIONS(np)=1     !To account for one version
        NUM_POINTS=NUM_POINTS+1
      ENDDO

      NUM_TRIANGLE_ATTRIBUTES=0
      ERROR_FLAG=0

      CALL triangulate(AREA_FIELD,ERROR_FLAG,FLAGS,IBT,
     '  IDO,INP,NAM,NBFM,NBJ,NCM,
     '  NEELEM,NEM,NHM,NIM,NJM,NKJE,NKM,NLM,NNM,NPF,NPM,NPNE,nr,
     '  NNT(nb),NUM_EDGES,NUM_HOLES,NPNODE(0,nr),NUM_NODE_VERSIONS,
     '  NUM_PARTITIONS,
     '  NUM_POINTS,0,NUM_SEGMENTS,
     '  NUM_STEINER_POINTS,NEELEM(0,nr),
     '  NUM_TRIANGLE_ATTRIBUTES,NVJE,NVJP,NVM,NXI,PLSG_EDGE_LIST,
     '  PLSG_EDGE_MARKER_LIST,PLSG_POINT_MARKER_LIST,PLSG_SEGMENT_LIST,
     '  PLSG_SEGMENT_MARKER_LIST,MAX_AREA,MIN_ANGLE,
     '  PLSG_ELEMENT_AREA,PLSG_HOLE_LIST,
     '  PLSG_PARTITION_LIST,POINT_ATTRIBUTE_LIST,SE,
     '  TRIANGLE_ATTRIBUTE_LIST,XA,XE,XP,ZA,ERROR_STRINGC)

      IF(ERROR_FLAG.NE.0) THEN
        CALL CSTRINGLEN(LENGTH,ERROR_STRINGC)
        CALL C2FSTRING(ERROR_STRINGC,LENGTH,ERROR_STRINGF)
        CALL ASSERT(.FALSE.,ERROR_STRINGF,ERROR,*9999)
      ENDIF

      NET(nr)=NEELEM(0,nr)       !highest element# in region nr
      NET(0)=NEELEM(0,nr)        !highest element# in all regions
      NEELEM(0,0) =NEELEM(0,nr)  !number elements in all regions
      NPT(nr)=NUM_POINTS      !highest node# in region nr
      NPT(0) =NUM_POINTS       !highest node# in all regions
      NPNODE(0,nr)=NUM_POINTS  !number nodes in region nr
      NPNODE(0,0) =NUM_POINTS  !number nodes in all regions
      CALL ASSERT(NPT(nr).LE.NP_R_M,'>>Increase NP_R_M',ERROR,*9999)
      CALL ASSERT(NET(nr).LE.NE_R_M,'>>Increase NE_R_M',ERROR,*9999)
      CALL_TRIANGLE=.TRUE.

C     Set up new nodes
C!!! Currently renumbers nodes from 1..the total number of nodes in nr
      DO nonode=1,NPNODE(0,nr)
        np=nonode
        NPNODE(nonode,nr)=np
        DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
          nj=NJ_LOC(NJL_GEOM,njj,nr)
          NVJP(nj,np)=1 !one version per geom nj per node
        ENDDO !nj
        DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          NKJ(nj,np)=NKJP(nj)
          NVJP(nj,np)=1 !one version per field nj per node
          DO nv=1,NVJP(nj,np)
            DO nk=1,NKJP(nj)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              XP(nk,nv,nj,np)=0.0d0
            ENDDO !nk
          ENDDO !nv
        ENDDO !njj
      ENDDO !nonode

C     Set up elements
C     NPNE for nb>1 is the same as for nb=1 
      DO nb1=2,NBT
        DO ne=1,NET(nr)
          DO nn=1,NNT(nb1)
            NPNE(nn,nb1,ne)=NPNE(nn,1,ne)
          ENDDO
        ENDDO
      ENDDO

      DO noelem=1,NEELEM(0,nr)
        ne=noelem
        NEELEM(noelem,nr)=ne 
        NRE(ne)=nr
        NW(ne,1,nx)=ELEMENT_TYPE
        DO nb1=1,NBFT
          DO njj1=1,3 !geometry/fibres/field
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              DO nn=1,NNT(nb1)
                NVJE(nn,nb1,nj,ne)=1
                NBJ(nj,ne)=nb
              ENDDO !nn
            ENDDO !njj2
          ENDDO !njj1
        ENDDO !nb
        DO njj1=1,3 !geometry/fibres/field
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            nb1=NBJ(nj,ne)
            DO nn=1,NNT(nb1)
              DO nk=1,NKT(nn,nb1)
                NKJE(nk,nn,nj,ne)=nk
              ENDDO !nk
            ENDDO !nn
          ENDDO !njj2
        ENDDO !njj1
      ENDDO

      DO nb1=1,NBFT
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ns=1,NST(nb1)+NAT(nb1)
            SE(ns,nb1,ne)=1.0D0
          ENDDO !ns
        ENDDO !noelem (ne)
      ENDDO !nb1

      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
      CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
      CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '  NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '  DL,SE,XP,ERROR,*9999)
      CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '  NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,NVJF,
     '  DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
      CALL GLOBALJ(NBJ,NEELEM,NKJ,NPLIST,NPNE,NPNODE,
     '  ERROR,*9999)


C     Set up dependent variable information
      IF(nx.ne.0) THEN
        ie=11 !plane stress
        nvar=0
        DO nhx=1,NLV(ie,NJT)
          NBH_LOCAL(nhx)=0
          nvar=nvar+1
          NHV(nvar,ie)=nhx
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NW(ne,1,nx)=ie
            NBH_LOCAL(nhx)=basis1
          ENDDO
        ENDDO

C       Set NHE(ne,nx) from NVE
C       Determine max# nhx's
        nhx_MAX=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          ie=NW(ne,1,nx)
          NHE(ne,nx)=NVE(ie)
          IF(NHE(ne,nx).GT.nhx_MAX) nhx_MAX=NHE(ne,nx)
        ENDDO !noelem
        
C       Set up NBH
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhx=1,NHE(ne,nx)
            nh=NH_LOC(nhx,nx)
            DO nc=1,2
              NBH(nh,nc,ne)=NBH_LOCAL(nhx)
            ENDDO !nc
          ENDDO !nhx
        ENDDO !noelem (ne)
      
C       Set NHP(np,nr,nx) from NHE(ne,nx).
C       (Needs to be done after nh_loc setup)
        CALL CALC_NHP(NBH,NEELEM,NHE,NHP(1,1,nx),NPNE,NPNODE,nr,nx,
     '    ERROR,*9999)
        CALL CALC_VERSIONS_DEP(NBH,NBJ,NEELEM,NHE,NHP(1,1,nx),
     '    NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NVJP,nx,ERROR,*9999)
C     Initialise all NKH in this region          
        CALL CALC_NKH(NBH,NEELEM,NHP(1,1,nx),NKH,NPNE,NPNODE,
     '    nr,NW(1,1,nx),nx,ERROR,*9999)
        CALL CALC_NY_MAPS_DEP(NBH,NEELEM,NHP(1,1,nx),NKH,
     '    NP_INTERFACE,NPNODE,NPNY,nr,NVHP,nx,
     '    NYNE,NYNP,NYNR,ERROR,*9999)
        CALL CALC_FACE_BASIS_DEP(IBT,NBH,NBHF,NEELEM,NFF,NHE,NNF,
     '    nr,nx,ERROR,*9999)
        WORK_PTR=0
        CALL ALLOCATE_MEMORY(NYT(1,1,nx)*NYT(2,1,nx),1,CHARTYPE,
     '    WORK_PTR,MEM_INIT,ERROR,*9999)
        CALC_SPARSITY=.TRUE.
        CALL CALC_SPARSE_SOLVE(NISC_GKM,NISR_GKM,ISC_GK,ISR_GK,LGE,
     '    NYT(1,1,nx),NYT(2,1,nx),NBH,1,NEELEM,NHE,NPNE,NPNY(0,1,0,nx),
     '    NRLIST,NVHE,nx,NYNE,NYNP,NYNR,NZ_GK_M,KTYP24,%VAL(WORK_PTR),
     '    FIX(1,1,nx),CALC_SPARSITY,ERROR,*9999)
        CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
        ASSEMBLE_GLOBAL(nr,nx)=.FALSE.
        CALL_EQUA=.TRUE.


        
C!!!  Set up initial conditions here.

C       Set up solve
        CALL GLOBALH(IBT,IDO,INP,NAN,NBH,NBJ,NELIST,NENP,NHE,NKHE,NKJE,
     '    NLL,NNB,NNF,NNL,NONY,NP_INTERFACE,NPF,NPL,NPNE,NPNY,nr,NRE,
     '    NVHE,NVHP,NVJE,NWP,nx,NXI,NYNE,NYNO,NYNP,NYNR(0,0,1,0,nx),
     '    NYQNR(0,0,1,0,nx),CONY,CYNO,SE,SP,XA,XE,XP,YP,
     '    FIX(1,1,nx),ERROR,*9999)
        CALL_SOLV=.TRUE.


C       Set up material
        DO il=1,ILT(ie,nr,nx)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CE(il,ne,nx)=CE(il,NEELEM(1,nr),nx)
          ENDDO
        ENDDO
        CALL_MATE=.TRUE.
      ENDIF


      CALL EXITS('REFINE_TRIANGLE_DYNAM')
      RETURN
 9999 CALL ERRORS('REFINE_TRIANGLE_DYNAM',ERROR)
      CALL EXITS('REFINE_TRIANGLE_DYNAM')
      RETURN 1
      END


C LKC 6-OCT-1999 Archived

      SUBROUTINE SENDFID(STRING,ERROR,*)
      
C#### Subroutine: SENDFID
C###  Description:
C###    SENDFID sends information about the fiducial markers from a 
C###    Motif front end, through the socket interface.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:pics00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('SENDFID',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','SENDFID',ERROR,*9999)
      ELSE
        IF(USE_SOCKET) THEN !Socket interface to front-end is being used
          CALL WRITE_SIGNAL(I2P,NUMSIGNALS,NUMSAMPLES,
     '      NUMSAMPLESPERSECOND,ERROR,*9999)
        ELSE     !Socket interface not being used
          !Cannot send info
        ENDIF
      ENDIF

      CALL EXITS('SENDFID')
      RETURN
 9999 CALL ERRORS('SENDFID',ERROR)
      CALL EXITS('SENDFID')
      RETURN 1
      END


C 18 August 1997
      SUBROUTINE UPGRID(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NENQ,
     '  NHE,NHP,NJE,NKE,NKH,NPF,NPNE,NPNODE,NPQ,NQE,NQGE,NQNE,NQXI,NQS,
     '  NQSCNB,NRLIST,NVHE,NVHP,NVJE,NW,NWQ,NXQ,NYNE,NYNP,
     '  CQ,DNUDXQ,DXDXIQ,FEXT,GCHQ,GUQ,PG,PROPQ,SE,
     '  XA,XE,XP,XQ,XQD,XQDRC,YP,YQ,ZA,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: UPGRID
C###  Description:
C###    UPGRID updates grid point parameters, including material 
C###    values (conductivity tensors) and potentials and source terms
C###    for the bidomain model.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),
     '  NENQ(0:8,NQM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NJE(NEM),
     '  NKE(NKM,NNM,NBFM,NEFM),NKH(NHM,NPM,NCM,0:NRM),NPF(15,NFM),
     '  NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),NPQ(NQM),
     '  NQGE(NGM,NEFM,NBM),NQNE(NEM,NQEM),NQS(NEM),NQSCNB(NQSCM),
     '  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEFM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NW(NEM,2),
     '  NWQ(6,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NQE(NSM,NBFM,NEFM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NVJE(NNM,NBFM,NJM,NEFM)
      REAL*8 CQ(NMM,NQM,NXM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     '  FEXT(NIFEXTM,NGM,NEM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PG(NSM,NUM,NGM,NBM),XA(NAM,NJM,NQM),
     '  PROPQ(3,3,4,2,NQM,NXM),PXI,
     '  SE(NSM,NBFM,NEFM),XQ(NJM,NQM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     '  XQD(NJM,NUM),XQDRC(NJM,NUM),
     '  YP(NYM,NIYM,NXM),YQ(NQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,i,ib,IFNTYP,II,IJ,IK,iL,INFO,
     '  IFROMC,j,jb,k,mLL,mq,mRR,nb,nb_extended,nbb,nc,ne,nee,neq,
     '  ng,ngBL,ngBR,ngTL,ngTR,NG_row,nii2,nij2,nik2,
     '  N3CO,ni,nii,nij,nik,NITB,nj,noelem,NU1(3),
     '  no_nrlist,NOQUES,nq,nq1,nq2,NQ_row,nr,nrr,nr_s,nr_t,ns,nTT,
     '  nx,nxc,nxielem,nx_s,nx_t,SCHEME
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  Current,DBM(3,3,3),DET,DIFFN,dNudXi(3,3),dPHIdX(3),
     '  dXidNu(3,3),d2PHIdX2(3,3),XI(3),DXIDXI(3),CHTOFF(3,3,3),
     '  Esac,Ext,GX,GY,PGG,PHIMQ(-1:1,-1:1,-1:1),GL(3,3),GU(3,3),
     '  SUM,SUM1,SUM2,Xi1,Xi2,XQ_TEMP(3),X3G(4,3)
      CHARACTER TYPE*22
      LOGICAL ALL_REGIONS,CBBREV,DIRECT,FILEIP
      DATA NU1/2,4,7/

      CALL ENTERS('UPGRID',*9999)

      nc=1 !temporary

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update grid geometry
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    This command calculates the geometric positions of the grid
C###    points and stores them in XQ. The grid is spaced evenly in
C###    material coordinates. The geometric positions are calculated
C###    by interpolation using the finite element basis functions used
C###    to describe the elements in which the grid points lie.

        OP_STRING(1)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid geometry deformed
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    This command changes the geometric positions of the grid
C###    points when the underlying finite element mesh has undergone
C###    a deformation. Because the grid points don't move in material
C###    space, they will deform with the finite elements.

        OP_STRING(1)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid metric
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    This command generates the metric tensor information
C###    associated with a grid. This information is necessary
C###    for activation solutions.

        OP_STRING(1)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid material/strain/source
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' material/strain/source'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid time
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' time'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid intracellular_bidomain
C###  Parameter:      <from_region #[2]>
C###  Parameter:      <from_class #[2]>
C###  Parameter:      <to_region #[1]>
C###  Parameter:      <to_class #[1]>
C###  Parameter:      <direct>
C###  Description: Updates the collocation array YQ from
C###    the extracellular potential.

        OP_STRING(1)=STRING(1:IEND)//' intracellular_bidomain'
        OP_STRING(2)=BLANK(1:IEND)//'<direct>'
        OP_STRING(3)=BLANK(1:IEND)//'<from_region #[2]>'
        OP_STRING(4)=BLANK(1:IEND)//'<from_class #[2]>'
        OP_STRING(5)=BLANK(1:IEND)//'<to_region #[1]>'
        OP_STRING(6)=BLANK(1:IEND)//'<to_class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid extracellular_bidomain
C###  Parameter:      <from_region #[1]>
C###  Parameter:      <from_class #[1>]
C###  Parameter:      <to_region #[2]>
C###  Parameter:      <to_class #[2]>
C###  Description: Updates the collocation array YQ from
C###    the 2nd derivative of intracellular potential.

        OP_STRING(1)=STRING(1:IEND)//' extracellular_bidomain'
        OP_STRING(2)=BLANK(1:IEND)//'<from_region #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<to_region #[2]>'
        OP_STRING(5)=BLANK(1:IEND)//'<to_class #[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','UPGRID',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'MATERIAL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='MATERIAL'
        ELSE IF(CBBREV(CO,'INTRACELLULAR_BIDOMAIN',1,noco+1,NTCO,N3CO))
     '    THEN
          TYPE='INTRACELLULAR_BIDOMAIN'
        ELSE IF(CBBREV(CO,'EXTRACELLULAR_BIDOMAIN',1,noco+1,NTCO,N3CO))
     '    THEN
          TYPE='EXTRACELLULAR_BIDOMAIN'
        ELSE IF(CBBREV(CO,'STRAIN',2,noco+1,NTCO,N3CO)) THEN
          TYPE='STRAIN'
        ELSE IF(CBBREV(CO,'TIME',1,noco+1,NTCO,N3CO)) THEN
          TYPE='TIME'
        ELSE IF(CBBREV(CO,'SOURCE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SOURCE'
        ELSE IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='GEOMETRY'
        ELSE IF(CBBREV(CO,'METRIC',2,noco+1,NTCO,N3CO)) THEN
          TYPE='METRIC'
        ENDIF

        IF(TYPE(1:22).EQ.'INTRACELLULAR_BIDOMAIN') THEN
          IF(CBBREV(CO,'FROM_REGION',6,noco+1,NTCO,N3CO)) THEN
            nr_s=IFROMC(CO(N3CO+1))
          ELSE
            nr_s=2
          ENDIF

          IF(CBBREV(CO,'FROM_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
      	  CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'Invalid source class',ERROR,*9999)

          CALL ASSERT(ITYP5(nr_s,nx_s).EQ.1.AND.ITYP2(nr_s,nx_s).EQ.5,
     '      'Source equation must be div(k.grad(u))=f',ERROR,*9999)

          IF(CBBREV(CO,'TO_REGION',4,noco+1,NTCO,N3CO)) THEN
            nr_t=IFROMC(CO(N3CO+1))
          ELSE
            nr_t=1
          ENDIF

          IF(CBBREV(CO,'TO_CLASS',4,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_t.GT.0,'>>Invalid target class',ERROR,*9999)
          CALL ASSERT(ITYP5(nr_t,nx_t).EQ.2.AND.ITYP2(nr_t,nx_t).EQ.9,
     '      'Target equation must be activation model',ERROR,*9999)
          
          IF(CBBREV(CO,'DIRECT',1,noco+1,NTCO,N3CO)) THEN
            DIRECT=.TRUE.
          ELSE
            DIRECT=.FALSE.
          ENDIF

        ELSE IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN
          IF(CBBREV(CO,'FROM_REGION',6,noco+1,NTCO,N3CO)) THEN
            nr_s=IFROMC(CO(N3CO+1))
          ELSE
            nr_s=1
          ENDIF

          IF(CBBREV(CO,'FROM_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
      	  CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'Invalid source class',ERROR,*9999)

          CALL ASSERT(ITYP5(nr_s,nx_s).EQ.2.AND.ITYP2(nr_s,nx_s).EQ.9,
     '      'Source equation must be activation model',ERROR,*9999)

          IF(CBBREV(CO,'TO_REGION',4,noco+1,NTCO,N3CO)) THEN
            nr_t=IFROMC(CO(N3CO+1))
          ELSE
            nr_t=2
          ENDIF

          IF(CBBREV(CO,'TO_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_t.GT.0,'>>Invalid target class',ERROR,*9999)
          CALL ASSERT(ITYP5(nr_t,nx_t).EQ.1.AND.ITYP2(nr_t,nx_t).EQ.5,
     '      'Target equation must be div(k.grad(u))=f',ERROR,*9999)

        ELSE IF((TYPE(1:8).EQ.'GEOMETRY').OR.(TYPE(1:6).EQ.
     '    'METRIC')) THEN
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
        ELSE
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,nxc,CO,ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '       ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(ITYP4(nr,nx).EQ.4,
     '        'Must be collocation solution',ERROR,*9999)
          ENDDO
        ENDIF !type

        IF(TYPE(1:6).EQ.'SOURCE') THEN
          NOQUES=0
          FILEIP=.FALSE.
          IOTYPE=1
          
          FORMAT='('' Enter function to use as source term [1]:'''//
     '      '/''   (1) -2(pi)^2 sin(pi x) sin(pi y)'''//
     '      '/''   (2) Unused'''//
     '      '/''   (3) Unused'''//
     '      '/''   (4) Unused'''//
     '      '/''   (5) Unused'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).GT.0) THEN
            IFNTYP=IDATA(1)
          ELSE
            IFNTYP=1
          ENDIF
        ENDIF !source/bidomain

        nb=1
        DO WHILE (nb.LE.NBT.AND.NBC(nb).NE.7)
          nb=nb+1
        ENDDO
! MLB 4/7/97 Removed temporarily, not required for new grid stuff
C        CALL ASSERT((NBC(nb).EQ.7),
C     '    'Extended basis function not defined',ERROR,*9999)
! AJP 19/7/96 Moved into above loop
c        DO no_nrlist=1,NRLIST(0)
c          nr=NRLIST(no_nrlist)
c          CALL ASSERT(ITYP4(nr,nx).EQ.4,
c     '      'Must be collocation solution',ERROR,*9999)
c        ENDDO
        NITB=NIT(nb)

C--type---------------------------------------------------------------

        IF(TYPE(1:8).EQ.'MATERIAL') THEN
C *** Compute Dij, the diffusion matrix, at each grid point in terms of
C *** xi coordinates.  This becomes a full matrix generated from the
C *** diagonal diffusion tensor in material coordinates.
C *** Similarly for the bidomain model, compute the conductivity
C *** tensors.

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)

C Find basis function for current region
            nxielem=NIT(NBJ(1,NEELEM(1,nr)))
            nb_extended=1
            DO nbb=1,NBT
              IF((nxielem.NE.NIT(nb_extended)).OR.(NBC(nb_extended)
     '          .NE.7)) THEN
                nb_extended=nb_extended+1
              ENDIF
            ENDDO
            CALL ASSERT((NBC(nb_extended).EQ.7),
     '        'Extended basis function not defined',ERROR,*9999)
            NITB=NIT(nb_extended)            
            
C Initialise entire arrays to start with
            DO i=1,3
              DO j=1,3
                dNudXi(i,j)=0.d0
                dXidNu(i,j)=0.d0
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO k=1,4
                    PROPQ(i,j,k,1,nq,nx)=0.d0
                    PROPQ(i,j,k,2,nq,nx)=0.d0
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            
            IF(DOP) THEN
C$            call mp_setlock()
              WRITE(OP_STRING,'('' C(I)ij, C(E)ij:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$            call mp_unsetlock()
            ENDIF
            DO nq=NQR(1,nr),NQR(2,nr)
C *** First compute dnu/dxi = dnu/dx * dx/dxi, where both of these come
C *** from DXQ, computed in DEGRID.
            
              DO i=1,NITB
                DO j=1,NITB
!new PJH 11Jan96 
                  dNudXi(i,j)=0.d0
                  DO nj=1,NJT
                    dNudXi(i,j)=dNudXi(i,j)+
     '                          DNUDXQ(i,nj,nq)*DXDXIQ(nj,j,nq)
                  ENDDO
                ENDDO
              ENDDO
C *** and dxi/dnu = 1/(dnu/dxi)
              IF(NITB.EQ.1) THEN
                IF(dNudXi(1,1).GT.1.0D-13) dXidNu(1,1)=1.d0/dNudXi(1,1)
              ELSE
                CALL INVERT(NITB,dNudXi,dXidNu,DET)
              ENDIF
           
              IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.5) THEN
C *** For div(k.grad(u))=f the conductivities are in CQ(2..4,nq,nx)
                DO i=1,NITB
                  DO j=1,NITB
                    PROPQ(i,j,1,1,nq,nx)=
     '                CQ(2,nq,nx)*dXidNu(i,1)*dNudXi(1,j)
     '               +CQ(3,nq,nx)*dXidNu(i,2)*dNudXi(2,j)
                    IF(NITB.EQ.3) THEN
                      PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)
     '                 +CQ(4,nq,nx)*dXidNu(i,3)*dNudXi(3,j)
                    ENDIF
                  ENDDO !j
                ENDDO !i

              ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
C *** For Bidomain
C ***   Compute and store C(I)ij(nq) and C(E)ij in PROPQ 
C ***     (already initialised)
C ***   Intracellular conductivities in CQ(3;4;5,nq,nx)
C ***   Extracellular conductivities in CQ(6;7;8,nq,nx)
              
                DO i=1,NITB
                  DO j=1,NITB
                    DO k=1,NITB
                      PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)+
     '                  CQ(k+2,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                      PROPQ(i,j,1,2,nq,nx)=PROPQ(i,j,1,2,nq,nx)+
     '                  CQ(k+5,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                    ENDDO
                  ENDDO !j
                ENDDO !i
              ENDIF !ityp4
           
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' nq: '',I5)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            ((PROPQ(i,j,1,1,nq,nx),j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            ((PROPQ(i,j,1,2,nq,nx),j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$              call mp_unsetlock()
              ENDIF
            ENDDO !nq
           
C *** Now compute Dij,k using a first order finite difference about
C *** each grid point.  This only needs to be computed for grid points
C *** that are internal (not on external bdy).
C *** Similarly for the bidomain model, the derivatives of the
C *** conductivity tensors.

            IF(DOP) THEN
C$            call mp_setlock()
              WRITE(OP_STRING,'('' C(I)ij,k, C(E)ij,k:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$            call mp_unsetlock()
            ENDIF
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(NWQ(1,nq,1).EQ.0) THEN !internal g.p.
                DO k=1,NITB
                  DO i=1,NITB
                    DO j=1,NITB
                      IF(NXQ(k,1,nq,1).GT.0) THEN
                        PROPQ(i,j,k+1,1,nq,nx)=
     '                               PROPQ(i,j,1,1,NXQ( k,1,nq,1),nx)
     '                              -PROPQ(i,j,1,1,NXQ(-k,1,nq,1),nx)
                        PROPQ(i,j,k+1,2,nq,nx)=
     '                               PROPQ(i,j,1,2,NXQ( k,1,nq,1),nx)
     '                              -PROPQ(i,j,1,2,NXQ(-k,1,nq,1),nx)
                      ELSE
                        PROPQ(i,j,k+1,1,nq,nx)=0.d0
                        PROPQ(i,j,k+1,2,nq,nx)=0.d0
                      ENDIF
                    ENDDO !j
                  ENDDO !i
                ENDDO !k
                IF(DOP) THEN
C$                call mp_setlock()
                  WRITE(OP_STRING,'('' nq: '',I5)') nq
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(9(F12.5))')
     '           (((PROPQ(i,j,k,1,nq,nx),k=2,nitb+1),j=1,nitb),i=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(9(F12.5))')
     '           (((PROPQ(i,j,k,2,nq,nx),k=2,nitb+1),j=1,nitb),i=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                call mp_unsetlock()
                ENDIF
              ENDIF !nwq
            ENDDO !nq
          ENDDO !no_nrlist

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:22).EQ.'INTRACELLULAR_BIDOMAIN') THEN
C Update the extracellular potentials at intracellular grid points 
C either from extracellular grid 
C or by interpolating from values at fem node points

          IF(ITYP4(nr_s,nx_s).EQ.4) THEN !source class is grid
            DO nq=1,NQT !Update extracellular potential
              YQ(nq,2,1,nx_t)=YQ(nq,1,1,nx_s)
            ENDDO
          
          ELSE IF(DIRECT) THEN         !source is grid mesh
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr_s,nx_s),
     '        NKH(1,1,1,nr_s),NPNODE,nr_s,NVHP(1,1,1,nr_s),nx_s,
     '        NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
            DO nq=1,NQT !Update extracellular potential
              YQ(nq,2,1,nx_t)=ZP(1,1,1,NPQ(nq),1)
            ENDDO
          
          ELSE                         !source class is fem mesh
            nb_extended=1
            DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
              nb_extended=nb_extended+1
            ENDDO
            CALL ASSERT((NBC(nb_extended).EQ.7),
     '        'Extended basis function not defined',ERROR,*9999)
          
! Check we have a 9x9(x9) grid point scheme
C GBS 3-July-96 Why????
C           NINE=.TRUE.
C           DO ni=1,NIT(nb_extended)
C             NINE=NINE.AND.(NGAP(ni,nb_extended).EQ.9)
C           ENDDO
C           CALL ASSERT(NINE,'Require 9x9(x9) grid',ERROR,*9999)
            
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr_s,nx_s),
     '        NKH(1,1,1,nr_s),NPNODE,nr_s,NVHP(1,1,1,nr_s),nx_s,
     '        NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
            DO noelem=1,NEELEM(0,nr_s)
              ne=NEELEM(noelem,nr_s)
              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx_s),NKE(1,1,1,ne),
     '          NPF,NPNE(1,1,ne),nr_s,NVHE(1,1,1,ne),NW(ne,1),nx_s,
     '          SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              DO ng=1,NGT(nb_extended)
                nq=NQGE(ng,ne,nb_extended)
                SUM=0.0d0
                DO ns=1,NST(nb_extended)
                  PGG=PG(ns,1,ng,nb_extended)
                  SUM=SUM+PGG*ZE(ns,1)
                ENDDO !ns
                YQ(nq,2,1,nx_t)=SUM !Update extracellular potential
              ENDDO !ng
            ENDDO !noelem
          ENDIF !ityp4
            
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN
C Update the RHS source term for the extracellular bidomain Poisson's
C eqn (class nx_t) by computing div(kgrad(phi_m)) at the extracellular 
C grid pts from a quad patch of intracellular grid pts on class nx_s.

C *** Loop over grid points
          DO nq=1,NQT
            DO nik=-1,1
              DO nij=-1,1
                DO nii=-1,1
                  PHIMQ(nii,nij,nik)=0.d0
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
          
C ***  Formulate a local quadratic element about nq defined by UMQ,
C ***  storing the value of phi(m) at each node point mq.
            IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
            IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                  IF(mq.GT.0) THEN
                    PHIMQ(nii,nij,nik)=YQ(mq,1,1,nx_s)
                  ELSE
                    PHIMQ(nii,nij,nik)=0.0d0 
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
            
C *** Compute PHI,k by 1st order finite differences about nq.
            dPHIdX(1)=PHIMQ(1,0,0)-PHIMQ(-1,0,0)
            dPHIdX(2)=PHIMQ(0,1,0)-PHIMQ(0,-1,0)
            dPHIdX(3)=PHIMQ(0,0,1)-PHIMQ(0,0,-1)

C *** Compute PHI,jk by 2nd order finite differences about nq
            DO ib=1,NITB
              DO jb=1,NITB
                IF(ib.EQ.jb) THEN
                  IF(ib.EQ.1) THEN !PHI,11
                    d2PHIdX2(ib,jb)=(PHIMQ(1,0,0)-
     '                2.0d0*PHIMQ(0,0,0)+PHIMQ(-1,0,0))*4.0d0
                  ELSE IF(ib.EQ.2) THEN !PHI,22
                    d2PHIdX2(ib,jb)=(PHIMQ(0,1,0)
     '                -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,-1,0))*4.0d0
                  ELSE IF(ib.EQ.3) THEN !PHI,33
                    d2PHIdX2(ib,jb)=(PHIMQ(0,0,1)
     '                -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,0,-1))*4.0d0
                  ENDIF               
                ELSE
                  IF(ib+jb.EQ.3) THEN !PHI,12 and PHI,21
                    d2PHIdX2(ib,jb)=PHIMQ(1,1,0)-PHIMQ(1,-1,0)
     '                -PHIMQ(-1,1,0)+PHIMQ(-1,-1,0)
                  ELSE IF(ib+jb.EQ.4) THEN !PHI,13 and PHI,31
                    d2PHIdX2(ib,jb)=PHIMQ(1,0,1)-PHIMQ(1,0,-1)
     '                -PHIMQ(-1,0,1)+PHIMQ(-1,0,-1)
                  ELSE IF(ib+jb.EQ.5) THEN !PHI,23 and PHI,32
                    d2PHIdX2(ib,jb)=PHIMQ(0,1,1)-PHIMQ(0,-1,1)
     '                -PHIMQ(0,1,-1)+PHIMQ(0,-1,-1)
                  ENDIF
                ENDIF !ib=jb
              ENDDO !jb
            ENDDO !ib

C *** Compute the diffusion in three parts.
            DIFFN=0.0d0
            DO nii=1,NITB
              DO nij=1,NITB           
                SUM1=0.0d0
                SUM2=0.0d0
                DO nik=1,NITB
                  SUM1=SUM1
     '                +PROPQ(nik,nii,nij+1,1,nq,nx_s)*dPHIdX(nik)
                  SUM2=SUM2
     '                +PROPQ(nik,nii,1,1,nq,nx_s)*d2PHIdX2(nij,nik)
                ENDDO
                DIFFN=DIFFN+(SUM1+SUM2)*GUQ(nii,nij,nq)
              ENDDO !nij
            ENDDO !nii
            DO nii=1,NITB
              SUM1=0.0d0
              DO nij=1,NITB
                SUM1=SUM1+PROPQ(nij,nii,1,1,nq,nx_s)*dPHIdX(nij)
              ENDDO
              DIFFN=DIFFN-SUM1*GCHQ(nii,nq)
            ENDDO !nii
            CQ(1,nq,nx_t)=-DIFFN !is source term for Poisson's eqn
          ENDDO !nqt

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:6).EQ.'STRAIN') THEN
            Esac=-20.d0 !(mV) reversal potential for SAC
          
          nb_extended=1
          DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
            nb_extended=nb_extended+1
          ENDDO
          CALL ASSERT((NBC(nb_extended).EQ.7),
     '      'Extended basis function not defined',ERROR,*9999)
          
C Generate inward Isac current from stretch if above threshold
C Do interior collocation pts first, then bdry pts
          
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(DOP) WRITE(*,'(/'' Interior pts: ne='',I5)') ne
          
              NQ_row=9
              NG_row=3
              
              DO j=2,NQ_row-1    !loop over collocation points
                DO i=2,NQ_row-1
                  nq=NQGE(i+(j-1)*NQ_row,ne,nb_extended) !is coll.n pt#
              
                  mLL=INT(REAL(i+1)/3.0)  !m(left   Gauss pt index)
                  nBB=INT(REAL(j+1)/3.0)  !n(bottom Gauss pt index)
                  mRR=min(mLL+1,NG_row)    !m(right  Gauss pt index)
                  nTT=min(nBB+1,NG_row)    !n(top    Gauss pt index)
              
                  ngBL=mLL+(nBB-1)*NG_row  !ng(bottom left  Gauss pt#)
                  ngBR=mRR+(nBB-1)*NG_row  !ng(bottom right Gauss pt#)
                  ngTL=mLL+(nTT-1)*NG_row  !ng(top left     Gauss pt#)
                  ngTR=mRR+(nTT-1)*NG_row  !ng(top right    Gauss pt#)
              
!                 WRITE(*,'(/'' i='',I1,'' j='',I1,'' nq='',I2)') i,j,nq
!                 WRITE(*,'('' mLL='',I1,'' nBB='',I1)') mLL,nBB
!                 WRITE(*,'('' ngBL='',I1,'' ngBR='',I1,'' ngTL='',I1,'
!    '              //''' ngTR='',I1)') ngBL,ngBR,ngTL,ngTR
              
                  iL=3*mLL-1  !locates bottom left Gauss pt index
                  jB=3*nBB-1  ! in the collocation grid
              
                  Xi1=DBLE(i-iL)/DBLE(NG_row) !how far ij colloc.n pt
                  Xi2=DBLE(j-jB)/DBLE(NG_row) ! is betw adj Gauss pts
              
                  Ext=(1.d0-Xi1)*(1.d0-Xi2)*FEXT(1,ngBL,ne) 
     '               +      Xi1 *(1.d0-Xi2)*FEXT(1,ngBR,ne)
     '               +(1.d0-Xi1)*      Xi2 *FEXT(1,ngTL,ne)
     '               +      Xi1 *      Xi2 *FEXT(1,ngTR,ne)
                  IF(DOP) THEN
C$                  call mp_setlock()
                    WRITE(*,'(/'' Xi1='',F5.2,'' Xi2='',F5.2)') Xi1,Xi2
                    WRITE(*,'('' Exten. ratio at nq='',I3,'
     '               //''' :'',F6.3)') nq,Ext
C$                  call mp_unsetlock()
                  ENDIF
          
C PJH 30jun96 - nx added. Check with Greg
                  Current=-(YQ(nq,1,1,nx)-Esac)*CQ(18,nq,nx)*(Ext-1.d0)
                  IF(Current.GT.100.d0) THEN !above threshold
                    YQ(nq,5,1,nx)=Current
                  ELSE
                    YQ(nq,5,1,nx)=0.d0
                  ENDIF
              
                ENDDO !i
              ENDDO !j
            ENDDO !noelem
          
C         Boundary collocation pts
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(DOP) WRITE(*,'(/'' Boundary pts: ne='',I5)') ne
          
C         RH bdry of element ne
              DO j=1,9 
                nq =NQGE(9*j,ne,nb_extended)
                nq1=NXQ(-1,1,nq,1)
                nq2=NXQ( 1,1,nq,1)
                IF(DOP) THEN
C$                call mp_setlock()
                  WRITE(*,'(/'' nq='',I3,'' nq1='',I3,'' nq2='',I3)')
     '              nq,nq1,nq2
C$                call mp_unsetlock()
                ENDIF
                IF(nq2.GT.0) THEN !RH point exists
                  Current=0.5d0*(YQ(nq1,5,1,nx)+YQ(nq2,5,1,nx))
                  IF(Current.GT.100.d0) THEN !above threshold
                    YQ(nq,5,1,nx)=Current
                  ELSE
                    YQ(nq,5,1,nx)=0.d0
                  ENDIF
                ENDIF
              ENDDO !j
          
C         TOP bdry of element ne
              DO i=1,9 
                nq =NQGE(72+i,ne,nb_extended)
                nq1=NXQ(-2,1,nq,1)
                nq2=NXQ( 2,1,nq,1)
                IF(DOP) THEN
C$                call mp_setlock()
                  WRITE(*,'(/'' nq='',I3,'' nq1='',I3,'' nq2='',I3)')
     '              nq,nq1,nq2
C$                call mp_unsetlock()
                ENDIF
                IF(nq2.GT.0) THEN !TOP point exists
                  Current=0.5d0*(YQ(nq1,5,1,nx)+YQ(nq2,5,1,nx))
                  IF(Current.GT.100.d0) THEN !above threshold
                    YQ(nq,5,1,nx)=Current
                  ELSE
                    YQ(nq,5,1,nx)=0.d0
                  ENDIF
                ENDIF
              ENDDO !i
            ENDDO !noelem
          ENDDO !no_nrlist
          
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:4).EQ.'TIME') THEN
C Update previous time step solution for multigrid transient heat eqtn 
          DO nq=1,NQT
            YQ(nq,6,1,nx)=YQ(nq,1,1,nx)
          ENDDO !nq
          
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:6).EQ.'SOURCE') THEN
C Calculate source term for collocation grid
          DO nq=1,NQT
            GX=XQ(1,nq)
            GY=XQ(2,nq)
            IF(IFNTYP.EQ.1) THEN
              CQ(1,nq,nx)=-2.d0*PI*PI*DSIN(PI*GX)*DSIN(PI*GY)
            ENDIF !ifntyp
          ENDDO !nq

C--type---------------------------------------------------------------
C*** Martin Buist, June 1997

        ELSE IF(TYPE(1:8).EQ.'GEOMETRY') THEN
! Calculate geometric coordinates of grid points
          CALL ASSERT(CALL_GRID,'>>No grid defined',ERROR,*9999)
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            DO nee=1,NEELEM(0,nr)
              ne=NEELEM(nee,nr)
              SCHEME=NQS(ne)
! Map global parametres to local element parameters
              CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          NQE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,
     '          XP,ERROR,*9999)
              II=MAX(1,NQXI(1,SCHEME))
              IJ=MAX(1,NQXI(2,SCHEME))
              IK=MAX(1,NQXI(3,SCHEME))
              DO i=1,3
                XQ_TEMP(i)=0.0d0
                XI(i)=0.0d0
              ENDDO !i
! Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME))+((nik-1)*
     '                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq)
! Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN 
! Local xi coordinates of grid point nq in element ne
                      IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                      IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                      IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
! Use basis function interpolation to get geometric position
                      nb=NQSCNB(SCHEME)
                      DO nj=1,NJT
C                        nb=NBJ(nj,ne)
                        XQ_TEMP(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,1,XI,XE(1,nj))
                      ENDDO !nj
! Change coordinates - all grid points in rectangular cartesian
                      CALL XZ(ITYP10(nr),XQ_TEMP,XQ(1,nq))
                    ENDIF
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
            ENDDO !element
          ENDDO !region
          IF(DOP) THEN
            DO nq=1,NQT
              WRITE(OP_STRING,'(''XQ'',3F10.4)') XQ(1,nq),XQ(2,nq),
     '          XQ(3,nq)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

C--type---------------------------------------------------------------
C*** Martin Buist, June 1997

        ELSE IF(TYPE(1:8).EQ.'METRIC') THEN
! Calculate metric tensor information for grid
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            DO i=1,3
              DO j=1,3
                GU(i,j)=0.0d0
                GL(i,j)=0.0d0
              ENDDO !j
            ENDDO !i
            DO nee=1,NEELEM(0,nr)
              ne=NEELEM(nee,nr)
              SCHEME=NQS(ne)
! Map global parametres to local element parameters
              CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          NQE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,
     '          XP,ERROR,*9999)
! dxi/dxi is the relationship between material space and the local 
!   quadratic element space. dxi(quad)/dxi(mate)
              DO ni=1,NQXI(0,SCHEME)
                DXIDXI(ni)=2.0d0/DBLE(NQXI(ni,SCHEME)-1)
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''dxidxi(ni)'',I3,F8.6)') ni,
     '              DXIDXI(ni)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !ni
              II=MAX(1,NQXI(1,SCHEME))
              IJ=MAX(1,NQXI(2,SCHEME))
              IK=MAX(1,NQXI(3,SCHEME))
              DO i=1,3
                XI(i)=0.0d0
              ENDDO !i
! Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME))+((nik-1)*
     '                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq)
! Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN 
! Local xi coordinates of grid point nq in element ne
                      IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                      IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                      IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)
                      CALL XEXW(IBT,IDO,INP,NAN,NBJ,nr,XE,XQD,XI,
     '                  ERROR,*9999)
! Gi (lower) = dx(nj)/dxi(ni) in direction nj
                      DO ni=1,NQXI(0,SCHEME)
                        CALL XZ_DERIV(ITYP10(nr),NU1(ni),XQD,XQDRC)
                        DO nj=1,NJT
                          DXDXIQ(nj,ni,nq)=XQDRC(nj,NU1(ni))*DXIDXI(ni)
                          IF(DOP) THEN
                            WRITE(OP_STRING,'(''dxdxiq'',3I3,F8.6)')
     '                        nj,ni,nq,DXDXIQ(nj,ni,nq)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDDO !nj
                      ENDDO !ni
! Calculate Gij (lower) and then Gij (upper) by inversion 
! Gij (lower) = dx(nj)/dxi(nii)*dx(nj)/dxi(nij)
                      DO nij2=1,NQXI(0,SCHEME)
                        DO nii2=1,NQXI(0,SCHEME)
                          GL(nii2,nij2)=0.0d0
                          DO nj=1,NJT
                            GL(nii2,nij2)=GL(nii2,nij2)+(DXDXIQ(nj,
     '                        nii2,nq)*DXDXIQ(nj,nij2,nq))
                          ENDDO !nj
                        ENDDO !nii2
                      ENDDO !nij2
! Gij (upper) = inv(Gij lower)
                      IF(NQXI(0,SCHEME).GT.1) THEN
                        CALL INVERT(NQXI(0,SCHEME),GL,GU,DET)
                      ELSE IF(GL(1,1).NE.0.0D0) THEN
                        GU(1,1)=1.0D0/GL(1,1)
                      ELSE
                        WRITE(OP_STRING,
     '                    '(''>Zero value for GL - cannot invert'')') 
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      ENDIF
! Write Gij (upper) at grid point nq into GUQ 
                      DO nij2=1,NQXI(0,SCHEME)
                        DO nii2=1,NQXI(0,SCHEME)
                          GUQ(nii2,nij2,nq)=GU(nii2,nij2)
                          IF(DOP) THEN
                            WRITE(OP_STRING,'(''guq'',3I3,F8.6)')
     '                        nii2,nij2,nq,GUQ(nii2,nij2,nq)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDDO !nij2
                      ENDDO !nii2
! Calculate the Christoffel symbol (CHTOFF)
                      CALL TOFFEL(NBJ(1,ne),NJE(ne),nr,CHTOFF,DBM,GU,
     '                  XQD,X3G,.FALSE.,ERROR,*9999)
! Calculate CHTOFF(k,i,j)*Gij (upper) = GCHQ(k)
                      DO nik2=1,NQXI(0,SCHEME)
                        SUM=0.0d0
                        DO nii2=1,NQXI(0,SCHEME)
                          DO nij2=1,NQXI(0,SCHEME)
                            SUM=SUM+(CHTOFF(nik2,nii2,nij2)*
     '                        GU(nii2,nij2))*DXIDXI(nii2)*DXIDXI(nij2)
                          ENDDO !nij2
                        ENDDO !nii2
                        GCHQ(nik,nq)=SUM*DXIDXI(nik)
                        IF(DOP) THEN 
                          WRITE(OP_STRING,'(''gchq'',2I6,F12.6)') 
     '                      nik,nq,GCHQ(nik,nq)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDDO !nik2
! Calculate dnu/dx (direction cosines of material coords)
                      CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '                  A_VECTOR,B_VECTOR,C_VECTOR,XE,XQD,XI,
     '                  ERROR,*9999)
                      DO nj=1,NJT
                        DNUDXQ(1,nj,nq)=A_VECTOR(nj)
                        IF(nj.GT.1) DNUDXQ(2,nj,nq)=B_VECTOR(nj)
                        IF(nj.GT.2) DNUDXQ(3,nj,nq)=C_VECTOR(nj)
                        IF(DOP) THEN
                          WRITE(OP_STRING,'(''dnu/dxq'',3F12.6)') 
     '                      DNUDXQ(1,nj,nq),DNUDXQ(2,nj,nq),
     '                      DNUDXQ(3,nj,nq)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDDO !nj
                    ENDIF
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
            ENDDO !element
          ENDDO !region

        ENDIF !type
      ENDIF

      CALL EXITS('UPGRID')
      RETURN
 9999 CALL ERRORS('UPGRID',ERROR)
      CALL EXITS('UPGRID')
      RETURN 1
      END


Module FE22
=========== 
C KAT 2001-12-14
      SUBROUTINE CAELEC(ISEG,ISELEC,NDDL,NDLT,NEELEM,STRING,ERROR,*)

C#### Subroutine: CAELEC
C###  Description:
C###    CAELEC cancels electrode segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISELEC(128),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),nd,nde,ne,ne1,noiw,nr,NTIW

      CALL ENTERS('CAELEC',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel electrodes;s
C###  Description:
C###    Cancel electrode segment. This command removes the
C###    electrodes from the screen but does not remove them
C###    from the problem.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAELEC',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,4,NTIW,noco,NTCO,CO,ERROR,*9999)

        nr=1
        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO ne1=1,NEELEM(0,nr)
              ne=NEELEM(ne1,nr)
              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)
                IF(ISELEC(nd).GT.0) THEN
                  CALL DELETE_SEGMENT(ISELEC(nd),ISEG,iw,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAELEC')
      RETURN
 9999 CALL ERRORS('CAELEC',ERROR)
      CALL EXITS('CAELEC')
      RETURN 1
      END


C KAT 2001-12-14
      SUBROUTINE HIELEC(ISEG,ISELEC,NDDL,NDLT,NEELEM,STRING,ERROR,*)

C#### Subroutine: HIELEC
C###  Description:
C###    HIELEC hides electrode segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISELEC(128),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),nd,nde,ne,ne1,noiw,nr,NTIW

      CALL ENTERS('HIELEC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide electrodes 
C###  Description:
C###    Hide electrodes on specified workstation.
C###  Parameter:    <on WS#[4]>
C###    Specify which window to hide the electrodes on. The
C###    default is to hide on window 4.

        OP_STRING(1)=STRING(1:IEND) //' <on WS#[4]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIELEC',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,4,NTIW,noco,NTCO,CO,ERROR,*9999)

        nr=1
        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO ne1=1,NEELEM(0,nr)
              ne=NEELEM(ne1,nr)
              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)
                IF(ISELEC(nd).GT.0) THEN
                  CALL VISIB(iw,ISEG,ISELEC(nd),'INVISIBLE',ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIELEC')
      RETURN
 9999 CALL ERRORS('HIELEC',ERROR)
      CALL EXITS('HIELEC')
      RETURN 1
      END


C KAT 2001-12-14
      SUBROUTINE SHELEC(ISEG,ISELEC,NDDL,NDLT,NEELEM,STRING,ERROR,*)

C#### Subroutine: SHELEC
C###  Description:
C###    SHELEC shows electrode segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISELEC(128),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),nd,nde,ne,ne1,nr,NTIW

      CALL ENTERS('SHELEC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show electrodes 
C###  Description:
C###    Make the electrodes visible on the specified workstation.
C###  Parameter:    <on WS#[4]>
C###    Specify which window to show the electrodes on. The
C###    default is to show electrodes on window 4.

        OP_STRING(1)=STRING(1:IEND) //' <on WS#[4]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHELEC',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        nr=1
        iw=1 !AJP 7/10/95 What should it be?
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO ne1=1,NEELEM(0,nr)
            ne=NEELEM(ne1,nr)
            DO nde=1,NDLT(ne)
              nd=NDDL(ne,nde)
              IF(ISEG(ISELEC(nd)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISELEC(nd),'VISIBLE',ERROR,*9999)
              ELSE IF(ISEG(ISELEC(nd)).EQ.0) THEN
                WRITE(OP_STRING,'('' >>Electrode number '',I4,'
     '            //''' is not defined on '',I1)') ne,iw
      		CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('SHELEC')
      RETURN
 9999 CALL ERRORS('SHELEC',ERROR)
      CALL EXITS('SHELEC')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE CAIMAG(ISEG,ISIMAG,STRING,ERROR,*)

C#### Subroutine: CAIMAG
C###  Description:
C###    CAIMAG cancels image segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISIMAG(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw

      CALL ENTERS('CAIMAG',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel image;s
C###  Description:
C###    Cancel image display.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAIMAG',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        DO iw=1,2
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISIMAG(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISIMAG(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('CAIMAG')
      RETURN
 9999 CALL ERRORS('CAIMAG',ERROR)
      CALL EXITS('CAIMAG')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE HIIMAG(ISEG,ISIMAG,STRING,ERROR,*)

C#### Subroutine: HIIMAG
C###  Description:
C###    HIIMAG hides image segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISIMAG(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('HIIMAG',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide image 
C###  Description:
C###    Hide image segments.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND)//' <on WS#s[1,2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIIMAG',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISIMAG(iw)).EQ.2) THEN
              CALL VISIB(iw,ISEG,ISIMAG(iw),'INVISIBLE',ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIIMAG')
      RETURN
 9999 CALL ERRORS('HIIMAG',ERROR)
      CALL EXITS('HIIMAG')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE SHIMAG(ISEG,ISIMAG,STRING,ERROR,*)

C#### Subroutine: SHIMAG
C###  Description:
C###    SHIMAG shows image segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISIMAG(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('SHIMAG',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C---------------------------------------------------------------------

C#### Command: FEM show image
C###  Description:
C###    Make the image visible.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s).
C###    The all keyword specifies all currently defined workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHIMAG',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISEG(ISIMAG(iw)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISIMAG(iw),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISIMAG(iw)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Image is not '
     '          //'defined on '',I1)') iw
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHIMAG')
      RETURN
 9999 CALL ERRORS('SHIMAG',ERROR)
      CALL EXITS('SHIMAG')
      RETURN 1
      END


      SUBROUTINE SHPMAR(ISEG,ISPMAR,STRING,ERROR,*)

C#### Subroutine: SHPMAR
C###  Description:
C###    SHPMAR shows polymarker segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISPMAR(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('SHPMAR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show polymarker 
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the 
C###    polymarker on. 
C###  Description:
C###    Make the polymarker segments visible on specified workstations.
        
        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHPMAR',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            IF(ISEG(ISPMAR(iw)).EQ.1) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              CALL VISIB(iw,ISEG,ISPMAR(iw),'VISIBLE',ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ELSE IF(ISEG(ISPMAR(iw)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>PMAR are not '
     '          //'defined on '',I1)') iw
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHPMAR')
      RETURN
 9999 CALL ERRORS('SHPMAR',ERROR)
      CALL EXITS('SHPMAR')
      RETURN 1
      END


C KAT 2001-12-17
      SUBROUTINE CATEXT(ISEG,ISTEXT,STRING,ERROR,*)

C#### Subroutine: CATEXT
C###  Description:
C###    CATEXT cancels text segments.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
!     Parameter List
      INTEGER ISTEXT(NWM),ISEG(*)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('CATEXT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel text;s 
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to cancel the 
C###    text on. 
C###  Description:
C###    Cancel text segments on specified workstations.

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
       
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CATEXT',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            IF(ISTEXT(iw).GT.0) THEN
              CALL DELETE_SEGMENT(ISTEXT(iw),ISEG,iw,ERROR,*9999)
            ENDIF
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
        NTTEXT=0
      ENDIF

      CALL EXITS('CATEXT')
      RETURN
 9999 CALL ERRORS('CATEXT',ERROR)
      CALL EXITS('CATEXT')
      RETURN 1
      END


C KAT 2001-12-17
      SUBROUTINE HITEXT(ISEG,ISTEXT,STRING,ERROR,*)

C#### Subroutine: HITEXT
C###  Description:
C###    HITEXY hides text segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISTEXT(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('HITEXT',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide text
C###  Parameter:    <on (all/WS#s)[all]> 
C###    Specify the workstation (GX window) to hide the 
C###    text on. 
C###  Description:
C###    Hide text segment on specified workstations.

        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
    
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HITEXT',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            IF(ISEG(ISTEXT(iw)).EQ.2) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              CALL VISIB(iw,ISEG,ISTEXT(iw),'INVISIBLE',ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HITEXT')
      RETURN
 9999 CALL ERRORS('HITEXT',ERROR)
      CALL EXITS('HITEXT')
      RETURN 1
      END


C KAT 2001-12-17
      SUBROUTINE SHTEXT(ISEG,ISTEXT,STRING,ERROR,*)

C#### Subroutine: SHTEXT
C###  Description:
C###    SHTEXT shows text segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISTEXT(NWM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),noiw,NTIW

      CALL ENTERS('SHTEXT',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show text
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to show the 
C###    text on. 
C###  Description:
C###    Make the text segments visible on the specified workstations.
  
        OP_STRING(1)=STRING(1:IEND) //' <on (all/WS#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

        OP_STRING(1)=STRING(1:IEND)
     '    //' <on WSS>[on all]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHTEXT',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            IF(ISEG(ISTEXT(iw)).EQ.1) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              CALL VISIB(iw,ISEG,ISTEXT(iw),'VISIBLE',ERROR,*9999)
              CALL DAWK(iw,1,ERROR,*9999)
            ELSE IF(ISEG(ISTEXT(iw)).EQ.0) THEN
              WRITE(OP_STRING,'('' >>Text is not '
     '          //'defined on '',I1)') iw
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SHTEXT')
      RETURN
 9999 CALL ERRORS('SHTEXT',ERROR)
      CALL EXITS('SHTEXT')
      RETURN 1
      END


C KAT 2001-12-14
      SUBROUTINE CATRAC(ISEG,ISSIGN,ISTRAC,STRING,ERROR,*)

C#### Subroutine: CATRAC
C###  Description:
C###    CATRAC cancels signal segments and window.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      
!     Parameter List
      INTEGER ISEG(*),ISSIGN(2,128),ISTRAC(10)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IEND,iw,nosign

      CALL ENTERS('CATRAC',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C---------------------------------------------------------------------

C#### Command: FEM cancel trace;s
C###  Description:
C###    Cancel signal segments and window.

        OP_STRING(1)=STRING(1:IEND)//';s'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
 
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CATRAC',ERROR,*9999)
      ELSE
        CALL CHECKQ('S',noco,1,CO,COQU,STRING,*1)
        iw=34
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          IF(ISSIGN(1,1).GT.0) THEN
            DO nosign=1,128
              CALL DELETE_SEGMENT(ISSIGN(1,nosign),ISEG,iw,ERROR,*9999)
              CALL DELETE_SEGMENT(ISSIGN(2,nosign),ISEG,iw,ERROR,*9999)
            ENDDO
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
          CALL CLOSE_WS(iw,ERROR,*9999)
          IWKS(iw)=0
        ENDIF
        iw=35
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO i=1,10
            IF(ISTRAC(i).GT.0) THEN
              CALL DELETE_SEGMENT(ISTRAC(i),ISEG,iw,ERROR,*9999)
            ENDIF
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
          CALL CLOSE_WS(iw,ERROR,*9999)
          IWKS(iw)=0
        ENDIF
c       IF(IWKS(36).GT.0) CALL INPUT_MODE(36,1,'CHOICE','REQUEST',
c    '    ERROR,*9999)
C        FIDUCALC=.FALSE.

      ENDIF

      CALL EXITS('CATRAC')
      RETURN
 9999 CALL ERRORS('CATRAC',ERROR)
      CALL EXITS('CATRAC')
      RETURN 1
      END


      SUBROUTINE CAUSER(STRING,ERROR,*)

C#### Subroutine: CAUSER
C###  Description:
C###    CAUSER cancels user defined names.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:user00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND

      CALL ENTERS('CAUSER',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel user
C###  Description:
C###    Cancel user segment.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAUSER',ERROR,*9999)
      ELSE
        NTUS=0
      ENDIF

      CALL EXITS('CAUSER')
      RETURN
 9999 CALL ERRORS('CAUSER',ERROR)
      CALL EXITS('CAUSER')
      RETURN 1
      END



Module FE23
=========== 


      SUBROUTINE DRAW(STRING,NOCO,CO,ERROR,*)

C**** Draws polyline on workstation.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER NOCO
      CHARACTER STRING*(MXCH),CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER I,IBEG,IEND,IFROMC,IW,NTITLV(20),NTLV,NTRL
      REAL*8 RL(3,100),RLX(100),RLY(100)
      LOGICAL ABBREV

      CALL ENTERS('DRAW',*9999)
      IF(CO(NOCO+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        OP_STRING(1)=STRING(IBEG:IEND)
     '    //' X_LIST Y_LIST'
     '    //' <on WS>[on 1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
      ELSE
        CALL PARSTR(CO(NOCO+1),20,NTLV,NTITLV,100,RLX,ERROR,*9999)
        NTRL=NTITLV(1)
        CALL PARSTR(CO(NOCO+2),20,NTLV,NTITLV,100,RLY,ERROR,*9999)
        IF(ABBREV(CO(NOCO+3),'ON',2)) THEN
          IW=IFROMC(CO(NOCO+4))
        ELSE
          IW=1
        ENDIF
        DO I=1,NTRL
          RL(1,I)=RLX(I)
          RL(2,I)=RLY(I)
        ENDDO
        CALL ACWK(IW,0,ERROR,*9999)
        CALL POLYLINE(2,40,2,RL,ERROR,*9999)
        CALL DAWK(IW,0,ERROR,*9999)
      ENDIF

 9998 CALL EXITS('DRAW')
      RETURN
 9999 CALL ERRORS('DRAW',ERROR)
      CALL EXITS('DRAW')
      RETURN 1
      END


      SUBROUTINE PLOT(NOCO,NTCO,NTCOQU,CO,COQU,ERROR,*)

C**** Plots vectors defined in CO and COQU.
C**** Valid commands:
C**** Plot line   list        on workstation
C**** Plot marker   "          "      "
C**** Plot line   list1;list2  "      "

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi10.cmn'
!     Parameter List
      INTEGER NOCO,NTCO,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IFROMC,IW,NORL,NTITLV(20),NTLV,NTRL
      REAL*8 RL1(100),RL2(100),RL3(3,100)
      LOGICAL ABBREV

      CALL ENTERS('PLOT',*9999)
      CALL PARSTR(CO(NOCO+2),20,NTLV,NTITLV,100,RL2,ERROR,*9999)
      NTRL=NTITLV(1)
      IF(NTCOQU(NOCO+2).EQ.0) THEN
        DO 100 NORL=1,NTRL
          RL1(NORL)=NORL/NTRL
 100    CONTINUE
      ELSE
        CALL PARSTR(COQU(NOCO+2,1),20,NTLV,NTITLV,100,RL1,ERROR,*9999)
      ENDIF
      IF(ABBREV(CO(NOCO+3),'ON',2)) THEN
        IW=IFROMC(CO(NOCO+4))
      ELSE
        IW=1
      ENDIF
      DO NORL=1,NTRL
        RL3(1,NORL)=RL1(NORL)
        RL3(2,NORL)=RL2(NORL)
      ENDDO

      CALL ACWK(IW,0,ERROR,*9999)
      IF(ABBREV(CO(NOCO+1),'LINE',1)) THEN
        CALL POLYLINE(1,IW,NTRL,RL3,ERROR,*9999)
      ELSE IF(ABBREV(CO(NOCO+1),'MARKER',1)) THEN
        CALL POLYMARKER(1,IW,NTRL,RL3,ERROR,*9999)
      ENDIF
      CALL DAWK(IW,0,ERROR,*9999)

 9998 CALL EXITS('PLOT')
      RETURN
 9999 CALL ERRORS('PLOT',ERROR)
      CALL EXITS('PLOT')
      RETURN 1
      END


      SUBROUTINE SG_ROW_NUMBER(ISEG,IS_ROW_NUMBER,CSEG,
     '  ERROR,*)

C**** Defines spreadsheet row segment.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:gen000.cmn'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),IS_ROW_NUMBER
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX,INDEX_OLD

      CALL ENTERS('SG_ROW_NUMBER',*9999)

      CALL OPEN_SEGMENT(IS_ROW_NUMBER,ISEG,50,'Row no.s ',INDEX,
     '  INDEX_OLD,1,1,CSEG,ERROR,*9999)

      CALL CLOSE_SEGMENT(IS_ROW_NUMBER,50,ERROR,*9999)

      CALL EXITS('SG_ROW_NUMBER')
      RETURN
 9999 CALL ERRORS('SG_ROW_NUMBER',ERROR)
      CALL EXITS('SG_ROW_NUMBER')
      RETURN 1
      END


Module FE24
===========

C Removed from
C#### Subroutine: DEDATA
C###  Description:
C###    DEDATA defines geometric data points.

C 9 Feb 01

C#### Command: FEM define data;c signal
C###  Description:
C###    Calculate the fidual marker times as data point values from
C###    the current signal (OBSELETE).
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';c signal'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C KAT 2001-12-13
C#### Command: FEM define data;c image
C###  Description:
C###    Calculate data points from the current image.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to associate data with.

        OP_STRING(1)=STRING(1:IEND)//';c image'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

          ELSE IF(CBBREV(CO,'IMAGE',2,noco+1,NTCO,N3CO)) THEN
            TYPE='IMAGE'

          ELSE IF(TYPE(1:5).EQ.'IMAGE') THEN
!         Reduces pixel data to those inside &/or near defined elements
!         (biquadratic elements only)
            CALL CALC_IMAGE_DATA(NBJ,NELIST,NKJE,NPF,NPNE,NRE,NVJE,
     '        SE,XA,XE,XP,ZD,ERROR,*9999) 

          ELSE IF(TYPE(1:6).EQ.'SIGNAL') THEN
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Transferring signal data'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO i=1,107
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Signal '',I3,'' is '',L5,'
     '            //''' with value '',F12.4)') i,ACCEPT(i),FIDMARK(i)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              ZD(4,i)=(FIDMARK(i)-DATUM)*1000.0D0
              IF(ACCEPT(i)) THEN
                WD(4,i)=1.0D0
              ELSE
                WD(4,i)=0.0D0
              ENDIF
            ENDDO

C KAT 2001-12-17 from after reading .pts map3d data file
                 CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//
     '            '.channels',STATUS,'SEQUEN','FORMATTED',
     '            132,ERROR,*9999)
C If local subdivion used then first number is not the total
C number of points.  AJP 10-3-94
C Could read in NT_SUB_DIV and compare with current value. However, this is
C only for files that have been created by cmiss.
                READ(IFILE,*) NT_CHANNELS
C                  IF(nt_channels.NE.NDT) THEN
C                    ERROR='>># Channels ne # points'
C                    GOTO 9999
C                  ENDIF
C                  DO i=1,NT_CHANNELS
                DO i=1,NT_CHANNELS
                  READ(IFILE,*) nd,CHANNELS(i)
                ENDDO
                CALL CLOSEF(IFILE,ERROR,*9999)
           


C Removed from
C#### Subroutine: DEFILE
C###  Description:
C###    DEFILE defines i/p file for another program.

C 3 Oct 2000
      
C#### Command: FEM define file;l/p/r/w vsaero
C###  Description:
C###    This command writes the appropriate type of file for VSAero.
C###  Parameter:      <basic_input/patch_geometry/wake_input/surface/
C###  Parameter:       streamlines/boundary_layer_input/
C###  Parameter:       off_body_velocity_scan/off_body_streamlines>

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w vsaero'
        OP_STRING(2)=BLANK(1:15)
     '    //'<basic_input/patch_geometry/wake_input'
        OP_STRING(3)=BLANK(1:15)//'/surface_streamlines'
        OP_STRING(4)=BLANK(1:15)//'/boundary_layer_input'
        OP_STRING(5)=BLANK(1:15)//'/off_body_velocity_scan'
        OP_STRING(6)=BLANK(1:15)//'/off_body_streamlines>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------
        ELSE IF(CBBREV(CO,'VSAERO',1,noco+1,NTCO,N3CO)) THEN
          PACKAGE ='VSAERO'
          IF(ABBREV(CO(noco+2),'BASIC_INPUT',1)) THEN
            FILE_EXT='IPVSA'
          ELSE IF(ABBREV(CO(noco+1),'PATCH_GEOMETRY',1)) THEN
            FILE_EXT='IPVSB'
          ELSE IF(ABBREV(CO(noco+1),'WAKE_INPUT',1)) THEN
            FILE_EXT='IPVSC'
          ELSE IF(ABBREV(CO(noco+1),'SURFACE_STREAMLINES',1)) THEN
            FILE_EXT='IPVSD'
          ELSE IF(ABBREV(CO(noco+1),'BOUNDARY_LAYER_INPUT',2)) THEN
            FILE_EXT='IPVSE'
          ELSE IF(ABBREV(CO(noco+1),'OFF_BODY_VELOCITY_SCAN',10)) THEN
            FILE_EXT='IPVSF'
          ELSE IF(ABBREV(CO(noco+1),'OFF_BODY_STREAMLINES',10)) THEN
            FILE_EXT='IPVSG'
          ENDIF
!--------------------------------- VSAERO ------------------------------------
          ELSE IF(PACKAGE(1:6).EQ.'VSAERO') THEN
            IF(FILE_EXT(1:5).EQ.'IPVSA') THEN
              CALL IPVSA(ERROR,*9999)
            ELSE IF(FILE_EXT(1:5).EQ.'IPVSB') THEN
              CALL IPVSB(ERROR,*9999)
            ELSE IF(FILE_EXT(1:5).EQ.'IPVSC') THEN
              CALL IPVSC(XP,ERROR,*9999)
            ELSE IF(FILE_EXT(1:5).EQ.'IPVSD') THEN
              CALL IPVSD(ERROR,*9999)
            ELSE IF(FILE_EXT(1:5).EQ.'IPVSE') THEN
              CALL IPVSE(ERROR,*9999)
            ELSE IF(FILE_EXT(1:5).EQ.'IPVSF') THEN
              CALL IPVSF(ERROR,*9999)
            ELSE IF(FILE_EXT(1:5).EQ.'IPVSG') THEN
              CALL IPVSG(ERROR,*9999)
            ENDIF

C KAT 2001-12-13 removed from DEGAUS

C#### Command: FEM define gauss;c/m image
C###  Description:
C###
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <element (#s/all)[all]>
C###    Specify the element numbers to used. The "all" keyword will
C###    use all currently defined elements in the given regions.
C###  Parameter:      <on WS#[1]>
C###    Specify the workstation (GX window) to draw the 
C###    (object) on. 
C###  Parameter:      <average #PIXELS[1]>

        OP_STRING(1)=STRING(1:IEND)//';c/m image'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<on WS#[1]>'
        OP_STRING(5)=BLANK(1:15)//'<average #PIXELS[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

        ELSE IF(CALCU.OR.MOUSE) THEN
          IF(CBBREV(CO,'IMAGE',1,noco+1,NTCO,N3CO)) THEN
            TYPE='IMAGE'
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            IF(CBBREV(CO,'AVERAGE',1,noco+1,NTCO,N3CO)) THEN
              NOPIXELS=IFROMC(CO(N3CO+1))
            ELSE
              NOPIXELS=1
            ENDIF
          ENDIF

        ELSE IF(CALCU) THEN
          IF(TYPE(1:5).EQ.'IMAGE') THEN
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              nr=NRE(ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Element '',I3)') ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA,XE,XP,ERROR,*9999)
              nb=NBJ(1,ne)
              DO ng=1,NGT(nb)
                CALL XEXG(NBJ(1,ne),ng,NRE(ne),PG,
     '            XE,XG,ERROR,*9999)
                IMAGE_X=1+NINT((XG(1,1)-DBLE(XMIN))*511.0D0/
     '            DBLE(XMAX-XMIN))
                IMAGE_Y=1+NINT((XG(2,1)-DBLE(YMIN))*511.0D0/
     '            DBLE(YMAX-YMIN))
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Image x,y is'',2I4)') IMAGE_X,
     '              IMAGE_Y
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                TOT=0.0D0
                DO i=MAX(IMAGE_X-(NOPIXELS-1),0),
     '            MIN(IMAGE_X+(NOPIXELS-1),IMGX)
                  DO j=MAX(IMAGE_Y-(NOPIXELS-1),0),
     '                 MIN(IMAGE_Y+(NOPIXELS-1),IMGY)
                    iw=1 !added for alpha PJH/5apr93
                    TOT=TOT+DBLE(I2P(i,j,iw))
                  ENDDO
                ENDDO
                YG(1,ng,ne)=TOT*RINTENSITY_MAX/DBLE((2*NOPIXELS-1)*
     '            (2*NOPIXELS-1))
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Density at gauss pt '',I3,'
     '              //''' is '',E12.3)') ng,YG(1,ng,ne)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
          ENDIF


C 18 August 1997
      SUBROUTINE DEGRID(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,
     '  NGAP,NHE,NHP,NHQ,NJE,
     '  NKE,NKH,NPF,NPNE,NPNODE,NQE,NQGE,NQNY,
     '  NRE,NRLIST,NVHE,NVHP,NVJE,
     '  NVQ,NW,NWQ,NXI,NXQ,NYNE,NYNP,NYNQ,NYNR,
     '  DNUDXQ,DXDXIQ,GCHQ,GUQ,PG,SE,XA,XE,XG,XGRC,XIG,
     '  XP,XQ,YP,YQ,ZA,ZE,ZG,ZP,STRING,ERROR,*)

C#### Subroutine: DEGRID
C###  Description:
C###    DEGRID defines finite difference grid.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:defn00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),NJE(NEM),
     '  NKE(NKM,NNM,NBFM,NEFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPF(15,NFM),NPNE(NNM,NBFM,NEFM),NPNODE(0:NP_R_M,0:NRM),
     '  NQE(NSM,NBFM,NEFM),NQGE(NGM,NEFM,NBM),NQNY(2,NYM,0:NRCM,NXM),
     '  NRE(NEM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEFM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEFM),
     '  NVQ(NQM,NAM),NW(NEM,2),NWQ(6,0:NQM,NAM),NXI(-NIM:NIM,0:4,0:NEM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNQ(NHM,NQM,0:NRCM,NXM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEFM),
     '  XA(NAM,NJM,NQM),XE(NSM,NJM),XG(NJM,NUM),XGRC(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  YP(NYM,NIYM,NXM),YQ(NQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,N3CO,no_nrlist,
     '  nr,nx,nxc,nx_d
      CHARACTER FILE*100,STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,DEFORMED,FILIO,GENER,MOUSE,
     '  NOINITIAL
          
      CALL ENTERS('DEGRID',*9999)
C      nc=1 !Temporary GBS 8-FEB-96

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)
        CALL TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define grid;c
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###    Define grid for equation "class"

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define grid;c deformed
C###  Parameter:      <from_class #[1]>
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###    Define grid using deformed coordinates specified by solutions
C###    from class "from_class".

        OP_STRING(1)=STRING(1:IEND)//';c deformed'
        OP_STRING(2)=BLANK(1:15)//'<from_class #[2]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM define grid;r/w;<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <noinitial>
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###    Read/write grid state information at current time. noinitial
C###    stops reading in of old initial conditions from file.

        OP_STRING(1)=STRING(1:IEND)//';r/w;'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<noinitial>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEGRID',ERROR,*9999)
      ELSE
c        CALL ASSERT(CALL_FIBR,'>>Fibres not defined',ERROR,*9999)
c        CALL ASSERT(CALL_ELFB,'>>Fibres elements not defined',
c     '    ERROR,*9999)
C NPS 15/6/97 commented ASSERTS which are not relivant for all 
C problems e.g 1d
!PJH 8Dec95 CALL ASSERT(CALL_EQUA,'>>Define equation first',ERROR,*9999)
        CALL ASSERT(NIQM.GE.6,'>>NIQM must be >= 6',ERROR,*9999)

        CALL PARSE_QUALIFIERS('CRW',noco,1,CO,COQU,CALCU,FILIO,
     '    GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,nxc,CO,ERROR,*9999)
      	CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'DEFORMED',2,noco+1,NTCO,N3CO)) THEN
          DEFORMED=.TRUE.
          IF(CBBREV(CO,'FROM_CLASS',2,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_d,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_d.NE.0,'Invalid class#',ERROR,*9999)
        ELSE
          DEFORMED=.FALSE.
          nx_d=0 !not used
        ENDIF

        IF(CBBREV(CO,'NOINITIAL',2,noco+1,NTCO,N3CO)) THEN
          NOINITIAL=.TRUE.
        ELSE
          NOINITIAL=.FALSE.
        ENDIF

        IF(CALCU) THEN
          ADD=.FALSE.
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF((ITYP4(nr,nx).EQ.4).OR.(ITYP4(nr,nx).EQ.3)) THEN 
C Collocation and Finite Difference regions only
              IF(DEFORMED) THEN
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx_d),NHP(1,nr,nx_d),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx_d,NYNE,NYNP,
     '            YP(1,1,nx_d),ZA,ZP,ERROR,*9999)
              ENDIF
              CALL IPGRID(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NGAP,NHE,
     '          NJE,NKE,NPF,NPNE,NQE,NQGE,nr,NRE,NVHE,NVJE,NVQ,NW,NWQ,
     '          nx_d,NXI,NXQ,DNUDXQ,DXDXIQ,GCHQ,GUQ,PG,SE,XA,XE,XG,
     '          XGRC,XIG,XP,XQ,YQ(1,1,1,nx),ZA,ZE,ZG,ZP,DEFORMED,
     '          NOINITIAL,'CALC',ERROR,*9999)
              CALL_GRID=.TRUE.
              ADD=.TRUE.
              IF (ITYP4(nr,nx).EQ.3) THEN
                CALL CALC_NY_GRID_DEP(NHQ,NQNY,nr,nx,NYNQ,
     '            NYNR(0,0,1,0,nx),ERROR,*9999)
              ENDIF
            ENDIF !Collocation region
          ENDDO !nr
        ELSE !read/write grid info at current time
          CALL ASSERT(CALL_GRID,'Grid points must be defined',ERROR,
     '      *9999)
          CALL TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.iogrid',STATUS,
     '      'SEQUEN','FORMATTED',160,ERROR,*9999)
          IF(IOTYPE.EQ.2) THEN ! read data file
            WRITE(OP_STRING,'('' Reading grid point information'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)            
          ELSE IF(IOTYPE.EQ.3) THEN ! write data file
            WRITE(OP_STRING,'('' Writing grid point information'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL IPGRID(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NGAP,NHE,
     '        NJE,NKE,NPF,NPNE,NQE,NQGE,nr,NRE,NVHE,NVJE,NVQ,NW,NWQ,
     '        nx_d,NXI,NXQ,DNUDXQ,DXDXIQ,GCHQ,GUQ,PG,SE,XA,XE,XG,
     '        XGRC,XIG,XP,XQ,YQ(1,1,1,nx),ZA,ZE,ZG,ZP,DEFORMED,
     '        NOINITIAL,'FILE',ERROR,*9999)
              IF (ITYP4(nr,nx).EQ.3) THEN
                CALL CALC_NY_GRID_DEP(NHQ,NQNY,nr,nx,NYNQ,
     '             NYNR(0,0,1,0,nx),ERROR,*9999)
              ENDIF
          ENDDO !no_nrlist
          CALL CLOSEF(IFILE,ERROR,*9999)          
        ENDIF !calc

      ENDIF

      CALL EXITS('DEGRID')
      RETURN
 9999 CALL ERRORS('DEGRID',ERROR)
      CALL EXITS('DEGRID')
      RETURN 1
      END


C from debase
        ELSE IF(MOUSE) THEN
          OPTION( 1)='1D_Linear'
          OPTION( 2)='1D_Quadratic'
          OPTION( 3)='1D_Cubic_Hermite'
          OPTION( 4)='2D_Bilinear'
          OPTION( 5)='2D_Biquadratic'
          OPTION( 6)='2D_Quadratic_Linear'
          OPTION( 7)='2D_Cubic_Linear'
          OPTION( 8)='2D_Bicubic_Hermite'
          OPTION( 9)='3D_Trilinear'
          OPTION(10)='3D_Cubic_Bilinear'
          OPTION(11)='3D_Bicubic_Linear'
          OPTION(12)='Exit'
c         CALL CHOICE('DEBASE',1,1,Input_Status,72,'EVENT',12,12,noch,
c    '      noco,9,CO,OPTION,STRING,0.01*XDISP,YDISP-0.45*XDISP,ERROR,
c    '      *9999)
          WRITE(OP_STRING,'(X,A)') '>>Choose basis type(s) then Exit'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CONTINUE=.TRUE.
          DO WHILE (CONTINUE)
c           CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,Input_Choice,
c    '        R4DATA,SDATA,ERROR,*9999)
            IF(CLASS(1:6).eq.'CHOICE') THEN
              IF(ID_WS.EQ.72) THEN
                CHOOSE=OPTION(Input_Choice)
                IF(.NOT.ABBREV(CHOOSE,'EXIT',2)) THEN
                  CALL ASSERT(NBT+1.LE.NBFM,'>>NBFMX too small',
     '              ERROR,*9999)
                  NBT=NBT+1
                  NBFT=NBFT+1
                  WRITE(OP_STRING,'('' Basis number '',I1,'
     '              //''' is '',A)')
     '              NBT,CHOOSE
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                IF(ABBREV(CHOOSE,'1D_LINEAR',9)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=1
                  NGT(NBT)=2
                  NKT(0,NBT)=1
                  NNT(NBT)=2
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=2
                  IBT(1,1,NBT)=1
                  IBT(2,1,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(1,nn,1,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  NBI(NBT)=1
                ELSE IF(ABBREV(CHOOSE,'1D_QUADRATIC',12)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=1
                  NGT(NBT)=3
                  NKT(0,NBT)=1
                  NNT(NBT)=3
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  IBT(1,1,NBT)=1
                  IBT(2,1,NBT)=2
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(1,nn,1,NBT)=1
                  ENDDO
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=3
                  NBI(NBT)=1
                ELSE IF(ABBREV(CHOOSE,'1D_CUBIC_HERMITE',16)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=1
                  NGT(NBT)=3
                  NKT(0,NBT)=2
                  NNT(NBT)=2
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  IBT(1,1,NBT)=2
                  IBT(2,1,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(2,nn,0,NBT)=2
                    IDO(1,nn,1,NBT)=1
                    IDO(2,nn,1,NBT)=2
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
C CPB 10/8/94 SWAPPING OVER NBI=4/5
C                  NBI(NBT)=4
                  NBI(NBT)=5
                ELSE IF(ABBREV(CHOOSE,'2D_BILINEAR',11)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=2
                  NGT(NBT)=4
                  NKT(0,NBT)=1
                  NNT(NBT)=4
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=2
                  NGAP(2,NBT)=2
                  IBT(1,1,NBT)=1
                  IBT(2,1,NBT)=1
                  IBT(1,2,NBT)=1
                  IBT(2,2,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=1
                  INP(4,1,NBT)=2
                  INP(1,2,NBT)=1
                  INP(2,2,NBT)=1
                  INP(3,2,NBT)=2
                  INP(4,2,NBT)=2
                  NBI(NBT)=1
                ELSE IF(ABBREV(CHOOSE,'2D_BIQUADRATIC',14)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=2
                  NGT(NBT)=9
                  NKT(0,NBT)=1
                  NNT(NBT)=9
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  NGAP(2,NBT)=3
                  IBT(1,1,NBT)=1 !for Lagrange
                  IBT(2,1,NBT)=2 !for quadratic
                  IBT(1,2,NBT)=1 !for Lagrange
                  IBT(2,2,NBT)=2 !for quadratic
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1 !for vertex 1
                  INP(2,1,NBT)=2 !for vertex 2
                  INP(3,1,NBT)=3 !for vertex 3
                  INP(4,1,NBT)=1 !for vertex 4
                  INP(5,1,NBT)=2 !for vertex 5
                  INP(6,1,NBT)=3 !for vertex 6
                  INP(7,1,NBT)=1 !for vertex 7
                  INP(8,1,NBT)=2 !for vertex 8
                  INP(9,1,NBT)=3 !for vertex 9
                  INP(1,2,NBT)=1 !for vertex 1
                  INP(2,2,NBT)=1 !for vertex 2
                  INP(3,2,NBT)=1 !for vertex 3
                  INP(4,2,NBT)=2 !for vertex 4
                  INP(5,2,NBT)=2 !for vertex 5
                  INP(6,2,NBT)=2 !for vertex 6
                  INP(7,2,NBT)=3 !for vertex 7
                  INP(8,2,NBT)=3 !for vertex 8
                  INP(9,2,NBT)=3 !for vertex 9
                  NBI(NBT)=1
                ELSE IF(ABBREV(CHOOSE,'2D_QUADRATIC_LINEAR',18)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=2
                  NGT(NBT)=6
                  NKT(0,NBT)=1
                  NNT(NBT)=6
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  NGAP(2,NBT)=2
                  IBT(1,1,NBT)=1 !for Lagrange
                  IBT(2,1,NBT)=2 !for quadratic
                  IBT(1,2,NBT)=1 !for Lagrange
                  IBT(2,2,NBT)=1 !for linear
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1 !for vertex 1
                  INP(2,1,NBT)=2 !for vertex 2
                  INP(3,1,NBT)=3 !for vertex 3
                  INP(4,1,NBT)=1 !for vertex 4
                  INP(5,1,NBT)=2 !for vertex 5
                  INP(6,1,NBT)=3 !for vertex 6
                  INP(1,2,NBT)=1 !for vertex 1
                  INP(2,2,NBT)=1 !for vertex 2
                  INP(3,2,NBT)=1 !for vertex 3
                  INP(4,2,NBT)=2 !for vertex 4
                  INP(5,2,NBT)=2 !for vertex 5
                  INP(6,2,NBT)=2 !for vertex 6
                  NBI(NBT)=1
                ELSE IF(ABBREV(CHOOSE,'2D_CUBIC_LINEAR',15)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=2
                  NGT(NBT)=6
                  NKT(0,NBT)=2
                  NNT(NBT)=4
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  NGAP(2,NBT)=2
                  IBT(1,1,NBT)=2
                  IBT(2,1,NBT)=1
                  IBT(1,2,NBT)=1
                  IBT(2,2,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(2,nn,0,NBT)=2
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                    IDO(2,nn,1,NBT)=2
                    IDO(2,nn,2,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=1
                  INP(4,1,NBT)=2
                  INP(1,2,NBT)=1
                  INP(2,2,NBT)=1
                  INP(3,2,NBT)=2
                  INP(4,2,NBT)=2
C CPB SWAPPING OVER NBI=4/5
C                  NBI(NBT)=4
                  NBI(NBT)=5
                ELSE IF(ABBREV(CHOOSE,'2D_BICUBIC_HERMITE',18)) THEN
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=2
                  NGT(NBT)=9
                  NKT(0,NBT)=4
                  NNT(NBT)=4
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  NGAP(2,NBT)=3
                  IBT(1,1,NBT)=2
                  IBT(2,1,NBT)=1
                  IBT(1,2,NBT)=2
                  IBT(2,2,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(2,nn,0,NBT)=2
                    IDO(3,nn,0,NBT)=4
                    IDO(4,nn,0,NBT)=6
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                    IDO(2,nn,1,NBT)=2
                    IDO(2,nn,2,NBT)=1
                    IDO(3,nn,1,NBT)=1
                    IDO(3,nn,2,NBT)=2
                    IDO(4,nn,1,NBT)=2
                    IDO(4,nn,2,NBT)=2
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=1
                  INP(4,1,NBT)=2
                  INP(1,2,NBT)=1
                  INP(2,2,NBT)=1
                  INP(3,2,NBT)=2
                  INP(4,2,NBT)=2
C CPB SWAPPING OVER NBI=4/5
C                  NBI(NBT)=4
                  NBI(NBT)=5
                ELSE IF(ABBREV(CHOOSE,'3D_TRILINEAR',12)) THEN
                  CALL ASSERT(NJT.EQ.3,
     '              '>>Not defined: you must have 3 global coordinates',
     '              ERROR,*9999)
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=3
                  NGT(NBT)=8
                  NKT(0,NBT)=1
                  NNT(NBT)=8
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=2
                  NGAP(2,NBT)=2
                  NGAP(3,NBT)=2
                  IBT(1,1,NBT)=1
                  IBT(2,1,NBT)=1
                  IBT(1,2,NBT)=1
                  IBT(2,2,NBT)=1
                  IBT(1,3,NBT)=1
                  IBT(2,3,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                    IDO(1,nn,3,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=1
                  INP(4,1,NBT)=2
                  INP(5,1,NBT)=1
                  INP(6,1,NBT)=2
                  INP(7,1,NBT)=1
                  INP(8,1,NBT)=2
                  INP(1,2,NBT)=1
                  INP(2,2,NBT)=1
                  INP(3,2,NBT)=2
                  INP(4,2,NBT)=2
                  INP(5,2,NBT)=1
                  INP(6,2,NBT)=1
                  INP(7,2,NBT)=2
                  INP(8,2,NBT)=2
                  INP(1,3,NBT)=1
                  INP(2,3,NBT)=1
                  INP(3,3,NBT)=1
                  INP(4,3,NBT)=1
                  INP(5,3,NBT)=2
                  INP(6,3,NBT)=2
                  INP(7,3,NBT)=2
                  INP(8,3,NBT)=2
                  NBI(NBT)=1
                ELSE IF(ABBREV(CHOOSE,'3D_CUBIC_BILINEAR',17)) THEN
                  CALL ASSERT(NJT.EQ.3,
     '              '>>Not defined: you must have 3 global coordinates',
     '              ERROR,*9999)
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=3
                  NGT(NBT)=12
                  NKT(0,NBT)=2
                  NNT(NBT)=8
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  NGAP(2,NBT)=2
                  NGAP(3,NBT)=2
                  IBT(1,1,NBT)=2
                  IBT(2,1,NBT)=1
                  IBT(1,2,NBT)=1
                  IBT(2,2,NBT)=1
                  IBT(1,3,NBT)=1
                  IBT(2,3,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(2,nn,0,NBT)=2
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                    IDO(1,nn,3,NBT)=1
                    IDO(2,nn,1,NBT)=2
                    IDO(2,nn,2,NBT)=1
                    IDO(2,nn,3,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=1
                  INP(4,1,NBT)=2
                  INP(5,1,NBT)=1
                  INP(6,1,NBT)=2
                  INP(7,1,NBT)=1
                  INP(8,1,NBT)=2
                  INP(1,2,NBT)=1
                  INP(2,2,NBT)=1
                  INP(3,2,NBT)=2
                  INP(4,2,NBT)=2
                  INP(5,2,NBT)=1
                  INP(6,2,NBT)=1
                  INP(7,2,NBT)=2
                  INP(8,2,NBT)=2
                  INP(1,3,NBT)=1
                  INP(2,3,NBT)=1
                  INP(3,3,NBT)=1
                  INP(4,3,NBT)=1
                  INP(5,3,NBT)=2
                  INP(6,3,NBT)=2
                  INP(7,3,NBT)=2
                  INP(8,3,NBT)=2
C CPB SWAPPING OVER NBI=4/5
C                  NBI(NBT)=4
                  NBI(NBT)=5
                ELSE IF(ABBREV(CHOOSE,'3D_BICUBIC_LINEAR',17)) THEN
                  CALL ASSERT(NJT.EQ.3,
     '              '>>Not defined: you must have 3 global coordinates',
     '              ERROR,*9999)
                  NBC(NBT)=1
                  NFBASE(1,NBT)=NBT
                  NIT(NBT)=3
                  NGT(NBT)=18
                  NKT(0,NBT)=4
                  NNT(NBT)=8
                  DO nn=1,NNT(NBT)
                    NKT(nn,NBT)=NKT(0,NBT)
                  ENDDO
                  NGAP(1,NBT)=3
                  NGAP(2,NBT)=3
                  NGAP(3,NBT)=2
                  IBT(1,1,NBT)=2
                  IBT(2,1,NBT)=1
                  IBT(1,2,NBT)=2
                  IBT(2,2,NBT)=1
                  IBT(1,3,NBT)=1
                  IBT(2,3,NBT)=1
                  DO nn=1,NNT(nb)
                    IDO(1,nn,0,NBT)=1
                    IDO(2,nn,0,NBT)=2
                    IDO(3,nn,0,NBT)=4
                    IDO(4,nn,0,NBT)=6
                    IDO(1,nn,1,NBT)=1
                    IDO(1,nn,2,NBT)=1
                    IDO(1,nn,3,NBT)=1
                    IDO(2,nn,1,NBT)=2
                    IDO(2,nn,2,NBT)=1
                    IDO(2,nn,3,NBT)=1
                    IDO(3,nn,1,NBT)=1
                    IDO(3,nn,2,NBT)=2
                    IDO(3,nn,3,NBT)=1
                    IDO(4,nn,1,NBT)=2
                    IDO(4,nn,2,NBT)=2
                    IDO(4,nn,3,NBT)=1
                  ENDDO !nn
                  INP(1,1,NBT)=1
                  INP(2,1,NBT)=2
                  INP(3,1,NBT)=1
                  INP(4,1,NBT)=2
                  INP(5,1,NBT)=1
                  INP(6,1,NBT)=2
                  INP(7,1,NBT)=1
                  INP(8,1,NBT)=2
                  INP(1,2,NBT)=1
                  INP(2,2,NBT)=1
                  INP(3,2,NBT)=2
                  INP(4,2,NBT)=2
                  INP(5,2,NBT)=1
                  INP(6,2,NBT)=1
                  INP(7,2,NBT)=2
                  INP(8,2,NBT)=2
                  INP(1,3,NBT)=1
                  INP(2,3,NBT)=1
                  INP(3,3,NBT)=1
                  INP(4,3,NBT)=1
                  INP(5,3,NBT)=2
                  INP(6,3,NBT)=2
                  INP(7,3,NBT)=2
                  INP(8,3,NBT)=2
C CPB SWAPPING OVER NBI=4/5
C                  NBI(NBT)=4
                  NBI(NBT)=5
                ENDIF
                IF(ABBREV(CHOOSE,'EXIT',2)) THEN
                  CONTINUE=.FALSE.
C                 set to request to eliminate choice menu prompt
c                 CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,*9999)
                ELSE
                  NAT(NBT)=0
                  NST(NBT)=NNT(NBT)*NKT(0,NBT)
                  NUT(NBT)=NIT(NBT)*NIT(NBT)+2
                  CALL GAUSS1(IBT(1,1,NBT),IDO(1,1,0,NBT),INP(1,1,NBT),
     '              NBT,NGAP(1,NBT),PG(1,1,1,NBT),WG(1,NBT),
     '              XIG(1,1,NBT),ERROR,*9999)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          CALL TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipbase',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          CALL IPBASE(IBT,IDO,INP,NAN,0,NDET,NGAP,NVJE,
     '      DET,PG,WG,XIG,ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)

C from decorn
        ELSE IF(MOUSE) THEN
          iw=IWK(1)
          CALL ASSERT(ISNONO(iw,NPNODE(1,1)).GT.0
     '            .OR.ISNONO(iw,NPNODE(1,2)).GT.0,
     '      '>>Node segments not defined',ERROR,*9999)
          CALL ACWK(iw,1,ERROR,*9999)
          DO nrr=1,NRT
            DO nonode=1,NPNODE(0,nrr)
              np=NPNODE(nonode,nrr)
              IF(ISEG(ISNONO(iw,np)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISNONO(iw,np),'VISIBLE',ERROR,*9999)
              ENDIF
              CALL DETECT(iw,ISEG,ISNONO(iw,np),'DETECTABLE',ERROR,
     '          *9999)
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
          CALL ACWK(iw,0,ERROR,*9999)
          WRITE(OP_STRING,'('' >>Pick node on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL PICK(iw,'EVENT',INSTAT,IPICKRET,IPICKID,ERROR,*9999)
          CONTINUE=.TRUE.
          DO WHILE (CONTINUE)
c           CALL EVENT(ID_WS,ID_Device,INSTAT,CLASS,IDATA,R4DATA,SDATA,
c    '        ERROR,*9999)
            IF(CLASS(1:4).EQ.'PICK') THEN
              IF(DOP) THEN
                WRITE(OP_STRING,'(X,A)') 'Input class is pick'
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(INSTAT.EQ.1) THEN
                np=IFROMC(CSEG(IDATA(1))(53:57))
                WRITE(OP_STRING,'('' Node '',I3)') np
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C **            Find the #elements that share the corner node np
                DO ne=1,NET(nr)
                  DO nn=1,NNT(NBJ(1,ne))
                    IF(np.EQ.NPNE(nn,NBJ(1,ne),ne)) THEN
C!!! needs rewriting
c                      NVHP(NVHP(nh,np,nc,nr),nh,np,nc,nr)=ne
c                      NVHP(nh,np,nc,nr)=NVHP(nh,np,nc,nr)+1
                    ENDIF
                  ENDDO
                ENDDO
c                NVHP(nh,np,nc,nr)=NVHP(nh,np,nc,nr)-1
c                CALL ASSERT(NVHP(nh,np,nc,nr).GT.2,
c     '            '>>This is not a corner node',ERROR,*9999)
              ELSE
C!!!            get rid of pick echos
c               CALL INPUT_MODE(iw,LD1,'PICK','REQUEST',ERROR,*9999)
                CONTINUE=.FALSE.
              ENDIF
            ENDIF
          ENDDO
          CALL DAWK(iw,0,ERROR,*9999)

          CALL ACWK(iw,1,ERROR,*9999)
          DO nrr=1,NRT
            DO nonode=1,NPNODE(0,nrr)
              np=NPNODE(nonode,nrr)
              CALL DETECT(iw,ISEG,ISNONO(iw,np),'UNDETECTABLE',ERROR,
     '          *9999)
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)

          CALL TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipcorn',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          CALL IPCORN(INP,NBJ,NJE,NKE,NPF,NPNE,NQE,
     '      nr,NVHP(1,1,1,nr),NVJE,PG,SE,XA,XE,XP,ERROR,*9999)
          CALL CLOSEF(IFILE,ERROR,*9999)        

C from deequa
        ELSE IF(MOUSE) THEN
          nr=NRLIST(1)
          OPTION( 1)='Heat_flow (static)'
          OPTION( 2)='Plane_stress'
          OPTION( 3)='Membrane_stress'
          OPTION( 4)='Heat_flow (dynamic)'
          OPTION( 5)='Advection_diffusion'
          OPTION( 6)='Return'
c         CALL CHOICE('FEM',1,1,INS,72,'EVENT',6,6,noch,noco,2,
c    '      CO,OPTION,STRING,XDISP,0.,ERROR,*9999)
          WRITE(OP_STRING,*) ' >>Choose equation type'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CONTINUE=.TRUE.
          DO WHILE (CONTINUE)
c           CALL EVENT(ID_WS,ID_Device,Input_Status,CLASS,Input_Choice,
c    '        R4DATA,SDATA,ERROR,*9999)
            IF(CLASS(1:6).EQ.'CHOICE') THEN
              IF(ID_WS.EQ.72) THEN
                CHOOSE=OPTION(Input_Choice)
                IF(.NOT.ABBREV(CHOOSE,'RETURN',2)) THEN
                  WRITE(OP_STRING,'('' Equation type is '',A)')CHOOSE
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
                DO il=1,NMM
                  ILP(il,1,nr,nx)=2  !material params defined by elements
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(il,ne,nx)=0.0d0
                  ENDDO
                ENDDO
                IF(ABBREV(CHOOSE,'HEAT_FLOW (STATIC)',15)) THEN
                  ITYP1(nr,nx)=3   !problem type is in FE30
                  ITYP2(nr,nx)=3   !Laplace equation
                  ITYP4(nr,nx)=1   !Galerkin technique
                  ITYP5(nr,nx)=1   !static problem
                  ITYP6(nr,nx)=1   !linear analysis
                  KTYP5=1          !initial solution is zero
                  KTYP6=0          !no boundary integral domains
                  KTYP7=1          !parameters constant wrt time
                  KTYP14=0         !no material params incremented
                  KTYP71=0         !no pressure loads
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    NHE(ne,nx)=1
                    NBH(NH_LOC(1,nx),1,ne)=1
                  ENDDO
                  DO nonode=1,NPNODE(0,nr)
                    np=NPNODE(nonode,nr)
                    NHP(np,nr,nx)=1
                    NKH(NH_LOC(1,nx),np,1,nr)=NKT(0,1)
                  ENDDO
                ELSE IF(ABBREV(CHOOSE,'PLANE_STRESS',10)) THEN
                  ITYP1(nr,nx)=4      ! problem type is in FE40
                  ITYP2(nr,nx)=1      ! linear elasticity problem
                  ITYP4(nr,nx)=1          ! Galerkin technique
                  ITYP5(nr,nx)=1          ! static problem
                  ITYP6(nr,nx)=1          ! linear analysis
                  KTYP5=1          ! initial solution is zero
                  KTYP6=0          ! no boundary integral domains
                  KTYP7=1          ! parameters constant wrt time
                  KTYP14=0         ! no material params incremented
                  KTYP71=0         ! no pressure loads
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    NHE(ne,nx)=1
                    NBH(NH_LOC(1,nx),1,ne)=1
                    NW(ne,1)=11
                  ENDDO
                  ETYP(11)=.TRUE.
                  DO nonode=1,NPNODE(0,nr)
                    np=NPNODE(nonode,nr)
                    NHP(np,nr,nx)=1
                    NKH(NH_LOC(1,nx),np,1,nr)=NKT(0,1)
                  ENDDO
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(1,ne,nx)=10.0D0
                    CE(2,ne,nx)=0.3D0
                    CE(3,ne,nx)=0.0D0
                    CE(4,ne,nx)=0.0D0
                  ENDDO
                ELSE IF(ABBREV(CHOOSE,'MEMBRANE_STRESS',10)) THEN
                ELSE IF(ABBREV(CHOOSE,'HEAT_FLOW (DYNAMIC)',15)) THEN
                  ITYP1(nr,nx)=3       ! problem type is in FE30
                  ITYP2(nr,nx)=3       ! advection-diffusion equation
                  ITYP4(nr,nx)=1       ! Galerkin technique
                  ITYP5(nr,nx)=2       ! time integration problem
                  KTYP22=1          ! linear algorithm for integration
                  KTYP23=2          ! automatic time stepping
                  THETA(1)=0.6667D0 ! integration parameter
                  T0=0.0D0          ! initial time
                  T1=1.0D0          ! final time
                  TINCR=0.001D0     ! initial time increment
                  FILE02=FILE00     ! history file is file00.HISTORY
                  LUMP=.FALSE.      ! no mass lumping
                  ITYP6(nr,nx)=1       ! linear analysis
                  KTYP5=1           ! initial solution is zero
                  KTYP6=0           ! no boundary integral domains
                  KTYP7=1           ! parameters constant wrt time
                  KTYP14=0          ! no material params incremented
                  KTYP71=0          ! no pressure loads
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    NHE(ne,nx)=1
                    NBH(NH_LOC(1,nx),1,ne)=1
                  ENDDO
                  DO nonode=1,NPNODE(0,nr)
                    np=NPNODE(nonode,nr)
                    NHP(np,nr,nx)=1
                    NKH(NH_LOC(1,nx),np,1,nr)=NKT(0,1)
                  ENDDO
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(1,ne,nx)=0.0D0
                    CE(2,ne,nx)=0.0D0
                    CE(3,ne,nx)=1.0D0
                    CE(4,ne,nx)=1.0D0
                    CE(5,ne,nx)=0.0D0
                    CE(6,ne,nx)=1.0D0
                    CE(7,ne,nx)=0.0D0
                  ENDDO
                ELSE IF(ABBREV(CHOOSE,'ADVECTION_DIFFUSION',15)) THEN
                  ITYP1(nr,nx)=3       ! problem type is in FE30
                  ITYP2(nr,nx)=3       ! advection-diffusion equation
                  ITYP4(nr,nx)=1       ! Galerkin technique
                  ITYP5(nr,nx)=2       ! time integration problem
                  KTYP22=1          ! linear algorithm for integration
                  KTYP23=1          ! constant time stepping
                  THETA(1)=0.6667D0 ! integration parameter
                  T0=0.0D0          ! initial time
                  T1=1.0D0          ! final time
                  TINCR=0.01D0      ! time increment
                  FILE02=FILE00     ! history file is file00.HISTORY
                  LUMP=.FALSE.      ! no mass lumping
                  ITYP6(nr,nx)=1       ! linear analysis
                  KTYP5=2           ! initial solution is read in
                  KTYP6=0           ! no boundary integral domains
                  KTYP7=1           ! parameters constant wrt time
                  KTYP14=0          ! no material params incremented
                  KTYP71=0          ! no pressure loads
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    NHE(ne,nx)=1
                    NBH(NH_LOC(1,nx),1,ne)=1
                  ENDDO
                  DO nonode=1,NPNODE(0,nr)
                    np=NPNODE(nonode,nr)
                    NHP(np,nr,nx)=1
                    NKH(NH_LOC(1,nx),np,1,nr)=NKT(0,1)
                  ENDDO
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    CE(1,ne,nx)=0.0D0
                    CE(2,ne,nx)=0.0D0
                    CE(3,ne,nx)=1.0D0
                    CE(4,ne,nx)=1.0D0
                    CE(5,ne,nx)=1.0D0
                    CE(6,ne,nx)=1.0D0
                    CE(7,ne,nx)=0.0D0
                  ENDDO
                ENDIF
                IWRIT1(nr,nx)=1
                IWRIT3(nr,nx)=1
                IWRIT4(nr,nx)=0
              ENDIF
              CONTINUE=.FALSE.
c             CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,*9999)
              DO nrc=0,2
                CALL ASSERT(NYNR(0,nrc,1,nr,nx).LE.NYM,
     '            '>>NYM too small',ERROR,*9999)
              ENDDO !nrc
            ENDIF !choice
          ENDDO !while continue

C     from subroutine DEWIND 

          IF(WINDOW_TYPE(1:3).EQ.'VWS') THEN
      	    !Display choice windows
      	    IF(NJT.LE.2) THEN
      	      IF(iw.EQ.1) THEN
      		OPTION( 1)='Archive..'
      		OPTION( 2)='Cancel..'
      		OPTION( 3)='Change..'
      		OPTION( 4)='Check..'
      		OPTION( 5)='Define..'
      		OPTION( 6)='Display..'
      		OPTION( 7)='Draw..'
      		OPTION( 8)='Hide/Show..'
      		OPTION( 9)='Label..'
      		OPTION(10)='List..'
      		OPTION(11)='Pick..'
      		OPTION(12)='Print..'
      		OPTION(13)='Read..'
      		OPTION(14)='Recall..'
      		OPTION(15)='Refine..'
      		OPTION(16)='Transform..'
      		OPTION(17)='Write..'
      		OPTION(18)='...'
      		OPTION(19)='Exit'
      		NTCH1=19
c     		CALL CHOICE('DETRAN',1,1,INSTAT,91,'EVENT',NTCH1,NTCH1,
c    '		  noch,noco,9,
c    '		  CO,OPTION,STRING,0.0,YDISP-0.50*XDISP,ERROR,*9999)
      	      ENDIF
      	    ELSE IF(NJT.EQ.3) THEN
      	      IF(iw.EQ.1) THEN
      		OPTION( 1)='Archive..'
      		OPTION( 2)='Cancel..'
      		OPTION( 3)='Change..'
      		OPTION( 4)='Check..'
      		OPTION( 5)='Define..'
      		OPTION( 6)='Display..'
      		OPTION( 7)='Draw..'
      		OPTION( 8)='Hide/Show..'
      		OPTION( 9)='Label..'
      		OPTION(10)='List..'
      		OPTION(11)='Pick..'
      		OPTION(12)='Print..'
      		OPTION(13)='Read..'
      		OPTION(14)='Recall..'
      		OPTION(15)='Refine..'
      		OPTION(16)='Transform..'
      		OPTION(17)='Write..'
      		OPTION(18)='...'
      		OPTION(19)='Exit'
      		NTCH1=19
c     		CALL CHOICE('DETRAN',1,1,INSTAT,91,'EVENT',NTCH1,NTCH1,
c    '		  noch,noco,9,
c    '		  CO,OPTION,STRING,0.0,0.50*DISP,ERROR,*9999)
      	      ELSE IF(iw.EQ.2) THEN
      		OPTION( 1)='Archive..'
      		OPTION( 2)='Cancel..'
      		OPTION( 3)='Change..'
      		OPTION( 4)='Define..'
      		OPTION( 5)='Draw..'
      		OPTION( 6)='Hide/Show..'
      		OPTION( 7)='Label..'
      		OPTION( 8)='Pick..'
      		OPTION( 9)='Print..'
      		OPTION(10)='Recall..'
      		OPTION(11)='Transform..'
      		OPTION(12)='Exit'
      		NTCH1=12
c     		CALL CHOICE('DETRAN',1,1,INSTAT,92,'EVENT',NTCH1,NTCH1,
c    '		  noch,noco,9,
c    '		  CO,OPTION,STRING,0.99*DISP,0.99*DISP,ERROR,*9999)
      	      ELSE IF(iw.EQ.3) THEN
      		OPTION( 1)='Archive..'
      		OPTION( 2)='Back plane'
      		OPTION( 3)='Cancel..'
      		OPTION( 4)='Define..'
      		OPTION( 5)='Draw..'
      		OPTION( 6)='Front plane'
      		OPTION( 7)='Hide/Show..'
      		OPTION( 8)='Label..'
      		OPTION( 9)='Pan'
      		OPTION(10)='Parallel'
      		OPTION(11)='Perspective'
      		OPTION(12)='Print..'
      		OPTION(13)='Proj ref pt'
      		OPTION(14)='Recall..'
      		OPTION(15)='Rescale'
      		OPTION(16)='Reset'
      		OPTION(17)='Rotate data'
      		OPTION(18)='Rotate view'
      		OPTION(19)='Save view'
      		OPTION(20)='Select view'
      		OPTION(21)='View ref pt'
      		OPTION(22)='View plane'
      		OPTION(23)='View dist'
      		OPTION(24)='View up'
      		OPTION(25)='Zoom'
      		OPTION(26)='Exit'
      		NTCH1=26
c     		CALL CHOICE('DETRAN',1,1,INSTAT,93,'EVENT',NTCH1,NTCH1,
c    '		  noch,noco,9,
c    '		  CO,OPTION,STRING,0.99*DISP,0.48*DISP,ERROR,*9999)
      	      ENDIF
      	    ENDIF

      	    IF(iw.EQ.4) THEN       !Hammer projection window
      	      OPTION( 1)='Archive..'
      	      OPTION( 2)='Cancel..'
      	      OPTION( 3)='Define..'
      	      OPTION( 4)='Draw..'
      	      OPTION( 5)='Hide/Show..'
      	      OPTION( 6)='Label..'
      	      OPTION( 7)='Pick..'
      	      OPTION( 8)='Print..'
      	      OPTION( 9)='Recall..'
      	      OPTION(10)='Transform..'
      	      OPTION(11)='Exit'
      	      NTCH1=11
c     	      CALL CHOICE('DETRAN',1,1,INSTAT,94,'EVENT',NTCH1,NTCH1,
c    '		noch,noco,9,CO,OPTION,STRING,0.20*DISP,0.25*DISP,ERROR,
c    '		*9999)

      	    ELSE IF(iw.EQ.35) THEN !signal trace window
      	      OPTION( 1)='Accept'
      	      OPTION( 2)='Reject'
      	      OPTION( 3)='Pick Signal'
      	      OPTION( 4)='Enter Signal'
      	      OPTION( 5)='Recalculate..'
      	      OPTION( 6)='Exit'
      	      NTCH1=6
c     	      CALL CHOICE('DITRAC',1,1,INSTAT,36,'EVENT',NTCH1,NTCH1,
c    '		noch,noco,9,CO,OPTION,STRING,0.0,YDISP-0.40*XDISP,ERROR,
c    '		*9999)

      	    ELSE IF(iw.EQ.95) THEN !epicardial mapping top level window
      	      OPTION( 1)='Read signals..'
      	      OPTION( 2)='Setup'
      	      OPTION( 3)='Display signals'
      	      OPTION( 4)='Write signals'
      	      OPTION( 5)='Display field'
      	      OPTION( 6)='Pick electrode'
      	      OPTION( 7)='Cancel..'
      	      OPTION( 8)='Change colour..'
      	      OPTION( 9)='Exit'
      	      NTCH1=9
c     	      CALL CHOICE('DITRAC',1,1,INSTAT,95,'EVENT',NTCH1,NTCH1,
c    '		noch,noco,9,CO,OPTION,STRING,0.0,0.30*XDISP,ERROR,*9999)
      	    ENDIF
          ENDIF


C KAT 2001-12-13
      SUBROUTINE DEIMAG(STRING,ERROR,*)

C#### Subroutine: DEIMAG
C###  Description:
C###    DEIMAG defines image.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:pics00.cmn'
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IFROMC,IEND,IEND1,IEND2,
     '  IOS,IUNIT,iw,N3CO
      REAL*8 RFROMC
      CHARACTER CHAR*3,FILE*100,STATUS*3
      LOGICAL CALCU,CBBREV,FILIO,GENER,MOUSE

      CALL ENTERS('DEIMAG',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define image;r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Read or write an image file FILENAME.img in the directory specified by PATH.
C###  Parameter:      <for WS_ARRAY#[1]>
C###    Specify an image number (for multiple images).
C###  Parameter:      <maximum INTENSITY#[1.0]>
C###    Specify a maximum intensity.
C###  Parameter:      <size X#[512] <by Y#[X#]>>
C###    Specify the image size in pixels.

        OP_STRING(1)=STRING(1:IEND)//';r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<for WS_ARRAY#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<maximum INTENSITY#[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<size X#[512] <by Y#[X#]>>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEIMAG',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('RW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN
          IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
            iw=IFROMC(CO(N3CO+1))
          ELSE
            iw=1
          ENDIF
          IF(CBBREV(CO,'MAXIMUM',1,noco+1,NTCO,N3CO)) THEN
            RINTENSITY_MAX=RFROMC(CO(N3CO+1))
          ELSE
            RINTENSITY_MAX=1.0D0
          ENDIF
          IF(CBBREV(CO,'SIZE',1,noco+1,NTCO,N3CO)) THEN
            IMGX=IFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
              IMGY=IFROMC(CO(N3CO+1))
            ELSE
              IMGY=IMGX
            ENDIF
          ELSE
            IMGX=512
            IMGY=512
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Reading image size '',I4,'' by '',I4)') IMGX,IMGY
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
c cpb 14/11/94 Temporarily commenting this out
C        CALL ASSERT(NXM.EQ.512,'>>Require NXM = 512 for image arrays',
C     '    ERROR,*9999)

        IF(FILIO) THEN
          IUNIT=8
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPEN(UNIT=IUNIT,FILE=FILE(IBEG:IEND)//'.IMG',STATUS=STATUS,
     '      FORM='UNFORMATTED',IOSTAT=IOS)
          IF(ios.NE.0) THEN
            WRITE(CHAR,'(I3)') IOS
            ERROR=' File open error: iostat='//CHAR
            GOTO 9999
          ENDIF
          IF(IOTYPE.EQ.2) THEN
c PJH 2/12/91 CALL IMAGE_READ(I2P(1,iw),IUNIT)
            CALL IMAGE_READ(I2P(1,1,iw),IUNIT)
          ELSE IF(IOTYPE.EQ.3) THEN
c PJH 2/12/91 CALL IMAGE_WRITE(I2P(1,iw),IUNIT)
            CALL IMAGE_WRITE(I2P(1,1,iw),IUNIT)
          ENDIF
          CALL CLOSEF(8,ERROR,*9999)
        ENDIF

        CALL_IMAG=.TRUE.
      ENDIF

      CALL EXITS('DEIMAG')
      RETURN
 9999 CALL ERRORS('DEIMAG',ERROR)
      CALL EXITS('DEIMAG')
      RETURN 1
      END


      SUBROUTINE DEMACR(STRING,ERROR,*)

C#### Subroutine: DEMACR
C###  Description:
C###    DEMACR defines macro commands.

C**** NT_MACRO_names(0) is total number of commands defined.
C**** NT_MACRO_names(macro_command) is  number of lines defined for 
C**** command.
C**** MACRO_names(macro_command) is command name for command.
C**** MACRO_COMMAND_buffer(nomacr,macro_command) is command string in 
C**** line nomacr of macro_command.
C**** MACRO_KEY_buffer(nomacr,macro_key),nomacr=1,NT_MACRO(macro_key) 
C**** are command lines associated with macro key number MACRO_KEY.
C**** MACRO_COMMAND_EXECUTE is set .TRUE. to cause parsing of macro.
C**** MACRO_KEY_EXECUTE is set .TRUE. to cause parsing of key macro.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IFROMC,IOSTAT,irec,LINE_number,
     '  macro_command,N3CO
      CHARACTER C1*255,FILE*100,STATUS*3,STRG*132,TYPE*7
      LOGICAL CALCU,CBBREV,CONTINUE,FILIO,GENER,MOUSE

      CALL ENTERS('DEMACR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)
        CALL TRIM(PATH00,IBEG2,IEND2)
C---------------------------------------------------------------------

C#### Command: define macro NAME <line (#/next)[next]> "COMMAND_STRING" 
C###  Description:
C###    

        OP_STRING(1)=' define macro NAME <line (#/next)[next]>'
     '    //' "COMMAND_STRING"'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command:  define macro key # as MACRO_NAME
C###  Description:
C###    
        OP_STRING(1)=' define macro key # as MACRO_NAME'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command:  define macro;l/r/w<;FILENAME[file]><;(PATH/example)[ ]> <(key/macro)[macro] #[1]> 
C###  Description:
C###    

        OP_STRING(1)=' define macro;l/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
     '    //' <(key/macro)[macro] #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEMACR',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS(' LRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO) THEN 
          IF(CBBREV(CO,'KEY',1,noco+1,NTCO,N3CO).AND.NTCO.EQ.4) THEN
            TYPE='KEY'
            macro_key=IFROMC(CO(N3CO+1))
          
          ELSE IF(CBBREV(CO,'MACRO',1,noco+1,NTCO,N3CO).AND.NTCO.EQ.4) 
     '      THEN
            TYPE='COMMAND'
            macro_command=IFROMC(CO(N3CO+1))

          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF

        ELSE IF(.NOT.FILIO) THEN 
          IF(CBBREV(CO,'KEY',1,noco+1,NTCO,N3CO)) THEN
            TYPE='KEY'
            macro_key=IFROMC(CO(N3CO+1))
            IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
C Need separate list of macro names for key #s
C At present the key # = macro command #
            ENDIF

          ELSE IF(NTCO.GT.noco+1) THEN !macro command
            TYPE='COMMAND'
            CALL TRIM(CO(noco+1),IBEG,IEND)
            CALL CUPPER(CO(noco+1),C1)
C            IF(CBBREV(MACRO_names,CUPPER(CO(noco+1)),IEND-IBEG+1,1,
            IF(CBBREV(MACRO_names,C1,IEND-IBEG+1,1,
     '        NT_MACRO_names(0),macro_command)) THEN
              !existing macro has ID macro_command
            ELSE !define new macro
              NT_MACRO_names(0)=NT_MACRO_names(0)+1
              macro_command=NT_MACRO_names(0)
              MACRO_names(macro_command)=CO(noco+1)(IBEG:IEND)
            ENDIF
            IF(DOP) THEN
              WRITE(*,'('' macro command ID = '',I1)') 
     '          macro_command
              WRITE(*,'('' #macro commands  = '',I1)') 
     '          NT_MACRO_names(0)
              WRITE(*,'(''        command   = '',A)') 
     '          MACRO_names(macro_command)
            ENDIF !dop

            IF(CBBREV(CO,'LINE',1,noco+1,NTCO,N3CO)) THEN
              LINE_number=IFROMC(CO(N3CO+1))
            ELSE
              LINE_number=NT_MACRO_names(macro_command)+1
            ENDIF
            IF(LINE_number.GT.NT_MACRO_names(macro_command)) THEN
              NT_MACRO_names(macro_command)=LINE_number
            ENDIF
            MACRO_COMMAND_buffer(LINE_number,macro_command)=CO(NTCO)

          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF !key or ntco>noco+1
        ENDIF

        IF(FILIO) THEN
          CALL TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.macro',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)

          IF(IOTYPE.EQ.1) THEN      !read macro from keyboard
            IF(TYPE(1:3).EQ.'KEY') THEN
              WRITE(OP_STRING,'('' Use DO key to define macro'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(TYPE(1:5).EQ.'MACRO') THEN
              CO(noco+1)='?'
              GO TO 1
            ENDIF

          ELSE IF(IOTYPE.EQ.2) THEN !read macro from file
            IF(TYPE(1:3).EQ.'KEY') THEN
              irec=1
              CONTINUE=.TRUE.
              DO WHILE(CONTINUE)
                READ(IFILE,FMT='(A)',REC=irec,IOSTAT=IOSTAT) STRG
                IF(IOSTAT.EQ.0) THEN
                  MACRO_KEY_buffer(irec,macro_key)=STRG
                  irec=irec+1
                ELSE IF(IOSTAT.EQ.36) THEN
                  CONTINUE=.FALSE.
                  CALL CLOSEF(IFILE,ERROR,*9999)
                  WRITE(OP_STRING,
     '              '('' >>Macro defined as key '',I1)') macro_key
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  NT_MACRO(macro_key)=irec-1  !#lines in macro
                ENDIF
              ENDDO !while continue

            ELSE IF(TYPE(1:7).EQ.'COMMAND') THEN
              irec=1
              CONTINUE=.TRUE.
              DO WHILE(CONTINUE)
                READ(IFILE,FMT='(A)',REC=irec,IOSTAT=IOSTAT) STRG
                IF(IOSTAT.EQ.0) THEN
                  MACRO_COMMAND_buffer(irec,macro_command)=STRG
                  irec=irec+1
                ELSE IF(IOSTAT.EQ.36) THEN
                  CONTINUE=.FALSE.
                  CALL CLOSEF(IFILE,ERROR,*9999)
                  WRITE(OP_STRING,
     '              '('' >>Macro defined as command '',I1)') 
     '              macro_command
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  NT_MACRO_names(macro_command)=irec-1 !#lines in macro
                ENDIF
              ENDDO !while continue
            ENDIF !type

          ELSE IF(IOTYPE.EQ.3) THEN !write macro to file
            IF(TYPE(1:3).EQ.'KEY') THEN
              DO irec=1,NT_MACRO(macro_key)
                STRG=MACRO_KEY_buffer(irec,macro_key)
                WRITE(IFILE,FMT='(A)',REC=irec) STRG
              ENDDO

            ELSE IF(TYPE(1:7).EQ.'COMMAND') THEN
              DO irec=1,NT_MACRO_names(macro_command)
                STRG=MACRO_COMMAND_buffer(irec,macro_command)
                WRITE(IFILE,FMT='(A)',REC=irec) STRG
              ENDDO
            ENDIF
            CALL CLOSEF(IFILE,ERROR,*9999)

          ENDIF !iotype
        ENDIF !filio
      ENDIF

      CALL EXITS('DEMACR')
      RETURN
 9999 CALL ERRORS('DEMACR',ERROR)
      CALL EXITS('DEMACR')
      RETURN 1
      END        


      SUBROUTINE DETEXT(NOCO,NTCO,NTCOQU,CO,COQU,STRING,ERROR,*)

C#### Subroutine: DETEXT
C###  Description: 
C**** Defines text.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:text00.cmn'
!     Parameter List
      INTEGER NOCO,NTCO,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER I,IBEG,IBEG1,IEND,IEND1,IOSTAT,
     '  INSTAT,IW,IWK(6),NOTEXT,NTIW
      LOGICAL CALCU,FILIO,GENER,MOUSE
      CHARACTER FILE*100,STATUS*3

      CALL ENTERS('DETEXT',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//']<;doc>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//' TEXT_STRING'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DETEXT',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS(' DLMPRW',NOCO,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)
        IF(MOUSE) CALL WS_LIST(IWK,0,NTIW,NOCO,NTCO,CO,ERROR,*9999)

        IF(MOUSE) THEN
          CALL ACWK(IW,0,ERROR,*9999)
          WRITE(OP_STRING,'('' >>Locate text on '',I1)') IW
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL LOCATOR(1,IW,INSTAT,'REQUEST',3,0.0D0,XWC_TEXT,0.0D0,
     '      YWC_TEXT,ERROR,*9999)
          CALL DAWK(IW,0,ERROR,*9999)
        ENDIF

        IF(FILIO) THEN
          CALL TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.iptext',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          DO I=1,20
            NOTEXT=NTTEXT+I
            CALL ASSERT(NOTEXT.LE.20,'>>Too many lines of text',
     '        ERROR,*9999)
            READ(IFILE,FMT='(A)',REC=I,IOSTAT=IOSTAT)
     '        TEXT(NOTEXT)(1:80)
            IF(IOSTAT.NE.0) GO TO 200
            IF(DOP) THEN
              WRITE(*,*) TEXT(NOTEXT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
 200      NTTEXT=NOTEXT
          CALL CLOSEF(IFILE,ERROR,*9999)

        ELSE
          NTTEXT=NTTEXT+1
          TEXT(NTTEXT)=CO(NOCO+1)
        ENDIF

      ENDIF

      CALL EXITS('DETEXT')
      RETURN
 9999 CALL ERRORS('DETEXT',ERROR)
      CALL EXITS('DETEXT')
      RETURN 1
      END


      SUBROUTINE DEVECT(noco,NTCOQU,CO,COQU,STRING,ERROR,*)

C#### Subroutine: DEVECT
C###  Description:
C###    Calculates vector stored in AVECT.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:coef00.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER noco,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,ICHAR,IEND,IEND1,INFO,ncvect,NOQUES
      CHARACTER CHAR*3,FILE*100,STATUS*3
      LOGICAL CALCU,FILEIP,FILIO,GENER,MOUSE

      CALL ENTERS('DEVECT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//']<;doc>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEVECT',ERROR,*9999)
      ELSE
        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        FILEIP=.FALSE.
        NOQUES=0

        IF(FILIO) THEN
          CALL TRIM(FILE,IBEG,IEND)
          CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipvect',STATUS,
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          IDEFLT(1)=NT_MATR
          WRITE(CHAR,'(I2)') NT_MATR
          CALL TRIM(CHAR,IBEG,IEND)
          FORMAT='($,'' Enter the vector length ['//CHAR(IBEG:IEND)//
     '      ']: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NT_VECT
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,1,50,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
c         CALL IINOUT(IOTYPE,IVDU,IFILE,FORMAT,1,IDATA,IDEFLT,1,50,
c    '      INFO,ERROR,*9999)
          IF(iotype.NE.3) NT_VECT=IDATA(1)
          FORMAT='($,'' Enter coefficients: '',10E11.4)'
          IF(IOTYPE.EQ.3) THEN
            DO ncvect=1,NT_VECT
              RDATA(ncvect)=AVECT(ncvect)
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '      IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '      INFO,ERROR,*9999)
c         CALL RINOUT(IOTYPE,IVDU,IFILE,FORMAT,NT_VECT,RDATA,RDEFLT,
c    '      -RMAX,RMAX,INFO,ERROR,*9999)
          IF(iotype.NE.3) THEN
            DO ncvect=1,NT_VECT
              AVECT(ncvect)=RDATA(ncvect)
            ENDDO
          ENDIF
          CALL CLOSEF(IFILE,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('DEVECT')
      RETURN
 9999 CALL ERRORS('DEVECT',ERROR)
      CALL EXITS('DEVECT')
      RETURN 1
      END


      SUBROUTINE DEXI(IBT,IDO,INP,NBJ,NEELEM,NELIST,
     '  NJE,NJP,NKE,noco,
     '  NPF,NPLIST,NPNE,NPNODE,NQE,NRLIST,NTCO,NTCOQU,NVJE,
     '  SE,XA,XE,XP,CO,COQU,STRING,ERROR,*)

C#### Subroutine: DEXI
C###  Description:
C###    Defines nodal parameters.
C###    NPNODE(0,nr) is the total number of nodes in region nr.
C###    NPNODE(nonode,nr), nonode=1..NPNODE(0,nr) are the node numbers.
C###    NPT(nr) is the highest node number.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,*),NEELEM(0:NE_R_M,0:*),NELIST(0:*),NJE(*),NJP(*),
     '  NKE(NKM,NNM,NBFM,*),noco,NPF(15,*),NPLIST(0:*),NPNE(NNM,NBFM,*),
     '  NPNODE(0:NP_R_M,0:*),NQE(NSM,NBFM,*),NRLIST(0:*),
     '  NTCO,NTCOQU(*),NVJE(NNM,NBFM,NJM,*)
      REAL*8 SE(NSM,NBFM,*),XA(NAM,NJM,*),XE(NSM,*),XP(NKM,NVM,NJM,*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,N3CO,nb,ne,ni,nj,nolist,
     '  nonode,np,nv
      REAL*8 X(3),XI(3)
      CHARACTER BASIS_TYPE*9,FILE*100,STATUS*3,TYPE*7
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,GENER,MOUSE

      CALL ENTERS('DEXI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME>['//FILE00(IBEG1:IEND1)//']<;doc>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(3)=BLANK(1:15)//'<node #s>[all]'
        OP_STRING(4)=BLANK(1:15)//'<in ELEMENTS>[all]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c image'
        OP_STRING(2)=BLANK(1:15)
     '    //'<linear/quadratic/cubic>[quadratic]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEXI',ERROR,*9999)
      ELSE
        nv=1 ! Temporary MPN 12-Nov-94
        CALL PARSE_QUALIFIERS('CDLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(CALCU) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          IF(CBBREV(CO,'IMAGE',1,noco+1,NTCO,N3CO)) THEN
            TYPE='IMAGE'
          ELSE
            TYPE='DEFAULT'
          ENDIF

          IF(TYPE(1:5).EQ.'IMAGE') THEN
            IF(CBBREV(CO,'LINEAR',1,noco+1,NTCO,N3CO)) THEN
              BASIS_TYPE='LINEAR'
            ELSE IF(CBBREV(CO,'QUADRATIC',1,noco+1,NTCO,N3CO)) THEN
              BASIS_TYPE='QUADRATIC'
            ELSE IF(CBBREV(CO,'CUBIC',1,noco+1,NTCO,N3CO)) THEN
              BASIS_TYPE='CUBIC'
            ENDIF

          ELSE IF(TYPE(1:7).EQ.'DEFAULT') THEN
            IF(CBBREV(CO,'NODES',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),ERROR,
     '          *9999)
            ELSE
              DO nonode=1,NPNODE(0,2)
                NPLIST(nonode)=NPNODE(nonode,2)
              ENDDO
              NPLIST(0)=NPNODE(0,2)
            ENDIF
          ENDIF
        ENDIF

        IF(FILIO) THEN

        ELSE IF(CALCU) THEN !Calculate Xi positions
          IF(TYPE(1:5).EQ.'IMAGE') THEN
c            CALL CALC_IMAGE_XI(BASIS_TYPE,ERROR,*9999)

          ELSE IF(TYPE(1:7).EQ.'DEFAULT') THEN !Calc Xi for nodes
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              DO nj=1,NJP(np)
                X(nj)=XP(1,nv,nj,np)  !are Xj coords of node np
              ENDDO
              CALL XCOORD(IBT,IDO,INP,NBJ,ne,NEELEM,
     '          NJE,NJP,NKE,NPF,NPNE,NPNODE,NQE,NVJE,
     '          SE,XA,XE,XI,XP,X,ERROR,*9999)
              nb=NBJ(1,ne)
              DO ni=1,NIT(nb)
                XP(2,nv,ni,np)=XI(ni) !are Xi coords of node np
              ENDDO
              WRITE(OP_STRING,'('' Node '',I3,'' lies in '
     '          //'element '',I3,'' at Xi coords: '',3E11.3)')
     '          np,ne,(XI(ni),ni=1,NIT(nb))
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

        ENDIF
      ENDIF

      CALL EXITS('DEXI')
      RETURN
 9999 CALL ERRORS('DEXI',ERROR)
      CALL EXITS('DEXI')
      RETURN 1
      END


      SUBROUTINE DEXI(IBT,IDO,INP,LD,LN,NBJ,
     '  NBH,NDDL,NDLT,NEELEM,NELIST,NHE,NJE,NKE,
     '  NLS_NDATA_CONT,NLS_NDATA_IMAG,noco,NPF,NPNE,NQE,NRE,NTCO,
     '  NTCOQU,NVHE,NVJE,NW,NXI,
     '  EDD,NLS_CON_PSI,NLS_CON_XI,
     '  SE,SQ,WD,XA,XE,XID,XP,ZA,ZD,ZDD,ZP,
     '  CO,COQU,STRING,ERROR,*)

C#### Subroutine: DEXI
C###  Description:
C###    DEXI defines data xi points.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:b10.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:bem000.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:data00.cmn'
      INCLUDE 'cmiss$reference:defn00.cmn'
      INCLUDE 'cmiss$reference:dft00.cmn'
      INCLUDE 'cmiss$reference:ditr00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:four00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:hist00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ioda00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:ktyp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
      INCLUDE 'cmiss$reference:map3d.cmn'
      INCLUDE 'cmiss$reference:moti00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:read00.cmn'
      INCLUDE 'cmiss$reference:scal01.cmn'
      INCLUDE 'cmiss$reference:time00.cmn'
      INCLUDE 'cmiss$reference:tree00.cmn'
!     Parameter List
C*** Pass NRLIST
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(*),LN(0:*),
     '  NBJ(NJM,NEM),NBH(NHM,NCM,NEM),NDDL(NEM,*),NDLT(*),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NE_R_M),NHE(NEM,NXM),NJE(NEM),
     '  NKE(NKM,NNM,NBFM,NEFM),
     '  NLS_NDATA_CONT(NDM),NLS_NDATA_IMAG(NDM),noco,
     '  NPF(15,NFM),NPNE(NNM,NBFM,NEFM),
     '  NQE(NSM,NBFM,NEFM),NRE(NEM),NTCO,
     '  NTCOQU(*),NVHE(NNM,NBFM,NHM,NEFM),
     '  NVJE(NNM,NBFM,NJM,NEFM),NW(NEM,2),NXI(-NIM:NIM,0:NEM)
      REAL*8 EDD(*),
     '  NLS_CON_PSI(16,2,1000,3,100),NLS_CON_XI(3,2,1000,100),
     '  SE(NSM,NBFM,NEFM),SQ(*),WD(NJM,*),
     '  XA(NAM,NJM,NQM),XE(NSM,NJM),XID(NIM,*),XP(NKM,NVM,NJM,NPM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,*),ZDD(NDM,*),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER DEFLIST(3),i,IBEG,IBEG1,icol,IEND,IEND1,

!!! ALISTAIR: check this: PJH 2Sept95
     '  ICON_X(200),ISURF_ELEM(1,10),

     '  INDEX_SURF(4),IT,ITMAX,
     '  l,n,nb,n1list,N3CO,nc,NCLUSTER,nd,ND0,ND1,ne,
     '  ni,NITB,nj,nj1,nj2,nk,NKTB,nn,nimag,
     '  noelem,nolist,np,nr,NRLIST(0:9),
     '  ns,ns1,ns2,ns3,ns4,nt,NTILDEF,NTIL,NTIMES,nv,nx
      REAL*8 A1,A2,AB(3),ALFA,B1,B2,BDOTC,BETA,BMAG(3),
     '  C1,C2,CMAG(3),CVEC(3,3),D1,D2,DELTA,DENOM1,DENOM2,DIFF,DIST,
     '  GAMA,PSI1,PXI,SLOPE,SQMAX,SQND,
     '  THETAMAX,THETAMIN,TMIN(2),TMAX(2),TREF,
     '  X0,X1,X2,XD(3),XI(3),XI3OFF,XI_3,Z1(4),Z2(3),Z3(3),Z4(3)
      CHARACTER CIW*1,FILE*100,STATUS*3,TYPE*20
      LOGICAL ALL_REGIONS,CALCU,CBBREV,CENTROID,
     '  CUBIC,DEFORM,EXCLUDE,EXTRAPOLATE,FILIO,FINISHED,
     '  FOUND,GENER,IN_ELEM,INLIST,MOUSE,NEW,ORTHOG,QUADRATIC,
     '  SET_XI_3,TAG2D
      EXTERNAL E04JBQ,MONIT

      CALL ENTERS('DEXI',*9999)

C ??? Why is nx used here ???
      nx=1 ! temporary cpb 22/11/94 
      nv=1 !temporary

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)
        WRITE(CIW,'(I1)') 2*NJT-3+IMAP

        OP_STRING(1)=STRING(1:IEND)//';r/w;FILENAME'
        OP_STRING(2)=BLANK(1:15)//'<total NIT>[NJT-1]'
        OP_STRING(3)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(4)=BLANK(1:15)//'<in ELEMENTS>[all]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C GMH 30/9/95 Allowing the user to force the projection process
C     to start from a specific position.
C     The documentation below is wrong - need to show that they
C     can specify - which data point(s)
C                 - start element number
C                 - start xi position
C                 - whether the search procedure can move to
C                     - different elements
C                     - move around in the same element
C                     - absolute value
        OP_STRING(1)=STRING(1:IEND)//';c specify'
        OP_STRING(2)=BLANK(1:15)//'<from DATA_POINT>[all]'
        OP_STRING(3)=BLANK(1:15)
     '    //'<linear/quadr/cubic/orthogonal>[linear]'
        OP_STRING(4)=BLANK(1:15)//'<extrapolate>'
        OP_STRING(5)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(6)=BLANK(1:15)//'<in ELEMENTS>[all]'
        OP_STRING(7)=BLANK(1:15)
     '    //'<undeform/deform DEFLIST[xyz]>[undeform]'
        OP_STRING(8)=BLANK(1:15)//'<max DIST[10000]>'
        OP_STRING(9)=BLANK(1:15)//'<accept XI3_OFFSET[0]>'
        OP_STRING(10)=BLANK(1:15)//'<xi_3 XI_3>[0.0]'
        OP_STRING(11)=BLANK(1:15)//'<element NUMBER[0]>'
        OP_STRING(12)=BLANK(1:15)//'<xi XIPOS[000]>'
        OP_STRING(13)=BLANK(1:15)//'<global/local/absolute>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'<new/old>[new]'
        OP_STRING(3)=BLANK(1:15)
     '    //'<linear/quadr/cubic/orthogonal>[linear]'
        OP_STRING(4)=BLANK(1:15)//'<extrapolate>'
        OP_STRING(5)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(6)=BLANK(1:15)//'<in ELEMENTS>[all]'
        OP_STRING(7)=BLANK(1:15)
     '    //'<undeform/deform DEFLIST[xyz]>[undeform]'
        OP_STRING(8)=BLANK(1:15)//'<max DIST[10000]>'
        OP_STRING(9)=BLANK(1:15)//'<accept XI3_OFFSET[0]>'
        OP_STRING(10)=BLANK(1:15)//'<xi_3 XI_3>[0.0]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'centroid'
        OP_STRING(3)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(4)=BLANK(1:15)//'<in ELEMENTS>[all]'
        OP_STRING(5)=BLANK(1:15)//'<cluster NCLUSTER>[72]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(3)=BLANK(1:15)//'<sub DXI1,DXI2>[sub 10,10]'
        OP_STRING(4)=BLANK(1:15)//'<in ELEMENTS>[all]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c'
        OP_STRING(2)=BLANK(1:15)//'tag2d/contour'
        OP_STRING(3)=BLANK(1:15)//'<max DIST[10000]>'
        OP_STRING(4)=BLANK(1:15)//'<surf INDEX>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';c beads'
        OP_STRING(2)=BLANK(1:15)//'<region #s/ALL>[1]'
        OP_STRING(3)=BLANK(1:15)//'<in ELEMENTS>[all]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEXI',ERROR,*9999)
      ELSE
        nc=1 ! Temporary MPN 12-Nov-94
        CALL PARSE_QUALIFIERS(' CDLMPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

c       IF(NTCO.GE.noco+1) THEN
          IF(CBBREV(CO,'BEADS',1,noco+1,NTCO,N3CO)) THEN
            TYPE='BEADS'
          ELSE !used to be: if(cbbrev(co,'xi',2,noco+1,ntco,n3co)) then
            TYPE='XIPOS'
            IF(CBBREV(CO,'LINEAR',1,noco+1,NTCO,N3CO)) THEN
              CUBIC=.FALSE.
              QUADRATIC=.FALSE.
              ORTHOG=.FALSE.
            ELSE IF(CBBREV(CO,'ORTHOGONAL',2,noco+1,NTCO,N3CO)) THEN
              CUBIC=.FALSE.
              QUADRATIC=.FALSE.
              ORTHOG=.TRUE.
            ELSE IF(CBBREV(CO,'CUBIC',2,noco+1,NTCO,N3CO)) THEN
              CUBIC=.TRUE.
              QUADRATIC=.FALSE.
              ORTHOG=.FALSE.
            ELSE IF(CBBREV(CO,'QUADRATIC',2,noco+1,NTCO,N3CO)) THEN
              QUADRATIC=.TRUE.
              CUBIC=.FALSE.
              ORTHOG=.FALSE.
            ELSE
              CUBIC=.FALSE.
              QUADRATIC=.FALSE.
              ORTHOG=.FALSE.
            ENDIF
            IF(CBBREV(CO,'NEW',1,noco+1,NTCO,N3CO)) THEN
              NEW=.TRUE.
            ELSE IF(CBBREV(CO,'OLD',2,noco+1,NTCO,N3CO)) THEN
              NEW=.FALSE.
            ELSE
              NEW=.TRUE.
            ENDIF
            IF(CBBREV(CO,'DEFORM',3,NOCO+1,NTCO,N3CO)) THEN
              DEFORM=.TRUE.
	      IF(NTCO.GE.N3CO+1)THEN
	        CALL PARSIL(CO(N3CO+1),3,NTILDEF,DEFLIST,ERROR,*9999)
	        IF(DOP)THEN
	  	  WRITE(*,*)' DEFLIST=',DEFLIST
                ENDIF
              ELSE
		NTILDEF=3
		DO i=1,3
		  DEFLIST(i)=i
                ENDDO
              ENDIF
            ENDIF
	    IF(CBBREV(CO,'MAX',3,NOCO+1,NTCO,N3CO)) THEN
	      CALL PARSRL(CO(N3CO+1),1,NTIL,SQMAX,ERROR,*9999)
	      SQMAX=SQMAX*SQMAX
            ELSE
	      SQMAX = 1.d4*1.d4
            ENDIF
	    IF(CBBREV(CO,'ACCEPT',3,NOCO+1,NTCO,N3CO)) THEN
	      CALL PARSRL(CO(N3CO+1),1,NTIL,XI3OFF,ERROR,*9999)
            ELSE
	      XI3OFF = 0.d0
            ENDIF
	    XI3MIN = -XI3OFF
	    XI3MAX = 1.d0 + XI3OFF
            IF(CBBREV(CO,'EXTRAPOLATE',2,noco+1,NTCO,N3CO)) THEN
              EXTRAPOLATE=.TRUE.
            ELSE
              EXTRAPOLATE=.FALSE.
            ENDIF
            IF(CBBREV(CO,'CENTROID',3,noco+1,NTCO,N3CO)) THEN
              CENTROID=.TRUE.
              IF(CBBREV(CO,'CLUSTER',3,noco+1,NTCO,N3CO)) THEN
                CALL PARSIL(CO(N3CO+1),1,NTIL,NCLUSTER,ERROR,*9999)
              ELSE
                NCLUSTER=72
              ENDIF
            ELSE
              CENTROID=.FALSE.
            ENDIF
!news       AAY 20 March 95 xi_3 keyword to specify what plane of xi_3
!           to project the data onto
!           also added "add" option to def da;c xi
            IF(CBBREV(CO,'XI_3',4,NOCO+1,NTCO,N3CO)) THEN
              SET_XI_3 = .TRUE.
              CALL PARSRL(CO(N3CO+1),1,NTIL,XI_3,ERROR,*9999)
            ELSE
              XI_3 = 0.0d0
              SET_XI_3 = .FALSE.
            ENDIF
            IF(ADD) THEN
              !calculate xi for NDTOLD to NDT
              ND0=NDTOLD+1
              ND1=NDT
            ELSE
              ND0=1
              ND1=NDT
            ENDIF
            TAG2D=.FALSE.              
            CONTOUR=.FALSE.              
          ENDIF !beads/xi
!newe
c       ENDIF !NTCO.gt.noco+1

        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(FILIO) THEN
C CPB 8/4/94 This needs to be generalised for NJ_LOC
          CALL TRIM(FILE,IBEG,IEND)

          IF(TYPE(1:5).EQ.'XIPOS') THEN
            IF(CBBREV(CO,'TOTAL',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),1,NTIL,NITB,ERROR,*9999)
            ELSE IF(NEELEM(1,1).NE.0.AND.NBJ(1,NEELEM(1,1)).NE.0)THEN
              NITB=NIT(NBJ(1,NEELEM(1,1)))
            ELSE
              NITB=NJT-1
              NITB=NITB+JTYP9+JTYP11
            ENDIF       
            NXIDEF=NITB
            CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipxi',STATUS,
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
            CALL IOXI(IFILE,LD,NITB,XID,ERROR,*9999)
            CALC_XI=.TRUE.
            IF(IOTYPE.EQ.2) THEN
              CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
            ENDIF
            CALL CLOSEF(IFILE,ERROR,*9999)
C           calc #points per frame (used in time history plots)
            TREF=XID(NITB,1)
            nd=1
            DO WHILE(TREF.EQ.XID(NITB,nd).and.nd.LE.NDT)
              nd=nd+1
            ENDDO
            NPPF=nd-1
          ENDIF !xipos
            
        ELSE IF(CALCU) THEN !Calculate data position info

            NXIDEF=NIT(NBJ(1,NEELEM(1,1)))
            L=0
            DO noelem=1,NEELEM(0,1) !to init LN (recalc.d by define fit)
              ne=NEELEM(noelem,1)
              L=L+1
              LN(L)=ne   
            ENDDO
            LN(0)=L

          IF(TYPE(1:5).EQ.'XIPOS') THEN
            CALL ASSERT(NET(1).GT.0,'>>no elements defined',ERROR,*9999)
            IF(ITYP10(1).EQ.1) THEN !rect. cart.
              NJ1=1
              NJ2=2
              IF(CENTROID) THEN
C**             find centroids and store in ZDD(nd,1..3)
                NTIMES=NDT/NCLUSTER
                DO nt=1,NTIMES
                  Z1(1)=0.0D0
                  Z1(2)=0.0D0
                  Z1(3)=0.0D0
                  DO nd=NCLUSTER*(nt-1)+1,NCLUSTER*nt
                    DO nj=1,NJT
                      Z1(nj)=Z1(nj)+ZD(nj,nd)
                    ENDDO
                  ENDDO
C                 centroid is
                  DO nj=1,NJT
                    Z1(nj)=Z1(nj)/NCLUSTER
                  ENDDO
                  DO nd=NCLUSTER*(nt-1)+1,NCLUSTER*nt
                    DO nj=1,NJT
                      ZDD(nd,nj)=Z1(nj)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF !centroid

            ELSE IF(ITYP10(1).GE.2) THEN
              IF(NJT.EQ.2) THEN
                NJ1=1
                NJ2=2
              ELSE IF(NJT.EQ.3) THEN
                NJ1=2
                NJ2=3
              ENDIF
            ENDIF
C
C CPB 20/1/93 - commenting out this if block so that the old data
C calculation routine is used.
C 
C            IF(ORTHOG.AND.(.NOT.NEW)) THEN           !AAY 2-JUL-90
C              !loop over data points then elements
C              !first find neighbouring elements
C              CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999)
C              !for each data point do
C              DO nd=1,NDT
C                IF(INLIST(LD(nd),NELIST(1),NELIST(0),n1list)) THEN
C                  IF(DOP) THEN
C                    WRITE(OP_STRING,
C     '                 '(//,'' *****>Data point '',I5)') nd
C                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  ENDIF
C                  !initialise ne and xi to starting position
C                  ne=LD(nd)
C                  XI(1)=XID(1,nd)
C                  XI(2)=XID(2,nd)
C                  FINISHED=.FALSE.
C                  DO WHILE(.NOT.FINISHED)
C                    IF(DOP) THEN
C                      WRITE(OP_STRING,
C     '                  '(/'' Element '',I4,'' Xi ='',2F10.4)')
C     '                  ne,XI(1),XI(2)
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF
C                    CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),
C     '                NPNE(1,1,ne),NPF(1,1),NQE(1,1,ne),
C     '                NRE(ne),NVJE(1,1,1,ne),
C     '                SE(1,1,ne),XA,XE,XP,ERROR,*9999)
C                    NITB=NIT(NBJ(NJ1,ne))
C                    ITMAX=20
C                    IF(NITB.EQ.1) THEN !1D elements
C                    ELSE IF(NITB.EQ.2) THEN !2D elements
C                      IF(ITYP10(1).EQ.1) THEN !rect cart coords
C                        CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
C     '                    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
C                      ENDIF
C                      IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0.AND
C     '                  .XI(2).GT.0.0D0.AND.XI(2).LE.1.0D0) THEN
C                        IF(IT.LT.ITMAX) THEN
C                          XID(1,nd)=XI(1)
C                          XID(2,nd)=XI(2)
C                          LD(nd)=ne
C                          SQ(nd)=SQND
C                          FINISHED=.TRUE.
C                          IF(DOP) THEN
C                            WRITE(OP_STRING,'('' Convergence reached''/'
C     '                        //''' nd='',I4,'' LD='',I4,'' XID='',
C     '                        2E12.6,'
C     '                        //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
C     '                        nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd),
C     '                        (ZD(nj,nd),nj=1,NJT)
C                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          ENDIF
C                        ELSE !still inside element but not converged
C                          WRITE(OP_STRING,'('' WARNING:'
C     '                      //' Convergence not reached in CLOS2X''/'
C     '                      //''' nd='',I4,'' LD='',I4,'' XID='',
C     '                      2E12.6,'
C     '                      //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
C     '                      nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd),
C     '                      (ZD(nj,nd),nj=1,NJT)
C                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                          FINISHED=.TRUE.
C                        ENDIF
C                      ELSE  !moved out of the element
CC                       save solution so far
C                        SQ(nd)=SQND
C                        XID(1,nd)=XI(1)
C                        XID(2,nd)=XI(2)
C                        LD(nd)=ne
CC                       find the element moved to
C                        IF(XI(1).LE.0.0D0) THEN
C                          XID(1,nd)=0.0D0
C                          XI(1)=1.0D0
C                          ne=NXI(-1,ne)
C                        ELSE IF(XI(1).GT.1.0D0) THEN
C                          XID(1,nd)=1.0D0
C                          XI(1)=0.0D0
C                          ne=NXI(1,ne)
C                        ENDIF
C                        IF(XI(2).LE.0.0D0) THEN
C                          XID(2,nd)=0.0D0
C                          XI(2)=1.0D0
C                          ne=NXI(-2,ne)
C                        ELSE IF(XI(2).GT.1.0D0) THEN
C                          XID(2,nd)=1.0D0
C                          XI(2)=0.0D0
C                          ne=NXI(2,ne)
C                        ENDIF
C                        FINISHED=(ne.EQ.0)
C                        IF(DOP) THEN
C                          WRITE(OP_STRING,
C     '                      '(/'' Element '',I4,'' Xi ='',
C     '                      2F10.4)') ne,XI(1),XI(2)
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                        ENDIF
C                      ENDIF !still in element
C                    ENDIF !2d element
C                  ENDDO !until finished
C                ENDIF !in element list
C              ENDDO  !for all data points
C
C            ELSE IF(CENTROID) THEN

            IF(CENTROID) THEN
              IF(NEW) THEN !assign approximate element locations
                DO nd=1,NDT !to initialize arrays
                  LD(nd)=0
                  SQ(nd)=0.0D0
                ENDDO
C               do for all elements
                DO nolist=1,NELIST(0)
                  ne=NELIST(nolist)
                  nb=NBJ(NJ1,ne)
                  NITB=NIT(NBJ(NJ1,ne))
                  IF(NITB.EQ.2) THEN          !2D elements
C                   do for all unassociated data points
                    DO nd=1,NDT
                      IF(LD(nd).EQ.0) THEN
C                       get vectors from centroid for all nodes of ne
                        DO nj=1,NJT
                          XD(nj)=ZD(nj,nd)-ZDD(nd,nj) !is proj vector
                          Z1(nj)=XP(1,nv,nj,NPNE(1,nb,ne))-ZDD(nd,nj)
                          Z2(nj)=XP(1,nv,nj,NPNE(2,nb,ne))-ZDD(nd,nj)
                          Z3(nj)=XP(1,nv,nj,NPNE(3,nb,ne))-ZDD(nd,nj)
                          Z4(nj)=XP(1,nv,nj,NPNE(4,nb,ne))-ZDD(nd,nj)
                        ENDDO
C                       convert to spherical polar coords
                        CALL ZX(3,XD,XD)
                        CALL ZX(3,Z1,Z1)
                        CALL ZX(3,Z2,Z2)
                        CALL ZX(3,Z3,Z3)
                        CALL ZX(3,Z4,Z4)
C                       do for both angles
                        DO IT=2,3
                          TMIN(IT-1)=DMIN1(Z1(IT),Z2(IT),Z3(IT),Z4(IT))
                          TMAX(IT-1)=DMAX1(Z1(IT),Z2(IT),Z3(IT),Z4(IT))
C                         assume element spans smallest angle between nodes
                          IF(TMAX(IT-1)-TMIN(IT-1).GT.
     '                      TMIN(IT-1)+2.0D0*PI-TMAX(IT-1)) THEN
                            IF(Z1(IT).LT.0.0D0) Z1(IT)=Z1(IT)+2.0D0*PI
                            IF(Z2(IT).LT.0.0D0) Z2(IT)=Z2(IT)+2.0D0*PI
                            IF(Z3(IT).LT.0.0D0) Z3(IT)=Z3(IT)+2.0D0*PI
                            IF(Z4(IT).LT.0.0D0) Z4(IT)=Z4(IT)+2.0D0*PI
                            IF(XD(IT).LT.0.0D0) XD(IT)=XD(IT)+2.0D0*PI
                            TMIN(IT-1)=DMIN1(Z1(IT),Z2(IT),Z3(IT),
     '                        Z4(IT))
                            TMAX(IT-1)=DMAX1(Z1(IT),Z2(IT),Z3(IT),
     '                        Z4(IT))
                          ENDIF
                        ENDDO
C                       IF(DOP)WRITE(*,*)' Data point',nd,'element',ne
C                       IF(DOP)WRITE(*,*)tmin(1),tmax(1),tmin(2),
C    '                    tmax(2),xd(2),xd(3)
                        IF(XD(2).GE.TMIN(1).AND.XD(2).LE.TMAX(1).AND.
     '                     XD(3).GE.TMIN(2).AND.XD(3).LE.TMAX(2)) THEN
C                         associate data point with element
                          LD(nd)=ne
                          XID(1,nd)=0.5D0
                          XID(2,nd)=0.5D0
                          IF(DOP) THEN
                            WRITE(OP_STRING,*)' Data point',nd,
     '                        ' is ','in element',ne
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !for all data points
                  ENDIF !2D element
                ENDDO !for all elements
              ENDIF !new
C             find neighbouring elements
              CALL NENXI(IBT,INP,NBJ,NEELEM,NPNE,NXI,ERROR,*9999)

C             copy data into common block for proj_funct and proj_hess
c moved to FE_ARCHIVE and EDATA00.cmn removed CPB 7Oct94
c             DO nb=1,NBT
c               DO ni=1,NIT(nb)
c                 DO II=1,2
c                   IBT2(II,ni,nb)=IBT(II,ni,nb)
c                 ENDDO
c                 DO nk=1,NKT(0,nb)
c                   IDO2(nk,ni,nb)=IDO(nk,1,ni,nb)
c                 ENDDO
c                 DO nn=1,NNT(nb)
c                   INP2(nn,ni,nb)=INP(nn,ni,nb)
c                 ENDDO
c               ENDDO
c               DO nk=1,NKT(0,nb)
c                 IDO2(nk,0,nb)=IDO(nk,1,0,nb)
c               ENDDO
c             ENDDO
c             DO ne=1,NEELEM(0,1)
c               DO nj=1,NJT
c                 NBJ2(nj,ne)=NBJ(nj,ne)
c               ENDDO
c             ENDDO

C             now loop over data points
              DO nd=1,NDT
                IF(LD(nd).EQ.0) THEN
                  LD(nd)=1
                  XID(1,nd)=0.5D0
                  XID(2,nd)=0.5D0
                ENDIF
                IF(INLIST(LD(nd),NELIST(1),NELIST(0),n1list)) THEN
                  IF(DOP) THEN
                    WRITE(OP_STRING,
     '                '(//,'' *****>Data point '',I5)') nd
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                  !initialise ne and xi to starting position
                  ne=LD(nd)
                  XI(1)=XID(1,nd)
                  XI(2)=XID(2,nd)
                  XI(3)=1.0D0
                  FINISHED=.FALSE.
                  DO WHILE(.NOT.FINISHED)
                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '(/'' Element '',I4,'' Xi ='',3F10.4)')
     '                  ne,XI(1),XI(2),XI(3)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),
     '                NPF(1,1),NPNE(1,1,ne),NQE(1,1,ne),
     '                NRE(ne),NVJE(1,1,1,ne),
     '                SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                    NITB=NIT(NBJ(NJ1,ne))
                    ITMAX=20
                    IF(NITB.EQ.1) THEN !1D elements
                    ELSE IF(NITB.EQ.2) THEN !2D elements
                      IF(ITYP10(1).EQ.1) THEN !rect cart coords
                        DO nj=1,NJT
                          Z1(nj)=ZD(nj,nd)
                          Z2(nj)=ZDD(nd,nj)
                        ENDDO
                        ERROR='>> PROJ21 Moved to FE_ARCHIVE'
                        GOTO 9999
C                        CALL PROJ21(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),ne,
C     '                    SQND,XE,XI,Z1,Z2,ERROR,*9999)
                      ENDIF
                      IF(DABS(SQND).LT.10D-1) THEN !converged
                        XID(1,nd)=XI(1)
                        XID(2,nd)=XI(2)
                        LD(nd)=ne
                        SQ(nd)=SQND
                        FINISHED=.TRUE.
                        IF(DOP) THEN
                          WRITE(OP_STRING,
     '                      '('' Convergence reached''/'
     '                      //''' nd='',I4,'' LD='',I4,'' XID='','
     '                      //'2E12.6,'//''' SQ='',E12.6,'' ZD='','
     '                      //'3(E12.6,1X))') 
     '                      nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd),
     '                      (ZD(nj,nd),nj=1,NJT)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ELSE IF(XI(1).GT.0.0D0.AND.XI(1).LT.1.0D0.AND
     '                  .XI(2).GT.0.0D0.AND.XI(2).LT.1.0D0) THEN
                        !still inside element and intersection not found
                        WRITE(OP_STRING,'('' WARNING:'
     '                    //' Convergence not reached in PROJ2X''/'
     '                    //''' nd='',I4,'' LD='',I4,'' XID='',2E12.6,'
     '                    //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
     '                    nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd),
     '                    (ZD(nj,nd),nj=1,NJT)
                        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        !save current, set ld(nd) to zero as flag
                        XID(1,nd)=XI(1)
                        XID(2,nd)=XI(2)
                        SQ(nd)=SQND
                        LD(nd)=0
                        FINISHED=.TRUE.
                      ENDIF
                      IF(.NOT.FINISHED) THEN !at a boundary
                        !shift to neighbour element & save solution so far
                        SQ(nd)=SQND
                        XID(1,nd)=XI(1)
                        XID(2,nd)=XI(2)
                        LD(nd)=ne
                        !find the element moved to
                        IF(XI(1).LE.0.0D0) THEN
                          XID(1,nd)=0.0D0
                          XI(1)=1.0D0
                          ne=NXI(-1,ne)
                        ELSE IF(XI(1).GE.1.0D0) THEN
                          XID(1,nd)=1.0D0
                          XI(1)=0.0D0
                          ne=NXI(1,ne)
                        ENDIF
                        IF(XI(2).LE.0.0D0) THEN
                          XID(2,nd)=0.0D0
                          XI(2)=1.0D0
                          ne=NXI(-2,ne)
                        ELSE IF(XI(2).GE.1.0D0) THEN
                          XID(2,nd)=1.0D0
                          XI(2)=0.0D0
                          ne=NXI(2,ne)
                        ENDIF
                        FINISHED=(ne.EQ.0)
                        IF(DOP) THEN 
                          WRITE(OP_STRING,'(/'' Move to Element '','
     '                      //'I4,'' Xi ='',2F10.4)') ne,XI(1),XI(2)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF !still in element
                    ENDIF !2d element
                  ENDDO !until finished
                ENDIF !in element list
              ENDDO  !for all data points

C
C CPB 22/1/93 Old data calculation routines
C
C            ELSE !redo this later to loop first over data points? AAY 2-JUL-90
C              DO nd=1,NDT !to initialize arrays
C                LD(nd)=0
C                SQ(nd)=0.0D0
C              ENDDO
C              DO nolist=1,NELIST(0)
C                ne=NELIST(nolist)
C                IF(DOP) THEN
C                  WRITE(OP_STRING,'(/'' Element '',I4)') ne
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                ENDIF
C                CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPNE(1,1,ne),
C     '            NPF(1,1),NQE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
C     '            SE(1,1,ne),XA,XE,XP,ERROR,*9999)
C                NITB=NIT(NBJ(NJ1,ne))
C                IF(NITB.EQ.1) THEN          !1D elements
C                  IF(ORTHOG) THEN           !use nonlinear Xi calculation
C                    DO nd=1,NDT
C                      IF(DOP) THEN
C                        WRITE(OP_STRING,'('' Data point '',I5)') nd
C                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                      ENDIF
C                      IF(NEW) THEN !start Xi from centre of element
C                        XI(1)=0.5D0
C                      ELSE         !start Xi from previous values
C                        XI(1)=XID(1,nd)
C                      ENDIF
C                      ITMAX=10
C                      IF(ITYP10(1).EQ.1) THEN      !rectanglar cartesian
C                        CALL CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
C     '                    SQND,XE,XI(1),ZD(1,nd),ERROR,*9999)
C                      ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
C                      ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
C                      ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
C                      ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
C                      ENDIF
C                      IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0) THEN
C                        IF(DOP) THEN
C                          WRITE(OP_STRING,
C     '                      '('' Xi= '',F7.3,'' in element'',
C     '                      I4,'' SQ='',E12.3)') XI(1),ne,SQND
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                        ENDIF
C                        IF(LD(nd).EQ.0.OR.SQND.LT.SQ(nd)) THEN
C                          XID(1,nd)=XI(1)
C                          LD(nd)=ne
C                          SQ(nd)=SQND
C                        ENDIF
C                      ENDIF
C                      IF(IT.GE.ITMAX) THEN
C                        WRITE(OP_STRING,'('' WARNING:'
C     '                    //' Convergence not reached in CLOS2X''/'
C     '                    //''' nd='',I4,'' LD='',I4,'' XID='',E12.6,'
C     '                    //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
C     '                    nd,LD(nd),XID(1,nd),SQ(nd),(ZD(nj,nd),nj=1,
C     '                    NJT)
C                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                      ENDIF
C                    ENDDO
C
C                  ELSE IF(.NOT.ORTHOG) THEN !use linear Xi calculation
C                    IF(ITYP10(1).EQ.2) THEN !cylindrical coords
C                      nb=NBJ(NJ2,ne) !is basis for theta
C                    ELSE
C                      nb=NBJ(NJ1,ne)
C                    ENDIF
C                    NKTB=NKT(0,nb)
C                    NS1=1
C                    NS2=1+NKTB
C                    A1=XE(NS2,NJ1)-XE(NS1,NJ1)
C                    IF(NJT.GT.2) THEN
C                      nb=NBJ(NJ2,ne)
C                      NKTB=NKT(0,nb)
C                      NS1=1
C                      NS2=1+NKTB
C                      A2=XE(NS2,NJ2)-XE(NS1,NJ2)
C                    ENDIF
C                    IF(DOP) THEN
C                      WRITE(OP_STRING,'(/'' Element '',I3)') ne
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                      WRITE(OP_STRING,'('' A1,A2='',2E11.3)') A1,A2
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF
C                    DO nd=1,NDT
C                      CALL ZX(ITYP10(1),ZD(1,nd),XD)
C                      IF(NJT.EQ.2.AND.ITYP10(1).EQ.1) THEN !2D rect.cart.
C                        D1=XD(NJ1)-XE(NS1,NJ1)
C                        D2=XD(NJ2)-XE(NS1,NJ2)
C                        IF(DABS(A2).LT.1.0D-6) THEN
C                          XI(1)=D1/A1
C                        ELSE IF(DABS(A1).LT.1.0D-6) THEN
C                          XI(1)=D2/A2
C                        ELSE
C                          SLOPE=A2/A1
C                          XI(1)=(SLOPE*D2+D1)/(SLOPE*A2+A1)
C                        ENDIF
C                        IF(XI(1).LE.1.0D0.AND.XI(1).GE.0.0D0) THEN
C                          X1=XE(NS1,NJ1)*(1.0D0-XI(1))+XE(NS2,NJ1)*
C     '                      XI(1)
C                          X2=XE(NS1,NJ2)*(1.0D0-XI(1))+XE(NS2,NJ2)*
C     '                      XI(1)
C                          DIST=(X1-XD(NJ1))**2+(X2-XD(NJ2))**2
C                          IF(LD(nd).EQ.0.OR.DIST.LT.EDD(nd)) THEN
C                            XID(1,nd)=XI(1)
C                            LD(nd)=ne
C                            EDD(nd)=DIST
C                          ENDIF
C                        ENDIF
C                        IF(DOP) THEN
C                          WRITE(OP_STRING,'('' ne='',I4,'' nd='',I4,
C     '                      '' Xi(1)='','//'E11.3,'' EDD(nd)='',
C     '                      E11.3)') ne,nd,XI(1),EDD(nd)
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                        ENDIF
C
C                      ELSE IF(NJT.EQ.2.AND.ITYP10(1).EQ.2) THEN !2D cyl.polar
C                        XI(1)=(XD(2)-XE(NS1,2))/(XE(NS2,2)-XE(NS1,2))
C                        IF(XI(1).LE.1.0D0.AND.XI(1).GE.0.0D0) THEN
C                          X1=XE(NS1,1)*(1.0D0-XI(1))+XE(NS2,1)*XI(1)
C                          DIST=(X1-XD(1))**2
C                          IF(LD(nd).EQ.0.OR.DIST.LT.EDD(nd)) THEN
C                            XID(1,nd)=XI(1)
C                            LD(nd)=ne
C                            EDD(nd)=DIST
C                          ENDIF
C                        ENDIF
C                        IF(DOP) THEN
C                          WRITE(OP_STRING,
C     '                      '('' ne='',I4,'' nd='',I4,'' Xi(1)='','
C     '                      //'E11.3,'' EDD(nd)='',
C     '                      E11.3)') ne,nd,XI(1),EDD(nd)
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                        ENDIF
C                      ELSE IF(NJT.EQ.3.OR.ITYP10(1).NE.1) THEN
C                      ENDIF
C                    ENDDO
C                  ENDIF
C
C                ELSE IF(NITB.GE.2) THEN     !2D elements
C                  IF(ORTHOG) THEN           !use nonlinear Xi calculation
C                    DO nd=1,NDT
C                      IF(DOP) THEN
C                        WRITE(OP_STRING,'('' Data point '',I5)') nd
C                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                      ENDIF
C                      IF(NEW) THEN !start Xi from centre of element
C                        XI(1)=0.5D0
C                        XI(2)=0.5D0
C                      ELSE         !start Xi from previous values
C                        XI(1)=XID(1,nd)
C                        XI(2)=XID(2,nd)
C                      ENDIF
C                      ITMAX=10
C                      IF(ITYP10(1).EQ.1) THEN !rectanglar cartesian
C                        CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
C     '                    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
C                      ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
C                      ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
C                      ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
C                        CALL CLOS24(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
C     '                    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
C                      ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
C                      ENDIF
C                      IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0.AND
C     '                  .XI(2).GT.0.0D0.AND.XI(2).LE.1.0D0) THEN
C                        IF(DOP) THEN
C                          WRITE(OP_STRING,
C     '                      '('' Xi: '',2F7.3,'' in element'',
C     '                      I4,'' SQ='',E12.3)') XI(1),XI(2),ne,SQND
C                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                        ENDIF
C                        IF(LD(nd).EQ.0.OR.SQND.LT.SQ(nd)) THEN
C                          XID(1,nd)=XI(1)
C                          XID(2,nd)=XI(2)
C                          LD(nd)=ne
C                          SQ(nd)=SQND
C                        ENDIF
C                      ENDIF
C                      IF(IT.GE.ITMAX) THEN
C                        WRITE(OP_STRING,'('' WARNING:'
C     '                    //' Convergence not reached in CLOS2X''/'
C     '                    //''' nd='',I4,'' ZD='',3(E12.6,1X)/'
C     '                    //''' LD='',I4,'' XID'',2(E12.6,1X),''SQ='',
C     '                    E12.6/)')nd,(ZD(nj,nd),nj=1,NJT),
C     '                    LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
C                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                      ENDIF
C                    ENDDO
C
C
C                  ELSE IF(CUBIC) THEN !use 1D nonlinear calculation
C                    nb=NBJ(NJ1,ne) !nb for first geometric coord
C                    NKTB=NKT(0,nb)
C                    DO nd=1,NDT
C                      IF(LD(nd).EQ.0) THEN
C                        CALL ZX(ITYP10(1),ZD(1,nd),XD)
CC ...                     Find XID([1,2],nd) for a point inside cubic hermite-
CC                         linear element
C                        XPXIF=ZD(1,nd)    !set x and y of data point
C                        YPXIF=ZD(2,nd)
C                        DO I=1,8
C                          XNODVAL(I,1)=XE(I,1)
C                          XNODVAL(I,2)=XE(I,2)
C                        ENDDO
C                        KTYP26=4   !inputs to e04jbf
C                        KTYP27=2
C                        NPXIF=1
C                        IPRINT=-1
C                        LOCSCH=.FALSE.
C                        INTYPE=0
C                        MAXCAL=240
C                        ETA=.5
C                        XTOL=1.0D-5
C                        STEPMX=2.0D0
C                        FEST=0.0D0
C                        DELXIF(1)=1.0D-10
C                        IBOUND=0
C                        BLXIF(1)=0.0000001D0
C                        BUXIF(1)=0.9999999D0
C                        XIF(1)=0.5D0
C                        LH=5   !can be 1
C                        LIW=5
C                        LW=18
C                        IFAIL=1
C                        CALL E04JBF(NPXIF,FUNCT1,MONIT,IPRINT,LOCSCH,
C     '                    INTYPE,E04JBQ,
C     '                    MAXCAL,ETA,XTOL,STEPMX,FEST,DELXIF,IBOUND,
C     '                    BLXIF,BUXIF,XIF,HESL,LH,HESD,ISTATE,FXIF,GXIF,
C     '                    IWXIF,LIW,WW,LW,IFAIL)
C                        XI(1)=XIF(1)
C
CC                       Find X_ENDO,X_EPI,Y_ENDO,Y_EPI using converged
CC                       value of XI1. Note most recent values of these
CC                       variables computed by NAG do not necessarily
CC                       correlate with converged xi1
C                        XI1=XI(1)
C                        S01=1.0D0-3.0D0*XI1*XI1+2.0D0*XI1*XI1*XI1  !evaluate hermite basis function
C                        S02=3.0D0*XI1*XI1-2.0D0*XI1*XI1*XI1
C                        S11=XI1-2.0D0*XI1*XI1+XI1*XI1*XI1
C                        S12=-(XI1*XI1-XI1*XI1*XI1)
C                        F01=XNODVAL(1,1)                  !Nodal values
C                        F02=XNODVAL(3,1)
C                        F11=XNODVAL(2,1)
C                        F12=XNODVAL(4,1)
C                        X_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
C                        F01=XNODVAL(1,2)
C                        F02=XNODVAL(3,2)
C                        F11=XNODVAL(2,2)
C                        F12=XNODVAL(4,2)
C                        Y_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
C                        F01=XNODVAL(5,1)
C                        F02=XNODVAL(7,1)
C                        F11=XNODVAL(6,1)
C                        F12=XNODVAL(8,1)
C                        X_EPI=S01*F01+S02*F02+S11*F11+S12*F12
C                        F01=XNODVAL(5,2)
C                        F02=XNODVAL(7,2)
C                        F11=XNODVAL(6,2)
C                        F12=XNODVAL(8,2)
C                        Y_EPI=S01*F01+S02*F02+S11*F11+S12*F12
C
CC ...                   Find value of objective function at converged XI(1)
C                        XDIFFERENCE=DABS(X_EPI-X_ENDO)
C                        IF(XDIFFERENCE.LT.1.0D-8) THEN     !singular case
C                          XNPB=X_EPI
C                          YNPB=YPXIF
C                        ELSE                                !Nonsingular case
C                          S1=X_ENDO-X_EPI
C                          S2=Y_ENDO-Y_EPI
C                          XNPB=(XPXIF*S1*S1+X_ENDO*S2*S2-S1*S2
C     '                      *(Y_ENDO-YPXIF))/(S1*S1+S2*S2)
C                          YNPB=Y_ENDO+S2*(XNPB-X_ENDO)/S1
C                        ENDIF
C                        FC=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2
C
CC ...                   flag points which are outside element
C                        IF(FC.GT.1.0D-14) XI(1)=99.0D0
C                        XDIFFERENCE=DABS(X_ENDO-X_EPI)
C                        IF(XDIFFERENCE.LT.1.0D-7) THEN
C                          XI(2)=(YPXIF-Y_ENDO)/(Y_EPI-Y_ENDO)
C                        ELSE
C                          XI(2)=(XPXIF-X_ENDO)/(X_EPI-X_ENDO)
C                        ENDIF
C
C                        IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0.AND.
C     '                     XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
C                          IF(NITB.EQ.2) THEN !2D basis
C                            IF(DOP) THEN
C                              WRITE(OP_STRING,'('' nd='',I4,'' Xi:'',
C     '                          3E11.3,'' (linear calc.)'')') nd,
C     '                          (XI(ni),ni=1,3)
C                              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                            ELSE
C                            ENDIF
C                            nb=NBJ(NJ1,ne)
C                            XID(1,nd)=XI(1)
C                            XID(2,nd)=XI(2)
C                            LD(nd)=ne !element no for data point nd
C                          ELSE IF(NITB.EQ.3) THEN !3D basis
C                          ENDIF
C                        ENDIF
C                      ENDIF !LD
C                    ENDDO !nd
C
C
C                  ELSE IF(QUADRATIC) THEN !use 2D nonlinear calculation
C                    nb=NBJ(NJ1,ne) !nb for first geometric coord
C                    NKTB=NKT(0,nb)
C                    DO nd=1,NDT
C                      IF(LD(nd).EQ.0) THEN
C                        CALL ZX(ITYP10(1),ZD(1,nd),XD)
CC ...                   Find XID([1,2],nd) for a point inside Lagrange-
CC                         quadratic element
C                        XPXIF=ZD(1,nd)            !set x and y of data point
C                        YPXIF=ZD(2,nd)
C
C                        DO I=1,9                  !Element Nodes
C                          XNODVAL(I,1)=XE(I,1)
C                          XNODVAL(I,2)=XE(I,2)
C                        ENDDO
C
C                        KTYP26=4   !inputs to e04jbf
C                        KTYP27=3
C                        NPXIF=2
C                        IPRINT=-1
C                        LOCSCH=.FALSE.
C                        INTYPE=0
C                        MAXCAL=560
C                        ETA=0.5D0
C                        XTOL=1.0D-5
C                        STEPMX=2.0D0
C                        FEST=0.0D0
C                        DELXIF(1)=1.0D-10
C                        DELXIF(2)=1.0D-10
C                        IBOUND=0
C                        BLXIF(1)=0.00001D0
C                        BUXIF(1)=0.99999D0
C                        BLXIF(2)=0.00001D0
C                        BUXIF(2)=0.99999D0
C                        XIF(1)=0.5D0
C                        XIF(2)=0.5D0
C                        LH=5   !can be 1
C                        LIW=5
C                        LW=18
C                        IFAIL=1
C
C                        CALL E04JBF(NPXIF,FUNCT1,MONIT,IPRINT,LOCSCH,
C     '                    INTYPE,E04JBQ,
C     '                    MAXCAL,ETA,XTOL,STEPMX,FEST,DELXIF,IBOUND,
C     '                    BLXIF,BUXIF,XIF,HESL,LH,HESD,ISTATE,FXIF,GXIF,
C     '                    IWXIF,LIW,WW,LW,IFAIL)
C                        XI(1)=XIF(1)          !XIF contain converged xi values
C                        XI(2)=XIF(2)
C
CC ...                   Recompute FC (which requires recomputing XNPB,YNPB)
C                        VL20=2.0D0*(XI(1)-0.5D0)*(XI(1)-1.0D0)
C                        VL21=-4.0D0*(XI(1)-1.0D0)*XI(1)
C                        VL22=2.0D0*XI(1)*(XI(1)-0.5D0)
C                        CL20=2.0D0*(XI(2)-0.5D0)*(XI(2)-1.0D0)
C                        CL21=-4.0D0*(XI(2)-1.0D0)*XI(2)
C                        CL22=2.0D0*XI(2)*(XI(2)-0.5D0)
C
C                        XNPB=VL20*CL20*XNODVAL(1,1)+VL21*CL20
C     '                                                  *XNODVAL(2,1)+
C     '                       VL22*CL20*XNODVAL(3,1)+VL20*CL21
C     '                                                  *XNODVAL(4,1)+
C     '                       VL21*CL21*XNODVAL(5,1)+VL22*CL21
C     '                                                  *XNODVAL(6,1)+
C     '                       VL20*CL22*XNODVAL(7,1)+VL21*CL22
C     '                                                  *XNODVAL(8,1)+
C     '                       VL22*CL22*XNODVAL(9,1)
C                        YNPB=VL20*CL20*XNODVAL(1,2)+VL21*CL20
C     '                                                  *XNODVAL(2,2)+
C     '                       VL22*CL20*XNODVAL(3,2)+VL20*CL21
C     '                                                  *XNODVAL(4,2)+
C     '                       VL21*CL21*XNODVAL(5,2)+VL22*CL21
C     '                                                  *XNODVAL(6,2)+
C     '                       VL20*CL22*XNODVAL(7,2)+VL21*CL22
C     '                                                  *XNODVAL(8,2)+
C     '                       VL22*CL22*XNODVAL(9,2)
C
C                        FC=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2  !value of obj. fn.
C
CC ...                   flag points which are outside element
C                         IF(FC.GT.1.0D-10) XI(1)=99.0D0
C                         IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0.AND.
C     '                      XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
C                          IF(NITB.EQ.2) THEN !2D basis
C                            IF(DOP) THEN
C                              WRITE(OP_STRING,'('' nd='',I4,'' Xi:'',
C     '                          3E11.3,'' (linear calc.)'')') nd,
C     '                          (XI(ni),ni=1,3)
C                              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                            ELSE
C                            ENDIF
C                            nb=NBJ(NJ1,ne)
C                            XID(1,nd)=XI(1)
C                            XID(2,nd)=XI(2)
C                            LD(nd)=ne !element no for data point nd
C                          ELSE IF(NITB.EQ.3) THEN !3D basis
C                          ENDIF
C                        ENDIF
C                      ENDIF !LD
C                    ENDDO !nd
C
C                  ELSE IF(.NOT.ORTHOG) THEN !use linear Xi calculation
C                    nb=NBJ(NJ1,ne) !nb for first geometric coord
C                    NKTB=NKT(0,nb)
C                    NS1=1
C                    NS2=1+NKTB
C                    NS3=1+2*NKTB
C                    NS4=1+3*NKTB
C                    A1=XE(NS2,NJ1)-XE(NS1,NJ1)
C                    B1=XE(NS3,NJ1)-XE(NS1,NJ1)
C                    C1=XE(NS1,NJ1)-XE(NS2,NJ1)-XE(NS3,NJ1)+XE(NS4,NJ1)
C                    nb=NBJ(NJ2,ne)
C                    NKTB=NKT(0,nb)
C                    NS1=1
C                    NS2=1+NKTB
C                    NS3=1+2*NKTB
C                    NS4=1+3*NKTB
C                    A2=XE(NS2,NJ2)-XE(NS1,NJ2)
C                    B2=XE(NS3,NJ2)-XE(NS1,NJ2)
C                    C2=XE(NS1,NJ2)-XE(NS2,NJ2)-XE(NS3,NJ2)+XE(NS4,NJ2)
C                    IF(DOP) THEN
C                      WRITE(OP_STRING,'(/'' Element '',I3)') ne
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                      WRITE(OP_STRING,'('' A1,B1,C1='',3E11.3)') A1,
C     '                  B1,C1
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                      WRITE(OP_STRING,'('' A2,B2,C2='',3E11.3)') A2,
C     '                  B2,C2
C                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                    ENDIF
C                    ALFA=A2*C1-A1*C2
C                    DO nd=1,NDT
C                      IF(LD(nd).EQ.0) THEN
C                        CALL ZX(ITYP10(1),ZD(1,nd),XD)
C                        EXCLUDE=.TRUE.
C                        IF(ITYP10(1).EQ.4) THEN !check if lies within element
C                          THETAMIN=DMIN1(XE(NS2,3),XE(NS4,3))
C                          THETAMAX=DMAX1(XE(NS1,3),XE(NS3,3))
C                          IF(XD(3).LT.THETAMAX.AND.
C     '                       XD(3).GT.THETAMIN) THEN
C                            EXCLUDE=.FALSE.
C                            D2=XD(NJ2)-XE(1,NJ2)
C                          ELSE IF((XD(3)+2.0D0*PI).LT.THETAMAX
C     '                       .AND.(XD(3)+2.0D0*PI).GT.THETAMIN) THEN
C                            EXCLUDE=.FALSE.
C                            D2=XD(NJ2)+2.0D0*PI-XE(1,NJ2)
C                          ELSE IF((XD(3)-2.0D0*PI).LT.THETAMAX
C     '                       .AND.(XD(3)-2.0D0*PI).GT.THETAMIN) THEN
C                            EXCLUDE=.FALSE.
C                            D2=XD(NJ2)-2.0D0*PI-XE(1,NJ2)
C                          ENDIF
C                        ELSE
C                          EXCLUDE=.FALSE.
C                          D2=XD(NJ2)-XE(1,NJ2)
C                        ENDIF
C                        IF(.NOT.EXCLUDE) THEN
C                          D1=XD(NJ1)-XE(1,NJ1)
C                          BETA=C2*D1-C1*D2+A2*B1-A1*B2
C                          GAMA=B2*D1-B1*D2
C                          DELTA=BETA*BETA-4.0D0*ALFA*GAMA
C                          IF(DOP) THEN
C                            WRITE(OP_STRING,'('' D1='',E11.3,'' D2='',
C     '                        E11.3,
C     '                        '' ALFA='',E11.3,'' BETA='',E11.3,
C     '                        '' GAMA='',E11.3,'' DELTA='',E11.3)')
C     '                        D1,D2,ALFA,BETA,GAMA,DELTA
C                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                          ENDIF
C                          IF(DABS(ALFA).GT.1D-6) THEN
C                            IF(DELTA.GE.0.0D0) THEN
C                              XI(1)=(-BETA+DSQRT(DELTA))/(2.0D0*ALFA)
C                              IF(XI(1).LT.0.0D0.OR.XI(1).GT.1.0D0) THEN
C                                XI(1)=(-BETA-DSQRT(DELTA))/(2.0D0*ALFA)
C                              ENDIF
C                            ENDIF
C                          ELSE IF(DABS(ALFA).LE.1D-6) THEN
C                            IF(DABS(BETA).GT.1D-6) THEN
C                              XI(1)=-GAMA/BETA
C                            ELSE
C                              XI(1)=-1.0D0
C                            ENDIF
C                          ENDIF
C                          IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0) THEN
C                            DENOM1=B1+C1*XI(1)
C                            IF(DABS(DENOM1).GT.1D-6) THEN
C                              XI(2)=(D1-A1*XI(1))/DENOM1
C                            ELSE
C                              DENOM2=B2+C2*XI(1)
C                              IF(DABS(DENOM2).GT.1D-6) THEN
C                                XI(2)=(D2-A2*XI(1))/DENOM2
C                              ELSE
C                                WRITE(OP_STRING,
C     '                            '('' Xi(2) cannot be defined for '
C     '                            //'nd='',I6,'' in element '',
C     '                            I5)') nd,ne
C                                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                              ENDIF
C                            ENDIF
C                            IF(XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
C                              IF(NITB.EQ.2) THEN
C                                IF(DOP) THEN
C                                  WRITE(OP_STRING,
C     '                              '('' nd='',I4,'' Xi:'',3E11.3,
C     '                              '' (linear calc.)'')') nd,(XI(ni),
C     '                              ni=1,3)
C                                  CALL WRITES(IODI,OP_STRING,ERROR,
C     '                              *9999)
C                                ENDIF
C                                nb=NBJ(NJ1,ne)
C                                IF(NKT(0,nb).GT.1) THEN !update Xi(2)
C                                  XI(2)=0.0D0
C                                  X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
C     '                              INP(1,1,nb),nb,1,XI,XE(1,NJ1))
C                                  XI(2)=1.0D0
C                                  X2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
C     '                              INP(1,1,nb),nb,1,XI,XE(1,NJ1))
C                                  IF(DABS(X2-X1).GT.1.0D-6) THEN !PJH 13APR91
C                                    XI(2)=(XD(NJ1)-X1)/(X2-X1)
C                                  ELSE
C                                    WRITE(OP_STRING,'('' Warning!!!!'',
C     '                                ''X1=X2'')')
C                                    CALL WRITES(IOOP,OP_STRING,ERROR,
C     '                                *9999)
C                                    XI(2)=1.0D0
C                                  ENDIF
C                                ENDIF
C                                XID(1,nd)=XI(1)
C                                XID(2,nd)=XI(2)
C                                LD(nd)=ne
C                              ELSE IF(NITB.EQ.3) THEN
C                                nb=NBJ(1,ne)
C                                XI(3)=0.0D0
C                                X0=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
C     '                            INP(1,1,nb),nb,1,XI,XE(1,1))
C                                XI(3)=1.0D0
C                                X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
C     '                            INP(1,1,nb),nb,1,XI,XE(1,1))
C                                DIFF=X1-X0
C                                IF(DABS(DIFF).GT.1D-6) THEN
C                                  XI(3)=(XD(1)-X0)/DIFF
C                                ELSE
C                                  XI(3)=0.0D0
C                                ENDIF
C                                IF(DOP) THEN
C                                  WRITE(OP_STRING,'('' X0='',E11.3,
C     '                              '' X1='',E11.3,'' XD(1)='',E11.3,
C     '                              '' Xi(3)='',E11.3)') X0,X1,
C     '                              XD(1),XI(3)
C                                  CALL WRITES(IODI,OP_STRING,ERROR,
C     '                              *9999)
C                                ENDIF
C                                IF(XI(3).GT.0.0D0.AND.XI(3).LE.1.0D0) THEN
C                                  LD(nd)=ne
C                                ELSE IF(EXTRAPOLATE.AND.
C     '                                
CXI(3).GT.1.0D0.AND.XI(3).LE.1.1D0) 
C     '                                 THEN
C      WRITE(OP_STRING,
C     '  '('' TAGGING: Found nd='',I6,'' with xi(3)='',F5.2)') nd,
C     '  XI(3)
C      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                                  IF(LD(nd).EQ.0) THEN
C                                    LD(nd)=-ne !to tag for 'extrapolate' below
C                                  ENDIF
C                                ENDIF
C                                XID(1,nd)=XI(1)
C                                XID(2,nd)=XI(2)
C                                XID(3,nd)=XI(3)
C                              ENDIF
C                              IF(DOP) THEN
C                                WRITE(OP_STRING,
C     '                            '('' nd='',I4,'' Xi:'',3E11.3)')
C     '                            nd,(XI(ni),ni=1,3)
C                                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                              ENDIF
C                            ENDIF !0<xi(2)<1
C                          ENDIF !0<xi(1)<1
C                        ENDIF !.NOT.EXCLUDE
C                      ENDIF !LD(nd)=0
C                    ENDDO !nd
C                  ENDIF !linear/orthog/cubic
C                ENDIF !nitb=1,2,3
C              ENDDO !nolist
C              IF(EXTRAPOLATE) THEN
C                DO nd=1,NDT
C                  IF(LD(nd).LT.0) THEN !data point currently not in any element
C                    LD(nd)=-LD(nd)
C                    XID(3,nd)=1.0D0
C                  ENDIF
C                ENDDO
C              ENDIF
C            ENDIF !orthog/centroid etc

!news AAY 9 Dec 94 
            ELSE IF(TAG2D) THEN
              !for all stripe data points do
              DO nd=STRIPE0,STRIPE1
                LD(nd)=0
                NT=NLS_NDATA_CONT(nd) !the tag identifier
                nimag =NLS_NDATA_IMAG(nd) !the image identifier
                !intersect the line segments with the image plane
                CALL CONT_PLANE_X(NBJ,NJE,NKE,NPNE,NPF,NQE,SE,XA,XE,ZP,
     '            NT,IMAGE_NORMAL(1,nimag),ZD(1,nd),LD(nd),XID(1,nd),
     '            SQMAX,ERROR,*9999)
              ENDDO !nd
            ELSE IF(CONTOUR) THEN
              NC=NTAG(0)+1 !use next tag data struct to store line segs
              CALL ASSERT(NC.LE.100,
     '          '>>Error: increase number of tags in /CONT01/',
     '          ERROR,*9999)
C             for all image planes do
              DO nimag=1,NIMAGE(0)
C               for all surfaces do
                DO nsurf=1,NTSURF
C!!! INDEX_SURF used before it is set
                  ns=INDEX_SURF(nsurf)
C                 intersect this surface with the image plane
c                 CALL SURF_PLANE_X(NBJ,NJE,NKE,NPNE,NPF,NQE,SE,
c    '              XA,XE,ZP,ns,NC,IMAGE_NORMAL(1,nimag),
c    '              IMAGE_POS(1,nimag),ERROR,*9999)
C!!! ICON_X used before it is set
C                 precalc the basis functions for the line segments
                  DO np=1,ICON_X(NC)
                    DO nj=1,NJT
C!!! ISURF_ELEM used before it is set
                      nb=NBJ(nj,ISURF_ELEM(1,ns))
                      DO n=1,2
                        DO nn=1,NNT(nb)
                          DO nk=1,NKT(0,nb)
                            NLS_CON_PSI(nk+(nn-1)*NKT(0,nb),n,np,nj,nc)
     '                        =PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
     '                        INP(1,1,nb),nb,1,nk,nn,
     '                        NLS_CON_XI(1,n,np,nc))
                          ENDDO !nk
                        ENDDO !nn
                      ENDDO  !n
                    ENDDO !nj
                  ENDDO !np
                  !for all data points on this surface and image do
                  DO nd=CONTOUR0,CONTOUR1
                    IF(NLS_NDATA_CONT(nd).EQ.ns
     '                .AND.NLS_NDATA_IMAG(nd).EQ.nimag) THEN
                      LD(nd)=0
!                     construct plane normal to image and contour
                      CALL CROSS(WD(1,nd),IMAGE_NORMAL(1,nimag),Z1)
                      !find the intersection with this plane 
                      CALL CONT_PLANE_X(NBJ,NJE,NKE,NPNE,NPF,NQE,SE,
     '                  XA,XE,ZP,NC,Z1,ZD(1,nd),LD(nd),XID(1,nd),
     '                  SQMAX,ERROR,*9999)
                    ENDIF
                  ENDDO !nd
                ENDDO !nsurf
              ENDDO !nimag
!newe
            ELSE !redo this later to 1st loop over data pts? AAY 2JUL90
C
C CPB 22/1/93 Reordering loops so that the outer loop is over the data
C points
C
!news         AAY initialize arrays to zero if new 23 March 95
              IF(NEW) THEN   !news AAY 15 Aug 91
                DO ND=ND0,ND1 !to initialize arrays
                  LD(ND)=0
                  SQ(ND)=0.d0
                ENDDO
              ENDIF       
!newe AAY 

!old	      DO nd=1,NDT
    	      DO nd=ND0,ND1  !new AAY 23 March 95 limits set with ADD keyword
      		IF(DOP) THEN
      		  WRITE(OP_STRING,'('' Data point '',I5)') nd
      		  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      		ENDIF
                FOUND=.FALSE.
 
C
C CPB 25/1/93 If the data projection has already been found previously, try
C and find the new projection in the same element
C
                IF(.NOT.NEW) THEN ! start Xi from previous values
	          IF(LD(nd).NE.0) THEN
                    IF(ORTHOG) THEN ! use nonlinear Xi calculation
      		      ne=LD(nd) ! start at the previous element
      		      IF(DOP) THEN
      			WRITE(OP_STRING,'(/'' Element '',I4)') ne
      			CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      		      ENDIF
                      IF(.NOT.DEFORM)THEN
      		        CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),
     '		          NPF(1,1),NPNE(1,1,ne),NQE(1,1,ne),
     '                    NRE(ne),NVJE(1,1,1,ne),
     '                    SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                      ELSE IF(DEFORM)THEN
c                       Note that deformed coords --> XE
                        CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '                    NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '                    NVHE(1,1,1,ne),NW(ne,1),nx,
     '                    SE(1,1,ne),ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
                      ENDIF
      		      NITB=NIT(NBJ(nj1,ne))
		      IF(NITB.EQ.1) THEN ! 1D elements
                        XI(1)=XID(1,nd)
      			ITMAX=10
      			IF(ITYP10(1).EQ.1) THEN      !rect cartesian
      			  CALL CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
     '			    SQND,XE,XI(1),ZD(1,nd),ERROR,*9999)
      			ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
      			ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
      			ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
      			ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
      			ENDIF
      			IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0) THEN
      			  IF(DOP) THEN
      			    WRITE(OP_STRING,
     '			      '('' Xi= '',F7.3,'' in element'','
     '			      //'I4,'' SQ='',E12.3)') XI(1),ne,SQND
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  ENDIF
      			  XID(1,nd)=XI(1)
      			  LD(nd)=ne
      			  SQ(nd)=SQND
                          FOUND=.TRUE.
      			ENDIF
      			IF(IT.GE.ITMAX) THEN
      			  WRITE(OP_STRING,'('' WARNING:'
     '			    //' Convergence not reached in CLOS2X''/'
     '			    //''' nd='',I4,'' LD='',I4,'' XID='',E12.6,'
     '			    //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
     '			    nd,LD(nd),XID(1,nd),SQ(nd),(ZD(nj,nd),nj=1,
     '			    NJT)
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      			ENDIF

		      ELSE IF(NITB.EQ.2) THEN !2D elements
                        XI(1)=XID(1,nd)
                        XI(2)=XID(2,nd)
      			ITMAX=10
      			IF(ITYP10(1).EQ.1) THEN !rectanglar cartesian
      			  CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
     '			    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
      			ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
      			ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
      			ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
      			  CALL CLOS24(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
     '			    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
      			ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
      			ENDIF
      			IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0.AND
     '			  .XI(2).GT.0.0D0.AND.XI(2).LE.1.0D0) THEN
      			  IF(DOP) THEN
      			    WRITE(OP_STRING,
     '			      '('' Xi: '',2F7.3,'' in element'',I4,'
     '                        //''' SQ='',E12.3)') XI(1),XI(2),ne,SQND
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  ENDIF
      			  XID(1,nd)=XI(1)
      			  XID(2,nd)=XI(2)
      			  LD(nd)=ne
      			  SQ(nd)=SQND
                          FOUND=.TRUE.
      			ENDIF
      			IF(IT.GE.ITMAX) THEN
      			  WRITE(OP_STRING,'('' WARNING:'
     '			    //' Convergence not reached in CLOS2X''/'
     '			    //''' nd='',I4,'' ZD='',3(E12.6,1X)/'
     '			    //''' LD='',I4,'' XID'',2(E12.6,1X),'
     '                      //'''SQ='',E12.6/)')nd,(ZD(nj,nd),nj=1,NJT),
     '			    LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
      			  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      			ENDIF ! it.ge.itmax

		      ELSE IF(NITB.EQ.3) THEN !3D elements
                        XI(1)=XID(1,ND)
                        XI(2)=XID(2,ND)
                        XI(3)=XID(3,ND)  !new AAY 30 May 1991
                        IN_ELEM = (LD(ND).EQ.NE) !new AAY 14 Aug 1991
                        ITMAX=15
                        IF(ITYP10(nr).EQ.1) THEN !rectanglar cartesian
			  IT=0
			  IF((NEW.AND.LD(ND).EQ.0).OR.
     '                      (.NOT.NEW.AND.IN_ELEM)) THEN
			    !not found yet
			    SQND = SQMAX
			    CALL CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ(1,NE),
     '                        SQND,XE,XI,ZD(1,ND),ERROR,*9999)
                            IF(SQND.GE.SQMAX) LD(nd)=0 !AAY 1 May 94
			  ELSE
			    XI(1)= -2.d0 !so it won't be put in xid
			  ENDIF
                        ENDIF !ityp10(nr)=1
      		      ENDIF !nitb=1,2,3

                    ENDIF !ld(nd).NE.0
	          ENDIF !orthog
	        ENDIF !not new
                IF(.NOT.FOUND) THEN !try and find the Xi point by
				    !looking at all the elements
      		  LD(nd)=0
      		  SQ(nd)=0.0D0
      		  XI(1)=0.5D0
      		  DO nolist=1,NELIST(0)
      		    ne=NELIST(nolist)
      		    IF(DOP) THEN
      		      WRITE(OP_STRING,'(/'' Element '',I4)') ne
      		      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      		    ENDIF
      		    CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),
     '                NPF(1,1),NPNE(1,1,ne),NQE(1,1,ne),
     '                NRE(ne),NVJE(1,1,1,ne),
     '                SE(1,1,ne),XA,XE,XP,ERROR,*9999)
      		    NITB=NIT(NBJ(NJ1,ne))
      		    IF(NITB.EQ.1) THEN !1D elements
      		      IF(ORTHOG) THEN  !use nonlinear Xi calculation
      			XI(1)=0.5D0 ! start at the centre of the element
      			ITMAX=10
      			IF(ITYP10(1).EQ.1) THEN      !rect cartesian
      			  CALL CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
     '			    SQND,XE,XI(1),ZD(1,nd),ERROR,*9999)
      			ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
      			ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
      			ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
      			ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
      			ENDIF
      			IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0) THEN
      			  IF(DOP) THEN
      			    WRITE(OP_STRING,
     '			      '('' Xi= '',F7.3,'' in element'','
     '			      //'I4,'' SQ='',E12.3)') XI(1),ne,SQND
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  ENDIF
      			  IF(LD(nd).EQ.0.OR.SQND.LT.SQ(nd)) THEN
      			    XID(1,nd)=XI(1)
      			    LD(nd)=ne
      			    SQ(nd)=SQND
      			  ENDIF
      			ENDIF
      			IF(IT.GE.ITMAX) THEN
      			  WRITE(OP_STRING,'('' WARNING:'
     '			    //' Convergence not reached in CLOS2X''/'
     '			    //''' nd='',I4,'' LD='',I4,'' XID='',E12.6,'
     '			    //''' SQ='',E12.6,'' ZD='',3(E12.6,1X))')
     '			    nd,LD(nd),XID(1,nd),SQ(nd),(ZD(nj,nd),nj=1,
     '			    NJT)
      			  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      			ENDIF

      		      ELSE IF(.NOT.ORTHOG) THEN !use linear Xi calcn
      			IF(ITYP10(1).EQ.2) THEN !cylindrical coords
      			  nb=NBJ(NJ2,ne) !is basis for theta
      			ELSE
      			  nb=NBJ(NJ1,ne)
      			ENDIF
      			NKTB=NKT(0,nb)
      			ns1=1
      			ns2=1+NKTB
      			A1=XE(ns2,NJ1)-XE(ns1,NJ1)
      			IF(NJT.GT.2) THEN
      			  nb=NBJ(NJ2,ne)
      			  NKTB=NKT(0,nb)
      			  ns1=1
      			  ns2=1+NKTB
      			  A2=XE(ns2,NJ2)-XE(ns1,NJ2)
      			ENDIF
      			IF(DOP) THEN
      			  WRITE(OP_STRING,'(/'' Element '',I3)') ne
      			  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  WRITE(OP_STRING,'('' A1,A2='',2E11.3)')
     '                      A1,A2
      			  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			ENDIF
      			CALL ZX(ITYP10(1),ZD(1,nd),XD)
      			IF(NJT.EQ.2.AND.ITYP10(1).EQ.1) THEN !2D rect.cart.
      			  D1=XD(NJ1)-XE(ns1,NJ1)
      			  D2=XD(NJ2)-XE(ns1,NJ2)
      			  IF(DABS(A2).LT.1.0D-6) THEN
      			    XI(1)=D1/A1
      			  ELSE IF(DABS(A1).LT.1.0D-6) THEN
      			    XI(1)=D2/A2
      			  ELSE
      			    SLOPE=A2/A1
      			    XI(1)=(SLOPE*D2+D1)/(SLOPE*A2+A1)
      			  ENDIF
      			  IF(XI(1).LE.1.0D0.AND.XI(1).GE.0.0D0) THEN
      			    X1=XE(ns1,NJ1)*(1.0D0-XI(1))+XE(ns2,NJ1)*
     '                        XI(1)
      			    X2=XE(ns1,NJ2)*(1.0D0-XI(1))+XE(ns2,NJ2)*
     '                        XI(1)
      			    DIST=(X1-XD(NJ1))**2+(X2-XD(NJ2))**2
      			    IF(LD(nd).EQ.0.OR.DIST.LT.EDD(nd)) THEN
      			      XID(1,nd)=XI(1)
      			      LD(nd)=ne
      			      EDD(nd)=DIST
      			    ENDIF
      			  ENDIF
      			  IF(DOP) THEN
      			    WRITE(OP_STRING,
     '                        '('' ne='',I4,'' nd='',I4,'
     '			      //''' Xi(1)='','//'E11.3,'' EDD(nd)='','
     '			      //'E11.3)') ne,nd,XI(1),EDD(nd)
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  ENDIF

      			ELSE IF(NJT.EQ.2.AND.ITYP10(1).EQ.2) THEN !2D cyl.polar
      			  XI(1)=(XD(2)-XE(ns1,2))/
     '				(XE(ns2,2)-XE(ns1,2))
      			  IF(XI(1).LE.1.0D0.AND.XI(1).GE.0.0D0) THEN
      			    X1=XE(ns1,1)*(1.0D0-XI(1))+XE(ns2,1)*XI(1)
      			    DIST=(X1-XD(1))**2
      			    IF(LD(nd).EQ.0.OR.DIST.LT.EDD(nd)) THEN
      			      XID(1,nd)=XI(1)
      			      LD(nd)=ne
      			      EDD(nd)=DIST
      			    ENDIF
      			  ENDIF
      			  IF(DOP) THEN
      			    WRITE(OP_STRING,
     '			      '('' ne='',I4,'' nd='',I4,'' Xi(1)='','
     '			      //'E11.3,'' EDD(nd)='',E11.3)')
     '                        ne,nd,XI(1),EDD(nd)
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  ENDIF
      			ELSE IF(NJT.EQ.3.OR.ITYP10(1).NE.1) THEN
      			ENDIF
      		      ENDIF

      		    ELSE IF(NITB.GE.2) THEN     !2D elements
      		      IF(ORTHOG) THEN           !use nonlinear Xi calculation
      			XI(1)=0.5D0 ! start at the centre of the element
      			XI(2)=0.5D0
      			ITMAX=10
      			IF(ITYP10(1).EQ.1) THEN !rectanglar cartesian
      			  CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
     '			    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
      			ELSE IF(ITYP10(1).EQ.2) THEN !cylindrical polar
      			ELSE IF(ITYP10(1).EQ.3) THEN !spherical polar
      			ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
      			  CALL CLOS24(IBT,IDO,INP,IT,ITMAX,NBJ(1,ne),
     '			    SQND,XE,XI,ZD(1,nd),ERROR,*9999)
      			ELSE IF(ITYP10(1).EQ.5) THEN !oblate spheroidal
      			ENDIF
      			IF(XI(1).GT.0.0D0.AND.XI(1).LE.1.0D0.AND
     '			  .XI(2).GT.0.0D0.AND.XI(2).LE.1.0D0) THEN
      			  IF(DOP) THEN
      			    WRITE(OP_STRING,
     '			      '('' Xi: '',2F7.3,'' in element'',I4,'
     '                        //''' SQ='',E12.3)') XI(1),XI(2),ne,SQND
      			    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  ENDIF
      			  IF(LD(nd).EQ.0.OR.SQND.LT.SQ(nd)) THEN
      			    XID(1,nd)=XI(1)
      			    XID(2,nd)=XI(2)
      			    LD(nd)=ne
      			    SQ(nd)=SQND
      			  ENDIF
      			ENDIF
      			IF(IT.GE.ITMAX) THEN
      			  WRITE(OP_STRING,'('' WARNING:'
     '			    //' Convergence not reached in CLOS2X''/'
     '			    //''' nd='',I4,'' ZD='',3(E12.6,1X)/'
     '			    //''' LD='',I4,'' XID'',2(E12.6,1X),'
     '                      //'''SQ='',E12.6/)')nd,(ZD(nj,nd),nj=1,NJT),
     '			    LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
      			  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      			ENDIF

      		      ELSE IF(CUBIC) THEN !use 1D nonlinear calculation
!old MPN 1-Nov-94: FUNCT1 is now in FE_ARCHIVE
c      			nb=NBJ(NJ1,ne) !nb for first geometric coord
c      			NKTB=NKT(0,nb)
c      			IF(LD(nd).EQ.0) THEN
c      			  CALL ZX(ITYP10(1),ZD(1,nd),XD)
cC ... 			    Find XID([1,2],nd) for a point inside cubic hermite-
cC     			    linear element
c      			  XPXIF=ZD(1,nd)    !set x and y of data point
c      			  YPXIF=ZD(2,nd)
c      			  DO I=1,8
c      			    XNODVAL(I,1)=XE(I,1)
c      			    XNODVAL(I,2)=XE(I,2)
c      			  ENDDO
c      			  KTYP26=4   !inputs to e04jbf
c      			  KTYP27=2
c      			  NPXIF=1
c      			  IPRINT=-1
c      			  LOCSCH=.FALSE.
c      			  INTYPE=0
c      			  MAXCAL=240
c      			  ETA=0.5D0
c      			  XTOL=1.0D-5
c      			  STEPMX=2.0D0
c      			  FEST=0.0D0
c      			  DELXIF(1)=1.0D-10
c      			  IBOUND=0
c      			  BLXIF(1)=0.0000001D0
c      			  BUXIF(1)=0.9999999D0
c      			  XIF(1)=0.5D0
c      			  LH=5   !can be 1
c      			  LIW=5
c      			  LW=18
c      			  IFAIL=1
c      			  CALL E04JBF(NPXIF,FUNCT1,MONIT,IPRINT,LOCSCH,
c     '			    INTYPE,E04JBQ,
c     '			    MAXCAL,ETA,XTOL,STEPMX,FEST,DELXIF,IBOUND,
c     '			    BLXIF,BUXIF,XIF,HESL,LH,HESD,ISTATE,FXIF,
c     '                      GXIF,IWXIF,LIW,WW,LW,IFAIL)
c      			  XI(1)=XIF(1)
c
cC     			  Find X_ENDO,X_EPI,Y_ENDO,Y_EPI using converged
cC     			  value of XI1. Note most recent values of these
cC     			  variables computed by NAG do not necessarily
cC     			  correlate with converged xi1
c      			  XI1=XI(1)
c      			  S01=1.0D0-3.0D0*XI1*XI1+2.0D0*XI1*XI1*XI1  !evaluate hermite basis function
c      			  S02=3.0D0*XI1*XI1-2.0D0*XI1*XI1*XI1
c      			  S11=XI1-2.0D0*XI1*XI1+XI1*XI1*XI1
c      			  S12=-(XI1*XI1-XI1*XI1*XI1)
c      			  F01=XNODVAL(1,1)                  !Nodal values
c      			  F02=XNODVAL(3,1)
c      			  F11=XNODVAL(2,1)
c      			  F12=XNODVAL(4,1)
c      			  X_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
c      			  F01=XNODVAL(1,2)
c      			  F02=XNODVAL(3,2)
c      			  F11=XNODVAL(2,2)
c      			  F12=XNODVAL(4,2)
c      			  Y_ENDO=S01*F01+S02*F02+S11*F11+S12*F12
c      			  F01=XNODVAL(5,1)
c      			  F02=XNODVAL(7,1)
c      			  F11=XNODVAL(6,1)
c      			  F12=XNODVAL(8,1)
c      			  X_EPI=S01*F01+S02*F02+S11*F11+S12*F12
c      			  F01=XNODVAL(5,2)
c      			  F02=XNODVAL(7,2)
c      			  F11=XNODVAL(6,2)
c      			  F12=XNODVAL(8,2)
c      			  Y_EPI=S01*F01+S02*F02+S11*F11+S12*F12
c
cC ... 			  Find value of objective function at converged XI(1)
c      			  XDIFFERENCE=DABS(X_EPI-X_ENDO)
c      			  IF(XDIFFERENCE.LT.1.0D-8) THEN     !singular case
c      			    XNPB=X_EPI
c      			    YNPB=YPXIF
c      			  ELSE                               !Nonsingular case
c      			    S1=X_ENDO-X_EPI
c      			    S2=Y_ENDO-Y_EPI
c      			    XNPB=(XPXIF*S1*S1+X_ENDO*S2*S2-S1*S2
c     '			      *(Y_ENDO-YPXIF))/(S1*S1+S2*S2)
c      			    YNPB=Y_ENDO+S2*(XNPB-X_ENDO)/S1
c      			  ENDIF
c      			  FC=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2
c
cC ... 			  flag points which are outside element
c      			  IF(FC.GT.1.0D-14) XI(1)=99.0D0
c      			  XDIFFERENCE=DABS(X_ENDO-X_EPI)
c      			  IF(XDIFFERENCE.LT.1.0D-7) THEN
c      			    XI(2)=(YPXIF-Y_ENDO)/(Y_EPI-Y_ENDO)
c      			  ELSE
c      			    XI(2)=(XPXIF-X_ENDO)/(X_EPI-X_ENDO)
c      			  ENDIF
c
c      			  IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0.AND.
c     '			     XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
c      			    IF(NITB.EQ.2) THEN !2D basis
c      			      IF(DOP) THEN
c      				WRITE(OP_STRING,
c     '                            '('' nd='',I4,'' Xi:'','
c     '				  //'3E11.3,'' (linear calc.)'')') nd,
c     '				  (XI(ni),ni=1,3)
c      				CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c      			      ELSE
c      			      ENDIF
c      			      nb=NBJ(NJ1,ne)
c      			      XID(1,nd)=XI(1)
c      			      XID(2,nd)=XI(2)
c      			      LD(nd)=ne !element no for data point nd
c      			    ELSE IF(NITB.EQ.3) THEN !3D basis
c      			    ENDIF
c      			  ENDIF
c      			ENDIF !LD

      		      ELSE IF(QUADRATIC) THEN !use 2D nonlinear calculation
!old MPN 1-Nov-94: FUNCT1 is now in FE_ARCHIVE
c      			nb=NBJ(NJ1,ne) !nb for first geometric coord
c      			NKTB=NKT(0,nb)
c      			IF(LD(nd).EQ.0) THEN
c      			  CALL ZX(ITYP10(1),ZD(1,nd),XD)
cC ... 			  Find XID([1,2],nd) for a point inside Lagrange-
cC     			    quadratic element
c      			  XPXIF=ZD(1,nd)            !set x and y of data point
c      			  YPXIF=ZD(2,nd)
c
c      			  DO I=1,9                  !Element Nodes
c      			    XNODVAL(I,1)=XE(I,1)
c      			    XNODVAL(I,2)=XE(I,2)
c      			  ENDDO
c
c      			  KTYP26=4   !inputs to e04jbf
c      			  KTYP27=3
c      			  NPXIF=2
c      			  IPRINT=-1
c      			  LOCSCH=.FALSE.
c      			  INTYPE=0
c      			  MAXCAL=560
c      			  ETA=0.5D0
c      			  XTOL=1.0D-5
c      			  STEPMX=2.0D0
c      			  FEST=0.0D0
c      			  DELXIF(1)=1.0D-10
c      			  DELXIF(2)=1.0D-10
c      			  IBOUND=0
c      			  BLXIF(1)=0.00001D0
c      			  BUXIF(1)=0.99999D0
c      			  BLXIF(2)=0.00001D0
c      			  BUXIF(2)=0.99999D0
c      			  XIF(1)=0.5D0
c      			  XIF(2)=0.5D0
c      			  LH=5   !can be 1
c      			  LIW=5
c      			  LW=18
c      			  IFAIL=1
c
c      			  CALL E04JBF(NPXIF,FUNCT1,MONIT,IPRINT,LOCSCH,
c     '			    INTYPE,E04JBQ,
c     '			    MAXCAL,ETA,XTOL,STEPMX,FEST,DELXIF,IBOUND,
c     '			    BLXIF,BUXIF,XIF,HESL,LH,HESD,ISTATE,FXIF,
c     '                      GXIF,IWXIF,LIW,WW,LW,IFAIL)
c      			  XI(1)=XIF(1)          !XIF contain converged xi values
c      			  XI(2)=XIF(2)
c
cC ... 			  Recompute FC (which requires recomputing XNPB,YNPB)
c      			  VL20=2.0D0*(XI(1)-0.5D0)*(XI(1)-1.0D0)
c      			  VL21=-4.0D0*(XI(1)-1.0D0)*XI(1)
c      			  VL22=2.0D0*XI(1)*(XI(1)-0.5D0)
c      			  CL20=2.0D0*(XI(2)-0.5D0)*(XI(2)-1.0D0)
c      			  CL21=-4.0D0*(XI(2)-1.0D0)*XI(2)
c      			  CL22=2.0D0*XI(2)*(XI(2)-0.5D0)
c
c      			  XNPB=VL20*CL20*XNODVAL(1,1)+VL21*CL20
c     '                                                *XNODVAL(2,1)+
c     '			       VL22*CL20*XNODVAL(3,1)+VL20*CL21
c     '						      *XNODVAL(4,1)+
c     '			       VL21*CL21*XNODVAL(5,1)+VL22*CL21
c     '						      *XNODVAL(6,1)+
c     '			       VL20*CL22*XNODVAL(7,1)+VL21*CL22
c     '						      *XNODVAL(8,1)+
c     '			       VL22*CL22*XNODVAL(9,1)
c      			  YNPB=VL20*CL20*XNODVAL(1,2)+VL21*CL20
c     '						      *XNODVAL(2,2)+
c     '			       VL22*CL20*XNODVAL(3,2)+VL20*CL21
c     '						      *XNODVAL(4,2)+
c     '			       VL21*CL21*XNODVAL(5,2)+VL22*CL21
c     '						      *XNODVAL(6,2)+
c     '			       VL20*CL22*XNODVAL(7,2)+VL21*CL22
c     '						      *XNODVAL(8,2)+
c     '			       VL22*CL22*XNODVAL(9,2)
c
c      			  FC=(XNPB-XPXIF)**2+(YNPB-YPXIF)**2  !value of obj. fn.
c
cC ... 			  flag points which are outside element
c      			   IF(FC.GT.1.0D-10) XI(1)=99.0D0
c      			   IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0.AND.
c     '			      XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
c      			    IF(NITB.EQ.2) THEN !2D basis
c      			      IF(DOP) THEN
c      				WRITE(OP_STRING,
c     '                            '('' nd='',I4,'' Xi:'','
c     '				  //'3E11.3,'' (linear calc.)'')') nd,
c     '				  (XI(ni),ni=1,3)
c      				CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c      			      ELSE
c      			      ENDIF
c      			      nb=NBJ(NJ1,ne)
c      			      XID(1,nd)=XI(1)
c      			      XID(2,nd)=XI(2)
c      			      LD(nd)=ne !element no for data point nd
c      			    ELSE IF(NITB.EQ.3) THEN !3D basis
c      			    ENDIF
c      			  ENDIF
c      			ENDIF !LD

      		      ELSE IF(.NOT.ORTHOG) THEN !use linear Xi calculation
      			nb=NBJ(NJ1,ne) !nb for first geometric coord
      			NKTB=NKT(0,nb)
      			ns1=1
      			ns2=1+NKTB
      			ns3=1+2*NKTB
      			ns4=1+3*NKTB
      			A1=XE(ns2,NJ1)-XE(ns1,NJ1)
      			B1=XE(ns3,NJ1)-XE(ns1,NJ1)
      			C1=XE(ns1,NJ1)-XE(ns2,NJ1)-XE(ns3,NJ1)+
     '                    XE(ns4,NJ1)
      			nb=NBJ(NJ2,ne)
      			NKTB=NKT(0,nb)
      			ns1=1
      			ns2=1+NKTB
      			ns3=1+2*NKTB
      			ns4=1+3*NKTB
      			A2=XE(ns2,NJ2)-XE(ns1,NJ2)
      			B2=XE(ns3,NJ2)-XE(ns1,NJ2)
      			C2=XE(ns1,NJ2)-XE(ns2,NJ2)-XE(ns3,NJ2)+
     '                    XE(ns4,NJ2)
      			IF(DOP) THEN
      			  WRITE(OP_STRING,'(/'' Element '',I3)') ne
      			  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  WRITE(OP_STRING,'('' A1,B1,C1='',3E11.3)')
     '			     A1,B1,C1
      			  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			  WRITE(OP_STRING,'('' A2,B2,C2='',3E11.3)')
     '			     A2,B2,C2
      			  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			ENDIF
      			ALFA=A2*C1-A1*C2
      			IF(LD(nd).EQ.0) THEN
      			  CALL ZX(ITYP10(1),ZD(1,nd),XD)
      			  EXCLUDE=.TRUE.
      			  IF(ITYP10(1).EQ.4) THEN !check if lies within element
      			    THETAMIN=DMIN1(XE(ns2,3),XE(ns4,3))
      			    THETAMAX=DMAX1(XE(ns1,3),XE(ns3,3))
      			    IF(XD(3).LT.THETAMAX.AND.
     '			       XD(3).GT.THETAMIN) THEN
      			      EXCLUDE=.FALSE.
      			      D2=XD(NJ2)-XE(1,NJ2)
      			    ELSE IF((XD(3)+2.0D0*PI).LT.THETAMAX
     '			       .AND.(XD(3)+2.0D0*PI).GT.THETAMIN) THEN
      			      EXCLUDE=.FALSE.
      			      D2=XD(NJ2)+2.0D0*PI-XE(1,NJ2)
      			    ELSE IF((XD(3)-2.0D0*PI).LT.THETAMAX
     '			       .AND.(XD(3)-2.0D0*PI).GT.THETAMIN) THEN
      			      EXCLUDE=.FALSE.
      			      D2=XD(NJ2)-2.0D0*PI-XE(1,NJ2)
      			    ENDIF
      			  ELSE
      			    EXCLUDE=.FALSE.
      			    D2=XD(NJ2)-XE(1,NJ2)
      			  ENDIF
      			  IF(.NOT.EXCLUDE) THEN
      			    D1=XD(NJ1)-XE(1,NJ1)
      			    BETA=C2*D1-C1*D2+A2*B1-A1*B2
      			    GAMA=B2*D1-B1*D2
      			    DELTA=BETA*BETA-4.0D0*ALFA*GAMA
      			    IF(DOP) THEN
      			      WRITE(OP_STRING,'('' D1='',E11.3,'
     ,                          //''' D2='',E11.3,'
     '				//''' ALFA='',E11.3,'' BETA='',E11.3,'
     '				//''' GAMA='',E11.3,'' DELTA='',E11.3)')
     '				D1,D2,ALFA,BETA,GAMA,DELTA
      			      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      			    ENDIF
      			    IF(DABS(ALFA).GT.1.0D-6) THEN
      			      IF(DELTA.GE.0.0D0) THEN
      				XI(1)=(-BETA+DSQRT(DELTA))/(2.0D0*ALFA)
      				IF(XI(1).LT.0.0D0.OR.
     '                            XI(1).GT.1.0D0) THEN
      				  XI(1)=(-BETA-DSQRT(DELTA))/
     '                              (2.0D0*ALFA)
      				ENDIF
      			      ENDIF
      			    ELSE IF(DABS(ALFA).LE.1.0D-6) THEN
      			      IF(DABS(BETA).GT.1.0D-6) THEN
      				XI(1)=-GAMA/BETA
      			      ELSE
      				XI(1)=-1.0D0
      			      ENDIF
      			    ENDIF
      			    IF(XI(1).GE.0.0D0.AND.XI(1).LE.1.0D0) THEN
      			      DENOM1=B1+C1*XI(1)
      			      IF(DABS(DENOM1).GT.1.0D-6) THEN
      				XI(2)=(D1-A1*XI(1))/DENOM1
      			      ELSE
      				DENOM2=B2+C2*XI(1)
      				IF(DABS(DENOM2).GT.1.0D-6) THEN
      				  XI(2)=(D2-A2*XI(1))/DENOM2
      				ELSE
      				  WRITE(OP_STRING,
     '				    '('' Xi(2) cannot be defined for '
     '				    //'nd='',I6,'' in element '',I5)')
     '                              nd,ne
      				  CALL WRITES(IOFI,OP_STRING,ERROR,
     '                              *9999)
      				ENDIF
      			      ENDIF
      			      IF(XI(2).GE.0.0D0.AND.XI(2).LE.1.0D0) THEN
      				IF(NITB.EQ.2) THEN
      				  IF(DOP) THEN
      				    WRITE(OP_STRING,
     '				      '('' nd='',I4,'' Xi:'',3E11.3,'
     '				      //''' (linear calc.)'')')
     '                                nd,(XI(ni),ni=1,3)
      				    CALL WRITES(IODI,OP_STRING,ERROR,
     '				      *9999)
      				  ENDIF
      				  nb=NBJ(NJ1,ne)
      				  IF(NKT(0,nb).GT.1) THEN !update Xi(2)
      				    XI(2)=0.0D0
      				    X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '				      INP(1,1,nb),nb,1,XI,XE(1,NJ1))
      				    XI(2)=1.0D0
      				    X2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '				      INP(1,1,nb),nb,1,XI,XE(1,NJ1))
      				    IF(DABS(X2-X1).GT.1.0D-6) THEN !PJH 13APR91
      				      XI(2)=(XD(NJ1)-X1)/(X2-X1)
      				    ELSE
      				      WRITE(OP_STRING,
     '                                  '('' Warning!!!!'',''X1=X2'')')
      				      CALL WRITES(IOOP,OP_STRING,ERROR,
     '					*9999)
      				      XI(2)=1.0D0
      				    ENDIF
      				  ENDIF
      				  XID(1,nd)=XI(1)
      				  XID(2,nd)=XI(2)
      				  LD(nd)=ne
      				ELSE IF(NITB.EQ.3) THEN
      				  nb=NBJ(1,ne)
      				  XI(3)=0.0D0
      				  X0=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '				    INP(1,1,nb),nb,1,XI,XE(1,1))
      				  XI(3)=1.0D0
      				  X1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '				    INP(1,1,nb),nb,1,XI,XE(1,1))
      				  DIFF=X1-X0
      				  IF(DABS(DIFF).GT.1.0D-6) THEN
      				    XI(3)=(XD(1)-X0)/DIFF
      				  ELSE
      				    XI(3)=0.0D0
      				  ENDIF
      				  IF(DOP) THEN
      				    WRITE(OP_STRING,
     '                                '('' X0='',E11.3,'
     '				      //''' X1='',E11.3,'' XD(1)='','
     '                                //'E11.3,'' Xi(3)='',E11.3)')
     '                                X0,X1,XD(1),XI(3)
      				    CALL WRITES(IODI,OP_STRING,ERROR,
     '				      *9999)
      				  ENDIF
      				  IF(XI(3).GT.0.0D0.AND.
     '                              XI(3).LE.1.0D0) THEN
      				    LD(nd)=ne
      				  ELSE IF(EXTRAPOLATE.AND.
     '			            XI(3).GT.1.0D0.AND.
     '                              XI(3).LE.1.1D0) THEN
      	                            WRITE(OP_STRING,
     '	      '('' TAGGING: Found nd='',I6,'' with xi(3)='',F5.2)') 
     '                              nd,XI(3)
      	                            CALL WRITES(IOOP,OP_STRING,
     '                              ERROR,*9999)
      				    IF(LD(nd).EQ.0) THEN
      				      LD(nd)=-ne !to tag for 'extrapolate' below
      				    ENDIF
      				  ENDIF
      				  XID(1,nd)=XI(1)
      				  XID(2,nd)=XI(2)
      				  XID(3,nd)=XI(3)
!news                             AAY use set xi_3 option if required 23 March 95
                                  IF(SET_XI_3) THEN
                                    XI(3)=XI_3
                                  ENDIF
!newe
      				ENDIF
      				IF(DOP) THEN
      				  WRITE(OP_STRING,
     '				    '('' nd='',I4,'' Xi:'',3E11.3)')
     '				    nd,(XI(ni),ni=1,3)
      				  CALL WRITES(IODI,OP_STRING,ERROR,
     '                              *9999)
      				ENDIF
      			      ENDIF !0<xi(2)<1
      			    ENDIF !0<xi(1)<1
      			  ENDIF !.NOT.EXCLUDE
      			ENDIF !LD(nd)=0
      		      ENDIF !linear/orthog/cubic
      		    ENDIF !nitb=1,2,3
      		  ENDDO !nolist
                ENDIF ! not found
              ENDDO ! nd=1,NDT loop
              IF(EXTRAPOLATE) THEN
                DO nd=1,NDT
                  IF(LD(nd).LT.0) THEN !data point currently not in any element
                    LD(nd)=-LD(nd)
                    XID(3,nd)=1.0D0
                  ENDIF
                ENDDO
              ENDIF
            ENDIF !orthog/centroid etc

            CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              CALL ASSERT(NDLT(ne).LE.NDEM,
     '          '>>NDEM too small',ERROR,*9999)
            ENDDO
            CALC_XI=.TRUE.

          ELSE IF(TYPE(1:5).EQ.'BEADS') THEN
            nb=NBJ(1,1)
            NITB=NIT(nb)
            NKTB=NKT(0,nb)
            IF(NET(1).GT.0.AND.NITB.EQ.3.AND.NJT.EQ.3.AND.NDT.GT.6) THEN
C ***         Calculate xi-coordinates of projections of bead positions
C ***         on to xi_3 lines of 3D elements using normalized
C ***         lengths in each column. It is assumed that data points
C ***         1,2,3 are epicardial beads numbered clockwise, that data
C ***         points 4,5,6 are bottom beads of cols 1,2,3, and that
C ***         elements numbers increase from epicardium to endocardium.
              DO icol=1,3
                DO ni=1,NITB
                  CVEC(ni,icol)=ZD(ni,icol+3)-ZD(ni,icol)
                ENDDO
                CMAG(icol)=0.0D0
                DO ni=1,NITB
                  CMAG(icol)=CMAG(icol)+CVEC(ni,icol)**2
                ENDDO
                CMAG(icol)=DSQRT(CMAG(icol))
              ENDDO
              XID(1,1)=0.0D0
              XID(2,1)=0.0D0
              XID(3,1)=0.0D0
              XID(1,2)=1.0D0
              XID(2,2)=0.0D0
              XID(3,2)=0.0D0
              XID(1,3)=0.5D0
              XID(2,3)=1.0D0
              XID(3,3)=0.0D0
              XID(1,4)=0.0D0
              XID(2,4)=0.0D0
              XID(3,4)=NET(1)
              XID(1,5)=1.0D0
              XID(2,5)=0.0D0
              XID(3,5)=NET(1)
              XID(1,6)=0.5D0
              XID(2,6)=1.0D0
              XID(3,6)=NET(1)
              DO nd=7,NDT
                DO icol=1,3
                  BMAG(icol)=0.0D0
                  BDOTC=0.0D0
                  DO ni=1,3
                    BMAG(icol)=BMAG(icol)+(ZD(ni,nd)-ZD(ni,icol))**2
                    BDOTC=BDOTC+(ZD(ni,nd)-ZD(ni,icol))*CVEC(ni,icol)
                  ENDDO
                  BMAG(icol)=DSQRT(BMAG(icol))
                  AB(icol)=DABS(BDOTC/(BMAG(icol)*CMAG(icol)))
                ENDDO
                IF(AB(1).GT.AB(2).AND.AB(1).GT.AB(3)) THEN
                  XID(1,nd)=0.0D0
                  XID(2,nd)=0.0D0
                  XID(3,nd)=BMAG(1)/CMAG(1)*NET(1)
                ELSE IF(AB(2).GT.AB(1).AND.AB(2).GT.AB(3)) THEN
                  XID(1,nd)=1.0D0
                  XID(2,nd)=0.0D0
                  XID(3,nd)=BMAG(2)/CMAG(2)*NET(1)
                ELSE IF(AB(3).GT.AB(1).AND.AB(3).GT.AB(2)) THEN
                  XID(1,nd)=0.5D0
                  XID(2,nd)=1.0D0
                  XID(3,nd)=BMAG(3)/CMAG(3)*NET(1)
                ENDIF
              ENDDO
              DO nd=1,NDT
                i=INT(XID(3,nd))
                IF(i.GE.0.AND.i.LE.(NET(1)-1)) THEN
                  LD(nd)=i+1
                  XID(3,nd)=XID(3,nd)-i
                ELSE IF(i.LT.0) THEN
                  LD(nd)=1
                  XID(3,nd)=0.0D0
                ELSE IF(i.GT.(NET(1)-1)) THEN
                  LD(nd)=i
                  XID(3,nd)=1.0D0
                ENDIF
                CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' Bead point '',I3,'' lies in '
     '              //'element '',I3,'' at Xi coords: '',3E11.3)')
     '              nd,LD(nd),(XID(ni,nd),ni=1,NITB)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ELSE IF(NET(1).EQ.0) THEN
              WRITE(OP_STRING,*) ' >>no elements defined'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(njt.NE.3.or.nitb.NE.3) THEN
              WRITE(OP_STRING,*) ' >>Bead option requires 3D elem.s'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(NDT.LT.6) THEN
              WRITE(OP_STRING,*) ' >>At least six beads are needed'
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

          ENDIF !xipos/beads

        ENDIF !filio/calc
      ENDIF

      CALL_DATA=.TRUE.
      CALL EXITS('DEXI')
      RETURN
 9999 CALL ERRORS('DEXI',ERROR)
      CALL EXITS('DEXI')
      RETURN 1
      END


C LC 25/2/97 removed section : 
C#### Subroutine: DEXI
C###  Description:
C###    DEXI defines data xi points.

! PJH 8Dec95  DO nd=STRIPE0,STRIPE1
!                LD(nd)=0
!                NT=NLS_NDATA_CONT(nd) !the tag identifier
!                nimag =NLS_NDATA_IMAG(nd) !the image identifier
!                !intersect the line segments with the image plane
!                CALL CONT_PLANE_X(NBJ,NJE,NKE,NPNE,NPF,NQE,SE,XA,XE,ZP,
!     '            NT,IMAGE_NORMAL(1,nimag),ZD(1,nd),LD(nd),XID(1,nd),
!     '            SQMAX,ERROR,*9999)
!              ENDDO !nd
!            ELSE IF(CONTOUR) THEN
C            ELSE !PJH temp 8Dec95
C              NC=NTAG(0)+1 !use next tag data struct to store line segs
C              CALL ASSERT(NC.LE.100,
C     '          '>>Error: increase number of tags in /CONT01/',
C     '          ERROR,*9999)
CC             for all image planes do
C              DO nimag=1,NIMAGE(0)
CC               for all surfaces do
C                DO nsurf=1,NTSURF
CC!!! INDEX_SURF used before it is set
C                  ns=INDEX_SURF(nsurf)
CC                 intersect this surface with the image plane
Cc                 CALL SURF_PLANE_X(NBJ,NJE,NKE,NPNE,NPF,NQE,SE,
Cc    '              XA,XE,ZP,ns,NC,IMAGE_NORMAL(1,nimag),
Cc    '              IMAGE_POS(1,nimag),ERROR,*9999)
CC!!! ICON_X used before it is set
CC                 precalc the basis functions for the line segments
C                  DO np=1,ICON_X(NC)
C                    DO nj=1,NJT
CC!!! ISURF_ELEM used before it is set
C                      nb=NBJ(nj,ISURF_ELEM(1,ns))
C                      DO n=1,2
C                        DO nn=1,NNT(nb)
C                          DO nk=1,NKT(0,nb)
C                            NLS_CON_PSI(nk+(nn-1)*NKT(0,nb),n,np,nj,nc)
C     '                        =PSI1(IBT(1,1,nb),IDO(1,1,0,nb),
C     '                        INP(1,1,nb),nb,1,nk,nn,
C     '                        NLS_CON_XI(1,n,np,nc))
C                          ENDDO !nk
C                        ENDDO !nn
C                      ENDDO  !n
C                    ENDDO !nj
C                  ENDDO !np
C                  !for all data points on this surface and image do
C! PJH 8Dec95      DO nd=CONTOUR0,CONTOUR1
C!                    IF(NLS_NDATA_CONT(nd).EQ.ns
C!     '                .AND.NLS_NDATA_IMAG(nd).EQ.nimag) THEN
C!                      LD(nd)=0
C!!                     construct plane normal to image and contour
C!                      CALL CROSS(WD(1,nd),IMAGE_NORMAL(1,nimag),Z1)
C!                      !find the intersection with this plane 
C!                      CALL CONT_PLANE_X(NBJ,NJE,NKE,NPNE,NPF,NQE,SE,
C!     '                  XA,XE,ZP,NC,Z1,ZD(1,nd),LD(nd),XID(1,nd),
C!     '                  SQMAX,ERROR,*9999)
C!                    ENDIF
C!                  ENDDO !nd
C                ENDDO !nsurf
C              ENDDO !nimag
C!newe


C KAT 2001-04-09
      SUBROUTINE DEPLANE(STRING,ERROR,*)
                                                       
C#### Subroutine: DEPLANE
C###  Description:
C###    DEPLANE defines image and tag planes.
    
      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:data00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERR,ibeg,IBEG1,IBEG2,iend,IEND1,IEND2,
     '  IOSTAT,IPFILE,N3CO,nd,nj
      CHARACTER FILE*80,STATUS*3,TYPE*20
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,GENER,MOUSE

      CALL ENTERS('DEPLANE',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN                           
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)                   
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define planes;r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines image and tag planes.
C###  Parameter:      <(stripes/images)[stripes]>

        OP_STRING(1)=STRING(1:IEND)//';r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<stripes/images[stripes]>' 
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('FE24','DOC','DEIMAG',ERROR,*9999)
      ELSE  
        CALL PARSE_QUALIFIERS('RW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(CBBREV(CO,'STRIPES',3,NOCO+1,NTCO,N3CO)) THEN
          TYPE='STRIPE'
        ELSE IF(CBBREV(CO,'IMAGES',3,NOCO+1,NTCO,N3CO)) THEN
          TYPE='IMAGE'
        ENDIF

        IF(FILIO) THEN
          ALL_REGIONS=.FALSE.
C LKC 24-APR-1998 Init IPFILE
          IPFILE=1

          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'plane',STATUS,
     '      ERR,ERROR,*9999)
          IF(IOTYPE.EQ.2) THEN
            IF(TYPE(1:6).EQ.'STRIPE') THEN
              !read in the tag planes
              NTAG(0)=0
              DO nd=1,NDM
                READ(UNIT=IFILE,FMT=*,IOSTAT=IOSTAT,END=9998)
     '            NTAG(nd),(TAG_POS(nj,nd),nj=1,3),
     '            (TAG_NORMAL(nj,nd),nj=1,3)
                NTAG(0)=NTAG(0)+1
              ENDDO !nd
            ELSE IF(TYPE(1:5).EQ.'IMAGE') THEN
              !read in the contour planes
              NIMAGE(0)=0
              DO nd=1,NDM
                READ(UNIT=IFILE,FMT=*,IOSTAT=IOSTAT,END=9998)
     '            NIMAGE(ND),(IMAGE_POS(nj,nd),nj=1,3),
     '            (IMAGE_NORMAL(nj,nd),nj=1,3)
                NIMAGE(0)=NIMAGE(0)+1
              ENDDO !nd
            ENDIF

          ELSE IF(IOTYPE.EQ.3) THEN
          ENDIF
!newe
          CALL CLOSEF(IFILE,ERROR,*9999)

        ENDIF

      ENDIF      
                     
 9998 CALL EXITS('DEPLANE')
      RETURN   
 9999 CALL ERRORS('DEPLANE',ERROR)
      CALL EXITS('DEPLANE')      
      RETURN 1       
      END  


C LKC 9-NOV-1998 Unused define potential
      SUBROUTINE DEPOTE(NDP,WD,ZD,STRING,ERROR,*)

C#### Subroutine: DEPOTE
C###  Description:
C###    DEPOTE defines potentials.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:file01.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER NDP(NDM)
      REAL*8 WD(NJM,NDM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER DATASET,IBEG,IBEG1,IBEG2,ICHAR,IEND,IEND1,IEND2,IFROMC,
     '  INFO,N3CO,nd,nj,NJTT,NOQUES,SERIES
      CHARACTER CHAR1*1,DATA_FORMAT*6,FILE*100,
     '  TITLE*80,TYPE*6,STATUS*3
      LOGICAL ABBREV,CALCU,CBBREV,CONTINUE,FILIO,GENER,MOUSE

      CALL ENTERS('DEPOTE',*9999)
      ICHAR=999


 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)
        CALL TRIM(PATH00,IBEG2,IEND2)

C#### Command: FEM define potential;l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <data_format (bifile/ipfile/map3d/emap)[ipfile]>
C###  Parameter:      <series #[1]>
C###  Parameter:      <dataset #[1]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//';l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
         OP_STRING(2)=BLANK(1:15)//
     '    '<data_format (bifile/ipfile/map3d/emap)[ipfile]>'
        OP_STRING(3)=BLANK(1:15)//'<series #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<dataset #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEPOTE',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'DATA_FORMAT',5,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'IPFILE',2)) THEN
            DATA_FORMAT='IPFILE'
          ELSE IF(ABBREV(CO(N3CO+1),'BIFILE',2)) THEN
            DATA_FORMAT='BIFILE'
          ELSE IF(ABBREV(CO(N3CO+1),'MAP3D',2)) THEN
            DATA_FORMAT='MAP3D'
          ELSE IF(ABBREV(CO(N3CO+1),'EMAP',2)) THEN
            DATA_FORMAT='EMAP'
          ELSE
            DATA_FORMAT='IPFILE'
          ENDIF
        ELSE
          DATA_FORMAT='IPFILE'
        ENDIF
        IF(DATA_FORMAT(1:6).EQ.'IPFILE') THEN
          CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,CALCU,FILIO,
     '      GENER,MOUSE,STATUS,ERROR,*1)
        ELSE
          CALL PARSE_QUALIFIERS('RW',noco,1,CO,COQU,CALCU,FILIO,
     '      GENER,MOUSE,STATUS,ERROR,*1)
        ENDIF

        FILEIP=.FALSE.
        NOQUES=0          
        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C MLB 19 Feb 1997 Both SERIES and DATASET are not used after they are
C set and are in the OMIT file. IF this is changed please remove them.
        IF(CBBREV(CO,'SERIES',1,noco+1,NTCO,N3CO)) THEN
          IF(DATA_FORMAT(1:6).EQ.'IPFILE') THEN
            IF(NTCO.GE.(N3CO+1)) THEN
              SERIES=IFROMC(CO(N3CO+1))
            ELSE
              SERIES=1
            ENDIF
          ELSE
            IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
              SERIES=0
            ELSE IF(NTCO.GE.(N3CO+1)) THEN
              SERIES=IFROMC(CO(N3CO+1))
            ELSE
              SERIES=1
            ENDIF
          ENDIF
        ELSE
          SERIES=1
        ENDIF
        IF(CBBREV(CO,'DATASET',5,noco+1,NTCO,N3CO)) THEN
          IF(DATA_FORMAT(1:6).EQ.'IPFILE') THEN
            IF(NTCO.GE.(N3CO+1)) THEN
              DATASET=IFROMC(CO(N3CO+1))
            ELSE
              DATASET=1
            ENDIF
          ELSE
            IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
              DATASET=0
            ELSE IF(NTCO.GE.(N3CO+1)) THEN
              DATASET=IFROMC(CO(N3CO+1))
            ELSE
              DATASET=1
            ENDIF
          ENDIF
        ELSE
          DATASET=1
        ENDIF

        IF(FILIO) THEN
          IF(DATA_FORMAT(1:6).EQ.'IPFILE') THEN
            CALL TRIM(FILE,IBEG,IEND)
            IF(IOTYPE.GT.1) THEN
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.ipdata',
     '          STATUS,'SEQUEN','FORMATTED',132,ERROR,*9999)
            ENDIF
            IF(IOTYPE.EQ.1) THEN
              nd=0
              CONTINUE=.TRUE.
              DO WHILE(CONTINUE)
      		nd=nd+1
                WRITE(CHAR1,'(I1)') NJT
      		FORMAT='(/$,'' Enter '//CHAR1(1:1)
     '            //' coords [Exit]: '')'
      		RDEFLT(1)=-1.0D6
      		CALL GINOUT(IOTYPE,5,IOIP,0,0,0,NOQUES,FILEIP,FORMAT,
     '		  NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     '		  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '		  ERROR,*9999)
      		IF(RDATA(1).GT.-1.0D5) THEN
      		  NDP(nd)=nd
      		  IF(ITYP10(1).GE.2) THEN
      		    RDATA(2)=RDATA(2)*PI/180.0D0
      		  ENDIF
      		  IF(ITYP10(1).GE.3) THEN
      		    RDATA(3)=RDATA(3)*PI/180.0D0
      		  ENDIF
      		  CALL XZ(ITYP10(1),RDATA,ZD(1,nd))
      		  DO nj=1,NJT
      		    WD(nj,nd)=1.0D0
      		  ENDDO
      		  RDEFLT(1)=0.0D0
      		  FORMAT='(/$,'' Enter potential value [0.0]: '')'
      		  CALL GINOUT(IOTYPE,5,IVDU,0,0,0,NOQUES,FILEIP,
     '		    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '		    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '		    RMAX,INFO,ERROR,*9999)
      		  ZD(NJT+1,nd)=RDATA(1)
      		  RDEFLT(1)=1.0D0
      		  FORMAT=
     '              '(/$,'' Enter potential value weight [1.0]: '')'
      		  CALL GINOUT(IOTYPE,5,IVDU,0,0,0,NOQUES,FILEIP,
     '		    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '		    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '		    RMAX,INFO,ERROR,*9999)
      		  WD(NJT+1,nd)=RDATA(1)
      		ELSE
      		  CONTINUE=.FALSE.
      		  NDT=nd-1
      		ENDIF
      	      ENDDO
      	    ELSE IF(IOTYPE.EQ.2) THEN !read potential data file
      	      NJTT=NJT+1
              TYPE='FIELD'
              WRITE(OP_STRING,'('' >>1 Potential data being read'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	      CALL IODATA('READ','RADIANS',TYPE,IFILE,NDP,NJTT,
     '          TITLE,WD,ZD,ERROR,*9999)
      	      WRITE(OP_STRING,*) TITLE
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	    ELSE IF(IOTYPE.EQ.3) THEN !write data file
      	      NJTT=NJT+1
              TITLE=' Potential data file'
              TYPE='FIELD'
      	      WRITE(OP_STRING,'('' >>1 Potential data being written'')')
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      	      CALL IODATA('WRITE','RADIANS',TYPE,IFILE,NDP,NJTT,
     '          TITLE,WD,ZD,ERROR,*9999)
      	    ENDIF
            CALL CLOSEF(IFILE,ERROR,*9999)
          ELSE IF(DATA_FORMAT(1:6).EQ.'BIFILE') THEN
            CALL TRIM(FILE,IBEG,IEND)
            IF(IOTYPE.GT.1) THEN
              CALL OPENF(IFILE,'DISK',FILE(IBEG:IEND)//'.bipote',
     '          STATUS,'SEQUEN','UNFORMATTED',132,ERROR,*9999)
            ENDIF
            CALL CLOSEF(IFILE,ERROR,*9999)
          ELSE IF(DATA_FORMAT(1:5).EQ.'MAP3D') THEN
      	    WRITE(OP_STRING,'('' >>Not implemented yet'')')
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(DATA_FORMAT(1:4).EQ.'EMAP') THEN
      	    WRITE(OP_STRING,'('' >>Not implemented yet'')')
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('DEPOTE')
      RETURN
 9999 CALL ERRORS('DEPOTE',ERROR)
      CALL EXITS('DEPOTE')
      RETURN 1
      END

C18-nov-1999
      SUBROUTINE DETIME(IBT,IDO,INP,NAN,NGAP,XE,STRING,ERROR,*)

C#### Subroutine: DETIME
C###  Description:
C###    Reads in a time variable from a .iptime file, and adds it to the
C###    available time variables. Time variables can be used to specify
C###    boundary conditions on nodes, elements, and grid points.
C *** Created : 9 August 1999 - David Nickerson

      IMPLICIT NONE

      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NGAP(NIM,NBM)
      REAL*8 XE(NSM,NJM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IPFILE,N3CO,VARIABLE
! not referenced
C     INTEGER nr,NTSL,nx,nxc
      REAL*8 START_TIME,END_TIME,GET_TV_VALUE_AT_TIME,VALUE,time
      LOGICAL ABBREV,ALL_REGIONS,CALCU,CBBREV,FILIO,GENER,MOUSE,
     '  FIRST_TIME,BACKUP
      LOGICAL START_TIME_SPECIFIED,END_TIME_SPECIFIED,DUMP
      CHARACTER FILE*100,STATUS*3

      CALL ENTERS('DETIME',*9999)

      START_TIME_SPECIFIED=.FALSE.
      END_TIME_SPECIFIED=.FALSE.

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)
        CALL TRIM(PATH00,IBEG2,IEND2)

C--------------------------------------------------------------------

C#### Command: FEM define time_variable;l/r/p/w/d<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:      <start_time #>
C###    Specify a start time to overwrite the value read from FILENAME
C###  Parameter:      <end_time #>
C###    Specify an end time to overwrite the value read from FILENAME
C###  Parameter:      <variable all/#[all]>
C###    Specify a variable to write, defaults to all variables
C###  Description:
C###    Reads in a time variable(s) from FILENAME.iptime and adds it 
C###    (them) to the variables available to be used to specify 
C###    boundary conditions on nodes, elements, and grid points.

        OP_STRING(1)=STRING(1:IEND)//';l/r/p/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<start_time #>'   
        OP_STRING(3)=BLANK(1:15)//'<end_time #>'
        OP_STRING(4)=BLANK(1:15)//'<variable all/#[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C--------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DETIME',ERROR,*9999)
      ELSE

        CALL PARSE_QUALIFIERS('LRPW',noco,1,CO,COQU,CALCU,FILIO,GENER,
     '    MOUSE,STATUS,ERROR,*1)
        
        IF(CBBREV(CO,'START_TIME',3,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSRE(CO(N3CO+1),START_TIME,ERROR,*1)
          START_TIME_SPECIFIED=.TRUE.
        ELSE       
          START_TIME_SPECIFIED=.FALSE.
        ENDIF

        IF(CBBREV(CO,'END_TIME',3,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSRE(CO(N3CO+1),END_TIME,ERROR,*1)
          END_TIME_SPECIFIED=.TRUE.
        ELSE       
          END_TIME_SPECIFIED=.FALSE.
        ENDIF

        IF(CBBREV(CO,'VARIABLE',3,NOCO+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
            VARIABLE=-1
          ELSE
            CALL PARSIN(CO(N3CO+1),VARIABLE,ERROR,*1)
            CALL ASSERT(VARIABLE.GT.0,'Variable needs to be > 0',
     '        ERROR,*1)
          ENDIF
        ELSE       
          VARIABLE=-1
        ENDIF

        IF(CBBREV(CO,'DUMP',3,NOCO+1,NTCO,N3CO)) THEN
          DUMP=.TRUE.
        ELSE       
          DUMP=.FALSE.
        ENDIF

        IF(FILIO) CALL CHECKF(2,NOCO,NTCOQU,CO,COQU,FILE,STRING,*1)

        IF(FILIO.AND..NOT.DUMP) THEN
          IPFILE=2
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          ALL_REGIONS=.FALSE. !when parse_regions not called
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'time',
     '        STATUS,ERROR,*9999)
            CALL IPTIME(IBT,IDO,INP,NAN,NGAP,VARIABLE,END_TIME,
     '        START_TIME,END_TIME_SPECIFIED,START_TIME_SPECIFIED,
     '        ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
c        ENDIF !filio

        ELSEIF(DUMP) THEN
          
          DO time=0.0d0,100.0d0,10.0d0
            VALUE=GET_TV_VALUE_AT_TIME(IBT,IDO,INP,VARIABLE,time,XE)
            WRITE(OP_STRING,'(D12.5,D12.5)') time,VALUE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ENDIF

      ENDIF

      CALL EXITS('DETIME')
      RETURN
 9999 CALL ERRORS('DETIME',ERROR)
      CALL EXITS('DETIME')
      RETURN 1
      END


Module FE25
=========== 

      SUBROUTINE SGBASE(INDEX,IBT,IDO,INP,ISBASE,ISEG,nb,CSEG,ERROR,*)

C#### Subroutine: SGBASE
C###  Description:
C###    SGBASE creates Basis segment ISBASE
C     This is an old routine. IW is not passed to it so take iw=1. AAY.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISBASE,ISEG(*),nb
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,INDEX_OLD,iw,nk,nn,ns
      REAL*8 POINTS(3,21),PXI,XE_LOCAL(8),XI(3),ZE_LOCAL(8)

      CALL ENTERS('SGBASE',*9999)
C MHT 23-03-00 iw=1, see above comment
      iw=1 
      CALL OPEN_SEGMENT(ISBASE,ISEG,iw,'BASE',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      IF(NIT(nb).EQ.1) THEN
        POINTS(1,1)=0.0D0
        POINTS(2,1)=0.0D0
        POINTS(1,2)=1.2D0
        POINTS(2,2)=0.0D0
        CALL POLYLINE(INDEX,1,2,POINTS,ERROR,*9999)
        POINTS(1,1)=0.0D0
        POINTS(2,1)=-0.5D0
        POINTS(1,2)=0.0D0
        POINTS(2,2)=1.2D0
        DO nn=1,NNT(nb)
          XE_LOCAL(1+(nn-1)*NKT(0,nb))=DBLE(nn-1)/DBLE(NNT(nb)-1)
          DO nk=2,NKT(0,nb)
            XE_LOCAL(nk+(nn-1)*NKT(0,nb))=1.0D0
          ENDDO
        ENDDO
        WRITE(OP_STRING,'('' XE_LOCAL: '',4E11.3)') 
     '    (XE_LOCAL(ns),ns=1,NST(nb))
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        DO i=1,10
          XI(1)=DBLE(i)/10.0D0
          POINTS(1,I)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '      XI,XE_LOCAL)
          WRITE(OP_STRING,'('' xi='',e11.3,'' points(1,i)='',e11.3)')
     '      xi(1),points(1,i)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          POINTS(2,i)=0.0D0
        ENDDO
        CALL POLYMARKER(INDEX,1,10,POINTS,ERROR,*9999)
        DO ns=1,NST(nb)
          ZE_LOCAL(ns)=0.0D0
        ENDDO
        DO nn=1,NNT(nb)
          DO nk=1,NKT(0,nb)
            ZE_LOCAL(nk+(nn-1)*NKT(0,nb))=1.0D0
            DO i=1,21
              XI(1)=DBLE(i-1)/20.0D0
              POINTS(1,i)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          1,XI,XE_LOCAL)
              POINTS(2,i)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          1,XI,ZE_LOCAL)
            ENDDO
            CALL POLYLINE(INDEX,1,21,POINTS,ERROR,*9999)
            ZE_LOCAL(nk+(nn-1)*NKT(0,nb))=0.0D0
          ENDDO
        ENDDO
      ELSE IF(NIT(nb).EQ.2) THEN
      ELSE IF(NIT(nb).EQ.3) THEN
      ENDIF
      CALL CLOSE_SEGMENT(ISBASE,iw,ERROR,*9999)

      CALL EXITS('SGBASE')
      RETURN
 9999 CALL ERRORS('SGBASE',ERROR)
      CALL EXITS('SGBASE')
      RETURN 1
      END


      SUBROUTINE SGFACE(INDEX,ISEG,ISFACE,ISFANO,iw,nf,NTDX,CSEG,
     '  POINTS,ERROR,*)

C#### Subroutine: SGFACE
C###  Description:
C###    SGFACE creates new face segment ISFACE(iw,nf) containing shading
C###    and ISFANO(iw,nf) containing face numbers.

C     I can't find any calls to this, AAY.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISFACE(NWM,NFM),ISFANO(NWM,NFM),iw,nf,NTDX
      REAL*8 POINTS(3,*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD,NMID
      CHARACTER CHAR*3

      CALL ENTERS('SGFACE',*9999)

      CALL OPEN_SEGMENT(ISFACE(iw,nf),ISEG,iw,'FACE',INDEX,INDEX_OLD,
     '  nf,1,CSEG,ERROR,*9999)
      CALL POLYLINE(INDEX,iw,NTDX,POINTS,ERROR,*9999)
      CALL CLOSE_SEGMENT(ISFACE(iw,nf),iw,ERROR,*9999)

      CALL OPEN_SEGMENT(ISFANO(iw,nf),ISEG,iw,'FACE',INDEX,INDEX_OLD,
     '  nf,1,CSEG,ERROR,*9999)
      WRITE(CHAR,'(I3)') nf
      CALL TRIM(CHAR,IBEG,IEND)
      IF(NTDX.EQ.2) THEN
        POINTS(1,1)=(POINTS(1,1)+POINTS(1,2))/2.0D0
        POINTS(2,1)=(POINTS(2,1)+POINTS(2,2))/2.0D0
        POINTS(3,1)=(POINTS(3,1)+POINTS(3,2))/2.0D0
      ELSE
        NMID=NTDX/2
        POINTS(1,1)=POINTS(1,NMID)
        POINTS(2,1)=POINTS(2,NMID)
        POINTS(3,1)=POINTS(3,NMID)
      ENDIF
      CALL TEXT(INDEX,iw,CHAR(IBEG:IEND),POINTS(1,1),ERROR,*9999)
      CALL CLOSE_SEGMENT(ISFANO(iw,nf),iw,ERROR,*9999)

      CALL EXITS('SGFACE')
      RETURN
 9999 CALL ERRORS('SGFACE',ERROR)
      CALL EXITS('SGFACE')
      RETURN 1
      END


      SUBROUTINE SGPLIN(INDEX,ISEG,ISPLIN,iw,NTPTS,CSEG,XLIST,YLIST,
     '  ERROR,*)

C#### Subroutine: SGPLIN
C###  Description:
C###    SGPLIN creates polyline segment ISPLIN.

      IMPLICIT NONE
!     Parameter List
      INTEGER INDEX,ISEG(*),ISPLIN,iw,NTPTS
      REAL*8 XLIST(*),YLIST(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD,MAXPTS,nopt
      PARAMETER (MAXPTS=500)  !note limit of 500 points
      REAL*8 PTS(3,MAXPTS)

      CALL ENTERS('SGPLIN',*9999)
      CALL OPEN_SEGMENT(ISPLIN,ISEG,iw,'PLIN',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      IF(NTPTS.GT.MAXPTS) THEN
        ERROR='Too many points in SGPLIN'
        GOTO 9999
      ENDIF
      DO nopt=1,NTPTS
        PTS(1,nopt)=XLIST(nopt)
        PTS(2,nopt)=YLIST(nopt)
      ENDDO
      CALL POLYLINE(INDEX,iw,NTPTS,PTS,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISPLIN,iw,ERROR,*9999)

      CALL EXITS('SGPLIN')
      RETURN
 9999 CALL ERRORS('SGPLIN',ERROR)
      CALL EXITS('SGPLIN')
      RETURN 1
      END


C KAT 2001-12-13 removed from SGPROF

      ELSE IF(TYPE(1:5).EQ.'IMAGE') THEN
        TITLE=' Image profile'
        NTCURV=1
        ND1=NSOBJE(3,NO_OBJECT) !is 1st nd in object
        ND2=NSOBJE(4,NO_OBJECT) !is 2nd nd in object
        nox=0
        DO nd=ND1,ND2
          IMAGE_X=1+NINT((ZD(1,nd)-DBLE(XMIN))*511.d0/DBLE(XMAX-XMIN))
          IMAGE_Y=1+NINT((ZD(2,nd)-DBLE(YMIN))*511.d0/DBLE(YMAX-YMIN))
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' At nd='',I5,''  Image x,y is'',2I4)')
     '        nd,IMAGE_X,IMAGE_Y
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          TOT=0.0D0
          DO j=MAX(IMAGE_Y-(NOPIXELS-1),0),MIN(IMAGE_Y+(NOPIXELS-1),
     '      IMGY)
            DO i=MAX(IMAGE_X-(NOPIXELS-1),0),MIN(IMAGE_X+(NOPIXELS-1),
     '        IMGX)
              TOT=TOT+DBLE(I2P((j-1)*IMGX+i,1,1))
            ENDDO
          ENDDO
          ZVAL=TOT*RINTENSITY_MAX/DBLE((2*NOPIXELS-1)*(2*NOPIXELS-1))
          IF(DOP) THEN
            WRITE(OP_STRING,'('' zval='',E12.3)') ZVAL
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          nox=nox+1
          IF(nox.EQ.1) THEN
            XPROF(nox)=0.0D0
          ELSE IF(nox.GT.1) THEN
            XPROF(nox)=XPROF(nox-1)+DSQRT((ZD(1,nd)-ZD(1,nd-1))**2
     '                                  +(ZD(2,nd)-ZD(2,nd-1))**2)
          ENDIF
          YPROF(nox,1)=ZVAL
        ENDDO
        NTX=nox



      SUBROUTINE SGSIGN(ICOL,INDEX,ISAMPLES,ISEG,ISIGNAL,ISSIGN,iw,CSEG,
     '  SIGHEADER,DATUM,FIDMARK,TSTART,TEND,ACCEPT,ERROR,*)

C#### Subroutine: SGSIGN
C###  Description:
C###    SGSIGN creates signal segment ISSIGN.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:pics00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
!     Parameter List
      INTEGER ICOL,INDEX(*),ISAMPLES,ISEG(*),ISIGNAL,ISSIGN(2),iw
      REAL*8 DATUM,FIDMARK,TEND,TSTART
      CHARACTER CSEG(*)*(*),ERROR*(*),SIGHEADER*(*)
      LOGICAL ACCEPT
!     Local Variables
      INTEGER ABSMAX,I,IBEG,IEND,IND,INDEX_OLD,ISTEP,J,K,NPOINTS
      REAL*8 DPOS,PL(3,2),PM(3,2048),PT(2),SAMPLES,XPOS,YPOS
      CHARACTER AXES*10

      CALL ENTERS('SGSIGN',*9999)

      SAMPLES=(TEND-TSTART)*DBLE(ISAMPLES)
      ISTEP=NINT(SAMPLES/200.0D0)
      IF(ISTEP.LT.1) ISTEP=1

      IF(ISIGNAL.LE.43) THEN      !Get Column (J) and Row (I)
        J=1
        I=ISIGNAL
      ELSE IF(ISIGNAL.LE.86) THEN
        J=2
        I=ISIGNAL-43
      ELSE
        J=3
        I=ISIGNAL-86
      ENDIF

      CALL OPEN_SEGMENT(ISSIGN(1),ISEG,iw,'SIGN',INDEX(1),INDEX_OLD,
     '  ISIGNAL,1,CSEG,ERROR,*9999)                 !Start of Axes Segment

      XPOS=DBLE(J)*10.0D0-8.5D0
      YPOS=435.0D0-DBLE(I)*10.0D0
      PL(1,1)=XPOS
      PL(2,1)=YPOS
      PL(1,2)=XPOS+8.5D0
      PL(2,2)=YPOS
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes
      PL(1,1)=XPOS
      PL(2,1)=YPOS+4.8D0
      PL(1,2)=XPOS
      PL(2,2)=YPOS-4.8D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes
      PT(1)=XPOS-0.3D0
      PT(2)=YPOS+0.3D0
      AXES=SIGHEADER
      CALL TRIM(AXES,IBEG,IEND)
      IF(FIDMARK.EQ.DATUM) THEN
        IND=8
      ELSE
        IND=5
      ENDIF
      CALL TEXT(INDEX(IND),iw,AXES(IBEG:IEND),PT,ERROR,*9999) !Label axes

      IF(DATUM.GE.TSTART.AND.DATUM.LE.TEND) THEN
        DPOS=(DATUM-TSTART)/(TEND-TSTART)
        PL(1,1)=XPOS+8.5D0*DPOS
        PL(2,1)=YPOS+5.0D0
        PL(1,2)=XPOS+8.5D0*DPOS
        PL(2,2)=YPOS-5.0D0
        CALL POLYLINE(INDEX(12),iw,2,PL,ERROR,*9999) !Draw Datum
      ENDIF
      IF(FIDMARK.GE.TSTART.AND.FIDMARK.LE.TEND.AND.ACCEPT) THEN
        PL(1,1)=XPOS+8.5D0*FIDMARK
        PL(2,1)=YPOS+5.0D0
        PL(1,2)=XPOS+8.5D0*FIDMARK
        PL(2,2)=YPOS-5.0D0
        CALL POLYLINE(INDEX(11),iw,2,PL,ERROR,*9999) !Draw Fiducial Marker
      ENDIF

      CALL CLOSE_SEGMENT(ISSIGN(1),iw,ERROR,*9999) !End of Axes Segment

      CALL OPEN_SEGMENT(ISSIGN(2),ISEG,iw,'SIGN',INDEX(ICOL),INDEX_OLD,
     '  ISIGNAL,1,CSEG,ERROR,*9999)                 !Start of Signal Segment

      NPOINTS=0
      ABSMAX=1000000
      DO K=1,INT(SAMPLES),ISTEP
        IF(ABS(I2P(ISIGNAL,K,1)).GT.ABSMAX) ABSMAX=ABS(I2P(ISIGNAL,K,1))
      ENDDO
      IF(DOP) THEN
        WRITE(OP_STRING,*) 'Absolute Maximum is ',ABSMAX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO K=1,INT(SAMPLES),ISTEP
        NPOINTS=NPOINTS+1
        PM(1,NPOINTS)=XPOS+DBLE(K)/SAMPLES*8.5D0
        PM(2,NPOINTS)=YPOS+4.8D0*DBLE(I2P(ISIGNAL,K,1))/DBLE(ABSMAX)
      ENDDO

      CALL POLYLINE(INDEX(ICOL),iw,NPOINTS,PM,ERROR,*9999) !Signal

      CALL CLOSE_SEGMENT(ISSIGN(2),iw,ERROR,*9999) !End of Signal Segment

      CALL EXITS('SGSIGN')
      RETURN
 9999 CALL ERRORS('SGSIGN',ERROR)
      CALL EXITS('SGSIGN')
      RETURN 1
      END


C GMH 18-4-96 I dont think this is used...
      SUBROUTINE SGTEXT(INDEX,ISEG,ISTEXT,iw,A_LIST,NTLIST,
     '  CSEG,STR,XWC,YWC,ERROR,*)

C#### Subroutine: SGTEXT
C###  Description:
C###    SGTEXT creates text segment ISTEXT.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:map000.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISTEXT,iw,A_LIST(*),NTLIST
      REAL*8 XWC,YWC
      CHARACTER CSEG(*)*(*),STR(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX_OLD,nolist,NOTEXT
      REAL*8 DINCR,DY,PT(3)

      CALL ENTERS('SGTEXT',*9999)
      CALL OPEN_SEGMENT(ISTEXT,ISEG,iw,'TEXT',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)
      IF(iw.LE.2) THEN
        DINCR=DBLE(YMAX-YMIN)/35.0D0
      ELSE IF(iw.EQ.4) THEN
        DINCR=0.04D0
      ENDIF

      DO nolist=1,NTLIST
        DY=(nolist-1)*DINCR
        NOTEXT=A_LIST(nolist)
        CALL TRIM(STR(NOTEXT),IBEG,IEND)
        PROJEC='RECTANGULAR'
        PT(1)=XWC
        PT(2)=YWC-DY
        CALL TEXT(INDEX,iw,STR(NOTEXT)(IBEG:IEND),PT,ERROR,*9999)
      ENDDO

      CALL CLOSE_SEGMENT(ISTEXT,iw,ERROR,*9999)

      CALL EXITS('SGTEXT')
      RETURN
 9999 CALL ERRORS('SGTEXT',ERROR)
      CALL EXITS('SGTEXT')
      RETURN 1
      END


      SUBROUTINE SGTRAC(ICOL,INDEX,ISAMPLES,ISEG,ISIGNAL,ISTRAC,iw,CSEG,
     '  SIGHEADER,DATUM,FIDMARK,TSTART,TEND,ERROR,*)

C#### Subroutine: SGTRAC
C###  Description:
C###    SGTRAC creates trace segment ISTRAC.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:pics00.cmn'
!     Parameter List
      INTEGER ICOL,INDEX(*),ISAMPLES,ISEG(*),ISIGNAL,ISTRAC(10),iw
      REAL*8 DATUM,FIDMARK,TEND,TSTART
      CHARACTER CSEG(*)*(*),ERROR*(*),SIGHEADER*(*)
!     Local Variables
      INTEGER IBEG,IEND,IEXPON,IFINISH,INDEX_OLD,ISTART,ISTEP,K,
     '  NPOINTS
      REAL*8 ABSMAX,ADIFF,AMAX,AMIN,DELAY,DPOS,FACTOR,FPOS,PL(3,2),
     '  PM(3,2048),PT(2),SAMPLES,XPOS,YPOS
      CHARACTER AXES*10,C1*10,CHARVAL*10

      CALL ENTERS('SGTRAC',*9999)

      CALL OPEN_SEGMENT(ISTRAC(1),ISEG,iw,'TRAC',INDEX(1),INDEX_OLD,
     '  10,1,CSEG,ERROR,*9999)                 !Start of Signal Segment

      SAMPLES=DBLE(ISAMPLES)
      ISTEP=5

      AMAX=1.0D-6
      AMIN=-1.0D-6
      DO K=1,ISAMPLES
        IF(I2P(ISIGNAL,K,1).GT.AMAX) AMAX=I2P(ISIGNAL,K,1)
        IF(I2P(ISIGNAL,K,1).LT.AMIN) AMIN=I2P(ISIGNAL,K,1)
      ENDDO
      ADIFF=AMAX-AMIN

      XPOS=0.9D0
      YPOS=9.95D0-1.9D0*AMAX/ADIFF
      PL(1,1)=XPOS
      PL(2,1)=YPOS
      PL(1,2)=XPOS+9.0D0
      PL(2,2)=YPOS
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Horizontal
      PT(1)=XPOS+0.2D0
      PT(2)=YPOS-0.3D0
      CALL TEXT(INDEX(7),iw,'0.0',PT,ERROR,*9999)
      PL(1,1)=XPOS+9.0D0/4.0D0
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0/4.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 0.25 Marker
      PL(1,1)=XPOS+9.0D0/2.0D0
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0/2.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 0.5 Marker
      PL(1,1)=XPOS+9.0D0/4.0D0*3.0D0
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0/4.0D0*3.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 0.75 Marker
      PL(1,1)=XPOS+9.0D0       
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 1.0 Marker
      PT(1)=XPOS+8.8D0
      PT(2)=YPOS-0.3D0
      CALL TEXT(INDEX(7),iw,'1.0',PT,ERROR,*9999)

      PL(1,1)=XPOS+9.0D0*DATUM
      PL(2,1)=8.3D0
      PL(1,2)=XPOS+9.0D0*DATUM
      PL(2,2)=9.7D0
      CALL POLYLINE(INDEX(12),iw,2,PL,ERROR,*9999) !Draw Datum
      PL(1,1)=XPOS+9.0D0*FIDMARK
      PL(2,1)=8.3D0
      PL(1,2)=XPOS+9.0D0*FIDMARK
      PL(2,2)=9.7D0
      CALL POLYLINE(INDEX(11),iw,2,PL,ERROR,*9999) !Draw Fiducial Marker

      PL(1,1)=XPOS
      PL(2,1)=9.95D0
      PL(1,2)=XPOS
      PL(2,2)=8.05D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Vertical
      PL(1,1)=XPOS-0.05D0
      PL(2,1)=9.95D0
      PL(1,2)=XPOS+0.05D0
      PL(2,2)=9.95D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Top Marker
      PL(1,1)=XPOS-0.05D0
      PL(2,1)=8.05D0
      PL(1,2)=XPOS+0.05D0
      PL(2,2)=8.05D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Bottom Marker

      PT(1)=XPOS-0.05D0
      PT(2)=YPOS+0.05D0
      AXES=SIGHEADER
      CALL TRIM(AXES,IBEG,IEND)
      CALL TEXT(INDEX(5),iw,AXES(IBEG:IEND),PT,ERROR,*9999) !Label axes #

      ABSMAX=DMAX1(AMAX,DABS(AMIN))
      IEXPON=INT(DLOG10(ABSMAX))
      FACTOR=10.0D0**IEXPON
      PT(1)=XPOS+1.6D0
      PT(2)=9.9D0
      WRITE(C1,'(I3)') IEXPON
      CALL TRIM(C1,IBEG,IEND)
      CHARVAL='(x1E'//C1(IBEG:IEND)//')'
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label exp

      PT(1)=XPOS
      PT(2)=9.9D0
      WRITE(CHARVAL,'(F5.2)') AMAX/FACTOR
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label max
      PT(1)=XPOS
      PT(2)=8.2D0
      WRITE(CHARVAL,'(F5.2)') AMIN/FACTOR
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label min

      NPOINTS=0
      DO K=1,INT(SAMPLES),ISTEP
        NPOINTS=NPOINTS+1
        PM(1,NPOINTS)=XPOS+DBLE(K)/SAMPLES*9.0D0
        PM(2,NPOINTS)=YPOS+1.9D0*DBLE(I2P(ISIGNAL,K,1))/ADIFF
      ENDDO
      CALL POLYLINE(INDEX(ICOL),iw,NPOINTS,PM,ERROR,*9999) !Signal

      PL(1,1)=0.0D0
      PL(2,1)=8.0D0
      PL(1,2)=10.0D0
      PL(2,2)=8.0D0
      CALL POLYLINE(INDEX(6),iw,2,PL,ERROR,*9999) !Partition line


      ISTART=INT(TSTART*SAMPLES)
      IFINISH=INT(TEND*SAMPLES)
      AMAX=1.0D-6
      AMIN=-1.0D-6
      DO K=ISTART+1,IFINISH
        IF(I2P(ISIGNAL,K,1).GT.AMAX) AMAX=I2P(ISIGNAL,K,1)
        IF(I2P(ISIGNAL,K,1).LT.AMIN) AMIN=I2P(ISIGNAL,K,1)
      ENDDO
      ADIFF=AMAX-AMIN
      ABSMAX=DMAX1(AMAX,DABS(AMIN))
      IEXPON=INT(DLOG10(ABSMAX))
      FACTOR=10.0D0**IEXPON

      XPOS=0.9D0
      YPOS=7.95D0-7.9D0*AMAX/ADIFF
      PL(1,1)=XPOS
      PL(2,1)=YPOS
      PL(1,2)=XPOS+9.0D0
      PL(2,2)=YPOS
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Horizontal
      PT(1)=XPOS+0.3D0
      PT(2)=YPOS-0.3D0
      WRITE(CHARVAL,'(F5.2)') TSTART
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label max
      PL(1,1)=XPOS+9.0D0/4.0D0
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0/4.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 0.25 Marker
      PL(1,1)=XPOS+9.0D0/2.0D0
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0/2.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 0.5 Marker
C CPB 21/3/94 is this + 9/(4*3) or 27/4 ???
      PL(1,1)=XPOS+9.0D0/4.0D0*3.0D0
      PL(2,1)=YPOS+0.1D0
C CPB 21/3/94 is this + 9/(4*3) or 27/4 ???
      PL(1,2)=XPOS+9.0D0/4.0D0*3.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 0.75 Marker
      PL(1,1)=XPOS+9.0D0
      PL(2,1)=YPOS+0.1D0
      PL(1,2)=XPOS+9.0D0
      PL(2,2)=YPOS-0.1D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - 1.0 Marker
      PT(1)=XPOS+8.7D0
      PT(2)=YPOS-0.3D0
      WRITE(CHARVAL,'(F5.2)') TEND
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999)

      IF(DATUM.GE.TSTART.AND.DATUM.LE.TEND) THEN
        DPOS=(DATUM-TSTART)/(TEND-TSTART)
        PL(1,1)=XPOS+9.0D0*DPOS
        PL(2,1)=0.3D0
        PL(1,2)=XPOS+9.0D0*DPOS
        PL(2,2)=7.7D0
        CALL POLYLINE(INDEX(12),iw,2,PL,ERROR,*9999) !Draw Datum
      ENDIF

      PL(1,1)=XPOS
      PL(2,1)=7.95D0
      PL(1,2)=XPOS
      PL(2,2)=0.05D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Vertical
      PL(1,1)=XPOS-0.05D0
      PL(2,1)=7.95D0
      PL(1,2)=XPOS+0.05D0
      PL(2,2)=7.95D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Top Marker
      PL(1,1)=XPOS-0.05D0
      PL(2,1)=0.05D0
      PL(1,2)=XPOS+0.05D0
      PL(2,2)=0.05D0
      CALL POLYLINE(INDEX(1),iw,2,PL,ERROR,*9999) !Draw Axes - Bottom Marker

      PT(1)=XPOS-0.05D0
      PT(2)=YPOS+0.05D0
      AXES=SIGHEADER
      CALL TRIM(AXES,IBEG,IEND)
      CALL TEXT(INDEX(5),iw,AXES(IBEG:IEND),PT,ERROR,*9999) !Label axes

      PT(1)=XPOS+1.6D0
      PT(2)=7.85D0
      WRITE(C1,'(I3)') IEXPON
      CALL TRIM(C1,IBEG,IEND)
      CHARVAL='(x1E'//C1(IBEG:IEND)//')'
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label exp

      PT(1)=XPOS
      PT(2)=7.9D0
      WRITE(CHARVAL,'(F5.2)') AMAX/FACTOR
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label max
      PT(1)=XPOS
      PT(2)=0.2D0
      WRITE(CHARVAL,'(F5.2)') AMIN/FACTOR
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,CHARVAL(IBEG:IEND),PT,ERROR,*9999) !Label min
      PT(1)=XPOS+8.0D0
      PT(2)=7.8D0
      DELAY=DBLE(INT((FIDMARK-DATUM)*2000.0D0))/2.0D0
      WRITE(CHARVAL,'(F6.1)') DELAY
      CALL TRIM(CHARVAL,IBEG,IEND)
      CALL TEXT(INDEX(7),iw,'Delay : '//CHARVAL(IBEG:IEND)//' ms',PT,
     '  ERROR,*9999)                                             !Label delay

      NPOINTS=0
      DO K=ISTART+1,IFINISH
        NPOINTS=NPOINTS+1
        PM(1,NPOINTS)=XPOS+DBLE(K-ISTART)/DBLE(IFINISH-ISTART)*9.0D0
        PM(2,NPOINTS)=YPOS+7.9D0*DBLE(I2P(ISIGNAL,K,1))/ADIFF
      ENDDO
      CALL POLYLINE(INDEX(ICOL),iw,NPOINTS,PM,ERROR,*9999) !Signal

      CALL CLOSE_SEGMENT(ISTRAC(1),iw,ERROR,*9999) !End of Signal Segment

      CALL OPEN_SEGMENT(ISTRAC(2),ISEG,iw,'FID_MARK',INDEX(11),
     '  INDEX_OLD,11,1,CSEG,ERROR,*9999)                 !Start of Fiducial Segment

      IF(FIDMARK.GE.TSTART.AND.FIDMARK.LE.TEND) THEN
        FPOS=(FIDMARK-TSTART)/(TEND-TSTART)
        PL(1,1)=XPOS+9.0D0*FPOS
        PL(2,1)=0.3D0
        PL(1,2)=XPOS+9.0D0*FPOS
        PL(2,2)=7.7D0
        CALL POLYLINE(INDEX(11),iw,2,PL,ERROR,*9999) !Draw Fiducial Marker
      ENDIF

      CALL CLOSE_SEGMENT(ISTRAC(2),iw,ERROR,*9999) !End of Fiducial Segment

      CALL EXITS('SGTRAC')
      RETURN
 9999 CALL ERRORS('SGTRAC',ERROR)
      CALL EXITS('SGTRAC')
      RETURN 1
      END


C LC 25/2/97 removed from :
C#### Subroutine: SGSURF
C###  Description:
C###    SGSURF creates element surface segment ISSURF.
C**** 22Aug89: Modified to interpolate surface points correctly in
C****          Lagrange/Hermite tensor-product elements, using ZEZW.
C****          Added DEFORM,NAN,NHE,NW,XG,ZG,ERROR to param list. ADMcC

C****   I don't think anyone wants to use this, delete it? AAY 18-12-90
CX      ELSE IF(SURFACE_TYPE(1:5).EQ.'LINES') THEN
CX        nb=NBJ(1,ne)
CX        IF(IBT(1,1,nb).NE.5) THEN
CX***       Lagrange/Hermite tensor product
CX          ISPL=1
CX        ELSE
CX***       B-spline tensor product
CX          ISPL=2
CX        ENDIF
CX
CX        IF(NIT(nb).EQ.3.AND.STATIC) THEN
CX***       Define nodes at interpolated position given by Xi_3
CX          NSTB2=NST(nb)/2
CX          DO nj=1,NJE(ne)
CX            DO ns=1,NSTB2
CX              XE(ns,nj)=(1.0D0-XI3)*XE(ns,nj)+XI3*XE(ns+NSTB2,nj)
CX            ENDDO
CX          ENDDO
CX***       Calculate control point array
CX          CALL XECP(IBT,NBJ(1,ne),NJE(ne),NUM_CP,CP,XE,ERROR,*9999)
CX***       Calculate refined control point array
CX          CALL CPCP2(ISPL,NQE(1,1,ne),NUM_CP,CP,NUM_CPSUB,CPSUB)
CX        ELSE IF(.NOT.STATIC) THEN
CX          XI(3)=XI3
CX          XI(NIT(NBJ(1,ne))+1)=TIME
CX          NUM_CPSUB(1)=7
CX          NUM_CPSUB(2)=7
CX          DO NI1=1,NUM_CPSUB(1)
CX            XI(1)=DBLE(NI1-1)/DBLE(NUM_CPSUB(1)-1)
CX            DO NI2=1,NUM_CPSUB(2)
CX              XI(2)=DBLE(NI2-1)/DBLE(NUM_CPSUB(2)-1)
CX              DO nj=1,NJT
CX                nb=NBJ(nj,ne)
CX                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
CX   '              nb,1,XI,XE(1,nj))
CX                nb=NBH(nj,nc,ne)
CX                Z2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
CX   '              nb,1,XI,ZE(1,nj))
CX              ENDDO
CX              CALL XZ(ITYP10(1),X,Z)
CX              DO nj=1,NJT
CX                IF(KTYP58(nr).EQ.2) THEN !displacement
CX                  CPSUB(nj,NI1,NI2)=Z(nj)+Z2(nj)
CX                ELSE
CX                  CPSUB(nj,NI1,NI2)=Z2(nj)
CX                ENDIF
CX              ENDDO
CX            ENDDO
CX          ENDDO
CX        ELSE IF(NIT(nb).EQ.2.AND.ISPL.EQ.1) THEN
CX          NUM_CPSUB(1)=9
CX          NUM_CPSUB(2)=9
CX          DO NI1=1,NUM_CPSUB(1)
CX            XI(1)=DBLE(NI1-1)/DBLE(NUM_CPSUB(1)-1)
CX            DO NI2=1,NUM_CPSUB(2)
CX              XI(2)=DBLE(NI2-1)/DBLE(NUM_CPSUB(2)-1)
CX              IF(DEFORM) THEN
CX                CALL ZEZW(0,IBT,IDO,INP,NAN,NBH(1,1,ne),NJE(ne),
CX   '               NRE(ne),nx,DXIX,ZE,ZG,XI,ERROR,*9999)
CX                DO nj=1,NJE(ne)
CX                  X(nj)=ZG(nj,1)
CX                ENDDO
CX              ELSE
CX                CALL XEXW(IBT,IDO,INP,NAN,NBJ(1,ne),NRE(ne),XE,XG,XI,
CX   '              ERROR,*9999)
CX                DO nj=1,NJE(ne)
CX                  X(nj)=XG(nj,1)
CX                ENDDO
CX              ENDIF
CX              CALL XZ(ITYP10(1),X,Z)
CX              DO nj=1,NJE(ne)
CX                CPSUB(nj,NI1,NI2)=Z(nj)
CX              ENDDO
CX            ENDDO
CX          ENDDO
CX        ELSE
CX***       Calculate control point array
CX          CALL XECP(IBT,NBJ(1,ne),NJE(ne),NPF(1,1),NUM_CP,CP,XE,
CX    '       ERROR,*9999)
CX***       Calculate refined control point array
CX          CALL CPCP2(ISPL,NQE(1,1,ne),NUM_CP,CP,NUM_CPSUB,CPSUB)
CX        ENDIF
CX
CX***     Draw lines in xi1 direction
CX        DO NI2=1,NUM_CPSUB(2)
CX          DO NI1=1,NUM_CPSUB(1)
CX            DO nj=1,NJT
CX              PXL(nj,NI1)=CPSUB(nj,NI1,NI2)
CX            ENDDO
CX            IF(ktyp8.NE.5) CALL XZ(ITYP10(1),PXL(1,NI1),PXL(1,NI1))
CX          ENDDO
CX          CALL POLYLINE(3,iw,NUM_CPSUB(1),PXL,ERROR,*9999)
CX        ENDDO
CX***     Draw lines in xi2 direction
CX        DO NI1=1,NUM_CPSUB(1)
CX          DO NI2=1,NUM_CPSUB(2)
CX            DO nj=1,NJT
CX              PXL(nj,NI2)=CPSUB(nj,NI1,NI2)
CX            ENDDO
CX            IF(ktyp8.NE.5) CALL XZ(ITYP10(1),PXL(1,NI2),PXL(1,NI2))
CX          ENDDO
CX          CALL POLYLINE(3,iw,NUM_CPSUB(2),PXL,ERROR,*9999)
CX        ENDDO


Module FE25
===========

C KAT 2001-12-14
      SUBROUTINE SGELEC(ISEG,ISELEC,iw,nd,CSEG,ZD,ERROR,*)

C#### Subroutine: SGELEC
C###  Description:
C###    SGELEC creates new electrode segment ISELEC(nd).

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:read00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISELEC,iw,nd
      REAL*8 ZD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX1,INDEX2,INDEX_OLD,INDEX_TEXT
      REAL*8 Z(3)
      CHARACTER CHAR*5

      CALL ENTERS('SGELEC',*9999)

      INDEX1=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')      !Accepted
      INDEX2=INDEX_TEXT(0,'WIDTH1','FONT1','WHITE')      !Rejected

      CALL OPEN_SEGMENT(ISELEC,ISEG,iw,'ELECTRODE',INDEX,INDEX_OLD,
     '  nd,1,CSEG,ERROR,*9999)
      CHAR=SIGHEADER(1,nd)
      CALL STRING_TRIM(CHAR,IBEG,IEND)
      CALL ZX(ITYP10(1),ZD(1,nd),Z) !transform to polar coords
      IF(ACCEPT(nd)) THEN
        CALL TEXT(INDEX1,iw,CHAR(IBEG:IEND),Z,ERROR,*9999)
      ELSE
        CALL TEXT(INDEX2,iw,CHAR(IBEG:IEND),Z,ERROR,*9999)
      ENDIF
      CALL CLOSE_SEGMENT(ISELEC,iw,ERROR,*9999)

      CALL EXITS('SGELEC')
      RETURN
 9999 CALL ERRORS('SGELEC',ERROR)
      CALL EXITS('SGELEC')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE SGIMAG(INDEX,ISEG,ISIMAG,iw,NDIMX,NDIMY,CSEG,XMIN_CA,
     '  XMAX_CA,YMIN_CA,YMAX_CA,ERROR,*)

C#### Subroutine: SGIMAG
C###  Description:
C###    SGIMAG creates image segment ISIMAG.

      IMPLICIT NONE
!     Parameter List
      INTEGER INDEX,ISEG(*),ISIMAG,iw,NDIMX,NDIMY
      REAL XMAX_CA,XMIN_CA,YMAX_CA,YMIN_CA
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER INDEX_OLD

      CALL ENTERS('SGIMAG',*9999)

      CALL OPEN_SEGMENT(ISIMAG,ISEG,iw,'IMAG',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)
      CALL CELL_ARRAY(iw,NDIMX,NDIMY,XMIN_CA,XMAX_CA,YMIN_CA,YMAX_CA,
     '  ERROR,*9999)
      CALL CLOSE_SEGMENT(ISIMAG,iw,ERROR,*9999)

      CALL EXITS('SGIMAG')
      RETURN
 9999 CALL ERRORS('SGIMAG',ERROR)
      CALL EXITS('SGIMAG')
      RETURN 1
      END


Module FE26
===========

C KAT 2001-12-13 removed from DIPROF

C#### Command: FEM display profile image
C###  Parameter:       <average #PIXELS[1]>
C###    Specify the average number of pixels to display.
C###  Parameter:       <with OBJECT_NAME[]>
C###    Specify and object name
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan) 
C###  Description:
C###    Displays an image of a profile

        OP_STRING(1)=STRING(1:IEND) //' image'
        OP_STRING(2)=BLANK(1:15)//'<average #PIXELS[1]>'
        OP_STRING(3)=BLANK(1:15)//'<with OBJECT_NAME['
     '                          //OBJECT(1:IEND3)//']>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

        ELSE IF(CBBREV(CO,'IMAGE',1,noco+1,NTCO,N3CO)) THEN
          TYPE='IMAGE'
          ID_TYPE=1

        IF(CBBREV(CO,'AVERAGE',1,noco+1,NTCO,N3CO)) THEN
          NOPIXELS=IFROMC(CO(N3CO+1))
        ELSE
          NOPIXELS=1
        ENDIF

        ELSE IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT'.
     '    OR.TYPE(1:5).EQ.'IMAGE') THEN



C LKC 8-JUN-1999 EXGEOM - map3d code
C LKC 8-JUN-1999 EXGEOM_NE - map3d code

      SUBROUTINE EXGEOM(IBT,IDO,INP,NBJ,NEELEM,NELIST,
     '  NKE,NPF,NPLIST,NPNE,NRE,NRLIST,NVJE,
     '  SE,STRING,XA,XE,XP,ERROR,*)

C#### Subroutine: EXGEOM
C###  Description:
C###    EXGEOM exports geometry data from finite element data base. 

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grou00.cmn'
      INCLUDE 'cmiss$reference:map3d.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NKE(NKM,NNM,NBFM,NEM),NPF(9,NFM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NRE(NEM),NRLIST(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ELGROUP,I,IBEG,IBEG1,IBEG2,ICHANNELS(1000,2),
     '  ICH,ICHT,IEND,IEND1,IEND2,LOCAL_NN,LOCAL_NP(100),N3CO,
     '  nb,ne,nj,N_NODE,N_CHANNEL,nn,np,nr,NOELEM,NOGREL,
     '  NOGRNO,NOGROUP,NONODE,NO_NRLIST,nv
      CHARACTER FILE*100,GROUP_NAME*30,OUTPUT*8,TYPE*6
      LOGICAL ALL_REGIONS,CBBREV,INLIST

      CALL ENTERS('EXGEOM',*9999)
      nv=1 ! Temporary
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM export geometry<;FILENAME[default]>
C###  Parameter:      <to map3d[map3d]>
C###  Parameter:      <number (ELEMENTS/GROUP/all)[all]>
C###    Defines the elements or group to include.
C###  Parameter:      <region (#s/all)[1]>
C###    Limit to elements of specified regions.
C###  Description:
C###    Export geometry data from finite element data base.

        OP_STRING(1)=STRING(1:IEND)
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<to map3d[map3d]>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<number (ELEMENTS/GROUP/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','EXELEM',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
        ELSE
          FILE=FILE00
        ENDIF
        CALL TRIM(FILE,IBEG,IEND)
       
        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'MAP3D',2,noco+1,NTCO,N3CO)) THEN
            OUTPUT='MAP3D'
          ELSE
            OUTPUT='MAP3D'
          ENDIF
        ELSE
          OUTPUT='MAP3D'
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'NUMBER',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,*9999)
          TYPE='NUMBER'
        ELSE IF(CBBREV(CO,'GROUP',2,noco+1,NTCO,N3CO)) THEN
          CALL TRIM(CO(N3CO+1),IBEG,IEND)
          GROUP_NAME=CO(N3CO+1)(IBEG:IEND)
          TYPE='GROUP'
c cpb 12/9/94 Need to find out the group numbers !!!!!!!!!!!!!
        ELSE
          DO NO_NRLIST=1,NRLIST(0)
            nr=NRLIST(NO_NRLIST)
            NELIST(0)=NEELEM(0,NR)
            DO NOELEM=1,NEELEM(0,NR)
              NELIST(NOELEM)=NEELEM(NOELEM,NR)
            ENDDO
          ENDDO
          TYPE='NUMBER'
        ENDIF

        IF(OUTPUT(1:5).EQ.'MAP3D') THEN

C CPB 13/2/94 Only have put in very simple element type checking

          CALL ASSERT(CALL_EXPO,'>>Call DEFINE EXPORT first',
     '      ERROR,*9999)
          nb=NBJ(1,NELIST(1))
          IF(NIT(NBJ(1,NB)).NE.2) THEN
            ERROR='>>Can only export geometry in map3d format for '
     '        //'2D elements'
            GOTO 9999
          ENDIF
          IF(.NOT.(IBT(1,1,NB).GE.2.AND.IBT(1,2,NB).GE.2)) THEN
            ERROR='>>Can only export geometry in map3d format for '
     '        //'BiHermite elements'
            GOTO 9999
          ENDIF

          IF(TYPE(1:6).EQ.'NUMBER') THEN
            NPLIST(0)=0
            DO NOELEM=1,NELIST(0)
              NE=NELIST(NOELEM)
              nb=NBJ(1,NE)
              DO nn=1,NNT(nb)
                NP=NPNE(nn,nb,NE)
                INLIST=.FALSE.
                DO NONODE=1,NPLIST(0)
                  IF(NP.EQ.NPLIST(NONODE)) INLIST=.TRUE.
                ENDDO
                IF(.NOT.INLIST) THEN
                  NPLIST(0)=NPLIST(0)+1
                  NPLIST(NPLIST(0))=NP
                ENDIF
              ENDDO
            ENDDO
          ENDIF

C 25/2/97 LC archived section : 
C Below is CPBs original .fac, .pts and .channels output.
C Altered AJP 9-3-94

C Write out nodal coordinates in order

C Begin new output. Divides each fe element locally NT_SUB_DIV times
C (NT_SUB_DIV is entered in DEEXPO).
C Note: The locally refined .channels file is only valid for the torso
C surface - the heart .channels file currently contains incorrect 
C information (but this file should not need to be used).

          CALL ASSERT(NT_SUB_DIV.LE.8,
     '      ' >> LOCAL_NP array too small',ERROR,*9999)
          !Array needs to be at least (NT_SUB_DIV+2)^2 
          CALL OPENF(IOFILE2,'DISK',FILE(IBEG:IEND)//'.pts','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          CALL OPENF(IOFILE3,'DISK',FILE(IBEG:IEND)//'.fac','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          CALL OPENF(IOFILE4,'DISK',FILE(IBEG:IEND)//'.channels',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)

          !Firstly write out all nodal coordinates in order to .pts
          !file and renumbered nodes to .channels file
          !AJP 14-5-94.  The first number in the .channels files is
          !the number of nodes in the map3d mesh.  If local subdivision
          !is used, then this number is not know apriori.  Thus
          !need to store all map3d nodes, and write them out after
          !they have all been calculated. This is stored in ICHANNELS
          !array (ICHANNELS(i,1)=map3d node number; ICHANNELS(i,2) is
          !the cmiss node number.
          ICHT=0 !Total number of channels
          IF(TYPE(1:6).EQ.'NUMBER') THEN
c            WRITE(IOFILE4,'(I3,'' Nodes and NT_SUB_DIV = '',I3)')
c     '        NPLIST(0),NT_SUB_DIV
            DO NONODE=1,NPLIST(0)
              WRITE(IOFILE2,*) (XP(1,nv,nj,NPLIST(NONODE)),NJ=1,NJT)
c              WRITE(IOFILE4,*) NONODE,NPLIST(NONODE)
              ICHT=ICHT+1
              ICHANNELS(ICHT,1)=NONODE
              ICHANNELS(ICHT,2)=NPLIST(NONODE)
            ENDDO
          ELSE IF(TYPE(1:5).EQ.'GROUP') THEN
            CALL TRIM(GROUP_NAME,IBEG1,IEND1)
            NOGRNO=0
            DO NOGROUP=1,NTGRNO
              CALL TRIM(LAGRNO(NOGROUP),IBEG2,IEND2)
              IF(GROUP_NAME(IBEG1:IEND1).EQ.
     '          LAGRNO(NOGROUP)(IBEG2:IEND2)) NOGRNO=NOGROUP
            ENDDO
            IF(NOGRNO.EQ.0) THEN
              ERROR='>>Group not found'
              GOTO 9999
            ENDIF
c            WRITE(IOFILE4,'(I3,'' Nodes (and local subdivision'')')
c     '         LIGRNO(0,NOGRNO)
            DO NONODE=1,LIGRNO(0,NOGRNO)
              WRITE(IOFILE2,*) (XP(1,nv,nj,LIGRNO(NONODE,NOGRNO)),
     '          NJ=1,NJT)
c              WRITE(IOFILE4,*) NONODE,LIGRNO(NONODE,NOGRNO)
              ICHT=ICHT+1
              ICHANNELS(ICHT,1)=NONODE
              ICHANNELS(ICHT,2)=LIGRNO(NONODE,NOGRNO)
            ENDDO
          ENDIF

          !Loop over each element and carry out local refinement.
          N_CHANNEL=NPLIST(0)
          N_NODE=NPLIST(NPLIST(0))
          LOCAL_NN=NPLIST(0) !Refined global node number start
          !Note: If NPLIST(0) < NPT(0) then the files can contain 
          !incorrect information.  However, one can't use NPT(0) 
          !above since files won't contain the correct information 
          !for Map3d.  The files for the torso surface will be correct
          !but only the geometry files for the heart surface will be
          !correct (and not the channels file).
          IF(TYPE(1:6).EQ.'NUMBER') THEN
            DO NOELEM=1,NELIST(0)
              NE=NELIST(NOELEM)
              CALL XPXE(NBJ(1,NE),NKE(1,1,1,NE),NPF(1,1),
     '          NPNE(1,1,NE),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,NE),XA,XE,XP,ERROR,*9999)
              nb=NBJ(1,NE)
              !Find the corresponding vertex number for the local 
              !nodes from NPLIST
              DO nn=1,NNT(nb)
                NP=NPNE(nn,nb,NE)
                I=1
                DO WHILE(NP.NE.NPLIST(I).AND.I.LE.NPLIST(0))
                  I=I+1
                ENDDO
                LOCAL_NP(nn)=I
              ENDDO

              CALL EXGEOM_NE(IBT,ICHANNELS,ICHT,IDO,INP,LOCAL_NN,
     '          LOCAL_NP,N_CHANNEL,N_NODE,nb,NBJ,ne,XE,ERROR,*9999)

            ENDDO !Element loop

          ELSE IF(TYPE(1:5).EQ.'GROUP') THEN

            CALL TRIM(GROUP_NAME,IBEG1,IEND1)
            NOGREL=0
            DO ELGROUP=1,NTGREL
              CALL TRIM(LAGREL(ELGROUP),IBEG2,IEND2)
              IF(GROUP_NAME(IBEG1:IEND1).EQ.
     '          LAGREL(ELGROUP)(IBEG2:IEND2)) NOGREL=ELGROUP
            ENDDO
            IF(NOGREL.EQ.0) THEN
              ERROR='>>Group not found'
              GOTO 9999
            ENDIF

            DO NOELEM=1,LIGREL(0,NOGREL)
              NE=LIGREL(NOELEM,NOGREL)
              CALL XPXE(NBJ(1,NE),NKE(1,1,1,NE),NPF(1,1),
     '          NPNE(1,1,NE),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,NE),XA,XE,XP,ERROR,*9999)
              nb=NBJ(1,NE)
              !Find the corresponding vertex number for the local 
              !nodes from LIGRNO
              DO nn=1,NNT(nb)
                NP=NPNE(nn,nb,NE)
                I=1
                DO WHILE(NP.NE.LIGRNO(I,NOGRNO).AND.
     '            I.LE.LIGRNO(0,NOGRNO))
                  I=I+1
                ENDDO
                LOCAL_NP(nn)=I
              ENDDO

              CALL EXGEOM_NE(IBT,ICHANNELS,ICHT,IDO,INP,LOCAL_NN,
     '          LOCAL_NP,N_CHANNEL,N_NODE,nb,NBJ,ne,XE,ERROR,*9999)
              

            ENDDO !Element loop
          ENDIF
C Write out the .channels file.
          CALL ASSERT(ICHT.LE.1000,
     '      ' >> ICHANNELS array too small',ERROR,*9999)
          !Array needs to be at least (ICHT,2) 
          WRITE(IOFILE4,'(I3,'' Nodes'')')
     '        ICHT
          !First line of .channels file can only contain the number of
          !nodes only
          DO ICH=1,ICHT 
            WRITE(IOFILE4,*) ICHANNELS(ICH,1),ICHANNELS(ICH,2)
          ENDDO
          CALL CLOSEF(IOFILE2,ERROR,*9999)
          CALL CLOSEF(IOFILE3,ERROR,*9999)
          CALL CLOSEF(IOFILE4,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('EXGEOM')
      RETURN
 9999 CALL ERRORS('EXGEOM',ERROR)
      CALL EXITS('EXGEOM')
      RETURN 1
      END


      SUBROUTINE EXGEOM_NE(IBT,ICHANNELS,ICHT,IDO,INP,LOCAL_NN,
     '  LOCAL_NP,N_CHANNEL,N_NODE,nb,NBJ,ne,XE,ERROR,*)

C#### Subroutine: EXGEOM_NE
C###  Description:
C###    EXGEOM_NE exports geometry data for element ne.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:map3d.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ICHANNELS(1000,*),ICHT,
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LOCAL_NN,LOCAL_NP(*),N_CHANNEL,N_NODE,nb,NBJ(NJM,NEM),ne
      REAL*8 XE(NSM,NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER DUMMY(3),I,IJ_END,ISTEP,J,NB2,nj,nn
      REAL*8 PXI,X(3),XI(2)

      CALL ENTERS('EXGEOM_NE',*9999)

      IF(IBT(1,1,NB).EQ.2.AND.IBT(1,2,NB).EQ.2) THEN ! Bicubic Hermite
        DUMMY(1)=LOCAL_NP(2)
        DUMMY(2)=LOCAL_NP(3)
        DUMMY(3)=LOCAL_NP(4)
        IJ_END=NT_SUB_DIV+2
        NN=1
        DO J=1,IJ_END
          DO I=1,IJ_END
            IF(.NOT.(I.EQ.1.AND.J.EQ.1))THEN
              NN=NN+1
              IF(I.EQ.IJ_END.AND.J.EQ.1)THEN
                LOCAL_NP(nn)=DUMMY(1)
              ELSEIF(I.EQ.1.AND.J.EQ.IJ_END)THEN
                LOCAL_NP(nn)=DUMMY(2)
              ELSEIF(I.EQ.IJ_END.AND.J.EQ.IJ_END)THEN
                LOCAL_NP(nn)=DUMMY(3)
              ELSE
                LOCAL_NN=LOCAL_NN+1
                LOCAL_NP(nn)=LOCAL_NN
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(DOP) THEN
          IJ_END=IJ_END*IJ_END
          WRITE(OP_STRING,*)' NE=',NE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,*)
     '      ' LOCAL_NP(nn)=',(LOCAL_NP(nn),NN=1,IJ_END)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        !Write out .fac information for this element
        ISTEP=NT_SUB_DIV+2
        DO J=0,NT_SUB_DIV
          DO I=0,NT_SUB_DIV
            WRITE(UNIT=IOFILE3,FMT=*) LOCAL_NP(J+I*ISTEP+1),
     '        LOCAL_NP(J+I*ISTEP+2),
     '        LOCAL_NP(J+NT_SUB_DIV+4+I*ISTEP)
            WRITE(UNIT=IOFILE3,FMT=*) LOCAL_NP(J+I*ISTEP+1),
     '        LOCAL_NP(J+NT_SUB_DIV+4+I*ISTEP),
     '        LOCAL_NP(J+NT_SUB_DIV+3+I*ISTEP)
            IF(DOP)THEN
              WRITE(OP_STRING,*) LOCAL_NP(J+I*ISTEP+1),
     '          LOCAL_NP(J+I*ISTEP+2),
     '          LOCAL_NP(J+NT_SUB_DIV+4+I*ISTEP)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,*) LOCAL_NP(J+I*ISTEP+1),
     '          LOCAL_NP(J+NT_SUB_DIV+4+I*ISTEP),
     '          LOCAL_NP(J+NT_SUB_DIV+3+I*ISTEP)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
        !Write out .pts and .channels information for this element
        DO J=0,NT_SUB_DIV+1
          DO I=0,NT_SUB_DIV+1
            IF(.NOT.((I.EQ.0.AND.J.EQ.0).OR.
     '         (I.EQ.0.AND.J.EQ.NT_SUB_DIV+1).OR.
     '         (I.EQ.NT_SUB_DIV+1.AND.J.EQ.0).OR.
     '         (I.EQ.NT_SUB_DIV+1.AND.J.EQ.NT_SUB_DIV+1)))THEN
              N_CHANNEL=N_CHANNEL+1
              N_NODE=N_NODE+1
              XI(1)=DBLE(I)/DBLE(NT_SUB_DIV+1)
              XI(2)=DBLE(J)/DBLE(NT_SUB_DIV+1)
              DO nj=1,NJT
                NB2=NBJ(nj,NE)
                X(nj)=PXI(IBT(1,1,NB2),IDO(1,1,0,NB2),
     '              INP(1,1,NB2),NB2,1,XI,XE(1,NJ))
              ENDDO
              WRITE(IOFILE2,*)(X(nj),NJ=1,NJT)
C              WRITE(IOFILE4,*)N_CHANNEL,N_NODE
              ICHT=ICHT+1
              ICHANNELS(ICHT,1)=N_CHANNEL
              ICHANNELS(ICHT,2)=N_NODE
              IF(DOP)THEN
                WRITE(OP_STRING,*)(X(nj),NJ=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,*)N_CHANNEL,N_NODE
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ELSE IF(IBT(1,2,NB).EQ.3.AND.IBT(1,2,NB).EQ.3) THEN 
                                                   ! Hermite-Simplex
        IF(NKT(1,NB).EQ.1) THEN ! Apex at node 1
          LOCAL_NP(NT_SUB_DIV+3)=LOCAL_NP(3)
          DO I=3,NT_SUB_DIV+2
            LOCAL_NN=LOCAL_NN+1
            LOCAL_NP(I)=LOCAL_NN
          ENDDO
          IF(DOP) THEN
            IJ_END=NT_SUB_DIV+3
            WRITE(OP_STRING,*)' NE=',NE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)
     '        ' LOCAL_NP(nn)=',(LOCAL_NP(nn),NN=1,IJ_END)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          !Write out .fac information for this element
          DO I=2,NT_SUB_DIV+2
            WRITE(UNIT=IOFILE3,FMT=*) LOCAL_NP(1),
     '          LOCAL_NP(I+1),LOCAL_NP(I)
            IF(DOP)THEN
              WRITE(OP_STRING,FMT=*) LOCAL_NP(1),
     '            LOCAL_NP(I+1),LOCAL_NP(I)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
          !Write out .pts and .channels information for this element
          DO I=1,NT_SUB_DIV
            N_CHANNEL=N_CHANNEL+1
            N_NODE=N_NODE+1
            XI(1)=DBLE(I)/DBLE(NT_SUB_DIV+1)
            XI(2)=1.0d0
            DO nj=1,NJT
              NB2=NBJ(nj,NE)
              X(nj)=PXI(IBT(1,1,NB2),IDO(1,1,0,NB2),
     '              INP(1,1,NB2),NB2,1,XI,XE(1,NJ))
            ENDDO
            WRITE(IOFILE2,*)(X(nj),NJ=1,NJT)
C            WRITE(IOFILE4,*)N_CHANNEL,N_NODE
            ICHT=ICHT+1
            ICHANNELS(ICHT,1)=N_CHANNEL
            ICHANNELS(ICHT,2)=N_NODE
            IF(DOP)THEN
              WRITE(OP_STRING,*)(X(nj),NJ=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,*)N_CHANNEL,N_NODE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ELSE IF(NKT(3,NB).EQ.1) THEN ! Apex at node 3
          LOCAL_NP(NT_SUB_DIV+2)=LOCAL_NP(2)
          LOCAL_NP(NT_SUB_DIV+3)=LOCAL_NP(3)
          DO I=2,NT_SUB_DIV+1
            LOCAL_NN=LOCAL_NN+1
            LOCAL_NP(I)=LOCAL_NN
          ENDDO
          IF(DOP) THEN
            IJ_END=NT_SUB_DIV+3
            WRITE(OP_STRING,*)' NE=',NE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)
     '        ' LOCAL_NP(nn)=',(LOCAL_NP(nn),NN=1,IJ_END)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          !Write out .fac information for this element
          DO I=1,NT_SUB_DIV+1
            WRITE(UNIT=IOFILE3,FMT=*) LOCAL_NP(I),
     '         LOCAL_NP(I+1),LOCAL_NP(NT_SUB_DIV+3)
            IF(DOP)THEN
              WRITE(OP_STRING,FMT=*) LOCAL_NP(I),
     '           LOCAL_NP(I+1),LOCAL_NP(NT_SUB_DIV+3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
          !Write out .pts and .channels information for this element
          DO I=1,NT_SUB_DIV
            N_CHANNEL=N_CHANNEL+1
            N_NODE=N_NODE+1
            XI(1)=DBLE(I)/DBLE(NT_SUB_DIV+1)
            XI(2)=0
            DO nj=1,NJT
              NB2=NBJ(nj,NE)
              X(nj)=PXI(IBT(1,1,NB2),IDO(1,1,0,NB2),
     '          INP(1,1,NB2),NB2,1,XI,XE(1,NJ))
            ENDDO
            WRITE(IOFILE2,*)(X(nj),NJ=1,NJT)
C            WRITE(IOFILE4,*)N_CHANNEL,N_NODE
            ICHT=ICHT+1
            ICHANNELS(ICHT,1)=N_CHANNEL
            ICHANNELS(ICHT,2)=N_NODE
            IF(DOP)THEN
              WRITE(OP_STRING,*)(X(nj),NJ=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,*)N_CHANNEL,N_NODE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDIF !Apex choice
      ENDIF !Element choice

      CALL EXITS('EXGEOM_NE')
      RETURN
 9999 CALL ERRORS('EXGEOM_NE',ERROR)
      CALL EXITS('EXGEOM_NE')
      RETURN 1
      END




      SUBROUTINE DIBASE(ISEG,CSEG,STRING,ERROR,*)

C#### Subroutine: DIBASE
C###  Description:
C###    DIBASE displays basis functions.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER I,IBEG,IBEG1,IDATA(2),ID_WS,IEND,
     '  IEND1,IFROMC,INDEX,INDEX_OLD,INSTAT,ISEGM,
     '  ISG1,ISG2,ISG3,ISG4,ISG5,ISG6,ISG7,ISG8,IPICK,
     '  N3CO,nb,NOCH
C     INTEGER ID_DEVICE,ID_STATUS,INS
      REAL R4DATA(2),VALUE,VALUE1,VALUE2,XWC
      REAL*8 PL(2,3),XI
      CHARACTER CHAR*6,CLASS*8
c     CHARACTER OPTION(50)*20,SDATA*10 
      LOGICAL CBBREV,CONTINUE,LBASE10,LBASE11,LBASE20,LBASE21,
     '  LBASE30

      CALL ENTERS('DIBASE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM display bases 
C###  Parameter:       <number BASIS[1]>
C###  Description:
C###    Display the basis function shape - 1D only at present.

        OP_STRING(1)=STRING(1:IEND)//' <number BASIS[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DIBASE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICSM.EQ.1,
     '    '>>Set USE_GRAPHICSM=1',ERROR,*9999)
        IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
          nb=IFROMC(CO(N3CO+1))
        ELSE
          nb=1
        ENDIF

        CALL ACWK(22,0,ERROR,*9999)
        IF(NIT(nb).EQ.1) THEN
          LBASE10=.FALSE.
          LBASE11=.FALSE.
          LBASE20=.FALSE.
          LBASE21=.FALSE.
          LBASE30=.FALSE.

C ***     Draw axes and title in segment ISG1
          CALL OPEN_SEGMENT(ISG1,ISEG,22,'base',INDEX,INDEX_OLD,
     '      1,1,CSEG,ERROR,*9999)
          PL(1,1)=0.0D0
          PL(2,1)=0.0D0
          PL(1,2)=1.0D0
          PL(2,2)=0.0D0
          CALL POLYLINE(1,22,2,PL,ERROR,*9999)
          PL(1,1)=0.0D0
          PL(2,1)=0.0D0
          PL(1,2)=0.0D0
          PL(2,2)=1.0D0
          CALL POLYLINE(1,22,2,PL,ERROR,*9999)
          DO I=1,3
            PL(1,1)=0.0D0
            PL(1,2)=0.02D0
            PL(2,1)=DBLE(I-1)/2.0D0
            PL(2,2)=PL(2,1)
            CALL POLYLINE(1,22,2,PL,ERROR,*9999)
            WRITE(CHAR,'(F6.1)') PL(2,1)
            CALL TRIM(CHAR,IBEG,IEND)
            PL(1,1)=-0.01D0
            CALL TEXT(1,22,CHAR(IBEG:IEND),PL(1,1),ERROR,*9999)
          ENDDO
          CALL ELLIPSE(0.007D0,0.02D0,0.0D0,0.0D0,10,ERROR,*9999)
          CALL ELLIPSE(0.007D0,0.02D0,1.0D0,0.0D0,10,ERROR,*9999)
          PL(1,1)=0.5D0
          PL(2,1)=1.0D0
          CALL TEXT(1,22,'Basis functions',PL(1,1),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISG1,22,ERROR,*9999)

C ***     Draw circles at nodes in segments ISG2 and ISG3
          CALL OPEN_SEGMENT(ISG2,ISEG,22,'base',INDEX,INDEX_OLD,
     '      2,1,CSEG,ERROR,*9999)
          CALL ELLIPSE(0.007D0,0.02D0,0.3334D0,0.0D0,10,ERROR,*9999)
          PL(1,1)=0.3334D0
          PL(2,1)=-0.5D0
          CALL TEXT(1,22,'Node 1',PL(1,1),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISG2,22,ERROR,*9999)

          CALL OPEN_SEGMENT(ISG3,ISEG,22,'base',INDEX,INDEX_OLD,
     '      3,1,CSEG,ERROR,*9999)
          CALL ELLIPSE(0.007D0,0.02D0,0.6667D0,0.0D0,10,ERROR,*9999)
          PL(1,1)=0.6667D0
          PL(2,1)=-0.05D0
          CALL TEXT(1,22,'Node 2',PL(1,1),ERROR,*9999)
          CALL CLOSE_SEGMENT(ISG3,22,ERROR,*9999)

C ***     Choose option
c         OPTION(1)='Set nodal params'
c         OPTION(2)='Zero nodal params'
c         OPTION(3)='Locate position'
c         OPTION(4)='Return'
c         CALL CHOICE('DIBASE',1,1,INS,72,'EVENT',4,4,NOCH,noco,1,
c    '      CO,OPTION,STRING,0.5*XDISP,YDISP-0.562*XDISP,ERROR,*9999)
          CONTINUE=.TRUE.
          DO WHILE (CONTINUE)
c           CALL EVENT(ID_WS,ID_DEVICE,ID_STATUS,CLASS,IDATA,R4DATA,
c    '        SDATA,ERROR,*9999)
            IF((CLASS(1:6).EQ.'CHOICE'.AND.ID_WS.EQ.72).OR
     '        .(CLASS(1:8).EQ.'VALUATOR'.AND.ID_WS.EQ.81))THEN
              NOCH=IDATA(1)
              IF(NOCH.EQ.1) THEN
                CALL DETECT(22,ISEG,ISG2,'DETECTABLE',ERROR,*9999)
                CALL DETECT(22,ISEG,ISG3,'DETECTABLE',ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' Pick nodes to include in interpolation '
     '            //'(use 2nd button to exit)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                INSTAT=1
                DO WHILE(INSTAT.EQ.1)
                  CALL PICK(22,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,
     '              *9999)
                  IF(INSTAT.EQ.1) THEN
                    IF(ISEGM.EQ.ISG2) THEN
                      WRITE(OP_STRING,*) 
     '                  ' Node 1 included in interpolation:'
     '                  //' set value'
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                     CALL VALUATOR('1',81,'REQUEST',5,0.,1.,0.,VALUE1,
c    '                  0.3334,-0.1,ERROR,*9999)

C ***                 Draw parameter height for node 1
                      IF(LBASE10) CALL DELETE_SEGMENT(ISG4,ISEG,22,
     '                  ERROR,*9999)

                      CALL OPEN_SEGMENT(ISG4,ISEG,22,'base',INDEX,
     '                  INDEX_OLD,4,1,CSEG,ERROR,*9999)
                      LBASE10=.TRUE.
                      PL(1,1)=0.3334D0
                      PL(1,2)=PL(1,1)
                      PL(2,1)=0.0D0
                      PL(2,2)=DBLE(VALUE1)
                      CALL POLYLINE(4,22,2,PL,ERROR,*9999)
                      WRITE(CHAR,'(F5.2)') DBLE(VALUE1)
                      PL(1,2)=-0.20d0
                      CALL TEXT(1,22,'u1='//CHAR(1:5),PL(1,1),ERROR,
     '                  *9999)
                      CALL CLOSE_SEGMENT(ISG4,22,ERROR,*9999)

C ***                 Draw basis function for node 1
                      IF(LBASE11) CALL DELETE_SEGMENT(ISG6,ISEG,22,
     '                  ERROR,*9999)

                      CALL OPEN_SEGMENT(ISG6,ISEG,22,'base',INDEX,
     '                  INDEX_OLD,6,1,CSEG,ERROR,*9999)
                      LBASE11=.TRUE.
                      PL(1,1)=0.0D0
                      PL(1,2)=0.3334D0
                      PL(1,3)=0.6667D0
                      PL(2,1)=0.0D0
                      PL(2,2)=DBLE(VALUE1)
                      PL(2,3)=0.0D0
                      CALL POLYLINE(3,22,2,PL,ERROR,*9999)
                      CALL CLOSE_SEGMENT(ISG6,22,ERROR,*9999)

                    ELSE IF(ISEGM.EQ.ISG3) THEN
                      WRITE(OP_STRING,*) 
     '                  ' Node 2 included in interpolation:'
     '                  //' set value'
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c                     CALL VALUATOR('2',81,'REQUEST',5,0.,1.,0.,VALUE2,
c    '                  0.6667,-0.1,ERROR,*9999)

C ***                 Draw parameter height for node 2
                      IF(LBASE20) CALL DELETE_SEGMENT(ISG5,ISEG,22,
     '                  ERROR,*9999)

                      CALL OPEN_SEGMENT(ISG5,ISEG,22,'base',INDEX,
     '                  INDEX_OLD,5,1,CSEG,ERROR,*9999)
                      LBASE20=.TRUE.
                      PL(1,1)=0.6667D0
                      PL(1,2)=PL(1,1)
                      PL(2,1)=0.0D0
                      PL(2,2)=DBLE(VALUE2)
                      CALL POLYLINE(4,22,2,PL,ERROR,*9999)
                      WRITE(CHAR,'(F5.2)') DBLE(VALUE1)
                      PL(1,2)=-0.20D0
                      CALL TEXT(1,22,'u2='//CHAR(1:5),PL(1,1),ERROR,
     '                  *9999)
                      CALL CLOSE_SEGMENT(ISG5,22,ERROR,*9999)

C ***                 Draw basis function for node 2
                      IF(LBASE21) CALL DELETE_SEGMENT(ISG7,ISEG,22,
     '                  ERROR,*9999)

                      CALL OPEN_SEGMENT(ISG7,ISEG,22,'base',INDEX,
     '                  INDEX_OLD,7,1,CSEG,ERROR,*9999)
                      LBASE21=.TRUE.
                      PL(1,1)=0.3334D0
                      PL(1,2)=0.6667D0
                      PL(1,3)=1.0D0
                      PL(2,1)=0.0D0
                      PL(2,2)=DBLE(VALUE2)
                      PL(2,3)=0.0D0
                      CALL POLYLINE(3,22,2,PL,ERROR,*9999)
                      CALL CLOSE_SEGMENT(ISG7,22,ERROR,*9999)

                    ENDIF
                  ENDIF
                ENDDO
              ELSE IF(NOCH.EQ.2) THEN
                WRITE(OP_STRING,
     '            '('' Pick nodes to exclude from interpolation '
     '            //'(use 2nd button to exit)'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                INSTAT=1
                DO WHILE(INSTAT.EQ.1)
                  CALL PICK(22,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,
     '              *9999)
                  IF(INSTAT.EQ.1) THEN
                    IF(ISEGM.EQ.ISG2) THEN
                      WRITE(OP_STRING,*) 
     '                  ' Node 1 excluded from interpolation'
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      IF(LBASE10) THEN
                        CALL DELETE_SEGMENT(ISG4,ISEG,22,ERROR,*9999)
                        LBASE10=.FALSE.
                      ENDIF
                      IF(LBASE11) THEN
                        CALL DELETE_SEGMENT(ISG6,ISEG,22,ERROR,*9999)
                        LBASE11=.FALSE.
                      ENDIF
                    ELSE IF(ISEGM.EQ.ISG3) THEN
                      WRITE(OP_STRING,*) 
     '                  ' Node 2 excluded from interpolation'
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                      IF(LBASE20) THEN
                        CALL DELETE_SEGMENT(ISG5,ISEG,22,ERROR,*9999)
                        LBASE20=.FALSE.
                      ENDIF
                      IF(LBASE21) THEN
                        CALL DELETE_SEGMENT(ISG7,ISEG,22,ERROR,*9999)
                        LBASE21=.FALSE.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ELSE IF(NOCH.EQ.3) THEN
                WRITE(OP_STRING,*)
     '            ' Set position to evaluate function'
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C               CALL LOCATOR(22,INSTAT,'REQUEST',4,0.0D0,XWC,0.0DD0,
C     '            YWC,ERROR,*9999)
                XWC=REAL(R4DATA(1))
                XI=3.0D0*(DBLE(XWC)-0.3334D0)
                IF(XI.LT.0.0D0) THEN
                  XI=0.0D0
                ELSE IF(XI.GT.1.0D0) THEN
                  XI=1.0D0
                ENDIF
                WRITE(OP_STRING,'('' Xi-position= '',F5.2)') XI
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                VALUE=REAL((1.0D0-XI)*DBLE(VALUE1)+XI*DBLE(VALUE2))

C ***           Draw height for interpolated function
                IF(LBASE30) CALL DELETE_SEGMENT(ISG8,ISEG,22,ERROR,
     '            *9999)
                CALL OPEN_SEGMENT(ISG8,ISEG,22,'base',INDEX,INDEX_OLD,
     '            8,1,CSEG,ERROR,*9999)
                LBASE30=.TRUE.
                PL(1,1)=DBLE(XWC)
                PL(1,2)=PL(1,1)
                PL(2,1)=0.0D0
                PL(2,2)=DBLE(VALUE)
                CALL POLYLINE(1,22,2,PL,ERROR,*9999)
                WRITE(CHAR,'(F5.2)') XI
                PL(1,2)=-0.45D0
                CALL TEXT(1,22,'Xi ='//CHAR(1:5),PL(1,1),ERROR,*9999)
                WRITE(CHAR,'(F5.2)') DBLE(VALUE)
                PL(1,2)=-0.60D0
                CALL TEXT(1,22,'u  ='//CHAR(1:5),PL(1,1),ERROR,*9999)
                CALL CLOSE_SEGMENT(ISG8,22,ERROR,*9999)

C ***           Redraw parameter height for node 1
                IF(LBASE10) CALL DELETE_SEGMENT(ISG4,ISEG,22,ERROR,
     '            *9999)
                CALL OPEN_SEGMENT(ISG4,ISEG,22,'base',INDEX,INDEX_OLD,
     '            4,1,CSEG,ERROR,*9999)
                LBASE10=.TRUE.
                PL(1,1)=0.3334D0
                PL(1,2)=PL(1,1)
                PL(2,1)=0.0D0
                PL(2,2)=DBLE(VALUE1)*(1.0D0-XI)
                CALL POLYLINE(4,22,2,PL,ERROR,*9999)
                CALL CLOSE_SEGMENT(ISG4,22,ERROR,*9999)

C ***           Redraw parameter height for node 2
                IF(LBASE20) CALL DELETE_SEGMENT(ISG5,ISEG,22,ERROR,
     '            *9999)
                CALL OPEN_SEGMENT(ISG5,ISEG,22,'base',INDEX,INDEX_OLD,
     '            5,1,CSEG,ERROR,*9999)
                LBASE20=.TRUE.
                PL(1,1)=0.6667D0
                PL(1,2)=PL(1,1)
                PL(2,1)=0.0D0
                PL(2,2)=DBLE(VALUE2)*XI
                CALL POLYLINE(4,22,2,PL,ERROR,*9999)
                CALL CLOSE_SEGMENT(ISG5,22,ERROR,*9999)

              ELSE IF(NOCH.EQ.4) THEN
c               CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,*9999)
                CONTINUE=.FALSE.
              ENDIF
            ENDIF
          ENDDO
        ELSE IF(NIT(nb).EQ.2) THEN
          WRITE(OP_STRING,
     '      '('' 2D basis functions not yet displayable'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ELSE IF(NIT(nb).EQ.3) THEN
          WRITE(OP_STRING,
     '      '('' 3D basis functions not yet displayable'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        CALL DAWK(22,0,ERROR,*9999)
      ENDIF

      CALL EXITS('DIBASE')
      RETURN
 9999 CALL ERRORS('DIBASE',ERROR)
      CALL EXITS('DIBASE')
      RETURN 1
      END


      SUBROUTINE DITRAC(ISEG,ISSIGN,ISTRAC,CSEG,STRING,ERROR,*)

C#### Subroutine: DITRAC
C###  Description:
C###    DITRAC displays trace.  If sockets are being used, then this 
C###    subroutine simply transfers the data to the Motif front-end.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cmis00.cmn'
      INCLUDE 'cmiss$reference:ditr00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:pics00.cmn'
      INCLUDE 'cmiss$reference:read00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSIGN(2,128),ISTRAC(10)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER I,IBEG,IBEG1,IBEG2,ICOL,IDATA(2),ID_WS,IEND,
     '  IEND1,IEND2,IFROMC,INDEX(20),INDEX_OLD,INDEX_POLYLINE,
     '  INDEX_POLYMARKER,INDEX_TEXT,INS,INST,INSTAT,IPICK,ISAMPLES,
     '  ISEG0,ISEG1,ISEG2,ISEG3,ISEG4,ISEG5,
     '  ISEGM,ISIG,IWS,IWT,N3CO,NAV,NOCH,NOCH1,NOCH2
c     INTEGER ID_DEVICE
c     REAL R4DATA(2)
      REAL*8 PL(3,2),POSN,RFROMC,TDIFF,TEMP,TSEND,TSSTART,TTBAND,
     '  XPOS,XWC,YPOS,YWC 
      CHARACTER ACCREJ*6,CLASS*10,C1*10,C2*10,
     '  OPTION(20)*80,OPTION2(5)*80
      LOGICAL CBBREV,CHOOSE,CONTINUE,REGEN,REGENALL,REMARK

      CALL ENTERS('DITRAC',*9999)
      IF(.NOT.FIDUCALC) THEN
        ISIGNAL=1
        ISIGLAST=0
        TTSTART=0.0D0
        TTEND=0.3D0
      ENDIF
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE02,IBEG1,IEND1)
        WRITE(C1,'(I3)') ISIGNAL
        CALL TRIM(C1,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM display trace 
C###  Parameter:       <number SIGNAL_NUMBER[1]>
C###  Parameter:       <from START_TIME[0.0]>
C###  Parameter:       <to END_TIME[1.0]>
C###  Parameter:       <average #[6]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)
     '    //' <number SIGNAL_NUMBER['//C1(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<from START_TIME[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<to END_TIME[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<average #[6]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM display trace pick
C###  Parameter:       <from START_TIME[0.0]>
C###  Parameter:       <to END_TIME[1.0]>
C###  Parameter:       <AVERAGE #[6]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND) //' pick'
        OP_STRING(2)=BLANK(1:15)//'<from START_TIME[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<to END_TIME[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<average #[6]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DITRAC',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICSM.EQ.1,
     '    '>>Set USE_GRAPHICSM=1',ERROR,*9999)
        IF(USE_SOCKET) THEN    !Socket interface to front-end is being used
          CALL WRITE_SIGNAL(I2P,NUMSIGNALS,NUMSAMPLES,
     '      NUMSAMPLESPERSECOND,ERROR,*9999)
        ELSE     !Socket interface not being used
          CHOOSE=.FALSE.
          IF(CBBREV(CO,'PICK',1,noco+1,NTCO,N3CO)) THEN
            CHOOSE=.TRUE.
          ENDIF

          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            TSSTART=RFROMC(CO(N3CO+1))
          ELSE
            TSSTART=0.0D0
          ENDIF

          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            TSEND=RFROMC(CO(N3CO+1))
          ELSE
            TSEND=1.0D0
          ENDIF

          IF(CBBREV(CO,'AVERAGE',1,noco+1,NTCO,N3CO)) THEN
            NAV=IFROMC(CO(N3CO+1))
          ELSE
            NAV=6
          ENDIF

          INDEX( 1)=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')   !Axes
          INDEX( 2)=INDEX_POLYLINE(0,'SOLID','WIDTH1','CYAN')   !Accepted signal
          INDEX( 3)=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')    !Rejected signal
          INDEX( 4)=INDEX_POLYLINE(0,'SOLID','WIDTH1','WHITE')  !Current signal
          INDEX( 5)=INDEX_TEXT(0,'WIDTH1','FONT1','YELLOW')     !Signal number
          INDEX( 6)=INDEX_POLYLINE(0,'DASHED','WIDTH1','BLACK') !Separator
          INDEX( 7)=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')      !Axes text
          INDEX( 8)=INDEX_TEXT(0,'WIDTH1','FONT1','WHITE')      !Pick,datum text
          INDEX( 9)=INDEX_POLYMARKER(0,'PLUS','SIZE1','LTBLUE') !Pick marker
          INDEX(10)=INDEX_POLYLINE(0,'SOLID','WIDTH1','LTBLUE') !Pick line
          INDEX(11)=INDEX_POLYLINE(0,'SOLID','WIDTH1','YELLOW') !Fiducial marker
          INDEX(12)=INDEX_POLYLINE(0,'SOLID','WIDTH1','WHITE')  !Datum line

          OPTION( 1)='Accept'
          OPTION( 2)='Reject'
          OPTION( 3)='Pick Signal'
          OPTION( 4)='Enter Signal'
          OPTION( 5)='Recalculate..'
          OPTION( 6)='Exit'
          ISAMPLES=2048

          IWS=34  !Signal window
          IWT=35  !Trace window

          IF(.NOT.FIDUCALC) THEN   !If signals not previously displayed

            WRITE(OP_STRING,'('' >>Calculating Fiducial Markers'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DATUM=1.0D0
            DO ISIG=1,128
              CALL CALC_FID(ISAMPLES,ISIG,NAV,FIDMARK(ISIG),TTSTART,
     '          TTEND,ACCEPT(ISIG),ERROR,*9999)
              IF(FIDMARK(ISIG).LT.DATUM.
     '          AND.ACCEPT(ISIG)) DATUM=FIDMARK(ISIG)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Fiducial marker '',I3,'' is '','
     '            //'F6.3)')ISIG,FIDMARK(ISIG)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            WRITE(OP_STRING,
     '        '('' >>>Datum calculated : '',F5.1,'' milliseconds.'')')
     '        DATUM*1000.0D0
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

            CALL ACWK(IWS,0,ERROR,*9999)
            DO ISIG=1,128
              IF(ISIG.EQ.ISIGNAL) THEN
                ICOL=4
              ELSE IF(ACCEPT(ISIG)) THEN
                ICOL=2
              ELSE !Reject
                ICOL=3
              ENDIF
              CALL SGSIGN(ICOL,INDEX,ISAMPLES,ISEG,ISIG,ISSIGN(1,ISIG),
     '          IWS,CSEG,SIGHEADER(1,ISIG),DATUM,FIDMARK(ISIG),TSSTART,
     '          TSEND,ACCEPT(ISIG),ERROR,*9999)
            ENDDO
            CALL DAWK(IWS,0,ERROR,*9999)
            FIDUCALC=.TRUE.
          ENDIF

          IF(CHOOSE) THEN ! Pick Signal with Mouse
            ISIGLAST=ISIGNAL
            CALL ACWK(IWS,1,ERROR,*9999)
            DO I=1,128
              CALL DETECT(IWS,ISEG,ISSIGN(2,I),'DETECTABLE',ERROR,*9999)
            ENDDO
            CALL DAWK(IWS,1,ERROR,*9999)
            CALL ACWK(IWS,0,ERROR,*9999)
            WRITE(OP_STRING,'('' >>Choose signal trace'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL PICK(IWS,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,*9999)
            ISIGNAL=IFROMC(CSEG(ISEGM)(53:57))
            CALL TRIM(SIGHEADER(1,ISIGNAL),IBEG,IEND)
            WRITE(OP_STRING,
     '        '('' >>>Signal #'',A5)') SIGHEADER(1,ISIGNAL)(IBEG:IEND)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            CALL DAWK(IWS,0,ERROR,*9999)
            CALL ACWK(IWS,1,ERROR,*9999)
            DO I=1,128
              CALL DETECT(IWS,ISEG,ISSIGN(2,I),'UNDETECTABLE',ERROR,
     '          *9999)
            ENDDO
            CALL DAWK(IWS,1,ERROR,*9999)
          ELSE
            IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
              CALL TRIM(CO(N3CO+1),IBEG2,IEND2)
              C2=CO(N3CO+1)(IBEG2:IEND2)
              CALL TRIM(C2,IBEG2,IEND2)
              ISIGNAL=1
              DO I=1,128
                C1=SIGHEADER(1,I)
                CALL TRIM(C1,IBEG1,IEND1)
                IF(C1(IBEG1:IEND1).EQ.C2(IBEG2:IEND2)) THEN
                  ISIGNAL=I
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          CALL ACWK(IWS,1,ERROR,*9999)

          IF(ISIGLAST.NE.ISIGNAL) THEN
            IF(ISIGLAST.GT.0) THEN
              IF(ACCEPT(ISIGLAST)) THEN
                ICOL=2
              ELSE
                ICOL=3
              ENDIF
              CALL SGSIGN(ICOL,INDEX,ISAMPLES,ISEG,ISIGLAST,
     '          ISSIGN(1,ISIGLAST),IWS,CSEG,SIGHEADER(1,ISIGLAST),
     '          DATUM,FIDMARK(ISIGLAST),TSSTART,TSEND,
     '          ACCEPT(ISIGLAST),ERROR,*9999)
                !Show last signal in correct colour
              ISIGLAST=ISIGNAL
            ENDIF
          ENDIF

          CALL SGSIGN(4,INDEX,ISAMPLES,ISEG,ISIGNAL,ISSIGN(1,ISIGNAL),
     '      IWS,CSEG,SIGHEADER(1,ISIGNAL),DATUM,FIDMARK(ISIGNAL),
     '      TSSTART,TSEND,ACCEPT(ISIGNAL),ERROR,*9999)
            !Show current signal in diff col
          CALL DAWK(IWS,1,ERROR,*9999)

          CALL ACWK(IWT,0,ERROR,*9999)  !Display trace window
          IF(ISIGNAL.GT.0) THEN
            IF(ACCEPT(ISIGNAL)) THEN
              ICOL=2
              ACCREJ='REJECT'
            ELSE
              ICOL=3
              ACCREJ='ACCEPT'
            ENDIF
            CALL SGTRAC(ICOL,INDEX,ISAMPLES,ISEG,ISIGNAL,ISTRAC,IWT,
     '        CSEG,SIGHEADER(1,ISIGNAL),DATUM,FIDMARK(ISIGNAL),TTSTART,
     '        TTEND,ERROR,*9999)
          ENDIF

          XPOS=0.90D0
          YPOS=8.05D0

          CALL PICK(IWT,'EVENT',INSTAT,ISEGM,IPICK,ERROR,*9999)

          CALL DETECT(IWT,ISEG,ISTRAC(2),'DETECTABLE',ERROR,*9999) ! Fid M.

          ISEG0=0
          CALL OPEN_SEGMENT(ISEG0,ISEG,IWT,'START_MARKER',1,INDEX_OLD,
     '      0,1,CSEG,ERROR,*9999)                ! Change TTSTART with marker
          PL(1,1)=XPOS+9.0D0*TTSTART
          PL(2,1)=YPOS
          CALL POLYMARKER(INDEX(9),IWT,1,PL,ERROR,*9999)
          CALL CLOSE_SEGMENT(ISEG0,IWT,ERROR,*9999)
          CALL DETECT(IWT,ISEG,ISEG0,'DETECTABLE',ERROR,*9999)

          ISEG1=0
          CALL OPEN_SEGMENT(ISEG1,ISEG,IWT,'END_MARKER',1,INDEX_OLD,
     '      1,1,CSEG,ERROR,*9999)                ! Change TTEND with marker
          PL(1,1)=XPOS+9.0D0*TTEND
          PL(2,1)=YPOS
          CALL POLYMARKER(INDEX(9),IWT,1,PL,ERROR,*9999)
          CALL CLOSE_SEGMENT(ISEG1,IWT,ERROR,*9999)
          CALL DETECT(IWT,ISEG,ISEG1,'DETECTABLE',ERROR,*9999)

          ISEG2=0
          CALL OPEN_SEGMENT(ISEG2,ISEG,IWT,'PREV_LABEL',1,INDEX_OLD,
     '      2,1,CSEG,ERROR,*9999)                ! Change to PREV with label
          PL(1,1)=2.0D0
          PL(2,1)=0.5D0
          CALL TEXT(INDEX(8),IWT,'PREV',PL,ERROR,*9999)
          CALL CLOSE_SEGMENT(ISEG2,IWT,ERROR,*9999)
          CALL DETECT(IWT,ISEG,ISEG2,'DETECTABLE',ERROR,*9999)

          ISEG3=0
          CALL OPEN_SEGMENT(ISEG3,ISEG,IWT,'NEXT_LABEL',1,INDEX_OLD,
     '      3,1,CSEG,ERROR,*9999)                ! Change to NEXT with label
          PL(1,1)=8.0D0
          PL(2,1)=0.5D0
          CALL TEXT(INDEX(8),IWT,'NEXT',PL,ERROR,*9999)
          CALL CLOSE_SEGMENT(ISEG3,IWT,ERROR,*9999)
          CALL DETECT(IWT,ISEG,ISEG3,'DETECTABLE',ERROR,*9999)

          ISEG4=0
          CALL OPEN_SEGMENT(ISEG4,ISEG,IWT,'ACCEPT/REJECT_LABEL',1,
     '      INDEX_OLD,4,1,CSEG,ERROR,*9999)                ! Accept/Reject with label
          PL(1,1)=8.0D0
          PL(2,1)=1.3D0
          CALL TEXT(INDEX(8),IWT,ACCREJ,PL,ERROR,*9999)
          CALL CLOSE_SEGMENT(ISEG4,IWT,ERROR,*9999)
          CALL DETECT(IWT,ISEG,ISEG4,'DETECTABLE',ERROR,*9999)

          ISEG5=0
          CALL OPEN_SEGMENT(ISEG5,ISEG,IWT,'BAND_MARKER',1,INDEX_OLD,
     '      5,1,CSEG,ERROR,*9999)                ! Move band with marker
          PL(1,1)=XPOS+9.0D0*(TTSTART+(TTEND-TTSTART)/2.0D0)
          PL(2,1)=YPOS+0.1D0
          PL(1,2)=XPOS+9.0D0*(TTSTART+(TTEND-TTSTART)/2.0D0)
          PL(2,2)=YPOS-0.1D0
          CALL POLYLINE(INDEX(10),IWT,2,PL,ERROR,*9999)
          CALL CLOSE_SEGMENT(ISEG5,IWT,ERROR,*9999)
          CALL DETECT(IWT,ISEG,ISEG5,'DETECTABLE',ERROR,*9999)

          CALL DAWK(IWT,0,ERROR,*9999)

          CALL ACWK(IWT,0,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' >>Modify window or choose from menu '')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          REGEN=.FALSE.
          REGENALL=.FALSE.
          REMARK=.FALSE.
          CONTINUE=.TRUE.
          DO WHILE(CONTINUE)
            ISIGLAST=ISIGNAL
c           CALL EVENT(ID_WS,ID_DEVICE,INSTAT,CLASS,IDATA,
c    '        R4DATA,SDATA,ERROR,*9999)
c           IF(DOP) THEN
c             WRITE(OP_STRING,*)
c    '          ' INPUT_CLASS=',CLASS,' INSTAT=',INSTAT
c             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c           ENDIF
            IF(CLASS(1:4).EQ.'PICK') THEN
              IF(INSTAT.EQ.1) THEN
                ISEGM=IDATA(1)
                IF(ISEGM.EQ.ISEG0.OR.ISEGM.EQ.ISEG1) THEN    !Start/end marker
                  IF(ISEGM.EQ.ISEG0) THEN
                    WRITE(OP_STRING,
     '                '('' >>Relocate time start marker'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ELSE
                    WRITE(OP_STRING,
     '                '('' >>Relocate time end marker'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                  CALL LOCATOR(INS,0.0D0,XWC,0.0D0,
     '              YWC,ERROR,*9999)
                  IF(INS.EQ.1) THEN
                    IF(XWC.LT.XPOS) THEN
                      XWC=XPOS
                    ELSE IF(XWC.GT.XPOS+9.0D0) THEN
                      XWC=XPOS+9.0D0
                    ENDIF
                    IF(ISEGM.EQ.ISEG0) THEN
                      TTSTART=(XWC-XPOS)/9.0D0
                    ELSE
                      TTEND=(XWC-XPOS)/9.0D0
                    ENDIF
                    IF(TTSTART.GT.TTEND) THEN
                      TEMP=TTEND
                      TTEND=TTSTART
                      TTSTART=TEMP
                    ENDIF
                    TDIFF=(TTEND-TTSTART)/2.0D0
                    TTBAND=TTSTART+TDIFF

                    REMARK=.TRUE.
                    REGEN=.TRUE.
                  ENDIF
                ELSE IF(ISEGM.EQ.ISEG2) THEN  !Previous signal
                  IF(ISIGNAL.EQ.1) THEN
                    WRITE(OP_STRING,'('' >>>Showing First Signal'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ELSE
                    ISIGNAL=ISIGNAL-1
                  ENDIF
                  REGEN=.TRUE.
                ELSE IF(ISEGM.EQ.ISEG3) THEN  !Next signal
                  IF(ISIGNAL.EQ.128) THEN
                    WRITE(OP_STRING,'('' >>>Showing Last Signal'')')
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ELSE
                    ISIGNAL=ISIGNAL+1
                  ENDIF
                  REGEN=.TRUE.
                ELSE IF(ISEGM.EQ.ISEG4) THEN  !Accept/reject signal
                  IF(ACCREJ(1:6).EQ.'ACCEPT') THEN
                    ACCEPT(ISIGNAL)=.TRUE.
                  ELSE
                    ACCEPT(ISIGNAL)=.FALSE.
                  ENDIF
                  REGEN=.TRUE.
                ELSE IF(ISEGM.EQ.ISEG5) THEN  !Move Band
                  WRITE(OP_STRING,'('' >>Relocate center of band'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(INS,0.0d0,XWC,0.0d0,
     '              YWC,ERROR,*9999)
                  IF(INS.EQ.1) THEN
                    IF(XWC.LT.XPOS) THEN
                      XWC=XPOS
                    ELSE IF(XWC.GT.XPOS+9.0D0) THEN
                      XWC=XPOS+9.0D0
                    ENDIF
                    TTBAND=(XWC-XPOS)/9.0D0
                    TDIFF=(TTEND-TTSTART)/2.0D0
                    IF((TTBAND-TDIFF).LT.0.0D0) TTBAND=TDIFF
                    IF((TTBAND-TDIFF).GT.1.0D0) TTBAND=1.0D0-TDIFF
                    TTSTART=TTBAND-TDIFF
                    TTEND=TTBAND+TDIFF
                    REGEN=.TRUE.
                    REMARK=.TRUE.
                  ENDIF
                ELSE IF(ISEGM.EQ.ISTRAC(2)) THEN  !Move Fiducial Marker
                  WRITE(OP_STRING,
     '              '('' >>Relocate Fiducial Marker'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL LOCATOR(INS,0.0D0,XWC,0.0D0,
     '              YWC,ERROR,*9999)
                  IF(INS.EQ.1) THEN
                    IF(XWC.LT.XPOS) THEN
                      XWC=XPOS
                    ELSE IF(XWC.GT.XPOS+9.0D0) THEN
                      XWC=XPOS+9.0D0
                    ENDIF
                    POSN=(XWC-XPOS)/9.0D0
                    FIDMARK(ISIGNAL)=POSN*(TTEND-TTSTART)+TTSTART
                  ENDIF
                  REGEN=.TRUE.
                ENDIF
              ELSE      !INSTAT.ne.1
              ENDIF
            ELSE IF(CLASS(1:6).EQ.'CHOICE') THEN
              IF(DOP) THEN
                WRITE(OP_STRING,*) 
     '            ' Input class is choice on ','id_ws=',ID_WS
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              NOCH=IDATA(1)
              IF(DOP) THEN
                WRITE(OP_STRING,*) ' Input_Choice=',NOCH
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              NOCH1=NOCH
              IF(OPTION(NOCH1)(1:6).EQ.'Accept') THEN
                ACCEPT(ISIGNAL)=.TRUE.
                REGEN=.TRUE.
              ELSE IF(OPTION(NOCH1)(1:6).EQ.'Reject') THEN
                ACCEPT(ISIGNAL)=.FALSE.
                REGEN=.TRUE.
              ELSE IF(OPTION(NOCH1)(1:13).EQ.'Pick Signal') THEN
                ! Rechoose signal from signal window (mouse)
                CALL DAWK(IWT,0,ERROR,*9999)
                CALL ACWK(IWS,1,ERROR,*9999)
                DO I=1,128
                  CALL DETECT(IWS,ISEG,ISSIGN(2,I),'DETECTABLE',ERROR,
     '              *9999)
                ENDDO
                CALL DAWK(IWS,1,ERROR,*9999)
                CALL ACWK(IWS,0,ERROR,*9999)
                WRITE(OP_STRING,'('' >>Choose signal trace'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL PICK(IWS,'REQUEST',INST,ISEGM,IPICK,ERROR,
     '            *9999)
                ISIGNAL=IFROMC(CSEG(ISEGM)(53:57))
                CALL TRIM(SIGHEADER(1,ISIGNAL),IBEG,IEND)
                WRITE(OP_STRING,'('' >>>Signal #'',A5)')
     '            SIGHEADER(1,ISIGNAL)(IBEG:IEND)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL DAWK(IWS,0,ERROR,*9999)
                CALL ACWK(IWS,1,ERROR,*9999)
                DO I=1,128
                  CALL DETECT(IWS,ISEG,ISSIGN(2,I),'UNDETECTABLE',
     '              ERROR,*9999)
                ENDDO
                CALL DAWK(IWS,1,ERROR,*9999)
                CALL ACWK(IWT,0,ERROR,*9999)
                REGEN=.TRUE.
              ELSE IF(OPTION(NOCH1)(1:14).EQ.'Enter Signal') THEN
                WRITE(OP_STRING,'($,'' Enter signal number ? '')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                READ(IOIP,'(A)') C2
                CALL TRIM(C2,IBEG2,IEND2)
                DO I=1,128
                  C1=SIGHEADER(1,I)
                  CALL TRIM(C1,IBEG1,IEND1)
                  IF(C1(IBEG1:IEND1).EQ.C2(IBEG2:IEND2)) THEN
                    ISIGNAL=I
                  ENDIF
                ENDDO
                IF(ISIGNAL.LT.1) ISIGNAL=1
                IF(ISIGNAL.GT.128) ISIGNAL=128
                REGEN=.TRUE.
              ELSE IF(OPTION(NOCH1)(1:15).EQ.'Recalculate..') THEN
                OPTION2(1)='Fiducial'
                OPTION2(2)='Datum'
                OPTION2(3)='Exit'
c               CALL CHOICE('DITRA1',1,1,INS,72,'REQUEST',3,3,NOCH2,
c    '            noco,9,
c    '            CO,OPTION2,STRING,0.1*XDISP,0.35*XDISP,ERROR,*9999)
                IF(OPTION2(NOCH2)(1:8).EQ.'Fiducial') THEN
                  WRITE(OP_STRING,
     '              '('' >>Calculating Fiducial Markers'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  DATUM=1.0d0
                  DO ISIG=1,128
                    CALL CALC_FID(ISAMPLES,ISIG,NAV,FIDMARK(ISIG),
     '                TTSTART,TTEND, ACCEPT(ISIG),ERROR,*9999)
                    IF(FIDMARK(ISIG).LT.DATUM.AND.ACCEPT(ISIG))
     '                DATUM=FIDMARK(ISIG)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' Fiducial marker '',I3,'
     '                  //''' is '',F6.3)') ISIG,FIDMARK(ISIG)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDDO
                  WRITE(OP_STRING,
     '              '('' >>>Datum calc.d : '',F6.2,'' ms'')')
     '              DATUM*1000.0d0
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  REGENALL=.TRUE.
                  REGEN=.TRUE.
                ELSE IF(OPTION2(NOCH2)(1:5).EQ.'Datum') THEN
                  DATUM=1.0d0
                  DO ISIG=1,128
                    IF(FIDMARK(ISIG).LT.DATUM.AND.ACCEPT(ISIG))
     '                DATUM=FIDMARK(ISIG)
                  ENDDO                                         
                  WRITE(OP_STRING,
     '              '('' >>>Datum calc.d : '',F6.2,'' ms'')')
     '              DATUM*1000.0d0
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  REGENALL=.TRUE.
                  REGEN=.TRUE.
                ELSE IF(OPTION2(NOCH2)(1:4).EQ.'Exit') THEN
c                  IF(IWKS(72).GT.0)
c    '               CALL INPUT_MODE(72,1,'CHOICE','REQUEST',ERROR,
c    '                 *9999)
                ENDIF
              ELSE IF(OPTION(NOCH1)(1:4).EQ.'Exit') THEN
c               CALL INPUT_MODE(IWT,LD1,'PICK','REQUEST',ERROR,*9999)
                CONTINUE=.FALSE.
              ENDIF
            ENDIF
            IF(REGENALL) THEN
              CALL DAWK(IWT,0,ERROR,*9999)
              CALL ACWK(IWS,0,ERROR,*9999)
              DO ISIG=1,128
                IF(ISIG.EQ.ISIGNAL) THEN
                  ICOL=4
                ELSE IF(ACCEPT(ISIG)) THEN
                  ICOL=2
                ELSE !Reject
                  ICOL=3
                ENDIF
                CALL SGSIGN(ICOL,INDEX,ISAMPLES,ISEG,ISIG,
     '            ISSIGN(1,ISIG),IWS,CSEG,SIGHEADER(1,ISIG),DATUM,
     '            FIDMARK(ISIG),TSSTART,TSEND,ACCEPT(ISIG),ERROR,*9999)
              ENDDO
              CALL DAWK(IWS,0,ERROR,*9999)
              CALL ACWK(IWT,0,ERROR,*9999)
c             CALL FLUSH_DEVICE_EVENTS(IWT,'ALL',ERROR,*9999)
              REGENALL=.FALSE.
            ENDIF
            IF(REGEN) THEN
              IF(ISIGLAST.NE.ISIGNAL) THEN
                CALL DAWK(IWT,0,ERROR,*9999)
                CALL ACWK(IWS,1,ERROR,*9999)
                IF(ACCEPT(ISIGLAST)) THEN
                  ICOL=2
                ELSE
                  ICOL=3
                ENDIF
                CALL SGSIGN(ICOL,INDEX,ISAMPLES,ISEG,ISIGLAST,
     '            ISSIGN(1,ISIGLAST),IWS,CSEG,SIGHEADER(1,ISIGLAST),
     '            DATUM,FIDMARK(ISIGLAST),
     '            TSSTART,TSEND,ACCEPT(ISIGLAST),ERROR,*9999)
                  !Show last signal

                CALL SGSIGN(4,INDEX,ISAMPLES,ISEG,ISIGNAL,
     '            ISSIGN(1,ISIGNAL),IWS,CSEG,SIGHEADER(1,ISIGNAL),
     '            DATUM,FIDMARK(ISIGNAL),TSSTART,TSEND,
     '            ACCEPT(ISIGNAL),ERROR,*9999)
                  !Show current signal
                CALL DAWK(IWS,1,ERROR,*9999)
                CALL ACWK(IWT,0,ERROR,*9999)  !Display trace window
              ENDIF
              IF(ACCEPT(ISIGNAL)) THEN
                ICOL=2
                ACCREJ='REJECT'
              ELSE
                ICOL=3
                ACCREJ='ACCEPT'
              ENDIF
              CALL SGTRAC(ICOL,INDEX,ISAMPLES,ISEG,ISIGNAL,ISTRAC,IWT,
     '          CSEG,SIGHEADER(1,ISIGNAL),DATUM,FIDMARK(ISIGNAL),
     '          TTSTART,TTEND,ERROR,*9999)
              CALL DETECT(IWT,ISEG,ISTRAC(2),'DETECTABLE',ERROR,*9999) ! Fid M.
              CALL OPEN_SEGMENT(ISEG4,ISEG,IWT,'ACCEPT/REJECT_LABEL',
     '          1,INDEX_OLD,4,1,CSEG,ERROR,*9999)                ! Accept/Reject with label
              PL(1,1)=8.0D0
              PL(2,1)=1.3D0
              CALL TEXT(INDEX(8),IWT,ACCREJ,PL,ERROR,*9999)
              CALL CLOSE_SEGMENT(ISEG4,IWT,ERROR,*9999)
              CALL DETECT(IWT,ISEG,ISEG4,'DETECTABLE',ERROR,*9999)
              REGEN=.FALSE.
              CONTINUE=.TRUE.
c             CALL FLUSH_DEVICE_EVENTS(IWT,'ALL',ERROR,*9999)
            ENDIF
            IF(REMARK) THEN
              CALL OPEN_SEGMENT(ISEG0,ISEG,IWT,'START_MARKER',1,
     '          INDEX_OLD,0,1,CSEG,ERROR,*9999)                     ! TSTART marker
              PL(1,1)=XPOS+9.0D0*TTSTART
              PL(2,1)=YPOS
              CALL POLYMARKER(INDEX(9),IWT,1,PL,ERROR,*9999)
              CALL CLOSE_SEGMENT(ISEG0,IWT,ERROR,*9999)
              CALL DETECT(IWT,ISEG,ISEG0,'DETECTABLE',ERROR,*9999)

              CALL OPEN_SEGMENT(ISEG1,ISEG,IWT,'END_MARKER',1,INDEX_OLD,
     '          1,1,CSEG,ERROR,*9999)                     ! TEND marker
              PL(1,1)=XPOS+9.0D0*TTEND
              PL(2,1)=YPOS
              CALL POLYMARKER(INDEX(9),IWT,1,PL,ERROR,*9999)
              CALL CLOSE_SEGMENT(ISEG1,IWT,ERROR,*9999)
              CALL DETECT(IWT,ISEG,ISEG1,'DETECTABLE',ERROR,*9999)

              CALL OPEN_SEGMENT(ISEG5,ISEG,IWT,'BAND_MARKER',1,
     '          INDEX_OLD,5,1,CSEG,ERROR,*9999)             ! Move band with marker
              PL(1,1)=XPOS+9.0D0*TTBAND
              PL(2,1)=YPOS+0.1D0
              PL(1,2)=XPOS+9.0D0*TTBAND
              PL(2,2)=YPOS-0.1D0
              CALL POLYLINE(INDEX(10),IWT,2,PL,ERROR,*9999)
              CALL CLOSE_SEGMENT(ISEG5,IWT,ERROR,*9999)
              CALL DETECT(IWT,ISEG,ISEG5,'DETECTABLE',ERROR,*9999)
c             CALL INPUT_MODE(IWT,LD1,'PICK','EVENT',ERROR,*9999)
              REMARK=.FALSE.
            ENDIF
          ENDDO
          CALL DELETE_SEGMENT(ISEG0,ISEG,IWT,ERROR,*9999) !Delete markers and
          CALL DELETE_SEGMENT(ISEG1,ISEG,IWT,ERROR,*9999) !pick text
          CALL DELETE_SEGMENT(ISEG2,ISEG,IWT,ERROR,*9999)
          CALL DELETE_SEGMENT(ISEG3,ISEG,IWT,ERROR,*9999)
          CALL DELETE_SEGMENT(ISEG4,ISEG,IWT,ERROR,*9999)
          CALL DELETE_SEGMENT(ISEG5,ISEG,IWT,ERROR,*9999)
          IF(ACCEPT(ISIGNAL)) THEN
            ICOL=2
          ELSE
            ICOL=3
          ENDIF
          CALL SGTRAC(ICOL,INDEX,ISAMPLES,ISEG,ISIGNAL,ISTRAC,IWT,CSEG,
     '      SIGHEADER(1,ISIGNAL),DATUM,FIDMARK(ISIGNAL),TTSTART,TTEND,
     '      ERROR,*9999)
          CALL DAWK(IWT,0,ERROR,*9999)
        ENDIF  !USE_SOCKET
      ENDIF

      CALL EXITS('DITRAC')
      RETURN
 9999 CALL ERRORS('DITRAC',ERROR)
      CALL EXITS('DITRAC')
      RETURN 1
      END


C 25/2/97 LC removed from : 

C#### Subroutine: EXELEM
C###  Description:
C###    EXELEM exports element data from finite element data base. 

C CBP 3/8/95 extending sectors
C?????????????????????????DB.  Only good for Lagrange
C                          IF(2.EQ.NIT(nb)) THEN
C                            NN_TOT=NNT(nb)+IBT(2,1,NB)
C                            NN_LAY=NNT(nb)
C                          ELSE
C                            IF(IBT(1,3,NB).EQ.2) THEN
C                              NN_TOT=NNT(nb)+2*IBT(2,1,NB)
C                              NN_LAY=NNT(nb)/2
C                            ELSE
C                              NN_TOT=NNT(nb)+(IBT(2,3,NB)+1)*IBT(2,1,NB)
C                              NN_LAY=NNT(nb)/(IBT(2,3,NB)+1)
C                            ENDIF
C                          ENDIF
C                          IF(DATAFILE) THEN
C                            WRITE(IFILE,'(5X,''#Nodes='',I2)') NN_TOT
C                          ELSE
C                            IF(FSKWRITE(NN_TOT,SK_LONG_INT,1,
C     '                        CONNID2).EQ.-1) GOTO 9999
C                          ENDIF
C                          WRITE(IFILE,'(5X,''#Nodes='',I2)') NN_TOT
C                          DO nn=1,NNT(nb)
C                            NP=NPNE(nn,nb,NE)
C                            NP_ELEMENT=ELEMENT_NODE_LIST(0)
C                            DO WHILE ((NP_ELEMENT.GT.0).AND.
C     '                        (NP.NE.ELEMENT_NODE_LIST(NP_ELEMENT)))
C                              NP_ELEMENT=NP_ELEMENT-1
C                            ENDDO
C                            IF(MOD(nn,NN_LAY).EQ.1) THEN
C                              N=IBT(2,1,NB)+1
C                            ELSE
C                              N=1
C                            ENDIF
C                            DO WHILE (N.GT.0)
C                              IF(DATAFILE) THEN
C                                WRITE(IFILE,'(5X,I2,''.  #Values='','
C     '                            //'I1)') NP_ELEMENT,NKT(nn,NB)
C                              ELSE
C                                IF(FSKWRITE(NP_ELEMENT,SK_LONG_INT,1,
C     '                            CONNID2).EQ.-1) GOTO 9999
C                                IF(FSKWRITE(NKT(nn,NB),SK_LONG_INT,1,
C     '                            CONNID2).EQ.-1) GOTO 9999
C                              ENDIF
C                              NK_TOT=0
C                              IF(DATAFILE) THEN
C                                WRITE(IFILE,'(7X,''Value indices:  '','
C     '                            //'12(1X,I3))')
C     '                            (NK_TOT+NKE(MK,nn,nb,NE),MK=1,
C     '                            NKT(nn,NB))
C                                WRITE(IFILE,'(7X,'
C     '                            //'''Scale factor indices:'','
C     '                            //'12(1X,I3))')
C     '                            (NS_TOT+NK,NK=1,NKT(nn,NB))
C                              ELSE
C                                DO mk=1,NKT(nn,NB)
C                                  IF(FSKWRITE(NK_TOT+NKE(mk,nn,nb,NE),
C     '                              SK_LONG_INT,1,CONNID2).EQ.-1)
C     '                              GOTO 9999
C                                ENDDO
C                                DO nk=1,NKT(nn,NB)
C                                  IF(FSKWRITE(NS_TOT+nk,SK_LONG_INT,1,
C     '                              CONNID2).EQ.-1) GOTO 9999
C                                ENDDO
C                              ENDIF
C                              NS_TOT=NS_TOT+NKT(nn,NB)
C                              N=N-1
C                            ENDDO
C                          ENDDO


C 25/2/97 LC removed section from : 

C#### Subroutine: EXELEM
C###  Description:
C###    EXELEM exports element data from finite element data base. 

C CPB 3/8/95 Sector elements are now the same as the other elements
C for the ns calculation
C                      IF(SPECIAL_BASIS_FLAG.EQ.3) THEN
C                        IF(2.EQ.NIT(nb)) THEN
C                          NN_LAY=NNT(nb)
C                        ELSE
C                          IF(IBT(1,3,NB).EQ.2) THEN
C                            NN_LAY=NNT(nb)/2
C                          ELSE
C                            NN_LAY=NNT(nb)/(IBT(2,3,NB)+1)
C                          ENDIF
C                        ENDIF
C                        NS_TOT=0
C                        DO nn=1,NNT(nb)
C                          IF(MOD(nn,NN_LAY).EQ.1) THEN
C                            N=IBT(2,1,NB)+1
C                          ELSE
C                            N=1
C                          ENDIF
C                          DO WHILE (N.GT.0)
C                            IF(DATAFILE) THEN
C                              WRITE(IFILE,'(4X,5(1X,E24.16))')
C     '                          (SE(nk,nb,NE),NK=NS_TOT+1,
C     '                          NS_TOT+NKT(nn,NB))
C                            ELSE
C                              IF(FSKWRITE(SE(NS_TOT+1,nb,NE),
C     '                          SK_DOUBLE_FLOAT,NKT(nn,NB),
C     '                          CONNID2).EQ.-1) GOTO 9999
C                            ENDIF
C                            N=N-1
C                          ENDDO
C                          NS_TOT=NS_TOT+NKT(nn,NB)
C                        ENDDO  
C                      ELSE


C 25/2/97 LC removed section from : 
C#### Subroutine: EXGEOM
C###  Description:
C###    EXGEOM exports geometry data from finite element data base. 


C Below is CPBs original .fac, .pts and .channels output.
C Altered AJP 9-3-94
C Write out nodal coordinates in order
C
C          CALL OPENF(IOFILE2,'DISK',FILE(IBEG:IEND)//'.pts','NEW',
C     '      'SEQUEN',132,ERROR,*9999)
C
C          IF(TYPE(1:6).EQ.'NUMBER') THEN
C            DO NONODE=1,NPLIST(0)
C              WRITE(IOFILE2,*) (XP(1,nv,nj,NPLIST(NONODE)),NJ=1,NJT)
C            ENDDO
C          ELSE IF(TYPE(1:5).EQ.'GROUP') THEN
C            CALL TRIM(GROUP_NAME,IBEG1,IEND1)
C            NOGRNO=0
C            DO NOGROUP=1,NTGRNO
C              CALL TRIM(LAGRNO(NOGROUP),IBEG2,IEND2)
C              IF(GROUP_NAME(IBEG1:IEND1).EQ.
C     '          LAGRNO(NOGROUP)(IBEG2:IEND2)) NOGRNO=NOGROUP
C            ENDDO
C            IF(NOGRNO.EQ.0) THEN
C              ERROR='>>Group not found'
C              GOTO 9999
C            ENDIF
C            DO NONODE=1,LIGRNO(0,NOGRNO)
C              WRITE(IOFILE2,*) (XP(1,nv,nj,LIGRNO(NONODE,NOGRNO)),
C     '          NJ=1,NJT)
C            ENDDO
C          ENDIF
C
C          CALL CLOSEF(IOFILE2,ERROR,*9999)
C
CC Write out the .channels file
C         
C          CALL OPENF(IOFILE2,'DISK',FILE(IBEG:IEND)//'.channels',
C     '      'NEW','SEQUEN',132,ERROR,*9999)
C
C          IF(TYPE(1:6).EQ.'NUMBER') THEN
C            WRITE(IOFILE2,'(I3,'' Nodes'')') NPLIST(0)
C            DO NONODE=1,NPLIST(0)
C              WRITE(IOFILE2,*) NONODE,NPLIST(NONODE)
C            ENDDO
C          ELSE IF(TYPE(1:5).EQ.'GROUP') THEN
C            WRITE(IOFILE2,'(I3,'' Nodes'')') LIGRNO(0,NOGRNO)
C            DO NONODE=1,LIGRNO(0,NOGRNO)
C              WRITE(IOFILE2,*) NONODE,LIGRNO(NONODE,NOGRNO)
C            ENDDO
C          ENDIF
C
C          CALL CLOSEF(IOFILE2,ERROR,*9999)
C
CC Write out the facets
C
C          CALL OPENF(IOFILE2,'DISK',FILE(IBEG:IEND)//'.fac','NEW',
C     '      'SEQUEN',132,ERROR,*9999)
C
C          IF(TYPE(1:6).EQ.'NUMBER') THEN
C            DO NOELEM=1,NELIST(0)
C              NE=NELIST(NOELEM)
C              nb=NBJ(1,NE)
C
CC Find the corresponding vertex number for the local nodes from NPLIST
C
C              DO nn=1,NNT(nb)
C                NP=NPNE(nn,nb,NE)
C                I=1
C                DO WHILE(NP.NE.NPLIST(I).AND.I.LE.NPLIST(0))
C                  I=I+1
C                ENDDO
C                LOCAL_NP(nn)=I
C              ENDDO
C              IF(IBT(1,1,NB).EQ.2.AND.IBT(1,2,NB).EQ.2) THEN ! Bicubic Hermite
C
CC Split the Bicubic Hermite quadralateral element into two triangluar elements
CC The first triangle has local NN vertices of 1 2 3 and the second triangle
CC has local NN verticies of 2 4 3
C
CC Write out the two triangular facets
C            
C              WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(1),LOCAL_NP(2),
C     '          LOCAL_NP(3)
C              WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(2),LOCAL_NP(4),
C     '          LOCAL_NP(3)
C
C              ELSE IF(IBT(1,2,NB).EQ.3.AND.IBT(1,2,NB).EQ.3) THEN 
C                                                           ! Hermite-Simplex
C
CC Write out the Hermite-Simplex element as one triangular facet
C
C                IF(NKT(1,NB).EQ.1) THEN ! Apex at node 1
C                  WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(1),
C     '              LOCAL_NP(3),LOCAL_NP(2)
C                ELSE IF(NKT(3,NB).EQ.1) THEN ! Apex at node 3
C                  WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(1),LOCAL_NP(2),
C     '              LOCAL_NP(3)
C                ENDIF
C              ENDIF
C            ENDDO
C          ELSE IF(TYPE(1:5).EQ.'GROUP') THEN
C            CALL TRIM(GROUP_NAME,IBEG1,IEND1)
C            NOGREL=0
C            DO NOGROUP=1,NTGREL
C              CALL TRIM(LAGREL(NOGROUP),IBEG2,IEND2)
C              IF(GROUP_NAME(IBEG1:IEND1).EQ.
C     '          LAGREL(NOGROUP)(IBEG2:IEND2)) NOGREL=NOGROUP
C            ENDDO
C            IF(NOGREL.EQ.0) THEN
C              ERROR='>>Group not found'
C              GOTO 9999
C            ENDIF
C            DO NOELEM=1,LIGREL(0,NOGREL)
C              NE=LIGREL(NOELEM,NOGREL)
C              nb=NBJ(1,NE)
C              DO nn=1,NNT(nb)
C                NP=NPNE(nn,nb,NE)
C                I=1
C                DO WHILE(NP.NE.LIGRNO(I,NOGRNO).AND.
C     '            I.LE.LIGRNO(0,NOGRNO))
C                  I=I+1
C                ENDDO
C                LOCAL_NP(nn)=I
C              ENDDO
C              IF(IBT(1,1,NB).EQ.2.AND.IBT(1,2,NB).EQ.2) THEN ! Bicubic Hermite
C
CC Split the Bicubic Hermite quadralateral element into two triangluar elements
CC The first triangle has local NN veritices of 1 2 3 and the second triangle
CC has local NN verticies of 2 4 3
C
CC Write out the two triangular facets
C            
C              WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(1),LOCAL_NP(2),
C     '          LOCAL_NP(3)
C              WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(2),LOCAL_NP(4),
C     '          LOCAL_NP(3)
C
C              ELSE IF(IBT(1,2,NB).EQ.3.AND.IBT(1,2,NB).EQ.3) THEN 
C                                                           ! Hermite-Simplex
C
CC Write out the Hermite-Simplex element as one triangular facet
C
C                IF(NKT(1,NB).EQ.1) THEN ! Apex at node 1
C                  WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(1),LOCAL_NP(3),
C     '              LOCAL_NP(2)
C                ELSE IF(NKT(3,NB).EQ.1) THEN ! Apex at node 3
C                  WRITE(UNIT=IOFILE2,FMT=*) LOCAL_NP(1),LOCAL_NP(2),
C     '              LOCAL_NP(3)
C                ENDIF
C              ENDIF
C            ENDDO
C          ENDIF
C          CALL CLOSEF(IOFILE2,ERROR,*9999)
C End old CPB output



Module FE27
===========

      SUBROUTINE LIFIBR(IDO,NBJ,NEELEM,NKJ,NPLIST,NPNODE,NRLIST,
     '  XA,XP,STRING,ERROR,*)

C#### Subroutine: LIFIBR
C###  Description:
C###    LIFIBR lists fibre information.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NKJ(NJM,NPM),
     '  NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM)
      REAL*8 XA(NAM,NJM,NQM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,no_nrlist,nr
      CHARACTER FILE*100
      LOGICAL ALL_REGIONS,OPFILE

      CALL ENTERS('LIFIBR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list fibre<;FILENAME>
C###  Parameter:    <nodes (GROUP/#s/all)[all]>
C###  Parameter:    <region (#s/all)[1]>
C###  Description:
C###    Lists fibres to screen or file FILENAME.opfibr if qualifier
C###    present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<nodes (GROUP/#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region #s/ALL>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIFIBR',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opfibr','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          CALL OPFIBR(IDO,NBJ,NEELEM,NKJ,NPNODE,nr,XA,XP,
     '      ERROR,*9999)
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIFIBR')
      RETURN
 9999 CALL ERRORS('LIFIBR',ERROR)
      CALL EXITS('LIFIBR')
      RETURN 1
      END
 

      SUBROUTINE LIMACR(STRING,ERROR,*)

C#### Subroutine: LIMACR
C###  Description:
C###    LIMACR lists macro commands.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:gtstr00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,
     '  macro_command,macrokey,nomacro,N3CO
      CHARACTER C1*255
      LOGICAL CBBREV

      CALL ENTERS('LIMACR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: list macro key NUMBER# 
C###  Description:
C###    Lists macro commands.
        OP_STRING(1)=' list macro key NUMBER#'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: list macro <(NAME/all)[all]> 
C###  Description:
C###    Lists macro commands.

        OP_STRING(1)=' list macro <(NAME/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','LIMACR',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'KEY',1,noco+1,NTCO,N3CO)) THEN
          IF(N3CO+1.LE.NTCO) THEN !list specified key
            macrokey=IFROMC(CO(N3CO+1))
            CALL ASSERT(macrokey.GE.1.AND.macrokey.LE.9,
     '        'key# must be between 1 & 9',ERROR,*9999)
            WRITE(OP_STRING,'('' Macro key number '',I1,'':'')')
     '        macrokey
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DO nomacro=1,NT_MACRO(macrokey)
              WRITE(OP_STRING,'('' Line number '',I2,'': '',A)')
     '          nomacro,MACRO_KEY_buffer(nomacro,macrokey)(1:80)
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO !nomacro

          ELSE !list all keys
            DO macrokey=1,NT_KEY
              WRITE(OP_STRING,'('' Macro key number '',I1,'':'')')
     '          macrokey
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nomacro=1,NT_MACRO(macrokey)
                WRITE(OP_STRING,'('' Line number '',I2,'': '',A)')
     '            nomacro,MACRO_KEY_buffer(nomacro,macrokey)(1:80)
      	        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO !nomacro
            ENDDO !macrokey
          ENDIF

        ELSE !named macro
          CALL CUPPER(CO(noco+1),C1)
C          IF(CBBREV(MACRO_names,CUPPER(CO(noco+1)),1,1,
          IF(CBBREV(MACRO_names,C1,1,1,
     '      NT_MACRO_names(0),macro_command)) THEN
            WRITE(OP_STRING,'('' Macro command number '',I1,'':'')')
     '        macro_command
      	    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            DO nomacro=1,NT_MACRO_names(macro_command)
              WRITE(OP_STRING,'('' Line number '',I2,'': '',A)')
     '          nomacro,
     '          MACRO_COMMAND_buffer(nomacro,macro_command)(1:80)
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDDO !nomacro

          ELSE !list all macros
            DO macro_command=1,NT_MACRO_names(0)
              WRITE(OP_STRING,'(/'' Macro command number '',I1,'': '','
     '          //'A)') macro_command,MACRO_names(macro_command)
      	      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              DO nomacro=1,NT_MACRO_names(macro_command)
                WRITE(OP_STRING,'('' Line number '',I2,'': '',A)')
     '            nomacro,
     '            MACRO_COMMAND_buffer(nomacro,macro_command)(1:80)
      	        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDDO !nomacro
            ENDDO !macro_command
          ENDIF
        ENDIF !key/name
      ENDIF

      CALL EXITS('LIMACR')
      RETURN
 9999 CALL ERRORS('LIMACR',ERROR)
      CALL EXITS('LIMACR')
      RETURN 1
      END        


      SUBROUTINE LIMATR(noco,NTCOQU,CO,COQU,STRING,ERROR,*)

C#### Subroutine: LIMATR
C###  Description:
C**** Lists matrix stored in AMATR.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:coef00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER noco,NTCOQU(*)
      CHARACTER CO(*)*(*),COQU(NXCO,*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LIMATR',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIMATR',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opmatr','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL OPMATR(ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIMATR')
      RETURN
 9999 CALL ERRORS('LIMATR',ERROR)
      CALL EXITS('LIMATR')
      RETURN 1
      END

C KAT 2001-12-14
      SUBROUTINE LITEXT(STRING,ERROR,*)

C#### Subroutine: LITEXT
C###  Description:
C###    LITEXT lists text.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LITEXT',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LITEXT',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optext','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        CALL OPTEXT(ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LITEXT')
      RETURN
 9999 CALL ERRORS('LITEXT',ERROR)
      CALL EXITS('LITEXT')
      RETURN 1
      END


C18-nov-1999
      SUBROUTINE LITIME(STRING,ERROR,*)

C#### Subroutine: LITIME
C###  Description:
C###    LITIME lists time variables
C *** Created 10 August 1999 - David Nickerson

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,VARIABLE,N3CO
      CHARACTER FILE*100
      LOGICAL OPFILE,CBBREV,ABBREV,NODAL_INFO

      CALL ENTERS('LITIME',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list time_variable<;FILENAME>
C###  Parameter:      <variable all/#[all]>
C###    Specify a variable to list, defaults to all variables
C###  Parameter:      <node_information>
C###    If this qualifier is present, the time variable(s) nodal 
C###    information is written, otherwise it is not.
C###  Description:
C###    Lists time variable(s) to the screen or to FILENAME.optime

C---------------------------------------------------------------------

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<variable all/#[all]>'
        OP_STRING(3)=BLANK(1:15)//'<nodal_information>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LITIME',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optime','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF
        IF(CBBREV(CO,'VARIABLE',3,NOCO+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
            VARIABLE=-1
          ELSE
            CALL PARSIN(CO(N3CO+1),VARIABLE,ERROR,*1)
            CALL ASSERT(VARIABLE.GT.0,'Variable needs to be > 0',
     '        ERROR,*1)
          ENDIF
        ELSE       
          VARIABLE=-1
        ENDIF
        IF(CBBREV(CO,'NODE_INFORMATION',4,NOCO+1,NTCO,N3CO)) THEN
          NODAL_INFO=.TRUE.
        ELSE
          NODAL_INFO=.FALSE.
        ENDIF

        CALL OPTIME(VARIABLE,NODAL_INFO,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LITIME')
      RETURN
 9999 CALL ERRORS('LITIME',ERROR)
      CALL EXITS('LITIME')
      RETURN 1
      END


C LKC archived 11-NOV-1998 unused

      SUBROUTINE LITRAN(STRING,ERROR,*)

C#### Subroutine: LITRAN
C###  Description:
C###    LITRAN lists Phigs transformation parameters.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER FILE*100
      LOGICAL OPFILE

      CALL ENTERS('LITRAN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list transformation<;FILENAME>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LITRAN',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.optran','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL OPTRAN(ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LITRAN')
      RETURN
 9999 CALL ERRORS('LITRAN',ERROR)
      CALL EXITS('LITRAN')
      RETURN 1
      END


      SUBROUTINE LIVSAE(STRING,ERROR,*)

C#### Subroutine: LIVSAE
C###  Description:
C###    LIVSAE lists vsaero parameters.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND
      CHARACTER TYPE*22

      CALL ENTERS('LIVSAE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list vsaero basic_input/patch_geometry/wake_input/
C###  Parameter:            surface_streamlines/off_body_velocity_scan/
C###  Parameter:            off_body_streamlines
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' basic_input/'
        OP_STRING(2)=BLANK(1:15) //'patch_geometry/'
        OP_STRING(3)=BLANK(1:15) //'wake_input/'
        OP_STRING(4)=BLANK(1:15) //'surface_streamlines/'
        OP_STRING(5)=BLANK(1:15) //'off_body_velocity_scan/'
        OP_STRING(6)=BLANK(1:15) //'off_body_streamlines'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIVSAE',ERROR,*9999)
      ELSE

        CALL TRIM(CO(noco+1),IBEG,IEND)
        TYPE=CO(noco+1)(IBEG:IEND)
        CALL OPVSAE(TYPE,ERROR,*9999)

      ENDIF

      CALL EXITS('LIVSAE')
      RETURN
 9999 CALL ERRORS('LIVSAE',ERROR)
      CALL EXITS('LIVSAE')
      RETURN 1
      END


Module FE28
=========== 

C KAT 2001-12-14
      SUBROUTINE DRELEC(ISEG,ISELEC,NDDL,NDLT,NEELEM,ZD,CSEG,
     '  STRING,ERROR,*)

C#### Subroutine: DRELEC
C###  Description:
C###    DRELEC defines electrodes on hammer map of electrical field.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER ISEG(*),ISELEC(128),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM)
      REAL*8 ZD(NJM,NDM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),nd,nde,ne,ne1,nr,NTIW

      CALL ENTERS('DRELEC',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw electrodes 
C###  Description:
C###    This command draws electrodes on hammer map of 
C###    electrical field.
C###  Parameter:    <on WS_ID#[4]>
C###    Specify the number of the window to draw the electrodes on.

        OP_STRING(1)=STRING(1:IEND)//' <on WS_ID#[4]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRELEC',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        nr=1  !may need generalizing
        CALL WS_LIST(IWK,4,NTIW,noco,NTCO,CO,ERROR,*9999)

        IW=IWK(1)
        CALL ACWK(iw,1,ERROR,*9999)
        DO ne1=1,NEELEM(0,nr)
          ne=NEELEM(ne1,nr)
          DO nde=1,NDLT(ne)
            nd=NDDL(ne,nde)
            CALL SGELEC(ISEG,ISELEC(nd),iw,nd,CSEG,ZD,ERROR,*9999)
          ENDDO
        ENDDO
        CALL DAWK(iw,1,ERROR,*9999)

      ENDIF

      CALL EXITS('DRELEC')
      RETURN
 9999 CALL ERRORS('DRELEC',ERROR)
      CALL EXITS('DRELEC')
      RETURN 1
      END


C KAT 2001-12-13
      SUBROUTINE DRIMAG(ISEG,ISIMAG,CSEG,STRING,ERROR,*)

C#### Subroutine: DRIMAG
C###  Description:
C###    DRIMAG draws image.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:pics01.cmn'
!     Parameter List
      INTEGER ISEG(*),ISIMAG(NWM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IFROMC,IEND,IEND1,INDEX,iw,IWK(6),
     '  N3CO,NDIMX,NDIMY,noiw,NTIW
      REAL XMAX_CA,XMIN_CA,YMAX_CA,YMIN_CA
      LOGICAL CBBREV

      CALL ENTERS('DRIMAG',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw image 
C###  Description:
C###    Draw an image.
C###  Parameter:    <on WS_ID#[1]>
C###    Specify the workstation on which to draw.
C###  Parameter:    <resolution X#[$size] <by Y#[X#]>>
C###    Specify the resolution.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<on WS_ID#[1]'
        OP_STRING(3)=BLANK(1:15)
     '    //'<resolution X#[$size] <by Y#[X#]>>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRIMAG',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'RESOLUTION',1,noco+1,NTCO,N3CO)) THEN
          NDIMX=IFROMC(CO(N3CO+1))
          IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
            NDIMY=IFROMC(CO(N3CO+1))
          ELSE
            NDIMY=NDIMX
          ENDIF
        ELSE
          NDIMX=IMGX
          NDIMY=IMGY
        ENDIF
        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Display image size '',I4,'' by '',I4)') NDIMX,NDIMY
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        CALL ASSERT(ndimx.NE.0,'>>Invalid X Resolution',ERROR,*9999)
        CALL ASSERT(ndimy.NE.0,'>>Invalid Y Resolution',ERROR,*9999)
c cpb 14/11/94 Temporarily commenting this out
C        CALL ASSERT(NXM.EQ.512,'>>Require NXM = 512 for image arrays',
C     '    ERROR,*9999)
      
        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          XMIN_CA=REAL(XMIN)
          XMAX_CA=REAL(XMAX)
          YMIN_CA=REAL(YMIN)
          YMAX_CA=REAL(YMAX)
          CALL SGIMAG(INDEX,ISEG,ISIMAG(iw),iw,NDIMX,NDIMY,CSEG,
     '      XMIN_CA,XMAX_CA,YMIN_CA,YMAX_CA,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('DRIMAG')
      RETURN
 9999 CALL ERRORS('DRIMAG',ERROR)
      CALL EXITS('DRIMAG')
      RETURN 1
      END


C LC 25/2/97 removed section from :

C#### Subroutine: DRNODS
C###  Description:
C###    DRNODS draws nodal parameters.
C**** NPNODE(0,nr) is the total number of nodes in region nr.
C**** NPNODE(nonode,nr), nonode=1..NPNODE(0,nr) are the node numbers.
C**** NPT(nr) is the highest node number.

c Is this code used? PJH 7Oct94
c       IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO)) THEN
c         TYPE='REFERENCE'
c       ELSE IF(CBBREV(CO,'DEPENDENT',3,noco+1,NTCO,N3CO)) THEN
c         TYPE='DEPENDENT'
c         NIY=IFROMC(CO(N3CO+1))
c       ELSE
c         TYPE='REFERENCE'
c       ENDIF
c
c       IF(CBBREV(CO,'INTERPOLATE',1,noco+1,NTCO,N3CO)) THEN
c         WRITE(OP_STRING,*) 
c    '      ' Note: elements should be defined before nodes'
c    '      //' when using this command'
c         CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
c         INTER=.TRUE.
c         IF(NTCO.GT.N3CO) THEN
c           DELIST=.TRUE.
c           CALL PARSIL(CO(N3CO+1),7,NTLIST,NDLIST,ERROR,*9999)
c         ELSE
c           DELIST=.FALSE.
c         ENDIF
c       ENDIF


      SUBROUTINE DRTEXT(ISEG,ISTEXT,NOCO,NTCO,
     '  CO,CSEG,STRING,ERROR,*)

C**** Draw text.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:back00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:cbwk01.cmn'
      INCLUDE 'cmiss$reference:colo00.cmn'
      INCLUDE 'cmiss$reference:file00.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:graf00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:jtyp00.cmn'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:ntsg00.cmn'
      INCLUDE 'cmiss$reference:text00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISTEXT(*),NOCO,NTCO
      CHARACTER CO(*)*(*),CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IEND,IEND1,INDEX,INDEX_TEXT,
     '  IW,IWK(6),LIST(20),N3CO,nolist,NTIW,NTLIST,NTRL
      REAL*8 RL(2)
      CHARACTER FILE*100
      LOGICAL CBBREV,MOUSE

      CALL ENTERS('DRTEXT',*9999)
 1    IF(CO(NOCO+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)
        CALL TRIM(FILE00,IBEG1,IEND1)

        OP_STRING(1)= BLANK(1:IEND)
        OP_STRING(2)=BLANK(1:IEND+1)//'<on WS_ID>[1]'
        OP_STRING(3)=BLANK(1:IEND+1)//'<locate/at VIEWPORT_COORDS>[0,0]'
        OP_STRING(4)=BLANK(1:IEND+1)
     '                             //'<all/number TEXT_list/last>[last]'
        OP_STRING(5)=BLANK(1:IEND+1)//'<rgb=RGB>[black]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(NOCO+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRTEXT',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,NOCO,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'LOCATE',2,NOCO+1,NTCO,N3CO)) THEN
          MOUSE=.TRUE.
        ELSE IF(CBBREV(CO,'AT',2,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),2,NTRL,RL,ERROR,*9999)
          XWC_TEXT=RL(1)
          YWC_TEXT=RL(2)
        ELSE
          XWC_TEXT=0.0D0
          YWC_TEXT=0.0D0
        ENDIF
        IF(CBBREV(CO,'ALL',2,NOCO+1,NTCO,N3CO)) THEN
          NTLIST=NTTEXT
          DO nolist=1,NTLIST
            LIST(nolist)=nolist
          ENDDO
        ELSE IF(CBBREV(CO,'NUMBER',2,NOCO+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),20,NTLIST,LIST,ERROR,*9999)
        ELSE
          NTLIST=1
          LIST(1)=NTTEXT
        ENDIF

        IF(CBBREV(CO,'RGB',3,NOCO+1,NTCO,N3CO)) THEN
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
        ENDIF
     
        CALL ACWK(IW,1,ERROR,*9999)
        CALL SGTEXT(INDEX,ISEG,ISTEXT(IW),IW,LIST,NTLIST,CSEG,TEXT,
     '    XWC_TEXT,YWC_TEXT,ERROR,*9999)
        CALL DAWK(IW,1,ERROR,*9999)

      ENDIF

      CALL EXITS('DRTEXT')
      RETURN
 9999 CALL ERRORS('DRTEXT',ERROR)
      CALL EXITS('DRTEXT')
      RETURN 1
      END


Module FE29
===========

C 30-3-98 MLB old UPGRID including GBS's old commands

      SUBROUTINE UPGRID(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NENQ,
     '  NHE,NKE,NPF,NPNE,NQGE,NQNE,NQXI,NQS,
     '  NRLIST,NVHE,NVJE,NW,NWQ,NXLIST,NXQ,
     '  CQ,CURVCORRECT,DNUDXQ,DXDXIQ,FEXT,GCHQ,GUQ,PG,PROPQ,SE,
     '  XA,XE,XP,XQ,XQD,XQDRC,YQ,ZA,ZE,ZG,ZP,STRING,ERROR,*)

C#### Subroutine: UPGRID
C###  Description:
C###    UPGRID updates grid point parameters, including material 
C###    values (conductivity tensors) and potentials and source terms
C###    for the bidomain model.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:b00.cmn'
      INCLUDE 'cmiss$reference:b01.cmn'
      INCLUDE 'cmiss$reference:call00.cmn'
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:grid00.cmn'
      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NHE(NEM,NXM),
     '  NKE(NKM,NNM,NBFM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NQGE(NGM,NEM,NBM),NQNE(NEM,NQEM),NQS(NEM),
     '  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NW(NEM,2),NWQ(6,0:NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  NVJE(NNM,NBFM,NJM,NEM),NXLIST(0:NXM)
      REAL*8 CQ(NMM,NQM,NXM),CURVCORRECT(2,2,NNM,NEM),DNUDXQ(3,3,NQM),
     '  DXDXIQ(3,3,NQM),FEXT(NIFEXTM,NGM,NEM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PG(NSM,NUM,NGM,NBM),XA(NAM,NJM,NQM),
     '  PROPQ(3,3,4,2,NQM,NXM),
     '  SE(NSM,NBFM,NEM),XQ(NJM,NQM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     '  XQD(NJM,NUM),XQDRC(NJM,NUM),
     '  YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,i,ib,IFNTYP,II,IJ,IK,iL,INFO,
     '  IFROMC,j,jb,k,mLL,mq,mRR,nb_extended,nbb,nc,ne,nee,neq,
     '  ngBL,ngBR,ngTL,ngTR,NG_row,ng_temp,nh,nhx,nii2,nij2,nik2,
     '  N3CO,ni,nii,nij,nij1,nik,NITB,nj,noelem,NU1(3),
     '  no_nrlist,NOQUES,nq,nq1,nq2,NQ_row,nr,nrr,nr_s,nr_t,nTT,nu,
     '  nx,nxc,nx_d,nx_s,nx_t,SCHEME
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  Current,DBM(3,3,3),DET,DIFFN,dNudXi(3,3),dPHIdX(3),
     '  dXidNu(3,3),d2PHIdX2(3,3),XI(3),DXIDXI(3),DXIX(3,3),
     '  CHTOFF(3,3,3),
     '  Esac,Ext,GX,GY,PHIMQ(-1:1,-1:1,-1:1),GL(3,3),GU(3,3),
     '  SUM,SUM1,SUM2,Xi1,Xi2,XQ_TEMP(3),X3G(4,3)
      CHARACTER TYPE*22
      LOGICAL ALL_REGIONS,CBBREV,DEFORMED,FILEIP
      DATA NU1/2,4,7/

      CALL ENTERS('UPGRID',*9999)

      nc=1 !temporary
      ICHAR=999

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update grid geometry
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###    This command calculates the geometric positions of the grid
C###    points and stores them in XQ. The grid is spaced evenly in
C###    material coordinates. The geometric positions are calculated
C###    by interpolation using the finite element basis functions used
C###    to describe the elements in which the grid points lie.

        OP_STRING(1)=STRING(1:IEND)//' geometry'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid geometry deformed
C###  Parameter:      <class #[1]>
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <from_class #[1]>
C###  Description:
C###    This command changes the geometric positions of the grid
C###    points when the underlying finite element mesh has undergone
C###    a deformation. Because the grid points don't move in material
C###    space, they will deform with the finite elements.

        OP_STRING(1)=STRING(1:IEND)//' geometry deformed'
        OP_STRING(2)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid metric
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###    This command generates the metric tensor information
C###    associated with a grid. This information is necessary
C###    for activation solutions.

        OP_STRING(1)=STRING(1:IEND)//' metric'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid metric deformed
C###  Parameter:      <class #[1]>
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <from_class #[1]>
C###  Description:
C###    This command updates the metric tensor information
C###    associated with a grid after a deformation has occured. The
C###    metric properties are updated for the deformed fibre angles.

        OP_STRING(1)=STRING(1:IEND)//' metric deformed'
        OP_STRING(2)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid material/strain/source
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' material/strain/source'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update grid time
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <class #[1]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' time'
        OP_STRING(2)=BLANK(1:IEND)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C obselete MLB 27 march 98
CC#### Command: FEM update grid intracellular_bidomain
CC###  Parameter:      <from_region #[2]>
CC###  Parameter:      <from_class #[2]>
CC###  Parameter:      <to_region #[1]>
CC###  Parameter:      <to_class #[1]>
CC###  Parameter:      <direct>
CC###  Description: 
CC###    Updates the collocation array YQ from
CC###    the extracellular potential.
C
C        OP_STRING(1)=STRING(1:IEND)//' intracellular_bidomain'
C        OP_STRING(2)=BLANK(1:IEND)//'<direct>'
C        OP_STRING(3)=BLANK(1:IEND)//'<from_region #[2]>'
C        OP_STRING(4)=BLANK(1:IEND)//'<from_class #[2]>'
C        OP_STRING(5)=BLANK(1:IEND)//'<to_region #[1]>'
C        OP_STRING(6)=BLANK(1:IEND)//'<to_class #[1]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C
C---------------------------------------------------------------------
C
CC#### Command: FEM update grid extracellular_bidomain
CC###  Parameter:      <from_region #[1]>
CC###  Parameter:      <from_class #[1>]
CC###  Parameter:      <to_region #[2]>
CC###  Parameter:      <to_class #[2]>
CC###  Description: 
CC###    Updates the collocation array YQ from
CC###    the 2nd derivative of intracellular potential.
C
C        OP_STRING(1)=STRING(1:IEND)//' extracellular_bidomain'
C        OP_STRING(2)=BLANK(1:IEND)//'<from_region #[1]>'
C        OP_STRING(3)=BLANK(1:IEND)//'<from_class #[1]>'
C        OP_STRING(4)=BLANK(1:IEND)//'<to_region #[2]>'
C        OP_STRING(5)=BLANK(1:IEND)//'<to_class #[2]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','UPGRID',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'MATERIAL',2,noco+1,NTCO,N3CO)) THEN
          TYPE='MATERIAL'
        ELSE IF(CBBREV(CO,'INTRACELLULAR_BIDOMAIN',1,noco+1,NTCO,N3CO))
     '    THEN
          TYPE='INTRACELLULAR_BIDOMAIN'
        ELSE IF(CBBREV(CO,'EXTRACELLULAR_BIDOMAIN',1,noco+1,NTCO,N3CO))
     '    THEN
          TYPE='EXTRACELLULAR_BIDOMAIN'
        ELSE IF(CBBREV(CO,'STRAIN',2,noco+1,NTCO,N3CO)) THEN
          TYPE='STRAIN'
        ELSE IF(CBBREV(CO,'TIME',1,noco+1,NTCO,N3CO)) THEN
          TYPE='TIME'
        ELSE IF(CBBREV(CO,'SOURCE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='SOURCE'
        ELSE IF(CBBREV(CO,'GEOMETRY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='GEOMETRY'
        ELSE IF(CBBREV(CO,'METRIC',2,noco+1,NTCO,N3CO)) THEN
          TYPE='METRIC'
        ENDIF

        IF(TYPE(1:22).EQ.'INTRACELLULAR_BIDOMAIN') THEN
          IF(CBBREV(CO,'FROM_REGION',6,noco+1,NTCO,N3CO)) THEN
            nr_s=IFROMC(CO(N3CO+1))
          ELSE
            nr_s=2
          ENDIF

          IF(CBBREV(CO,'FROM_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
      	  CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'Invalid source class',ERROR,*9999)

          CALL ASSERT(ITYP5(nr_s,nx_s).EQ.1.AND.ITYP2(nr_s,nx_s).EQ.5,
     '      'Source equation must be div(k.grad(u))=f',ERROR,*9999)

          IF(CBBREV(CO,'TO_REGION',4,noco+1,NTCO,N3CO)) THEN
            nr_t=IFROMC(CO(N3CO+1))
          ELSE
            nr_t=1
          ENDIF

          IF(CBBREV(CO,'TO_CLASS',4,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_t.GT.0,'>>Invalid target class',ERROR,*9999)
          CALL ASSERT(ITYP5(nr_t,nx_t).EQ.2.AND.ITYP2(nr_t,nx_t).EQ.9,
     '      'Target equation must be activation model',ERROR,*9999)
          
C old MLB 18 August 1997
C          IF(CBBREV(CO,'DIRECT',1,noco+1,NTCO,N3CO)) THEN
C            DIRECT=.TRUE.
C          ELSE
C            DIRECT=.FALSE.
C          ENDIF

        ELSE IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN
          IF(CBBREV(CO,'FROM_REGION',6,noco+1,NTCO,N3CO)) THEN
            nr_s=IFROMC(CO(N3CO+1))
          ELSE
            nr_s=1
          ENDIF

          IF(CBBREV(CO,'FROM_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
      	  CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'Invalid source class',ERROR,*9999)

          CALL ASSERT(ITYP5(nr_s,nx_s).EQ.2.AND.ITYP2(nr_s,nx_s).EQ.9,
     '      'Source equation must be activation model',ERROR,*9999)

          IF(CBBREV(CO,'TO_REGION',4,noco+1,NTCO,N3CO)) THEN
            nr_t=IFROMC(CO(N3CO+1))
          ELSE
            nr_t=2
          ENDIF

          IF(CBBREV(CO,'TO_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_t.GT.0,'>>Invalid target class',ERROR,*9999)
          CALL ASSERT(ITYP5(nr_t,nx_t).EQ.1.AND.ITYP2(nr_t,nx_t).EQ.5,
     '      'Target equation must be div(k.grad(u))=f',ERROR,*9999)

        ELSE IF((TYPE(1:8).EQ.'GEOMETRY').OR.(TYPE(1:6).EQ.
     '    'METRIC')) THEN
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          IF(CBBREV(CO,'DEFORMED',3,noco+1,NTCO,N3CO)) THEN
            DEFORMED=.TRUE.
            IF(CBBREV(CO,'FROM_CLASS',2,noco+1,NTCO,N3CO)) THEN
              nxc=IFROMC(CO(N3CO+1))
            ELSE
              nxc=1
            ENDIF
            CALL NX_LOC(NX_INQUIRE,nxc,nx_d,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx_d.NE.0,'Invalid class#',ERROR,*9999)
          ELSE
            DEFORMED=.FALSE.
            nx_d=0 !not used
          ENDIF 

        ELSE
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '       ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(ITYP4(nr,nx).EQ.4,
     '        'Must be collocation solution',ERROR,*9999)
          ENDDO
        ENDIF !type

        IF(TYPE(1:6).EQ.'SOURCE') THEN
          NOQUES=0
          FILEIP=.FALSE.
          IOTYPE=1
          
          FORMAT='('' Enter function to use as source term [1]:'''//
     '      '/''   (1) -2(pi)^2 sin(pi x) sin(pi y)'''//
     '      '/''   (2) Unused'''//
     '      '/''   (3) Unused'''//
     '      '/''   (4) Unused'''//
     '      '/''   (5) Unused'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).GT.0) THEN
            IFNTYP=IDATA(1)
          ELSE
            IFNTYP=1
          ENDIF
        ENDIF !source/bidomain

C old MLB 18 August 1997
C        nb=1
C        DO WHILE (nb.LE.NBT.AND.NBC(nb).NE.7)
C          nb=nb+1
C        ENDDO
C        CALL ASSERT((NBC(nb).EQ.7),
C     '    'Extended basis function not defined',ERROR,*9999)
! AJP 19/7/96 Moved into above loop
c        DO no_nrlist=1,NRLIST(0)
c          nr=NRLIST(no_nrlist)
c          CALL ASSERT(ITYP4(nr,nx).EQ.4,
c     '      'Must be collocation solution',ERROR,*9999)
c        ENDDO
C        NITB=NIT(nb)

C--type---------------------------------------------------------------

        IF(TYPE(1:8).EQ.'MATERIAL') THEN
C *** Compute Dij, the diffusion matrix, at each grid point in terms of
C *** xi coordinates.  This becomes a full matrix generated from the
C *** diagonal diffusion tensor in material coordinates.
C *** Similarly for the bidomain model, compute the conductivity
C *** tensors.

          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)

C old MLB 18 August 1997
C Find basis function for current region
C            nxielem=NIT(NBJ(1,NEELEM(1,nr)))
C            nb_extended=1
C            DO nbb=1,NBT
C              IF((nxielem.NE.NIT(nb_extended)).OR.(NBC(nb_extended)
C     '          .NE.7)) THEN
C                nb_extended=nb_extended+1
C              ENDIF
C            ENDDO
C            CALL ASSERT((NBC(nb_extended).EQ.7),
C     '        'Extended basis function not defined',ERROR,*9999)
C            NITB=NIT(nb_extended)            
 
            NITB=NQXI(0,NQS(NEELEM(1,nr)))
           
C Initialise entire arrays to start with
            DO i=1,3
              DO j=1,3
                dNudXi(i,j)=0.d0
                dXidNu(i,j)=0.d0
                DO nq=NQR(1,nr),NQR(2,nr)
                  DO k=1,4
                    PROPQ(i,j,k,1,nq,nx)=0.d0
                    PROPQ(i,j,k,2,nq,nx)=0.d0
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            
            IF(DOP) THEN
C$            call mp_setlock()
              WRITE(OP_STRING,'('' C(I)ij, C(E)ij:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$            call mp_unsetlock()
            ENDIF
            DO nq=NQR(1,nr),NQR(2,nr)
C *** First compute dnu/dxi = dnu/dx * dx/dxi, where both of these come
C *** from DXQ, computed in DEGRID.
            
              DO i=1,NITB
                DO j=1,NITB
!new PJH 11Jan96 
                  dNudXi(i,j)=0.d0
                  DO nj=1,NJT
                    dNudXi(i,j)=dNudXi(i,j)+
     '                          DNUDXQ(i,nj,nq)*DXDXIQ(nj,j,nq)
                  ENDDO
                ENDDO
              ENDDO
C *** and dxi/dnu = 1/(dnu/dxi)
              IF(NITB.EQ.1) THEN
                IF(dNudXi(1,1).GT.ZERO_TOL) dXidNu(1,1)=1.d0/dNudXi(1,1)
              ELSE
                CALL INVERT(NITB,dNudXi,dXidNu,DET)
              ENDIF
           
              IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.5) THEN
C *** For div(k.grad(u))=f the conductivities are in CQ(2..4,nq,nx)
                DO i=1,NITB
                  DO j=1,NITB
                    PROPQ(i,j,1,1,nq,nx)=
     '                CQ(2,nq,nx)*dXidNu(i,1)*dNudXi(1,j)
     '               +CQ(3,nq,nx)*dXidNu(i,2)*dNudXi(2,j)
                    IF(NITB.EQ.3) THEN
                      PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)
     '                 +CQ(4,nq,nx)*dXidNu(i,3)*dNudXi(3,j)
                    ENDIF
                  ENDDO !j
                ENDDO !i

              ELSE IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9) THEN
C *** For Bidomain
C ***   Compute and store C(I)ij(nq) and C(E)ij in PROPQ 
C ***     (already initialised)
C ***   Intracellular conductivities in CQ(3;4;5,nq,nx)
C ***   Extracellular conductivities in CQ(6;7;8,nq,nx)
              
                DO i=1,NITB
                  DO j=1,NITB
                    DO k=1,NITB
                      PROPQ(i,j,1,1,nq,nx)=PROPQ(i,j,1,1,nq,nx)+
     '                  CQ(k+2,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                      PROPQ(i,j,1,2,nq,nx)=PROPQ(i,j,1,2,nq,nx)+
     '                  CQ(k+5,nq,nx)*dXidNu(i,k)*dNudXi(k,j)
                    ENDDO
                  ENDDO !j
                ENDDO !i
              ENDIF !ityp4
           
              IF(DOP) THEN
C$              call mp_setlock()
                WRITE(OP_STRING,'('' nq: '',I5)') nq
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            ((PROPQ(i,j,1,1,nq,nx),j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(9(F12.5))')
     '            ((PROPQ(i,j,1,2,nq,nx),j=1,nitb),i=1,nitb)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$              call mp_unsetlock()
              ENDIF
            ENDDO !nq
           
C *** Now compute Dij,k using a first order finite difference about
C *** each grid point.  This only needs to be computed for grid points
C *** that are internal (not on external bdy).
C *** Similarly for the bidomain model, the derivatives of the
C *** conductivity tensors.

            IF(DOP) THEN
C$            call mp_setlock()
              WRITE(OP_STRING,'('' C(I)ij,k, C(E)ij,k:'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$            call mp_unsetlock()
            ENDIF
            DO nq=NQR(1,nr),NQR(2,nr)
              IF(NWQ(1,nq,1).EQ.0) THEN !internal g.p.
                DO k=1,NITB
                  DO i=1,NITB
                    DO j=1,NITB
                      IF(NXQ(k,1,nq,1).GT.0) THEN
                        PROPQ(i,j,k+1,1,nq,nx)=
     '                               PROPQ(i,j,1,1,NXQ( k,1,nq,1),nx)
     '                              -PROPQ(i,j,1,1,NXQ(-k,1,nq,1),nx)
                        PROPQ(i,j,k+1,2,nq,nx)=
     '                               PROPQ(i,j,1,2,NXQ( k,1,nq,1),nx)
     '                              -PROPQ(i,j,1,2,NXQ(-k,1,nq,1),nx)
                      ELSE
                        PROPQ(i,j,k+1,1,nq,nx)=0.d0
                        PROPQ(i,j,k+1,2,nq,nx)=0.d0
                      ENDIF
                    ENDDO !j
                  ENDDO !i
                ENDDO !k
                IF(DOP) THEN
C$                call mp_setlock()
                  WRITE(OP_STRING,'('' nq: '',I5)') nq
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(9(F12.5))')
     '           (((PROPQ(i,j,k,1,nq,nx),k=2,nitb+1),j=1,nitb),i=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(9(F12.5))')
     '           (((PROPQ(i,j,k,2,nq,nx),k=2,nitb+1),j=1,nitb),i=1,nitb)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C$                call mp_unsetlock()
                ENDIF
              ENDIF !nwq
            ENDDO !nq
          ENDDO !no_nrlist

          UP_GRID_MATERIAL=.TRUE.

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:22).EQ.'INTRACELLULAR_BIDOMAIN') THEN

          CALL ASSERT(.FALSE.,'This code is no longer supported',
     '      ERROR,*9999)

C MLB 18 August 1997
C This block of code has been commented out because it was written
C to solve the bidomain equations when it was thought that the 
C transmembrane and extracellular potentials had different time and
C space scales. Bidomain problems can now be solved by selecting
C the implicit grid solver.

C Update the extracellular potentials at intracellular grid points 
C either from extracellular grid 
C or by interpolating from values at fem node points
C
C          IF(ITYP4(nr_s,nx_s).EQ.4) THEN !source class is grid
C            DO nq=1,NQT !Update extracellular potential
C              YQ(nq,2,1,nx_t)=YQ(nq,1,1,nx_s)
C            ENDDO
C          
C          ELSE IF(DIRECT) THEN         !source is grid mesh
C            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr_s,nx_s),
C     '        NKH(1,1,1,nr_s),NPNODE,nr_s,NVHP(1,1,1,nr_s),nx_s,
C     '        NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C            DO nq=1,NQT !Update extracellular potential
C              YQ(nq,2,1,nx_t)=ZP(1,1,1,NPQ(nq),1)
C            ENDDO
C          
C          ELSE                         !source class is fem mesh
C            nb_extended=1
C            DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
C              nb_extended=nb_extended+1
C            ENDDO
C            CALL ASSERT((NBC(nb_extended).EQ.7),
C     '        'Extended basis function not defined',ERROR,*9999)
C
C! Check we have a 9x9(x9) grid point scheme
C GBS 3-July-96 Why????
C           NINE=.TRUE.
C           DO ni=1,NIT(nb_extended)
C             NINE=NINE.AND.(NGAP(ni,nb_extended).EQ.9)
C           ENDDO
C           CALL ASSERT(NINE,'Require 9x9(x9) grid',ERROR,*9999)
C            
C            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr_s,nx_s),
C     '        NKH(1,1,1,nr_s),NPNODE,nr_s,NVHP(1,1,1,nr_s),nx_s,
C     '        NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C            DO noelem=1,NEELEM(0,nr_s)
C              ne=NEELEM(noelem,nr_s)
C              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx_s),NKE(1,1,1,ne),
C     '          NPF,NPNE(1,1,ne),nr_s,NVHE(1,1,1,ne),NW(ne,1),nx_s,
C     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
C     '          ERROR,*9999)
C              DO ng=1,NGT(nb_extended)
C                nq=NQGE(ng,ne,nb_extended)
C                SUM=0.0d0
C                DO ns=1,NST(nb_extended)
C                  PGG=PG(ns,1,ng,nb_extended)
C                  SUM=SUM+PGG*ZE(ns,1)
C                ENDDO !ns
C                YQ(nq,2,1,nx_t)=SUM !Update extracellular potential
C              ENDDO !ng
C            ENDDO !noelem
C          ENDIF !ityp4
            
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN
C Update the RHS source term for the extracellular bidomain Poisson's
C eqn (class nx_t) by computing div(kgrad(phi_m)) at the extracellular 
C grid pts from a quad patch of intracellular grid pts on class nx_s.

C *** Loop over grid points
          DO nq=1,NQT
            DO nik=-1,1
              DO nij=-1,1
                DO nii=-1,1
                  PHIMQ(nii,nij,nik)=0.d0
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
          
C ***  Formulate a local quadratic element about nq defined by UMQ,
C ***  storing the value of phi(m) at each node point mq.
            IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
            IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                  IF(mq.GT.0) THEN
                    PHIMQ(nii,nij,nik)=YQ(mq,1,1,nx_s)
                  ELSE
                    PHIMQ(nii,nij,nik)=0.0d0 
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
            
C *** Compute PHI,k by 1st order finite differences about nq.
            dPHIdX(1)=PHIMQ(1,0,0)-PHIMQ(-1,0,0)
            dPHIdX(2)=PHIMQ(0,1,0)-PHIMQ(0,-1,0)
            dPHIdX(3)=PHIMQ(0,0,1)-PHIMQ(0,0,-1)

C *** Compute PHI,jk by 2nd order finite differences about nq
            DO ib=1,NITB
              DO jb=1,NITB
                IF(ib.EQ.jb) THEN
                  IF(ib.EQ.1) THEN !PHI,11
                    d2PHIdX2(ib,jb)=(PHIMQ(1,0,0)-
     '                2.0d0*PHIMQ(0,0,0)+PHIMQ(-1,0,0))*4.0d0
                  ELSE IF(ib.EQ.2) THEN !PHI,22
                    d2PHIdX2(ib,jb)=(PHIMQ(0,1,0)
     '                -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,-1,0))*4.0d0
                  ELSE IF(ib.EQ.3) THEN !PHI,33
                    d2PHIdX2(ib,jb)=(PHIMQ(0,0,1)
     '                -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,0,-1))*4.0d0
                  ENDIF               
                ELSE
                  IF(ib+jb.EQ.3) THEN !PHI,12 and PHI,21
                    d2PHIdX2(ib,jb)=PHIMQ(1,1,0)-PHIMQ(1,-1,0)
     '                -PHIMQ(-1,1,0)+PHIMQ(-1,-1,0)
                  ELSE IF(ib+jb.EQ.4) THEN !PHI,13 and PHI,31
                    d2PHIdX2(ib,jb)=PHIMQ(1,0,1)-PHIMQ(1,0,-1)
     '                -PHIMQ(-1,0,1)+PHIMQ(-1,0,-1)
                  ELSE IF(ib+jb.EQ.5) THEN !PHI,23 and PHI,32
                    d2PHIdX2(ib,jb)=PHIMQ(0,1,1)-PHIMQ(0,-1,1)
     '                -PHIMQ(0,1,-1)+PHIMQ(0,-1,-1)
                  ENDIF
                ENDIF !ib=jb
              ENDDO !jb
            ENDDO !ib

C *** Compute the diffusion in three parts.
            DIFFN=0.0d0
            DO nii=1,NITB
              DO nij=1,NITB           
                SUM1=0.0d0
                SUM2=0.0d0
                DO nik=1,NITB
                  SUM1=SUM1
     '                +PROPQ(nik,nii,nij+1,1,nq,nx_s)*dPHIdX(nik)
                  SUM2=SUM2
     '                +PROPQ(nik,nii,1,1,nq,nx_s)*d2PHIdX2(nij,nik)
                ENDDO
                DIFFN=DIFFN+(SUM1+SUM2)*GUQ(nii,nij,nq)
              ENDDO !nij
            ENDDO !nii
            DO nii=1,NITB
              SUM1=0.0d0
              DO nij=1,NITB
                SUM1=SUM1+PROPQ(nij,nii,1,1,nq,nx_s)*dPHIdX(nij)
              ENDDO
              DIFFN=DIFFN-SUM1*GCHQ(nii,nq)
            ENDDO !nii
            CQ(1,nq,nx_t)=-DIFFN !is source term for Poisson's eqn
          ENDDO !nqt

C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:6).EQ.'STRAIN') THEN

          WRITE(OP_STRING,'(''Code not currently supported, '
     '      //'results are unreliable'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)

            Esac=-20.d0 !(mV) reversal potential for SAC
 
C old MLB 18 August 1997         
          nb_extended=1
          DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
            nb_extended=nb_extended+1
          ENDDO
          CALL ASSERT((NBC(nb_extended).EQ.7),
     '      'Extended basis function not defined',ERROR,*9999)
          
C Generate inward Isac current from stretch if above threshold
C Do interior collocation pts first, then bdry pts
          
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(DOP) WRITE(*,'(/'' Interior pts: ne='',I5)') ne
          
              NQ_row=9
              NG_row=3
              
              DO j=2,NQ_row-1    !loop over collocation points
                DO i=2,NQ_row-1
                  nq=NQGE(i+(j-1)*NQ_row,ne,nb_extended) !is coll.n pt#
              
                  mLL=INT(REAL(i+1)/3.0)  !m(left   Gauss pt index)
                  nBB=INT(REAL(j+1)/3.0)  !n(bottom Gauss pt index)
                  mRR=min(mLL+1,NG_row)    !m(right  Gauss pt index)
                  nTT=min(nBB+1,NG_row)    !n(top    Gauss pt index)
              
                  ngBL=mLL+(nBB-1)*NG_row  !ng(bottom left  Gauss pt#)
                  ngBR=mRR+(nBB-1)*NG_row  !ng(bottom right Gauss pt#)
                  ngTL=mLL+(nTT-1)*NG_row  !ng(top left     Gauss pt#)
                  ngTR=mRR+(nTT-1)*NG_row  !ng(top right    Gauss pt#)
              
!                 WRITE(*,'(/'' i='',I1,'' j='',I1,'' nq='',I2)') i,j,nq
!                 WRITE(*,'('' mLL='',I1,'' nBB='',I1)') mLL,nBB
!                 WRITE(*,'('' ngBL='',I1,'' ngBR='',I1,'' ngTL='',I1,'
!    '              //''' ngTR='',I1)') ngBL,ngBR,ngTL,ngTR
              
                  iL=3*mLL-1  !locates bottom left Gauss pt index
                  jB=3*nBB-1  ! in the collocation grid
              
                  Xi1=DBLE(i-iL)/DBLE(NG_row) !how far ij colloc.n pt
                  Xi2=DBLE(j-jB)/DBLE(NG_row) ! is betw adj Gauss pts
              
                  Ext=(1.d0-Xi1)*(1.d0-Xi2)*FEXT(1,ngBL,ne) 
     '               +      Xi1 *(1.d0-Xi2)*FEXT(1,ngBR,ne)
     '               +(1.d0-Xi1)*      Xi2 *FEXT(1,ngTL,ne)
     '               +      Xi1 *      Xi2 *FEXT(1,ngTR,ne)
                  IF(DOP) THEN
C$                  call mp_setlock()
                    WRITE(*,'(/'' Xi1='',F5.2,'' Xi2='',F5.2)') Xi1,Xi2
                    WRITE(*,'('' Exten. ratio at nq='',I3,'
     '               //''' :'',F6.3)') nq,Ext
C$                  call mp_unsetlock()
                  ENDIF
          
C PJH 30jun96 - nx added. Check with Greg
                  Current=-(YQ(nq,1,1,nx)-Esac)*CQ(18,nq,nx)*(Ext-1.d0)
                  IF(Current.GT.100.d0) THEN !above threshold
                    YQ(nq,5,1,nx)=Current
                  ELSE
                    YQ(nq,5,1,nx)=0.d0
                  ENDIF
              
                ENDDO !i
              ENDDO !j
            ENDDO !noelem
          
C         Boundary collocation pts
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(DOP) WRITE(*,'(/'' Boundary pts: ne='',I5)') ne
          
C         RH bdry of element ne
              DO j=1,9 
                nq =NQGE(9*j,ne,nb_extended)
                nq1=NXQ(-1,1,nq,1)
                nq2=NXQ( 1,1,nq,1)
                IF(DOP) THEN
C$                call mp_setlock()
                  WRITE(*,'(/'' nq='',I3,'' nq1='',I3,'' nq2='',I3)')
     '              nq,nq1,nq2
C$                call mp_unsetlock()
                ENDIF
                IF(nq2.GT.0) THEN !RH point exists
                  Current=0.5d0*(YQ(nq1,5,1,nx)+YQ(nq2,5,1,nx))
                  IF(Current.GT.100.d0) THEN !above threshold
                    YQ(nq,5,1,nx)=Current
                  ELSE
                    YQ(nq,5,1,nx)=0.d0
                  ENDIF
                ENDIF
              ENDDO !j
          
C         TOP bdry of element ne
              DO i=1,9 
                nq =NQGE(72+i,ne,nb_extended)
                nq1=NXQ(-2,1,nq,1)
                nq2=NXQ( 2,1,nq,1)
                IF(DOP) THEN
C$                call mp_setlock()
                  WRITE(*,'(/'' nq='',I3,'' nq1='',I3,'' nq2='',I3)')
     '              nq,nq1,nq2
C$                call mp_unsetlock()
                ENDIF
                IF(nq2.GT.0) THEN !TOP point exists
                  Current=0.5d0*(YQ(nq1,5,1,nx)+YQ(nq2,5,1,nx))
                  IF(Current.GT.100.d0) THEN !above threshold
                    YQ(nq,5,1,nx)=Current
                  ELSE
                    YQ(nq,5,1,nx)=0.d0
                  ENDIF
                ENDIF
              ENDDO !i
            ENDDO !noelem
          ENDDO !no_nrlist
          
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:4).EQ.'TIME') THEN
C Update previous time step solution for multigrid transient heat eqtn 
          DO nq=1,NQT
            YQ(nq,6,1,nx)=YQ(nq,1,1,nx)
          ENDDO !nq
          
C--type---------------------------------------------------------------

        ELSE IF(TYPE(1:6).EQ.'SOURCE') THEN
C Calculate source term for collocation grid
          DO nq=1,NQT
            GX=XQ(1,nq)
            GY=XQ(2,nq)
            IF(IFNTYP.EQ.1) THEN
              CQ(1,nq,nx)=-2.d0*PI*PI*DSIN(PI*GX)*DSIN(PI*GY)
            ENDIF !ifntyp
          ENDDO !nq

C--type---------------------------------------------------------------
C*** Martin Buist, June 1997

        ELSE IF(TYPE(1:8).EQ.'GEOMETRY') THEN
! Calculate geometric coordinates of grid points
          CALL ASSERT(CALL_GRID,'>>No grid defined',ERROR,*9999)
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)

! Initialise XQ
            DO nq=NQR(1,nr),NQR(2,nr)
              DO nj=1,NJM
                XQ(nj,nq)=0.0d0
              ENDDO
            ENDDO

            DO nee=1,NEELEM(0,nr)
              ne=NEELEM(nee,nr)
              SCHEME=NQS(ne)

! Map global parameters to local element parameters
              CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,
     '          XP,ERROR,*9999)
              IF(DEFORMED) THEN
                CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx_d),NKE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '            NW(ne,1),nx_d,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '            ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              ENDIF

              II=MAX(1,NQXI(1,SCHEME))
              IJ=1
              IK=1
              IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME))
              IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))
              DO i=1,3
                XQ_TEMP(i)=0.0d0
                XI(i)=0.0d0
              ENDDO !i
! Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME))
                    IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)*
     '                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq)
! Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN 
! Local xi coordinates of grid point nq in element ne
                      IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                      IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                      IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)

! Use basis function interpolation to get geometric position

C This is old, changing to use XEXW which is more consistent with 
C GBS's IPGRID
C                      nb=NQSCNB(SCHEME)
C                      DO nj=1,NJT
C                        XQ_TEMP(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
C     '                    INP(1,1,nb),nb,1,XI,XE(1,nj))
C                      ENDDO !nj

                      IF(DEFORMED) THEN
                        CALL ZEZW(0,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '                    NHE(ne,nx_d),nr,nx_d,DXIX,ZE,ZG,XI,
     '                    ERROR,*9999)
                        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                          nj=NJ_LOC(NJL_GEOM,nhx,nr)
                          nh=NH_LOC(nhx,nx_d)
                          DO nu=1,NUT(NBH(nh,nc,ne))
                            XQD(nj,nu)=ZG(nhx,nu)
                          ENDDO !nu
                        ENDDO !nh/nj
                      ELSE
                        CALL XEXW(IBT,IDO,INP,NAN,NBJ(1,ne),
     '                    nr,XE,XQD,XI,ERROR,*9999)
                      ENDIF

                      DO nj=1,NJT
                        XQ_TEMP(nj)=XQD(nj,1)
                      ENDDO !nj

! Change coordinates - all grid points in rectangular cartesian
                      CALL XZ(ITYP10(nr),XQ_TEMP,XQ(1,nq))
                    ENDIF
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
            ENDDO !element
          ENDDO !region
          IF(DOP) THEN
            DO nq=1,NQT
              WRITE(OP_STRING,'(''XQ'',I6,3F10.4)') nq,XQ(1,nq),
     '          XQ(2,nq),XQ(3,nq)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

C--type---------------------------------------------------------------
C*** Martin Buist, June 1997

        ELSE IF(TYPE(1:8).EQ.'METRIC') THEN
! Calculate metric tensor information for grid
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            DO i=1,3
              DO j=1,3
                GU(i,j)=0.0d0
                GL(i,j)=0.0d0
              ENDDO !j
            ENDDO !i
            DO nee=1,NEELEM(0,nr)
              ne=NEELEM(nee,nr)
              SCHEME=NQS(ne)
! Map global parametres to local element parameters
              CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA,XE,
     '          XP,ERROR,*9999)
              IF(DEFORMED) THEN
                CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx_d),NKE(1,1,1,ne),
     '            NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '            NW(ne,1),nx_d,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '            ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
              ENDIF

! dxi/dxi is the relationship between material space and the local 
!   quadratic element space. dxi(quad)/dxi(mate)
              DO ni=1,NQXI(0,SCHEME)
                DXIDXI(ni)=2.0d0/DBLE(NQXI(ni,SCHEME)-1)
                IF(DOP) THEN
                  WRITE(OP_STRING,'(''dxidxi(ni)'',I3,F12.6)') ni,
     '              DXIDXI(ni)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !ni
              II=MAX(1,NQXI(1,SCHEME))
              IJ=1
              IK=1
              IF(NQXI(0,SCHEME).GT.1) IJ=MAX(1,NQXI(2,SCHEME))
              IF(NQXI(0,SCHEME).GT.2) IK=MAX(1,NQXI(3,SCHEME))
              DO i=1,3
                XI(i)=0.0d0
              ENDDO !i
! Loop over the grid points in each element
              DO nik=1,IK
                DO nij=1,IJ
                  DO nii=1,II
                    neq=nii+((nij-1)*NQXI(1,SCHEME))
                    IF(NQXI(0,SCHEME).GT.1) neq=neq+((nik-1)*
     '                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,neq)
! Evaluate each grid point only once
                    IF(NENQ(1,nq).EQ.ne) THEN 
! Local xi coordinates of grid point nq in element ne
                      IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                      IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                      IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)

                      IF(DEFORMED) THEN
                        CALL ZEZW(0,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '                    NHE(ne,nx_d),nr,nx_d,DXIX,ZE,ZG,XI,
     '                    ERROR,*9999)
                        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                          nj=NJ_LOC(NJL_GEOM,nhx,nr)
                          nh=NH_LOC(nhx,nx_d)
                          DO nu=1,NUT(NBH(nh,nc,ne))
                            XQD(nj,nu)=ZG(nhx,nu)
                          ENDDO !nu
                        ENDDO !nh/nj
                      ELSE
                        CALL XEXW(IBT,IDO,INP,NAN,NBJ(1,ne),nr,
     '                    XE,XQD,XI,ERROR,*9999)
                      ENDIF

! Gi (lower) = dx(nj)/dxi(ni) in direction nj
                      DO ni=1,NQXI(0,SCHEME)
                        CALL XZ_DERIV(ITYP10(nr),NU1(ni),XQD,XQDRC)
                        DO nj=1,NJT
                          DXDXIQ(nj,ni,nq)=XQDRC(nj,NU1(ni))*DXIDXI(ni)
                          IF(DOP) THEN
                            WRITE(OP_STRING,'(''dxdxiq'',3I6,F12.6)')
     '                        nj,ni,nq,DXDXIQ(nj,ni,nq)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDDO !nj
                      ENDDO !ni
! Calculate Gij (lower) and then Gij (upper) by inversion 
! Gij (lower) = dx(nj)/dxi(nii)*dx(nj)/dxi(nij)
                      DO nij2=1,NQXI(0,SCHEME)
                        DO nii2=1,NQXI(0,SCHEME)
                          GL(nii2,nij2)=0.0d0
                          DO nj=1,NJT
                            GL(nii2,nij2)=GL(nii2,nij2)+(DXDXIQ(nj,
     '                        nii2,nq)*DXDXIQ(nj,nij2,nq))
                          ENDDO !nj
                        ENDDO !nii2
                      ENDDO !nij2
                      IF(DABS(GL(1,1)+GL(2,1)+GL(3,1)).LT.ZERO_TOL) THEN
                        DO nij2=1,NQXI(0,SCHEME)
                          DO nii2=1,NQXI(0,SCHEME)
                            GUQ(nii2,nij2,nq)=0.0d0
                          ENDDO
                          GCHQ(nij2,nq)=0.0d0
                        ENDDO
                        DO nij1=1,NJT
                          DO nij2=1,NJT
                            DNUDXQ(nij1,nij2,nq)=0.0d0
                          ENDDO !nij2
                        ENDDO !nij1
                      ELSE
! Gij (upper) = inv(Gij lower)
                        IF(NQXI(0,SCHEME).GT.1) THEN
                          CALL INVERT(NQXI(0,SCHEME),GL,GU,DET)
                        ELSE IF(GL(1,1).NE.0.0D0) THEN
                          GU(1,1)=1.0D0/GL(1,1)
                        ELSE
                          WRITE(OP_STRING,
     '                      '(''>Zero value for GL - cannot invert'')') 
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ENDIF
! Write Gij (upper) at grid point nq into GUQ 
                        DO nij2=1,NQXI(0,SCHEME)
                          DO nii2=1,NQXI(0,SCHEME)
                            GUQ(nii2,nij2,nq)=GU(nii2,nij2)
                            IF(DOP) THEN
                              WRITE(OP_STRING,'(''guq'',3I6,F12.6)')
     '                          nii2,nij2,nq,GUQ(nii2,nij2,nq)
                              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                            ENDIF
                          ENDDO !nij2
                        ENDDO !nii2
! Calculate the Christoffel symbol (CHTOFF)
                        CALL TOFFEL(NBJ(1,ne),nr,CHTOFF,DBM,GU,
     '                    XQD,X3G,.FALSE.,ERROR,*9999)
! Calculate CHTOFF(k,i,j)*Gij (upper) = GCHQ(k)
                        DO nik2=1,NQXI(0,SCHEME)
                          SUM=0.0d0
                          DO nii2=1,NQXI(0,SCHEME)
                            DO nij2=1,NQXI(0,SCHEME)
                              SUM=SUM+(CHTOFF(nik2,nii2,nij2)*
     '                          GU(nii2,nij2))*DXIDXI(nii2)*DXIDXI(nij2)
                            ENDDO !nij2
                          ENDDO !nii2
                          GCHQ(nik2,nq)=SUM*DXIDXI(nik2)
                          IF(DOP) THEN 
                            WRITE(OP_STRING,'(''gchq'',2I6,F12.6)') 
     '                        nik2,nq,GCHQ(nik2,nq)
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                        ENDDO !nik2
                        IF(CALL_FIBR) THEN
! Calculate dnu/dx (direction cosines of material coords)
                          IF(DEFORMED) THEN
                            ng_temp=0 ! to compute at xi position
                            CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,NBH,
     '                        NBJ(1,ne),ng_temp,NHE(ne,nx_d),nr,nx_d,
     '                        A_VECTOR,B_VECTOR,C_VECTOR,PG,XE,XQD,XI,
     '                        ZE,ZG,ERROR,*9999)
                          ELSE
                            CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne),
     '                        nr,A_VECTOR,B_VECTOR,C_VECTOR,XE,XQD,XI,
     '                        ERROR,*9999)
                          ENDIF

                          DO nj=1,NJT
                            DNUDXQ(1,nj,nq)=A_VECTOR(nj)
                            IF(NJT.GT.1) DNUDXQ(2,nj,nq)=B_VECTOR(nj)
                            IF(NJT.GT.2) DNUDXQ(3,nj,nq)=C_VECTOR(nj)
                          ENDDO !nj
                        ELSE
                          DO nij1=1,NJT
                            DO nij2=1,NJT
                              IF(nij1.EQ.nij2) THEN
                                DNUDXQ(nij1,nij2,nq)=1.0d0
                              ELSE
                                DNUDXQ(nij1,nij2,nq)=0.0d0
                              ENDIF
                            ENDDO !nij2
                          ENDDO !nij1
                        ENDIF
                      ENDIF
                      DO nj=1,NJT
                        IF(DOP) THEN
                          WRITE(OP_STRING,'(''dnu/dxq'',3F12.6)') 
     '                      DNUDXQ(1,nj,nq),DNUDXQ(2,nj,nq),
     '                      DNUDXQ(3,nj,nq)
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
            ENDDO !element
          ENDDO !region

          UP_GRID_TENSOR=.TRUE.

        ENDIF !type
      ENDIF

      CALL EXITS('UPGRID')
      RETURN
 9999 CALL ERRORS('UPGRID',ERROR)
      CALL EXITS('UPGRID')
      RETURN 1
      END

C 30-3-98 MLB removed part of upgaus

C---------------------------------------------------------------------

C#### Command: FEM update gauss extracellular_bidomain
C###  Parameter:      <region (#s/all)[1]>
C###  Parameter:      <basis #[1]>
C###  Parameter:      <from_class #[1]>
C###  Parameter:      <to_class #[2]>
C###  Description: 
C###    Computes div(sigma_i grad V_m) i.e. the rhs for 
C###    the poisson equation.

        OP_STRING(1)=STRING(1:IEND)//' extracellular_bidomain'
        OP_STRING(2)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<from_class #[1]>'
        OP_STRING(5)=BLANK(1:IEND)//'<to_class #[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

        ELSE IF(CBBREV(CO,'EXTRACELLULAR_BIDOMAIN',1,noco+1,NTCO,N3CO))
     '    THEN
          TYPE='EXTRACELLULAR_BIDOMAIN'


        ELSE IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN
          IF(CBBREV(CO,'FROM_CLASS',1,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
      	  CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_s.NE.0,'>>Invalid source class',ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(ITYP5(nr,nx_s).EQ.2.AND.ITYP2(nr,nx_s).EQ.9,
     '        '>>Source equation must be activation model',ERROR,*9999)
          ENDDO

          ELSE IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN

            CALL ASSERT(.FALSE.,'This code is no longer supported',
     '        ERROR,*9999)

C MLB 18 August 1997
C This block of code has been commented out because it was written
C to solve the bidomain equations when it was thought that the 
C transmembrane and extracellular potentials had different time and
C space scales. Bidomain problems can now be solved by selecting
C the implicit grid solver.


C   Update the right hand side source term for the bidomain poissons
C   equations by computing div kgrad phi_m at the gauss points
C   in each element from the grid point scheme.
C   Now no longer requires certain#s of ng or nq GBS 3-7-96
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              nb=NBJ(1,ne)
C              DO ng=1,NGT(nb)
C                nq=NQGE(ng,ne,nb) !closest grid point
C                DO nik=-1,1
C                  DO nij=-1,1
C                    DO nii=-1,1
C                      PHIMQ(nii,nij,nik)=0.d0
C                    ENDDO !nii
C                  ENDDO !nij
C                ENDDO !nik
C
C***  Formulate a local quadratic element about nq defined by UMQ,
C***  storing the value of phi(m) at each node point mq.
C                NITB=NIT(nb)
C                IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
C                IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
C                DO nik=-IK,IK
C                  DO nij=-IJ,IJ
C                    DO nii=-1,1
C                      mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
C                      IF(mq.GT.0) THEN
C                        PHIMQ(nii,nij,nik)=YQ(mq,1,1,nx_s)
C                      ELSE
C                        PHIMQ(nii,nij,nik)=0.0d0
C                      ENDIF
C                    ENDDO !nii
C                  ENDDO !nij
C                ENDDO !nik
C *** Compute PHI,k by first order finite differences about nq.
C                dPHIdX(1)=PHIMQ(1,0,0)-PHIMQ(-1,0,0)
C                dPHIdX(2)=PHIMQ(0,1,0)-PHIMQ(0,-1,0)
C                dPHIdX(3)=PHIMQ(0,0,1)-PHIMQ(0,0,-1)
C *** Compute PHI,jk by second order f.d. about nq
C                DO ib=1,NITB
C                  DO jb=1,NITB
C                    IF(ib.EQ.jb) THEN
C                      IF(ib.EQ.1) THEN !PHI,11
C                        d2PHIdX2(ib,jb)=(PHIMQ(1,0,0)-
C     '                    2.0d0*PHIMQ(0,0,0)+PHIMQ(-1,0,0))*4.0d0
C                      ELSE IF(ib.EQ.2) THEN !PHI,22
C                        d2PHIdX2(ib,jb)=(PHIMQ(0,1,0)
C     '                    -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,-1,0))*4.0d0
C                      ELSE IF(ib.EQ.3) THEN !PHI,33
C                        d2PHIdX2(ib,jb)=(PHIMQ(0,0,1)
C     '                    -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,0,-1))*4.0d0
C                      ENDIF               
C                    ELSE
C                      IF(ib+jb.EQ.3) THEN !PHI,12 and PHI,21
C                        d2PHIdX2(ib,jb)=PHIMQ(1,1,0)-PHIMQ(1,-1,0)
C     '                    -PHIMQ(-1,1,0)+PHIMQ(-1,-1,0)
C                      ELSE IF(ib+jb.EQ.4) THEN !PHI,13 and PHI,31
C                        d2PHIdX2(ib,jb)=PHIMQ(1,0,1)-PHIMQ(1,0,-1)
C     '                    -PHIMQ(-1,0,1)+PHIMQ(-1,0,-1)
C                      ELSE IF(ib+jb.EQ.5) THEN !PHI,23 and PHI,32
C                        d2PHIdX2(ib,jb)=PHIMQ(0,1,1)-PHIMQ(0,-1,1)
C     '                    -PHIMQ(0,1,-1)+PHIMQ(0,-1,-1)
C                      ENDIF
C                    ENDIF
C                  ENDDO !jb
C                ENDDO !ib
! Compute the diffusion in three parts.
C                DIFFN=0.0d0
C                DO nii=1,NITB
C                  DO nij=1,NITB           
C                    SUM1=0.0d0
C                    SUM2=0.0d0
C                    DO nik=1,NITB
C                      SUM1=SUM1+PROPQ(nik,nii,nij+1,1,nq,nx_s)
C     '                         *dPHIdX(nik)
C                      SUM2=SUM2+PROPQ(nik,nii,1,1,nq,nx_s)
C     '                         *d2PHIdX2(nij,nik)
C                    ENDDO
C                    DIFFN=DIFFN+(SUM1+SUM2)*GUQ(nii,nij,nq)
C                  ENDDO
C                ENDDO
C                DO nii=1,NITB
C                  SUM1=0.0d0
C                  DO nij=1,NITB
C                    SUM1=SUM1+PROPQ(nij,nii,1,1,nq,nx_s)*dPHIdX(nij)
C                  ENDDO
C                  DIFFN=DIFFN-SUM1*GCHQ(nii,nq)
C                ENDDO
C OLD AJP/GBS Should be -ve (GBS Thesis, 3.3.2)   YG(ng,1,ne)=DIFFN
C Try again - I think it should be +ve            YG(ng,1,ne)=-DIFFN
C                YG(ng,1,ne)=DIFFN
C              ENDDO !ng
C            ENDDO !ne

C MLB 30-3-98 archiving upsour

      SUBROUTINE UPSOUR(NBJ,NEELEM,NPQ,NXQ,
     '  CP,GCHQ,GUQ,PROPQ,YQ,STRING,ERROR,*)

C#### Subroutine: UPSOUR
C###  Description:
C###    UPSOUR updates the source terms of the extracellular bidomain 
C###    from the intracellular solution.

      IMPLICIT NONE
      INCLUDE 'cmiss$reference:cbdi02.cmn'
      INCLUDE 'cmiss$reference:cbdi10.cmn'
      INCLUDE 'cmiss$reference:geom00.cmn'
      INCLUDE 'cmiss$reference:ityp00.cmn'
      INCLUDE 'cmiss$reference:loc00.inc'
      INCLUDE 'cmiss$reference:mxch.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NPQ(NQM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 CP(NMM,NPM,NXM),GCHQ(3,NQM),GUQ(3,3,NQM),
     '  PROPQ(3,3,4,2,NQM,NXM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,
     '  ib,IJ,IK,jb,mq,nb,nii,nij,nik,NITB,N3CO,
     '  nq,nr_s,nr_t,nxc,nx_s,nx_t
      REAL*8 DIFFN,dPHIdX(3),d2PHIdX2(3,3),PHIMQ(-1:1,-1:1,-1:1),
     '  SUM1,SUM2
      CHARACTER TYPE*22
      LOGICAL CBBREV

      CALL ENTERS('UPSOUR',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update source extracellular_bidomain
C###  Parameter:      <from_region #[1]>
C###  Parameter:      <from_class #[1]>
C###  Parameter:      <to_region #[2]>
C###  Parameter:      <to_class #[2]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' extracellular_bidomain'
        OP_STRING(2)=BLANK(1:IEND)//'<from_region #[1]>'
        OP_STRING(3)=BLANK(1:IEND)//'<from_class #[1]>'
        OP_STRING(4)=BLANK(1:IEND)//'<to_region #[2]>'
        OP_STRING(5)=BLANK(1:IEND)//'<to_class #[2]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','UPSOUR',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'EXTRACELLULAR_BIDOMAIN',2,noco+1,NTCO,N3CO)) THEN
          TYPE='EXTRACELLULAR_BIDOMAIN'
        ENDIF

        IF(CBBREV(CO,'FROM_REGION',6,noco+1,NTCO,N3CO)) THEN
          nr_s=IFROMC(CO(N3CO+1))
        ELSE
          nr_s=1
        ENDIF

        IF(CBBREV(CO,'FROM_CLASS',6,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=1
        ENDIF
      	CALL NX_LOC(NX_INQUIRE,nxc,nx_s,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_s.NE.0,'Invalid source class',ERROR,*9999)
        CALL ASSERT(ITYP5(nr_s,nx_s).EQ.2.AND.ITYP2(nr_s,nx_s).EQ.9,
     '    'Source equation must be activation model',ERROR,*9999)

        IF(CBBREV(CO,'TO_REGION',4,noco+1,NTCO,N3CO)) THEN
          nr_t=IFROMC(CO(N3CO+1))
        ELSE
          nr_t=2
        ENDIF

        IF(CBBREV(CO,'TO_CLASS',4,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=2
        ENDIF
      	CALL NX_LOC(NX_INQUIRE,nxc,nx_t,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_t.NE.0,'Invalid target class',ERROR,*9999)

        CALL ASSERT(ITYP5(nr_t,nx_t).EQ.1.AND.ITYP2(nr_t,nx_t).EQ.5,
     '    'Target equation must be div(k.grad(u))=f',ERROR,*9999)

C--type---------------------------------------------------------------

        IF(TYPE(1:22).EQ.'EXTRACELLULAR_BIDOMAIN') THEN

!   Update RHS source term for the bidomain Poisson eqtn (nx_t on nr_t)
!   from div(kgrad(phi_m)) at the intracellular grid pts (nx_s on nr_s)
          nb=NBJ(1,NEELEM(1,nr_s)) !basis# for 1st elem of source region
          DO nq=1,NQT !loop over all grid pts
            DO nik=-1,1
              DO nij=-1,1
                DO nii=-1,1
                  PHIMQ(nii,nij,nik)=0.d0
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik

!      Formulate a local quadratic element about nq defined by UMQ,
!      storing the value of phi(m) at each node point mq.
            NITB=NIT(nb)
            IK=MAX(0,NITB-2) !zero for 1,2D, one for 3D
            IJ=MIN(NITB-1,1) !zero for 1D, one for 2,3D
            DO nik=-IK,IK
              DO nij=-IJ,IJ
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                  IF(mq.GT.0) THEN
                    PHIMQ(nii,nij,nik)=YQ(mq,1,1,nx_s)
                  ELSE
                    PHIMQ(nii,nij,nik)=0.0d0
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik

!      Compute PHI,k by first order finite differences about nq.
            dPHIdX(1)=PHIMQ(1,0,0)-PHIMQ(-1,0,0)
            dPHIdX(2)=PHIMQ(0,1,0)-PHIMQ(0,-1,0)
            dPHIdX(3)=PHIMQ(0,0,1)-PHIMQ(0,0,-1)

!      Compute PHI,jk by second order f.d. about nq
            DO ib=1,NITB
              DO jb=1,NITB
                IF(ib.EQ.jb) THEN
                  IF(ib.EQ.1) THEN !PHI,11
                    d2PHIdX2(ib,jb)=(PHIMQ(1,0,0)-
     '                2.0d0*PHIMQ(0,0,0)+PHIMQ(-1,0,0))*4.0d0
                  ELSE IF(ib.EQ.2) THEN !PHI,22
                    d2PHIdX2(ib,jb)=(PHIMQ(0,1,0)
     '                -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,-1,0))*4.0d0
                  ELSE IF(ib.EQ.3) THEN !PHI,33
                    d2PHIdX2(ib,jb)=(PHIMQ(0,0,1)
     '                -2.0d0*PHIMQ(0,0,0)+PHIMQ(0,0,-1))*4.0d0
                  ENDIF               
                ELSE
                  IF(ib+jb.EQ.3) THEN !PHI,12 and PHI,21
                    d2PHIdX2(ib,jb)=PHIMQ(1,1,0)-PHIMQ(1,-1,0)
     '                -PHIMQ(-1,1,0)+PHIMQ(-1,-1,0)
                  ELSE IF(ib+jb.EQ.4) THEN !PHI,13 and PHI,31
                    d2PHIdX2(ib,jb)=PHIMQ(1,0,1)-PHIMQ(1,0,-1)
     '                -PHIMQ(-1,0,1)+PHIMQ(-1,0,-1)
                  ELSE IF(ib+jb.EQ.5) THEN !PHI,23 and PHI,32
                    d2PHIdX2(ib,jb)=PHIMQ(0,1,1)-PHIMQ(0,-1,1)
     '                -PHIMQ(0,1,-1)+PHIMQ(0,-1,-1)
                  ENDIF
                ENDIF
              ENDDO !jb
            ENDDO !ib

!      Compute the diffusion in three parts.
            DIFFN=0.0d0
            DO nii=1,NITB
              DO nij=1,NITB           
                SUM1=0.0d0
                SUM2=0.0d0
                DO nik=1,NITB
                  SUM1=SUM1+PROPQ(nik,nii,nij+1,1,nq,nx_s)
     '                     *dPHIdX(nik)
                  SUM2=SUM2+PROPQ(nik,nii,1,1,nq,nx_s)
     '                     *d2PHIdX2(nij,nik)
                ENDDO
                DIFFN=DIFFN+(SUM1+SUM2)*GUQ(nii,nij,nq)
              ENDDO !nij
            ENDDO !nii
            DO nii=1,NITB
              SUM1=0.0d0
              DO nij=1,NITB
                SUM1=SUM1+PROPQ(nij,nii,1,1,nq,nx_s)*dPHIdX(nij)
              ENDDO
              DIFFN=DIFFN-SUM1*GCHQ(nii,nq)
            ENDDO !nii
            CP(1,NPQ(nq),nx_t)=-DIFFN !is extracellular source term
          ENDDO !nq

C--type---------------------------------------------------------------

        ENDIF !type
      ENDIF

      CALL EXITS('UPSOUR')
      RETURN
 9999 CALL ERRORS('UPSOUR',ERROR)
      CALL EXITS('UPSOUR')
      RETURN 1
      END


