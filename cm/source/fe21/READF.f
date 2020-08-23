      SUBROUTINE READF(IBT,IDO,INP,ISC_GKK,
     '  ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,ISR_GKK,ITHRES,
     '  LD_NP,MAP_ART_VEIN,NAN,NBH,NBJ,NBJF,NDET,NEELEM,NELIST,NENP,NEL,
     &  NFF,NGAP,NHE,NHP,NHQ,NKJE,NKEF,NKH,NKJ,NLF,NLL,NNF,NNL,NONY,NPF,
     &  NPL,NPNE,NPNODE,NPNY,NP_INTERFACE,NQNE,NQNY,NQS,NQXI,NRE,NRLIST,
     '  NRLIST2,NVHE,NVHP,
     '  NVJE,NVJP,NW,NWQ,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '  CE,CONY,CP,DET,DL,FEXT,GKK,GR,GRR,MFI,PE,PF,PG,
     '  PHI,PHI_H,PHI_H_EXACT,SE,T_BH,T_BH_INV,THRES,WG,XA,XIG,XP,YP,YQ,
     '  YQS,ZCROSSING,STRING,FIX,FIXP,ERROR,*)

C#### Subroutine: READF
C###  Description:
C###    READF reads .COM, HISTORY, .IOD, MATRIX, .XP, or .SIGNAL files.

C LKC 5-NOV-97 unused ,NBHF
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GKK(NISC_GKKM,NXM),ISR_GKK(NISR_GKKM,NXM),
     '  ISIZE_MFI(3,NSSM),ISIZE_PHI(2),ISIZE_PHIH(2),ISIZE_TBH(2),
     '  ITHRES(3,NGM,NEM),LD_NP(NDM),MAP_ART_VEIN(0:NDM,NRM),
     &  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDET(NBFM,0:NNM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NEL(0:NELM,NLM),NFF(6,NEM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJ(NJM,NPM),NLF(4,NFM),
     '  NLL(12,NEM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NP_INTERFACE(0:NPM,0:3),
     '  NQNE(NEQM,NQEM),NQNY(2,NYQM,0:NRCM,NXM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,3,NXM),
     '  NWQ(8,0:NQM,NAM),NXLIST(0:NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CP(NMM,NPM,NXM),DET(NBFM,0:NNM,NGM,6),DL(3,NLM),
     '  FEXT(NIFEXTM,NGM,NEM),
     '  GKK(NZ_GKK_M,NXM),GR(NYROWM),
     '  GRR(NOM),MFI(NDM,NTSM,3,NSSM),
     '  PE(2,NEM),PF(2,NEM),PG(NSM,NUM,NGM,NBM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PHI_H_EXACT(NY_TRANSFER_M,NTSM),SE(NSM,NBFM,NEM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M),THRES(3,NGM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),
     '  ZCROSSING(NY_TRANSFER_M,NTSM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXP(2,NEM)

!     Local Variables
      INTEGER CLOCAT,I,IBEG,IBEG1,IEND,IEND1,IFROMC,IPFILE,IPOS,
     '  IUNIT,N3CO,na,NIQLIST(0:99),NIQSLIST(0:500),NIYLIST(0:99),
     '  niy,niq,niqs,nj,nk,nonode,nonr,nr,nss,NTVARLIST,nx,nxc
      REAL*8 RFROMC,TIME,YPMIN(16),YPMAX(16)
      CHARACTER FILE*100,FILEFORMAT*6,INOUTTYPE*10,
     '  TYPE*7
      CHARACTER*10 VARLIST(3)
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,END,ENDFILE,YPDATA,YQDATA,
     '  YQSDATA


C LKC 5-NOV-97 unused INTEGER NBHF(NHM,NCM,NFM),

C *** DPN 19 October 1999 - Adding ability to copy matrices from 1-D
C ***   strip to 2-D sheet.

      INTEGER*4 YQ_1D_PTR,YQS_1D_PTR
      INTEGER nee,ne,nq_1d_start,nq2,nq1,nq_1d,nq
      LOGICAL LAST_POINT_EXTERNAL,ONE2TWO

C#### Comment: CELL MATRIX FILES
C###  Description:
C###    <html>
C###    <p>It is now possible to take a YQ/YQS matrix file from a 1-D
C###    strip simulation and copy it onto a 2-D sheet of elements.
C###    This is done by copying the 1-D strip onto each row of grid
C###    points in the 2-D sheet.</p>
C###    <p>Maybe extend to:
C###    <ul>
C###      <li>Copy from 1-D/2-D to 3-D ??</li>
C###      <li>Be able to specify which xi direction is the rows for
C###          the copy - currently assumes xi(1)</li>
C###      <li>others....</li>
C###    </ul></p>
C###    <html>

      CALL ENTERS('READF',*9999)
 1    IF(CO(noco+1).EQ.'?'.OR.CO(noco+2).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        IF(CO(noco+1).NE.'?') THEN !? is second command option
C         delete first command option from string
          CALL STRING_TRIM(CO(noco+1),IBEG1,IEND1)
          IPOS=CLOCAT(CO(noco+1)(IBEG1:IEND1),STRING)
          STRING=STRING(1:IPOS-1)//STRING(IPOS+IEND1-IBEG1+1:IEND)
        ENDIF
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM read commands/iod/xp;<FILENAME[$current]><;example>
C###  Description:
C###    Read information from file FILENAME.  FILENAME.com is command
C###    file.  FILENAME.iod is input/oputput dump file of all variables
C###    and arrays.  FILENAME.xp reads a file for image analysis.

        OP_STRING(1)=STRING(1:IEND)
     '    //' commands/iod/xp;<FILENAME['//FILE00(IBEG1:IEND1)
     '    //']><;example>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM read history;<FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Reads history information from a history file FILENAME
C###    (extension .iphist for ascii files and .binhis for binary
C###    files) in the directory specified by PATH into the YP and
C###    YQ arrays.
C###  Parameter: <unit #[20]>
C###    Specify the unit number of the history file to read. This unit
C###    number must have been previously opened.
C###  Parameter: <variables VARLIST[yp,yq,yqs]>
C###    Specify which variable types to read. The available variable
C###    types are: yp, yq and yqs.
C###  Parameter: <niylist (#s/all)[all]>
C###    Specify the list of niy numbers of the yp data to store in the
C###    history file. The all value specifies all currently defined niy
C###    values. Any niy numbers contained in the file but are not
C###    specified will be skipped.
C###  Parameter: <niqlist (#s/all)[all]>
C###    Specify the list of niq numbers of the yq data to store in the
C###    history file. The all value specifies all currently defined niq
C###    values. Any niq numbers contained in the file but are not
C###    specified will be skipped.
C###  Parameter: <niqslist (#s/all)[all]>
C###    Specify the list of niqs numbers of the yqs data to store in the
C###    history file. The all value specifies all currently defined niqs
C###    values. Any niqs numbers contained in the file but are not
C###    specified will be skipped.
C###  Parameter: <time #[0.0]>
C###    Specify the time instance to read.
C###  Parameter: <(ascii/binary)[ascii]>
C###    Specify whether the history file to read from is an ascii or
C###    binary file.
C###  Parameter: <region (#s/all)[all]>
C###    Specify the history file region numbers to read. The all value
C###    specifies all currently defined regions. Any region numbers
C###    contained in the file but are not specified will be skipped.
C###  Parameter: <intoregion (#s/same)[same]>
C###    Specify the actual region numbers to read the history region
C###    numbers into. The same value indicates that the history region
C###    numbers (as specified by the region parameter) will be read
C###    into the same actual region numbers.
C###  Parameter: <class #[1]>
C###    Specify the class number (of solve type) to read the history
C###    data into.

        OP_STRING(1)=STRING(1:IEND)
     '    //' history;<FILENAME['//FILE00(IBEG1:IEND1)
     '    //']><;example>'
        OP_STRING(2)=BLANK(1:15)//'<unit #[20]>'
        OP_STRING(3)=BLANK(1:15)//'<variables VARLIST[yp,yq,yqs]>'
        OP_STRING(4)=BLANK(1:15)//'<niylist (#s/all)[all]>'
        OP_STRING(5)=BLANK(1:15)//'<niqlist (#s/all)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<niqslist (#s/all)[all]>'
        OP_STRING(7)=BLANK(1:15)//'<time #[0.0]>'
        OP_STRING(8)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[all]>'
        OP_STRING(10)=BLANK(1:15)//'<intoregion (#s/same)[same]>'
        OP_STRING(11)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM read matrix;<FILENAME[$current]><;(PATH/example)[$current]>
C###  Parameter:    <(sparsity/array)[both]>
C###    Specifies the type of matrix that is being read
C###  Parameter:    <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:    <region #[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:    <one2two>
C###    If present, attempts to copy a matrix from a 1-D strip into
C###    a 2-D sheet.
C###  Parameter:    <elements  (#s/all)[all]>
C###    If the <one2two> qualifier is present, specifies the elements to
C###    copy the 1-D strip into.
C###  Parameter:    <nss  #[1]>
C###    Specify the signal set to read the matrix into
C###  Description:
C###    Reads matrix from a file with a specified format

        OP_STRING(1)=STRING(1:IEND)//' matrix;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<(sparsity/array)[both]>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<one2two>'
        OP_STRING(7)=BLANK(1:15)//'<elements  (#s/all)[all]>'
        OP_STRING(8)=BLANK(1:15)//'<nss  (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','READF',ERROR,*9999)
      ELSE

C news MPN 6-Jul-95: changing "read file;com" to "read com;file" etc
        IF(ABBREV(CO(noco+1),'COMMANDS',2)) THEN
          TYPE='COM'
        ELSE IF(ABBREV(CO(noco+1),'HISTORY',2)) THEN
          TYPE='HISTORY'
        ELSE IF(ABBREV(CO(noco+1),'IOD',2)) THEN
          TYPE='IOD'
        ELSE IF(ABBREV(CO(noco+1),'MATRIX',2)) THEN
          TYPE='MATRIX'
        ELSE IF(ABBREV(CO(noco+1),'XP',2)) THEN
          TYPE='XP'
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        CALL CHECKF(1,noco+1,NTCOQU,CO,COQU,FILE,STRING,*1)

        CALL STRING_TRIM(FILE,IBEG,IEND)

        IOTYPE=2

        IF(TYPE(1:3).EQ.'COM') THEN
C KAT 25/2/00: Opening and closing COM_UNIT in the same routine.
C          COM_UNIT=COM_UNIT+1
C          NEST=NEST+1
C          CALL OPENF(COM_UNIT,'DISK',FILE(IBEG:IEND)//'.com',
C     '      'OLD','DIRECT',
C     '      'FORMATTED',132,ERROR,*9998)
          IF(FILE(IEND-3:IEND).NE.'.com') THEN
            FILE(IEND+1:)='.com'
            IEND=IEND+4
          ENDIF
          END=.FALSE.
          CALL READCOM(FILE(IBEG:IEND),END,ERROR,*9999)
C          CALL CLOSEF(COM_UNIT,ERROR,*9999)
C          COM_UNIT=COM_UNIT-1
C          NEST=NEST-1

        ELSE IF(TYPE(1:7).EQ.'HISTORY') THEN

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'UNIT',1,noco+1,NTCO,N3CO)) THEN
            IUNIT=IFROMC(CO(N3CO+1))
          ELSE
            IUNIT=20
          ENDIF

          IF(CBBREV(CO,'VARIABLES',1,noco+1,NTCO,N3CO)) THEN
            YPDATA=.FALSE.
            YQDATA=.FALSE.
            YQSDATA=.FALSE.
            CALL PARSSL(CO(N3CO+1),3,NTVARLIST,VARLIST,ERROR,*9999)
            DO i=1,NTVARLIST
              IF(ABBREV(VARLIST(i),'YP',2)) THEN
                YPDATA=.TRUE.
              ELSE IF(ABBREV(VARLIST(i),'YQS',3)) THEN
                YQSDATA=.TRUE.
              ELSE IF(ABBREV(VARLIST(i),'YQ',2)) THEN
                YQDATA=.TRUE.
              ENDIF
            ENDDO
          ELSE
            YPDATA=.TRUE.
            YQDATA=.TRUE.
            YQSDATA=.TRUE.
          ENDIF

          IF(YPDATA) THEN
            IF(CBBREV(CO,'NIYLIST',3,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),99,NIYLIST(0),NIYLIST(1),
     '          ERROR,*9999)
            ELSE
              CALL ASSERT(NIYM.LE.99,
     '          '>>Increase dimension of NIYLIST in READF',
     '          ERROR,*9999)
              NIYLIST(0)=NIYM
              DO niy=1,NIYM
                NIYLIST(niy)=niy
              ENDDO
            ENDIF
          ENDIF
          IF(YQDATA) THEN
            IF(CBBREV(CO,'NIQLIST',4,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),99,NIQLIST(0),NIQLIST(1),
     '          ERROR,*9999)
            ELSE
              CALL ASSERT(NIQM.LE.99,
     '          '>>Increase dimension of NIQLIST in READF',
     '          ERROR,*9999)
              NIQLIST(0)=NIQM
              DO niq=1,NIQM
                NIQLIST(niq)=niq
              ENDDO
            ENDIF
          ENDIF
          IF(YQSDATA) THEN
            IF(CBBREV(CO,'NIQSLIST',4,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),99,NIQSLIST(0),NIQSLIST(1),
     '          ERROR,*9999)
            ELSE
              CALL ASSERT(NIQM.LE.500,
     '          '>>Increase dimension of NIQSLIST in READF',
     '          ERROR,*9999)
              NIQSLIST(0)=NIQSM
              DO niqs=1,NIQSM
                NIQSLIST(niqs)=niqs
              ENDDO
            ENDIF
          ENDIF

          IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
            TIME=RFROMC(CO(N3CO+1))
          ELSE
            TIME=0.0d0
          ENDIF

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'INTOREGIONS',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '        ERROR,*9999)
          ELSE
            NRLIST2(0)=NRLIST(0)
            DO nonr=1,NRLIST(0)
              NRLIST2(nonr)=NRLIST(nonr)
            ENDDO
          ENDIF

          ENDFILE=.FALSE.
          na=1

          CALL IOHIST(IUNIT,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '      NUMTIMEDATA_HIST(IUNIT),nx,NYNR(0,0,1,0,nx),
     '      NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,
     '      YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,FILE,
     '      'TIME_DATA',ENDFILE,.FALSE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

          IF(ENDFILE) THEN
            ERROR='>>Time not found'
            GOTO 9999
          ENDIF

        ELSE IF(TYPE(1:3).EQ.'IOD') THEN
          IUNIT=71
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.iod','OLD',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          READ(UNIT=IUNIT,FMT='(37X,I2)') IPFILE
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' File version number is '',I2)')
     '        IPFILE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          READ(UNIT=IUNIT,FMT='(10X,A)') HEADING
          WRITE(OP_STRING,'('' File heading: '',A)') HEADING
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          READ(UNIT=IUNIT,FMT='(1X)')
          IF(ipfile.ne.3) THEN      !old file
            WRITE(OP_STRING,'('' Cannot read iod version '',I2)')
     '        IPFILE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE                      !new file
            CALL IOGEOM('READ',IUNIT,IBT,IDO,INP,ITHRES,NAN,NBH,
     '        NBJ,NBJF,NDET,NEELEM,NEL,NFF,NGAP,NHE,NHP,
     '        NKJE,NKEF,NKH,NKJ,NLF,NLL,NNF,NNL,NONY,NPF,NPL,NPNE,
     '        NPNODE,NPNY,NP_INTERFACE,NRE,NVHE,NVHP,
     '        NVJE,NVJP,NW,NXI,NYNE,NYNO,NYNP,
     '        CE,CONY,CP,DET,DL,FEXT,PE,PF,PG,SE,THRES,
     '        WG,XA,XIG,XP,YP,FIX,FIXP,ERROR,*9999)
          ENDIF
          CALL CLOSEF(IUNIT,ERROR,*9999)
          CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
C 2001-02-07 KAT: Can't see why NXI needs recalculating after we have
C         just read it in.  We haven't read in NNB, so we can't
C         recalculate it anyway.
C          CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)

        ELSE IF(TYPE(1:6).EQ.'MATRIX') THEN

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)

C LKC 17-APR-2002 Change this to a warning... for matrices like PHI,MFI
C            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     '        ERROR,*9999)
          IF(nx.LE.0) THEN
            OP_STRING(1)='>>WARNING: No nx defined for this solve class'
            OP_STRING(2)='>>WARNING: Forcing nx to 1'
            CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            nx=1

C LKC can't read into these arrays if nx isn't valid anyway.
            IF(ISC_GK_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISC_GKM*USE_SPARSE,0,INTTYPE,
     '          ISC_GK_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(ISR_GK_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISR_GKM*USE_SPARSE,0,INTTYPE,
     '          ISR_GK_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(GK_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NZ_GK_M,1,DPTYPE,GK_PTR(nx),MEM_INIT,
     '          ERROR,*9999)
            ENDIF
            IF(ISC_GQ_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISC_GQM*USE_SPARSE,0,INTTYPE,
     '          ISC_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(ISR_GQ_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NISR_GQM*USE_SPARSE,0,INTTYPE,
     '          ISR_GQ_PTR(nx),MEM_INIT,ERROR,*9999)
            ENDIF
            IF(GQ_PTR(nx).EQ.0) THEN
              CALL ALLOCATE_MEMORY(NZ_GQ_M,1,DPTYPE,GQ_PTR(nx),MEM_INIT,
     '          ERROR,*9999)
            ENDIF
          ENDIF

          IF(CBBREV(CO,'SPARSITY',2,noco+1,NTCO,N3CO)) THEN
            INOUTTYPE='SPARSITY'
          ELSE IF(CBBREV(CO,'ARRAY',2,noco+1,NTCO,N3CO)) THEN
            INOUTTYPE='ARRAYS'
          ELSE
            INOUTTYPE='ALL_INOUT'
          ENDIF

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'NSS',1,noco+1,NTCO,N3CO)) THEN
            nss=IFROMC(CO(N3CO+1))
          ELSE
            nss=1
          ENDIF


C *** DPN 19 October 1999 - Add ability to copy a matrix from a 1D
C ***   strip to a 2D sheet.
          IF(CBBREV(CO,'ONE2TWO',3,NOCO+1,NTCO,N3CO)) THEN
            ONE2TWO=.TRUE.
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
            !check that all are 2-d elements
            DO nee=1,NELIST(0)
              ne=NELIST(nee)
              CALL ASSERT(NQXI(0,NQS(ne)).EQ.2,'Elements must be 2-D',
     '          ERROR,*9999)
            ENDDO
          ELSE
            ONE2TWO=.FALSE.
          ENDIF

          IUNIT=20
          nr=NRLIST(1)
          IF(ONE2TWO) THEN !want to go from 1-d strip to 2-d sheet
            !Allocate memory for the 1-d strip solution
            YQ_1D_PTR=0
            YQS_1D_PTR=0
            CALL ALLOCATE_MEMORY(NYQM*NIQM,1,DPTYPE,YQ_1D_PTR,
     '        MEM_INIT,ERROR,*9999)
            CALL ALLOCATE_MEMORY(NIQSM*NQM,1,DPTYPE,YQS_1D_PTR,
     '        MEM_INIT,ERROR,*9999)

            !Read in the 1-D matrices
            CALL IOMATR('READ ',%VAL(ISC_GK_PTR(nx)),
     '        ISC_GKK(1,nx),%VAL(ISC_GQ_PTR(nx)),
     '        ISIZE_MFI(1,nss),ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '        %VAL(ISR_GK_PTR(nx)),ISR_GKK(1,nx),
     '        %VAL(ISR_GQ_PTR(nx)),IUNIT,LD_NP,MAP_ART_VEIN,nr,nx,
     &        %VAL(GK_PTR(nx)),GKK(1,nx),%VAL(GQ_PTR(nx)),GR,GRR,MFI(1,
     &        1,1,nss),PHI,PHI_H,PHI_H_EXACT,T_BH,T_BH_INV,YP,
     &        %VAL(YQ_1D_PTR),%VAL(YQS_1D_PTR),ZCROSSING,
     '        FILEFORMAT,FILE00,INOUTTYPE,ERROR,*9999)

            !Copy the 1-D strip into the grid point rows
            nq_1d_start=1
            LAST_POINT_EXTERNAL=.FALSE.
            DO nee=1,NELIST(0) !loop through all specified elements
              ne=NELIST(nee)
              IF(NWQ(1,NQNE(ne,1),NQS(ne)).GT.0) THEN
                !the first grid point in the element is on an external
                !boundary, so move back to the start of the 1-D strip.

                !Also need to check that the last grid point in the
                !previous element was on an external border so make sure
                !that the end of the first row of elements has been
                !reached, since the first grid pt. will always be on an
                !external boundary in the first row of elements
                IF(LAST_POINT_EXTERNAL) THEN
                  nq_1d_start=1
                ENDIF
              ENDIF
              DO nq2=1,NQXI(2,NQS(ne))
                !loop through all the "rows" of 1-d strips of grid pts.
                nq_1d=nq_1d_start
                DO nq1=1,NQXI(1,NQS(ne))
                  !loop through each of the grid points in the "row"
                  !get the current grid point number
                  nq=NQNE(ne,nq1+(nq2-1)*NQXI(1,NQS(ne)))
                  !copy over the values - need a subroutine since
                  !yq_1d and yqs_1d are dynamically allocated??
                  CALL ONE_2_TWO(nq_1d,nq,nx,%VAL(YQ_1D_PTR),
     '              %VAL(YQS_1D_PTR),YQ,YQS,ERROR,*9999)
                  !increment 1-d strip counter
                  nq_1d=nq_1d+1
                ENDDO !nq1=1,NQXI(1,NQS(ne))
              ENDDO !nq2=1,NQXI(2,NQS(ne))
              !set the starting grid point in 1-d strip for next element
              nq_1d_start=nq_1d_start+NQXI(1,NQS(ne))-1
              !check if the last grid point in the element is external
              IF(NWQ(1,NQNE(ne,NQXI(1,NQS(ne))*NQXI(2,NQS(ne))),
     '          NQS(ne)).GT.0) THEN
                LAST_POINT_EXTERNAL=.TRUE.
              ELSE
                LAST_POINT_EXTERNAL=.FALSE.
              ENDIF
            ENDDO !nee=1,NELIST(0)

            !Deallocate memory
            CALL FREE_MEMORY(YQ_1D_PTR,ERROR,*9999)
            CALL FREE_MEMORY(YQS_1D_PTR,ERROR,*9999)

          ELSE !just want to read in the matrices

            CALL IOMATR('READ ',%VAL(ISC_GK_PTR(nx)),
     '        ISC_GKK(1,nx),%VAL(ISC_GQ_PTR(nx)),
     '        ISIZE_MFI(1,nss),ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '        %VAL(ISR_GK_PTR(nx)),ISR_GKK(1,nx),
     '        %VAL(ISR_GQ_PTR(nx)),IUNIT,LD_NP,MAP_ART_VEIN,nr,nx,
     &        %VAL(GK_PTR(nx)),GKK(1,nx),%VAL(GQ_PTR(nx)),GR,GRR,
     &        MFI(1,1,1,nss),PHI,PHI_H,PHI_H_EXACT,T_BH,T_BH_INV,YP,
     &        YQ(1,1,1,nx),YQS,ZCROSSING,FILEFORMAT,FILE00,INOUTTYPE,
     &        ERROR,*9999)

          ENDIF !ONE2TWO

        ELSE IF(TYPE(1:2).EQ.'XP') THEN
          CALL OPENF(7,'DISK',FILE(IBEG:IEND)//'.xp','OLD','SEQUEN',
     '      'FORMATTED',132,ERROR,*9999)
          READ(7,'(5E13.5)')(((XP(nk,1,nj,NPNODE(nonode,1)),nk=1,4),
     '      nj=1,3),nonode=1,NPNODE(0,1))
          CALL CLOSEF(7,ERROR,*9999)

        ENDIF
      ENDIF

      CALL EXITS('READF')
      RETURN

C 9998 COM_UNIT=COM_UNIT-1
C      NEST=NEST-1
 9999 IF(TYPE(1:3).EQ.'IOD') THEN
        CLOSE(UNIT=71)
      ELSE IF(TYPE(1:2).EQ.'XP') THEN
        CLOSE(UNIT=7)
      ENDIF
      CALL ERRORS('READF',ERROR)
      CALL EXITS('READF')
      RETURN 1
      END


