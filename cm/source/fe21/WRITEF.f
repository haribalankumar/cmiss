      SUBROUTINE WRITEF(IBT,IDO,INP,ISC_GKK,
     '  ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,ISR_GKK,ITHRES,
     '  LD_NP,MAP_ART_VEIN,NAN,NBH,NBJ,NBJF,NDDL,NDET,NDLT,NDP,NEELEM,
     &  NEL,NELIST,NFF,NGAP,NHE,NHP,NHQ,NKEF,NKH,NKHE,NKJ,NKJE,NLF,
     '  NLL,NNF,NNL,NONY,NPF,NPL,NPNE,NPNODE,NPNY,
     '  NP_INTERFACE,NQNY,NRE,NRLIST,NRLIST2,NVHE,NVHP,NVJE,
     '  NVJP,NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '  CE,CG,CONY,CP,CURVCORRECT,DET,DL,FEXT,GKK,GR,GRR,MFI,PE,PF,
     '  PG,PHI,PHI_H,PHI_H_EXACT,SE,T_BH,T_BH_INV,THRES,
     '  WD,WDL,WG,XA,XE,XG,XID,XIDL,XIG,XP,YP,YQ,YQS,ZA,Z_CONT,
     '  ZCROSSING,
     '  ZD,ZD2,ZDL,ZE,ZG,ZP,STRING,FIX,FIXP,ERROR,*)

C#### Subroutine: WRITEF
C###  Description:
C###    WRITEF writes .DAT, .IOD, MATRIX, .OPG, .OPS, or .TRACE file.
C###    CPB 14/2/93 Writes map3d .data files

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'dtran00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'head00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'matr00.cmn'
      INCLUDE 'matr00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GKK(NISC_GKKM,NXM),ISIZE_MFI(3,NSSM),
     '  ISIZE_PHI(2),ISIZE_PHIH(2),ISIZE_TBH(2),ISR_GKK(NISR_GKKM,NXM),
     '  ITHRES(3,NGM,NEM),LD_NP(NDM),MAP_ART_VEIN(0:NDM,NRM),
     &  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDDL(NEM,NDEM),NDET(NBFM,0:NNM),NDLT(NEM),
     '  NDP(NDM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),NELIST(0:NEM),
     '  NFF(6,NEM),NGAP(NIM,NBM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NHQ(NRM,NXM),NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLF(4,NFM),NLL(12,NEM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NP_INTERFACE(0:NPM,0:3),
     '  NQNY(2,NYQM,0:NRCM,NXM),NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),
     '  CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),DET(NBFM,0:NNM,NGM,6),
     '  DL(3,NLM),FEXT(NIFEXTM,NGM,NEM),
     '  GKK(NZ_GKK_M,NXM),GR(NYROWM),
     '  GRR(NOM),MFI(NDM,NTSM,3,NSSM),
     '  PE(2,NEM),PF(2,NEM),PG(NSM,NUM,NGM,NBM),
     '  PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PHI_H_EXACT(NY_TRANSFER_M,NTSM),SE(NSM,NBFM,NEM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  THRES(3,NGM,NEM),WD(NJM,NDM),
     '  WDL(NHM,NDEM),WG(NGM,NBM),XA(NAM,NJM),XE(NSM,NJM),XG(NJM,NUM),
     '  XID(NIM,NDM),XIDL(NIM,NDEM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),Z_CONT(NDM,2,67),
     '  ZCROSSING(NY_TRANSFER_M,NTSM),
     '  ZD(NJM,NDM),ZD2(NJM,NDM),ZDL(NHM,NDEM),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXP(2,NEM)

!     Local Variables
      INTEGER CLOCAT,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,IPFILE,
     '  IPOS,IUNIT,matr,N,na,N3CO,nb,nd,nde,ne,ni,nj,noelem,nolist,
     '  nomatr,norl,nr,nss,NTITLV(2),NTLV,NTOT,NTRL1,NTVARLIST,
     '  nx,nxc,SERIES
C      INTEGER FCreateFile,FSetText,FSetTimeSeriesIndex,
C     '  FSetTimeSeriesLabel,FSetTimeSeriesFormat,
C     '  FSetTimeSeriesSpecs,FSetTimeSeriesData,FCloseFile
      REAL*8 D,PXI,RL(200,4),X(3),Z(3),Z2(3),
     '  EG(3,3),RFROMC,TIME,W(4),YPMAX(1),YPMIN(1)
      CHARACTER CHAR*1,FILE*100,FILEFORMAT*6,
     '  INOUTTYPE*10,TITLE*80,TYPE*7
      CHARACTER*10 VARLIST(3)
C      CHARACTER HEADERTEXT1*256,HEADERTEXT2*80
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,DEFORM,ENDFILE,FOUND,PROJ,
     '  STATIC,YPDATA,YQDATA,YQSDATA,ZEROS

C 17-JAN-98 191297
C!!! LKC for AJP the PHI(3) variable has been renamed to PHIL(3)


      CALL ENTERS('WRITEF',*9999)

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

C#### Command: FEM write trace;<FILENAME[$current]>
C###  Parameter:         <transform>
C###  Parameter:         <in (ELEMENT#s/all)[all]>
C###  Parameter:         <(static/dynamic)[static]>
C###  Parameter:         <(deform/undeform)[undeform]>
C###  Parameter:         <zeros>
C###  Parameter:         <project>
C###  Description:
C###    Write FILENAME.trace, containing data trace positions (imaging
C###    work).

        OP_STRING(1)=STRING(1:IEND)//' trace;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<transform>'
        OP_STRING(3)=BLANK(1:15)//'<in (ELEMENT#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<(static/dynamic)[static]>'
        OP_STRING(5)=BLANK(1:15)//'<(deform/undeform)[undeform]>'
        OP_STRING(6)=BLANK(1:15)//'<zeros>'  !use if only overwrite NJs with zero weight AAY
        OP_STRING(7)=BLANK(1:15)//'<project>'!use if project predicted point onto data plane AAYMar93
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C KAT 2001-04-09 Not working/used
CC#### Command: FEM write kinematics;<FILENAME<$current]>
CC###  Parameter:     <in (ELEMENT#s/all)[all]>
CC###    Specify for which elements to write out the kinematics
CC###  Description:
CC###    Write out a kinematics data file

C        OP_STRING(1)=STRING(1:IEND)//' kinematics;<FILENAME['
C     '    //FILE00(IBEG1:IEND1)//']>'
C        OP_STRING(2)=BLANK(1:15)//'<in (ELEMENT#s/all)[all]>'
C        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM write curv;<FILENAME<$current]>
C###  Parameter:     <in (ELEMENT#s/all)[all]>
C###    Specify for which elements to write out the kinematics
C###  Description:
C###    Write out a curvature data file

        OP_STRING(1)=STRING(1:IEND)//' curv;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<in (ELEMENT#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM write iod;<FILENAME<$current]>
C###  Description:
C###    Write FILENAME;iod, containing finite element arrays.

        OP_STRING(1)=STRING(1:IEND)//' iod;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM write dat;<FILENAME[$current]>
C###  Parameter:     <ARRAY_1><,ARRAY_2><,ARRAY_3>
C###  Description:
C###    Write specified arrays to FILENAME.dat.

        OP_STRING(1)=STRING(1:IEND)//' dat;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<ARRAY_1><,ARRAY_2><,ARRAY_3>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

C#### Command: FEM write history
C###  Description:
C###    Write history information at a given time to a file specified
C###    by the unit number.
C###  Parameter: <time #[0.0]>
C###    Specify the time at which to write the history data.
C###  Parameter: <unit #[20]>
C###    Specify the unit number of the history file to write to.
C###  Parameter: <variables VARLIST[yp,yq,yqs]>
C###    Specify the variable types to be written to the history file.
C###    The list of available variable types are: yp, yq and yqs.
C###  Parameter: <(ascii/binary)[ascii]>
C###    Specify whether the history file to close is an ascii or
C###    binary file.
C###  Parameter: <class #[1]>
C###    Specify the class number (of solve type) to write the history
C###    data from.


        OP_STRING(1)=STRING(1:IEND)//' history'
        OP_STRING(2)=BLANK(1:15)//'<unit #[20]>'
        OP_STRING(3)=BLANK(1:15)//'<variables VARLIST[yp,yq,yqs]>'
        OP_STRING(4)=BLANK(1:15)//'<time #[0.0]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM write matrix;<FILENAME[$current]>
C###  Parameter:     <matrices NAMELIST[GK]>
C###  Parameter:     <(sparsity/array)[both]>
C###  Parameter:     <(ascii/binary)[ascii]>
C###  Parameter:     <region #[1]>
C###  Parameter:     <using (fit/solve)[solve]>
C###  Parameter:     <class #[1]>
C###  Parameter:     <nss #[1]>
C###    The signal set to write out.
C###  Description:
C###    Write matrix information to a file.

        OP_STRING(1)=STRING(1:IEND)//' matrix;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<matrices NAMELIST[GK]>'
        OP_STRING(3)=BLANK(1:15)//'<(sparsity/array)[both]>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<nss #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM write map3d;<FILENAME[$current]>
C###  Parameter:     <series #[1]>
C###  Description:
        OP_STRING(1)=STRING(1:IEND)//' map3d;<FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<series #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','WRITEF',ERROR,*9999)
      ELSE

C news MPN 6-Jul-95: changing "write file;com" to "write com;file" etc
        IF(ABBREV(CO(noco+1),'DAT',3)) THEN
          TYPE='DAT'
        ELSE IF(ABBREV(CO(noco+1),'GKS',3)) THEN
          TYPE='GKS'
        ELSE IF(ABBREV(CO(noco+1),'HISTORY',3)) THEN
          TYPE='HISTORY'
        ELSE IF(ABBREV(CO(noco+1),'IOD',3)) THEN
          TYPE='IOD'
        ELSE IF(ABBREV(CO(noco+1),'MATRIX',2)) THEN
          TYPE='MATRIX'
        ELSE IF(ABBREV(CO(noco+1),'OPG',3)) THEN
          TYPE='OPG'
        ELSE IF(ABBREV(CO(noco+1),'OPS',3)) THEN
          TYPE='OPS'
        ELSE IF(ABBREV(CO(noco+1),'TRACE',2)) THEN
          TYPE='TRACE'
c MHT archived 18-Feb-99
c       ELSE IF(ABBREV(CO(noco+1),'SIGNAL',2)) THEN
c         TYPE='SIGNAL'
        ELSE IF(ABBREV(CO(noco+1),'MAP3D',2)) THEN
          TYPE='MAP3D'
          IF(CBBREV(CO,'SERIES',2,noco+1,NTCO,N3CO)) THEN
            SERIES=IFROMC(CO(N3CO+1))
          ELSE
            SERIES=1
          ENDIF
C GMH 13/2/97 Avoid unused error (should be deleted?)
          SERIES=SERIES
C KAT 2001-04-09 Not working/used
C        ELSE IF(ABBREV(COQU(noco+1,1),'KINEMATICS',2)) THEN
C          TYPE='KINEMA'
        ELSE IF(ABBREV(COQU(noco+1,1),'CURV',4)) THEN
          TYPE='CURV'
        ELSE
          ERROR='Unknown file type to write out'
          GOTO 9999
        ENDIF

C LKC 17-MAR-2003 Signal set, currently only for MFI        
        IF(CBBREV(CO,'NSS',1,noco+1,NTCO,N3CO)) THEN
          nss=IFROMC(CO(N3CO+1))
        ELSE
          nss=1
        ENDIF


C news MPN 12-Apr-96: Check handles defaults
        CALL CHECKF(1,noco+1,NTCOQU,CO,COQU,FILE,STRING,*1)
        CALL STRING_TRIM(FILE,IBEG,IEND)

        IF(TYPE(1:3).EQ.'IOD') THEN
C CPB 23/11/92 changed unit number for use with writes
          IUNIT=IOFILE2
C****
C**** Put back to previous open command.
C**** Apparantly changed for silicon graphics version, however under VMS
C**** new version of iod file is not created if an iod file with the
C**** same name already exists.
C**** MPN 27/4/92
C****
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.iod','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          IPFILE=3 !is input file version number on 26-May-1993
          WRITE(IUNIT,'(A,I2)') ' CMISS Version '//CMISS
     '      //' IOD File Version ',IPFILE
          WRITE(IUNIT,'(A)') ' Heading: '//HEADING
          WRITE(IUNIT,'(1X)')
          CALL IOGEOM('WRITE',IUNIT,IBT,IDO,INP,ITHRES,NAN,NBH,NBJ,
     '      NBJF,NDET,NEELEM,NEL,NFF,NGAP,NHE,NHP,
     '      NKJE,NKEF,NKH,NKJ,NLF,NLL,NNF,NNL,NONY,NPF,NPL,NPNE,
     '      NPNODE,NPNY,NP_INTERFACE,NRE,NVHE,NVHP,NVJE,NVJP,NW,
     '      NXI,NYNE,NYNO,NYNP,
     '      CE,CONY,CP,DET,DL,FEXT,PE,PF,PG,SE,THRES,
     '      WG,XA,XIG,XP,YP,FIX,FIXP,ERROR,*9999)
          CALL CLOSEF(IUNIT,ERROR,*9999)

        ELSE IF(TYPE(1:3).EQ.'DAT') THEN
          TITLE(1:9)='Data file'
C****
C**** Put back to previous open command.
C**** Apparantly changed for silicon graphics version, however under VMS
C**** new version of dat file is not created if a dat file with the
C**** same name already exists.
C**** MPN 27/4/92
C****
          IUNIT=IOFILE2
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.dat','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
!old      CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.dat','UNKNOWN',
!old '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          IF(noco+1.EQ.NTCO) THEN
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C            CALL IODATA('WRITE','RADIANS','GEOMETRY',IUNIT,NDP,
C     '        NJ_LOC(0,0,0),TITLE,WD,ZD,ERROR,*9999)
            CALL IODATA('WRITE','RADIANS','GEOMETRY',IUNIT,NDP,
     &        NJ_LOC(0,0,0),TITLE,WD,Z_CONT,ZD,ERROR,*9999)
          ELSE
            WRITE(IUNIT,*) TITLE(1:9)
            CALL PARSTR(CO(noco+2),2,NTLV,NTITLV,200,RL(1,1),
     '        ERROR,*9999)
            NTRL1=NTITLV(1)
            IF(NTCOQU(noco+2).GE.1) THEN
              CALL PARSTR(COQU(noco+2,1),2,NTLV,NTITLV,200,RL(1,2),
     '          ERROR,*9999)
            ENDIF
            IF(NTCOQU(noco+2).GE.2) THEN
              CALL PARSTR(COQU(noco+2,2),2,NTLV,NTITLV,200,RL(1,3),
     '          ERROR,*9999)
            ENDIF
            IF(NTCOQU(noco+2).GE.3) THEN
              CALL PARSTR(COQU(noco+2,3),2,NTLV,NTITLV,200,RL(1,4),
     '          ERROR,*9999)
            ENDIF
            NTOT=NTCOQU(noco+2)+1
            WRITE(CHAR,'(I1)') 2*NTOT
            FORMAT='(I4,'//CHAR//'E11.3,'' ,0,0,0,0'')'
            DO N=1,NTOT
              W(N)=1.0d0
            ENDDO
            DO norl=1,NTRL1
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' RL('',I4,'',n): '',4E11.3)')
     '            norl,(RL(N,norl),N=1,NTOT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(IUNIT,FORMAT) norl,(RL(norl,N),N=1,NTOT),
     '            (W(N),N=1,NTOT)
CC$              call mp_unsetlock()
              ENDIF
            ENDDO
          ENDIF
C!!! Why 9?  Same unit as IOOUT for `set output'!
          CALL CLOSEF(9,ERROR,*9999)


        ELSE IF(TYPE(1:7).EQ.'HISTORY') THEN

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)

          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)

          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'TIME',1,noco+1,NTCO,N3CO)) THEN
            TIME=RFROMC(CO(N3CO+1))
          ELSE
            TIME=0.0d0
          ENDIF

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

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF
          na=1

          CALL IOHIST(IUNIT,na,NHQ(1,nx),NIQLIST_HIST(0,IUNIT),
     '      NIQSLIST_HIST(0,IUNIT),NIYLIST_HIST(0,IUNIT),
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_HIST(0,IUNIT),
     '      NRLIST2,NUMTIMEDATA_HIST(IUNIT),nx,NYNR(0,0,1,0,nx),
     '      NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,
     '      YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,FILE,
     '      'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

        ELSE IF(TYPE(1:6).EQ.'MATRIX') THEN

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

C LKC 23-JUN-1999 generalisation for different types
C
C          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C     '      ERROR,*9999)

          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '          ERROR,*9999)
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ENDIF
          ELSE
C LKC 27-DEC-2002 Shouldn't have to have a nx defined to do this
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
C            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
C           '      ERROR,*9999)
            IF(nx.EQ.0) THEN
              nx=1 !set to prevent 0th components of arrays being passed
            ENDIF
          ENDIF

          IF(CBBREV(CO,'MATRICES',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSSL(CO(N3CO+1),NUMMATRMX,NUMMATRNAMES,MATRNAMELIST,
     '        ERROR,*9999)
            DO nomatr=1,NUMMATRNAMES
C              MATRNAMELIST(nomatr)=CUPPER(MATRNAMELIST(nomatr))
              CALL CUPPER(MATRNAMELIST(nomatr),MATRNAMELIST(nomatr))
            ENDDO !nomatr
          ELSE
            NUMMATRNAMES=1
            MATRNAMELIST(1)='GK         '
          ENDIF

          MATRLIST(0)=0
          VECTLIST(0)=0
          DO nomatr=1,NUMMATRNAMES
            CALL STRING_TRIM(MATRNAMELIST(nomatr),IBEG1,IEND1)
            matr=1
            FOUND=.FALSE.
            DO WHILE(.NOT.FOUND.AND.matr.LE.NUMMATRMX)
              CALL STRING_TRIM(MATRNAME(matr),IBEG2,IEND2)
              IF(MATRNAMELIST(nomatr)(IBEG1:IEND1).EQ.
     '          MATRNAME(matr)(IBEG2:IEND2)) THEN
                FOUND=.TRUE.
              ELSE
                matr=matr+1
              ENDIF
            ENDDO
            IF(FOUND) THEN
              IF(MATRTYPE(matr).EQ.0) THEN !is a vector
                VECTLIST(0)=VECTLIST(0)+1
                VECTLIST(VECTLIST(0))=matr
              ELSE IF(MATRTYPE(matr).EQ.1) THEN !is a matrix
                MATRLIST(0)=MATRLIST(0)+1
                MATRLIST(MATRLIST(0))=matr
              ENDIF
            ELSE
              WRITE(ERROR,'('' >>Matrix '',A,'' is unknown'')')
     '          MATRNAMELIST(nomatr)(IBEG1:IEND1)
              GOTO 9999
            ENDIF
          ENDDO !matr

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

          IUNIT=20
          nr=NRLIST(1)

C LKC 17-MAR-2003 Now command line specified
C          nss=1          
          CALL IOMATR('WRITE',%VAL(ISC_GK_PTR(nx)),
     '      ISC_GKK(1,nx),%VAL(ISC_GQ_PTR(nx)),ISIZE_MFI(1,nss),
     '      ISIZE_PHI,ISIZE_PHIH,
     '      ISIZE_TBH,%VAL(ISR_GK_PTR(nx)),ISR_GKK(1,nx),
     '      %VAL(ISR_GQ_PTR(nx)),IUNIT,LD_NP,MAP_ART_VEIN,nr,nx,
     &      %VAL(GK_PTR(nx)),GKK(1,nx),%VAL(GQ_PTR(nx)),GR,GRR,
     &      MFI(1,1,1,nss),PHI,PHI_H,PHI_H_EXACT,T_BH,T_BH_INV,
     '      YP,YQ(1,1,1,nx),YQS,ZCROSSING,FILEFORMAT,FILE00,
     '      INOUTTYPE,ERROR,*9999)

        ELSE IF(TYPE(1:3).EQ.'OPG') THEN
C 25/2/97 LC archived section  : Not used (i hope). GBS 27-OCT-1994
        ELSE IF(TYPE(1:3).EQ.'OPS') THEN
C LC 25/2/97 archived section  : Not used. PJH and AJP  26-4-93
        ELSE IF(TYPE(1:5).EQ.'TRACE') THEN
!news     AAY 21 Jun 91  need to zero zd2
          DO nd=1,NDT
            DO nj=1,NJT
              ZD2(nj,nd)=0.0d0
            ENDDO
          ENDDO
!newe
          IF(CBBREV(CO,'IN',1,noco,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),
     '        ERROR,*9999)
          ELSE
            NELIST(0)=0
            DO nr=1,NRT
              DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
                NELIST(noelem)=NEELEM(noelem,nr)
              ENDDO
              NELIST(0)=NELIST(0)+NEELEM(0,nr)
            ENDDO
          ENDIF
          nr=1 !Needs fixing
          nx=1 !temporary
          STATIC=CBBREV(CO,'STATIC',2,noco,NTCO,N3CO)
          ZEROS=CBBREV(CO,'ZEROS',3,noco,NTCO,N3CO)   !new AAY 30Aug95
          DEFORM=CBBREV(CO,'DEFORM',3,noco,NTCO,N3CO) !new AAY 30Aug95
          PROJ=CBBREV(CO,'PROJECT',3,noco,NTCO,N3CO)  !new AAY 30Aug95
          IF(.NOT.STATIC)THEN !motion field in YP
            CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '        YP(1,1,nx),ZA,ZP,ERROR,*9999)
          ENDIF
          TITLE(1:10)='Trace file'
C****
C**** Put back to previous open command.
C**** Apparantly changed for silicon graphics version, however under VMS
C**** new version of ipdata file is not created if an ipdata file with
C**** the same name already exists.
C**** MPN 27/4/92
C****
!old      CALL OPENF(9,'DISK',FILE(IBEG:IEND)//'.ipdata','UNKNOWN',
!old '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          IUNIT=IOFILE2
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.ipdata','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            nr=NRE(ne)
            ni=NIT(NBJ(1,ne))
!news AAY 30 May 1991
            IF(DEFORM)THEN
c           Note use of XE instead of ZE
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),
     '          NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '          ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
            ELSE
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA,XE,XP,ERROR,*9999)
            ENDIF
!newe
            CALL ZDZDL(0,NBH(1,1,ne),NBJ(1,ne),NDDL,NDLT(ne),ne,1,WD,
     '        WDL,XID,XIDL,ZD,ZDL,ERROR,*9999)

            IF(.NOT.STATIC)THEN !motion field in ZE
              CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),
     '          NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '          ZE,ZP,ERROR,*9999)
            ENDIF

            DO nde=1,NDLT(ne)
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XIDL(1,NDE),XE(1,nj))
                IF(.NOT.STATIC)THEN !motion
                  nb=NBH(nj,1,ne)
                  Z2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XIDL(1,NDE),ZE(1,nj))
                ENDIF
              ENDDO !nj
              CALL XZ(ITYP10(nr),X,Z)

              DO nj=1,NJT
                !new AAY leave nj alone if ZEROS and weight not zero
                IF(.NOT.ZEROS.OR.WD(nj,NDDL(ne,NDE)).EQ.0.0d0)THEN
                  ZD2(nj,NDDL(ne,NDE))=Z(nj)
                ELSE
                  ZD2(nj,NDDL(ne,NDE))=ZD(nj,NDDL(ne,NDE))
                ENDIF
              ENDDO !nj

!news         AAY 10 Mar 93 project onto plane for tag stripe fits
              ND=NDDL(NE,NDE)
              IF(PROJ) THEN
                D=0.d0
                DO nj=1,NJT
c                 Note: stripe normals in ZD(NJT+nj,nd)
                  D=D+(ZD2(nj,nd)-ZD(nj,nd))*ZD(NJT+nj,nd)
                ENDDO
                DO nj=1,NJT
                  ZD2(nj,nd)=ZD2(nj,nd)-D*ZD(NJT+nj,nd)
                ENDDO
              ENDIF !proj
!newe
            ENDDO !nde
          ENDDO !nolist

          IF(ABBREV(CO(noco+1),'TRANSFORM',1)) THEN
            DO nd=1,NDT
              CALL ZZ(ZD2(1,nd),ZD2(1,nd),DATRAN)
            ENDDO
          ENDIF
C*** 10/10/08 JHC Pass in Z_CONT to IODATA.f
C          CALL IODATA('WRITE','RADIANS','GEOMETRY',IUNIT,NDP,NJT,
C     '      TITLE,WD,ZD2,ERROR,*9999)
          CALL IODATA('WRITE','RADIANS','GEOMETRY',IUNIT,NDP,NJT,
     &      TITLE,WD,Z_CONT,ZD2,ERROR,*9999)
          CALL CLOSEF(IUNIT,ERROR,*9999)

C KAT 2001-04-09 Not working/used
!news   AAY 30 Apr 92  kinematic data file
C        ELSE IF(TYPE(1:6).EQ.'KINEMA') THEN
C          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.KINE','UNKNOWN',
C     '      'SEQUEN','FORMATTED',360,ERROR,*9999)
C          NITB=NIT(NBJ(1,1))
C          IF(NITB.EQ.3)THEN
C            WRITE(IUNIT,'(A)')'  nd        Undef pos              '
C     '        //'Def pos'
C     '        //'      w1  w2  w3  ne  xi1      xi2      xi3'
C     '        //'   alpha     r1      r2      r3'
C     '        //'     e11     e12     e13     e22     e23     e33'
C     '        //'    pst1    pst2    pst3     phi1     phi2     phi3'
C          ELSE IF(NITB.EQ.2)THEN
C            WRITE(IUNIT,'(A)')'  nd        Undef pos              '
C     '        //'Def pos'
C     '        //'      w1  w2  w3  ne  xi1      xi2'
C     '        //'   alpha'
C     '        //'     e11     e12     e22'
C     '        //'    pst1    pst2    phi1'
C          ENDIF

C          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
C     '      ERROR,*9999)

C          IF(NITB.EQ.3)FORMAT='(I5,6F7.2,3F4.1,I4,16F8.3,3F8.3)'
C          IF(NITB.EQ.2)FORMAT='(I5,4F7.2,2F4.1,I4,10F8.3)'
C          DO nolist=1,NELIST(0)
C            ne=NELIST(nolist)
C            nr=NRE(ne)
C            NITB=NIT(NBJ(1,ne))
C            CALL XPXE(NBJ(1,ne),NKE(1,1,1,ne),NPF(1,1),
C     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
C         '        SE(1,1,ne),XA,XE,XP,ERROR,*9999)
C            CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKE(1,1,1,ne),
C     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
C     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
C     '        ERROR,*9999)
C            CALL ZDZDL(0,NBH(1,1,ne),NBJ(1,ne),NDDL,NDLT(ne),ne,1,WD,
C     '        WDL,XID,XIDL,ZD,ZDL,ERROR,*9999)
C            DO NDE=1,NDLT(NE)
C              ND=NDDL(NE,NDE)
C              DO nj=1,NJT
C                nb=NBJ(nj,NE)
C                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,1,XIDL(1,NDE),XE(1,nj))
C                Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C     '            nb,1,XIDL(1,NDE),ZE(1,nj))
C              ENDDO
C              CALL ZEEX50('Reference',IBT,IDO,INP,NAN,
C     '          NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nx),
C     '          NPNE(1,1,ne),nr,nx,
C     '          DZDX,CE(1,ne,nx),CG,CP(1,1,nx),EG,PG,PHI,PST,
C     '          R,RGX,RI1,RI2,RI3,RM,U,
C     '          XE,XG,XIDL(1,nde),ZE,ZG,ERROR,*9999)
C              IF(NJT.EQ.3) CALL RALPHA3D(R,ALPHA,AXIS)
C              IF(NJT.EQ.2) CALL RALPHA2D(R,ALPHA,AXIS)
C              ALPHA=ALPHA*180.0d0/PI
CC             swap axis direction to point in +ve alpha direction
C              IF(ALPHA.LT.0.0)THEN
C                AXIS(1) = -AXIS(1)
C                AXIS(2) = -AXIS(2)
C                AXIS(3) = -AXIS(3)
C                ALPHA = -ALPHA
C              ENDIF
C              CALL ZEEX50('Fibre',IBT,IDO,INP,NAN,
C     '          NBH(1,1,ne),NBJ(1,ne),0,NHE(ne,nx),
C     '          NPNE(1,1,ne),nr,nx,
C     '          DZDX,CE(1,ne,nx),CG,CP(1,1,nx),EG,PG,PHI,PST,
C     '          R,RGX,RI1,RI2,RI3,RM,U,
C     '          XE,XG,XIDL(1,nde),ZE,ZG,ERROR,*9999)


CC LKC for AJP 191297
CC !!! PHI(3) changed to PHIL(3) in 2 places below

C              IF(NITB.EQ.3)THEN
C                WRITE(IUNIT,FMT=FORMAT)ND,(X(nj),nj=1,3),(Z(nj),nj=1,3),
C     '            (WD(nj,ND),nj=1,3),ne,(XIDL(NI,NDE),NI=1,3),ALPHA,
C     '            (AXIS(nj),nj=1,3),EG(1,1),EG(1,2),EG(1,3),EG(2,2),
C     '            EG(2,3),EG(3,3),(PST(nj),nj=1,3),(PHIL(nj),nj=1,3)
C              ELSE IF(NITB.EQ.2)THEN
C                WRITE(IUNIT,FMT=FORMAT)ND,(X(nj),nj=1,2),(Z(nj),nj=1,2),
C     '            (WD(nj,ND),nj=1,2),ne,(XIDL(NI,NDE),NI=1,2),ALPHA,
C     '            EG(1,1),EG(1,2),EG(2,2),
C     '            (PST(nj),nj=1,2),(PHIL(nj),nj=1,2)
C              ENDIF
C            ENDDO
C          ENDDO
C          CALL CLOSEF(IUNIT,ERROR,*9999)
!newe
!news   AAY 30 Apr 94  curvature data file
        ELSE IF(TYPE(1:4).EQ.'CURV') THEN
          CALL OPENF(IUNIT,'DISK',FILE(IBEG:IEND)//'.CURV','UNKNOWN',
     '      'SEQUEN','FORMATTED',360,ERROR,*9999)
          WRITE(IUNIT,'(A)')'  nd        Curv pos              Rect pos'
     '      //'      ne  xi1      xi2      xi3'
     '      //'      p1      p2'
     '      //'      p11     p22'

          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)

          FORMAT='(I5,6F7.2,I4,7F9.4)'
           DO NOLIST=1,NELIST(0)
            ne=NELIST(NOLIST)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA,XE,XP,ERROR,*9999)
            CALL ZDZDL(0,NBH(1,1,ne),NBJ(1,ne),NDDL,NDLT(ne),ne,1,WD,
     '        WDL,XID,XIDL,ZD,ZDL,ERROR,*9999)
            DO NDE=1,NDLT(ne)
              ND=NDDL(ne,NDE)
              DO nj=1,NJT
                nb=NBJ(nj,ne)
                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XIDL(1,NDE),XE(1,nj))
              ENDDO
              CALL XZ(ITYP10(nr),X,Z)
              CALL XECURV(IBT,IDO,INP,NBJ(1,ne),
     '          EG,Z2,XE,XIDL(1,NDE),ERROR,*9999)
              WRITE(IUNIT,FMT=FORMAT) ND,(X(nj),nj=1,3),(Z(nj),nj=1,3),
     '          ne,(XIDL(NI,NDE),NI=1,3),Z2(1),Z2(2),EG(1,1),EG(2,2)
            ENDDO
          ENDDO
          CALL CLOSEF(IUNIT,ERROR,*9999)
!newe
        ELSE IF(TYPE(1:6).EQ.'MAP3D') THEN
C          ERRORLEVEL=1
C          ERR=FCreateFile(FILE(IBEG:IEND)//'.data',DALBSPM,
C     '      ERRORLEVEL,THEFILE)
C          IF(err.ne.0) THEN
C            ERROR='>>Error creating the data file'
C            GOTO 9999
C          ENDIF
C          HEADERTEXT='Mapping data created with CMISS'
C          ERR=FSetText(THEFILE,HEADERTEXT)
C          IF(err.ne.0) THEN
C            ERROR='>>Error setting the header text'
C            GOTO 9999
C          ENDIF
C          ERR=FSetTimeSeriesIndex(THEFILE,SERIES)
C          IF(ERR.LT.0) THEN
C            ERROR='>>Error setting the series index'
C            GOTO 9999
C          ENDIF
C          HEADERTEXT='Series #'//CFROMI(SERIES,'(I1)')
C          ERR=FSetTimeSeriesLabel(THEFILE,HEADERTEXT)
C          IF(ERR.LT.0) THEN
C            ERROR='>>Error setting the time series header text'
C            GOTO 9999
C          ENDIF
C          ERR=FSetTimeSeriesFormat(THEFILE,DALBSPM)
C          IF(ERR.LT.0) THEN
C            ERROR='>>Error setting the time series format'
C            GOTO 9999
C          ENDIF
C          ERR=FSetTimeSeriesSpecs(THEFILE,NT_LEADS,NT_DATASETS)
C          IF(ERR.LT.0) THEN
C            ERROR='>>Error setting the time series specifications'
C            GOTO 9999
C          ENDIF
C          ERR=FSetTimeSeriesData(THEFILE,R4DATA2)
C          IF(ERR.LT.0) THEN
C            ERROR='>>Error setting the time series data'
C            GOTO 9999
C          ENDIF
C          ERR=FCloseFile(THEFILE)
C          IF(ERR.LT.0) THEN
C            ERROR='>>Error closing the data file'
C            GOTO 9999
C          ENDIF
        ELSE
C LKC 10-JUL-2003 Shifted the error checking into the command parsing
C   code and put out an explicit error
C         CO(noco+1)='?'
C         GO TO 1
          ERROR='>> Unknown file type: Update Code'
          GOTO 9999

        ENDIF
      ENDIF

      CALL EXITS('WRITEF')
      RETURN
 9999 CALL ERRORS('WRITEF',ERROR)
      CALL EXITS('WRITEF')
      RETURN 1
      END
