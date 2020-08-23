      SUBROUTINE DIHIST(IBT,IDO,INP,ISEG,ISHIST,NBH,NBJ,NEELEM,NHE,
     '  NHP,NHQ,NKH,NKHE,NKJE,NPF,NPLIST,NPNE,NPNODE,NPNY,NQNY,NRE,
     '  NRLIST,NRLIST2,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,NYNR,NYQNR,
     '  CURVCORRECT,SE,XA,XE,XID,XP,YP,YQ,YQS,ZA,ZE,ZP,CSEG,STRING,
     '  ERROR,*)

C#### Subroutine: DIHIST
C###  Description:
C###    DIHIST displays time history of nodal variables.

C**** Note: ISHIST(0) is segment containing axes and labels.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISEG(*),ISHIST(0:NPM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPLIST(0:NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NQNY(2,NYQM,0:NRCM,NXM),NRE(NEM),
     '  NRLIST(0:NRM),NRLIST2(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,
     '  INDEX_DATA,INDEX_AXES,INDEX_POLYLINE,N3CO,na,nb,nc,
     '  ne,NIQLIST(0:1),NIQSLIST(0:1),NIYLIST(0:16),niy,nj,nohist,
     '  NOLIST,no_niylist,nonr,nonrlist,no_time,np,nr,NTRL,
     '  NUMTIMEDATA,nx,nxc
      REAL*8 PXI,RFROMC,TIME,TMAX,TMIN,X(3),XHIST(5000),YHIST(5000),
     '  YMAX,YMIN,YPMAX(16),YPMIN(16),Z(3),Z2(3),RL(5)
      CHARACTER FILE*255,FILE2*255,ERROR1*255,FILEFORMAT*6
      LOGICAL ABBREV,ALL_REGIONS,BINARYFILE,CBBREV,
     '  ENDFILE,OUTFILE,STATIC,YQDATA,YQSDATA,YPDATA

      SAVE YMAX

      CALL ENTERS('DIHIST',*9999)
        CALL STRING_TRIM(FILE02,IBEG1,IEND1)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM display history<;FILENAME>
C###  Parameter:       <nodes (#s/all)[all]>
C###    The number of the nodes to be displayed
C###  Parameter:       <scale FACTOR#[1.0]>
C###    Specified the scale factor
C###  Parameter:       <rgb=RGB[blue]>
C###    Define colour (eg red,green,blue,cyan)
C###  Parameter:       <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:       <nc #[1]>
C###    Specify the variable to display
C###  Parameter:       <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:       <intoregion (#/same)[same]>
C###    Reads the solutions in to different regions
C###  Parameter:       <outfile OUTFILE>
C###    Specify the ascii data file to write the history to.
C###  Parameter:       <using (fit/solve)[solve]>
C###    Specify the nx type to use.
C###  Parameter:       <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Display the value of the solution variable vs time at specified
C###    nodes with the specified scale factor.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME['
     '    //FILE02(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<nodes (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<scale FACTOR#[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[blue]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<nc #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<intoregion (#/same)[same]>'
        OP_STRING(9)=BLANK(1:15)//'<outfile OUTFILE>'
        OP_STRING(10)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(11)=BLANK(1:15)//'<class #[1]>'

        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM display history dynamic
C###  Parameter:       <at ELEM#,XI1#,XI2#[1,0.0,0.0]>
C###  Parameter:       <scale FACTOR#[1.0]>
C###  Parameter:       <coord NJ#[1]>
C###  Parameter:       <rgb=RGB[blue]>
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' dynamic'
        OP_STRING(2)=BLANK(1:15)//'<at ELEM#,XI1#,XI2#[1,0.0,0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<scale FACTOR#[1.0]>'
        OP_STRING(4)=BLANK(1:15)//'<coord NJ#[1]>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[blue]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe26','doc','DIHIST',ERROR,*9999)
      ELSE
        nx=0
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '        ERROR,*9999)
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDIF
        ELSE
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ENDIF

        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        STATIC=.NOT.CBBREV(CO,'DYNAMIC',2,noco+1,NTCO,N3CO)

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX_DATA=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX_DATA=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
        ENDIF
        INDEX_AXES=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')

        CALL ACWK(10,0,ERROR,*9999)

        IF(STATIC)THEN
          IF(NTCOQU(noco).EQ.1) THEN
            FILE=COQU(noco,1)
          ELSE
            FILE=FILE02
          ENDIF
          CALL STRING_TRIM(FILE,IBEG,IEND)

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
C Moved above - before nxlist is used!
C          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)

          IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
            FACTOR=RFROMC(CO(N3CO+1))
          ELSE
            FACTOR=1.0d0
          ENDIF

          IF(CBBREV(CO,'BINARY',1,noco+1,NTCO,N3CO)) THEN
            BINARYFILE=.TRUE.
            FILEFORMAT='BINARY'
          ELSE
            BINARYFILE=.FALSE.
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
            nc=IFROMC(CO(N3CO+1))
          ELSE
            nc=1
          ENDIF

          IF(CBBREV(CO,'INTOREGION',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '        ERROR,*9999)
          ELSE
            NRLIST2(0)=NRLIST(0)
            DO nonr=1,NRLIST(0)
              NRLIST2(nonr)=NRLIST(nonr)
            ENDDO
          ENDIF

          IF(CBBREV(CO,'OUTFILE',3,noco+1,NTCO,N3CO)) THEN
            OUTFILE=.TRUE.
            CALL STRING_TRIM(CO(N3CO+1),IBEG2,IEND2)
            FILE2=CO(N3CO+1)(IBEG2:IEND2)
            IOFI=IOFILE3
            CALL OPENF(IOFI,'DISK',FILE2(IBEG2:IEND2)//'.dihist','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
C LKC 26-NOV-1998 Adding else
          ELSE
            OUTFILE=.FALSE.
          ENDIF

C***      Determine the maximum and minimum values
C Why was this set after parse regions? MLB
C          NRLIST(0)=1
          nr=NRLIST(1)
          NIYLIST(0)=1
          NIYLIST(1)=1
          NIQLIST(0)=0
          NIQSLIST(0)=0
          YPDATA=.TRUE.
          YQDATA=.FALSE.
          YQSDATA=.FALSE.
          na=1
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,FILE,'OPEN',
     '      ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          YMAX=-RMAX
          YMIN=RMAX
          TMIN=0.0d0
          IF(BINARYFILE) THEN
            DO no_time=1,NUMTIMEDATA

              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,
     '          'READ',FILEFORMAT,FILE,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              DO no_niylist=1,NIYLIST(0)
                niy=NIYLIST(no_niylist)
                CALL YPZP(niy,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),
     '            nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                DO nolist=1,NPLIST(0)
                  np=NPLIST(nolist)
                  IF(ZP(1,1,1,np,nc).GT.YMAX) YMAX=ZP(1,1,1,np,nc)
                  IF(ZP(1,1,1,np,nc).LE.YMIN) YMIN=ZP(1,1,1,np,nc)
                ENDDO !nolist
              ENDDO !no_niylist
            ENDDO !no_time
          ELSE
            ENDFILE=.FALSE.
            NUMTIMEDATA=0
            DO WHILE(.NOT.ENDFILE)

              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP,YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,
     '          FILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '          ERROR,*9999)
              IF(.NOT.ENDFILE) THEN
                NUMTIMEDATA=NUMTIMEDATA+1
                DO no_niylist=1,NIYLIST(0)
                  niy=NIYLIST(no_niylist)
                  CALL YPZP(niy,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '              NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),
     '              nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                  DO nolist=1,NPLIST(0)
                    np=NPLIST(nolist)
                    IF(ZP(1,1,1,np,nc).GT.YMAX) YMAX=ZP(1,1,1,np,nc)
                    IF(ZP(1,1,1,np,nc).LE.YMIN) YMIN=ZP(1,1,1,np,nc)
                  ENDDO !nolist
                ENDDO !no_niylist (niy)
              ENDIF
            ENDDO
          ENDIF
          TMAX=TIME
          CALL ASSERT(NUMTIMEDATA.LE.5000,'>>Too many points',ERROR,
     '      *9999)
          NTHIST=NUMTIMEDATA
C*** Loop over each node in the node list displaying the history
          DO nolist=1,NPLIST(0)
            np=NPLIST(nolist)
C*** Reset the file to the begining of the time data

            CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '        nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '        YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,FILE,
     '        'RESET',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

            DO nohist=1,NUMTIMEDATA

              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,
     '          'READ',FILEFORMAT,FILE,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)

              DO nonrlist=1,NRLIST(0)
                nr=NRLIST(nonrlist)
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),
     '            nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
              ENDDO !nonrlist

              XHIST(nohist)=TIME
              YHIST(nohist)=ZP(1,1,1,np,nc)*FACTOR

            ENDDO !nohist

C*** Plot history for node np
            CALL SGHIST(INDEX_DATA,ISEG,ISHIST(np),np,NTHIST,CSEG,TMAX,
     '        TMIN,XHIST,YHIST,YMAX,YMIN,ERROR,*9999)

            IF(OUTFILE) THEN
              DO nohist=1,NUMTIMEDATA
                WRITE(IOFI,*)XHIST(nohist),YHIST(nohist)
              ENDDO
            ENDIF

          ENDDO !nolist

          IF(OUTFILE) THEN
            CALL CLOSEF(IOFI,ERROR,*9999)
          ENDIF

C*** Draw labels
          CALL SGHIST(INDEX_AXES,ISEG,ISHIST(0),0,0,CSEG,TMAX,TMIN,
     '      XHIST,YHIST,YMAX,YMIN,ERROR,*9999)

C*** Close history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,FILE,' ',
     '      ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

        ELSE !plot time behaviour

C **      Draw time plots
          CALL YPZP(MOTION_IY,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '      NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),4,NTRL,RL,ERROR,*9999)
            NE=INT(RL(1))
          ELSE
            NE=1
            RL(2)=0.0D0
            RL(3)=0.0D0
          ENDIF
          IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
            FACTOR=RFROMC(CO(N3CO+1))
          ELSE
            FACTOR=1.0D0
          ENDIF
          IF(CBBREV(CO,'COORD',3,noco+1,NTCO,N3CO)) THEN
            NJHIST=IFROMC(CO(N3CO+1))
          ELSE
            NJHIST=1
          ENDIF
          T1=XID(NXIDEF,NDT) !take final time as that of last data pt
          YMIN=1000000.0d0
          YMAX=-1000000.0d0
          CALL XPXE(NBJ(1,NE),NKJE(1,1,1,NE),NPF(1,1),
     '      NPNE(1,1,NE),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,NE),XA(1,1,ne),XE,XP,
     '      ERROR,*9999)
          CALL ZPZE(NBH(1,1,NE),nc,NHE(ne,nx),NKHE(1,1,1,NE),
     '      NPF(1,1),NPNE(1,1,NE),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,NE),ZA(1,1,1,NE),
     '      ZE,ZP,ERROR,*9999)
          DO nj=1,NJT
            nb=NBJ(nj,NE)
            X(nj)=PXI(IBT(1,1,NB),IDO(1,1,0,NB),INP(1,1,NB),
     '        NB,1,RL(2),XE(1,NJ))
          ENDDO
          CALL XZ(ITYP10(nr),X,Z)
          NTHIST=100          !plot 100 points
          DO NOHIST=1,NTHIST
            XHIST(NOHIST)=DBLE(NOHIST-1)*T1/DBLE(NTHIST-1)  !time
            RL(NIT(NBJ(NJHIST,NE))+2)=XHIST(NOHIST)
            NJ=NJHIST
            nb=NBH(nj,1,NE)
            Z2(nj)=PXI(IBT(1,1,NB),IDO(1,1,0,NB),INP(1,1,NB),
     '        NB,1,RL(2),ZE(1,NJ))
            IF(KTYP58(nr).EQ.2)THEN !displacement
              Z2(nj)=Z(nj)+Z2(nj)
            ELSE
              Z2(nj)=Z2(nj)
            ENDIF
            YHIST(NOHIST)=Z2(nj)                              !position
            IF(DOP) THEN
              WRITE(OP_STRING,*)' z,z2= ',Z(NJHIST),Z2(NJHIST)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(YHIST(NOHIST).GT.YMAX) YMAX=YHIST(NOHIST)
            IF(YHIST(NOHIST).LT.YMIN) YMIN=YHIST(NOHIST)
            XHIST(NOHIST)=XHIST(NOHIST)
          ENDDO
          TMIN=0.0d0
          TMAX=T1
          YMAX=YMAX+0.01D0*DMAX1(YMAX-YMIN,1.0D0)*FACTOR
          YMIN=YMIN-0.01D0*DMAX1(YMAX-YMIN,1.0D0)*FACTOR
          IF(DOP)THEN
            WRITE(OP_STRING,*)' xhist   yhist'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO NOHIST=1,NTHIST
                WRITE(OP_STRING,*) XHIST(NOHIST),YHIST(NOHIST)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
          YDMIN=YMIN
          YDMAX=YMAX
          CALL SGHIST(INDEX_DATA,ISEG,ISHIST(ne),ne,NTHIST,CSEG,TMAX,
     '      TMIN,XHIST,YHIST,YMAX,YMIN,ERROR,*9999)

C **      Draw labels
          CALL SGHIST(INDEX_AXES,ISEG,ISHIST(0),0,0,CSEG,TMAX,TMIN,
     '      XHIST,YHIST,YMAX,YMIN,ERROR,*9999)

        ENDIF
        CALL DAWK(10,0,ERROR,*9999)
      ENDIF

      CALL EXITS('DIHIST')
      RETURN
 9999 IF(nx.GT.0) THEN
        CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '    NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA,
     '    nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP,YPMAX,YPMIN,
     '    YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,FILE,' ',ENDFILE,.TRUE.,
     '    YPDATA,YQDATA,YQSDATA,ERROR1,*9998)
      ENDIF
 9998 IF(OUTFILE) CALL CLOSEF(IOFI,ERROR1,*1111)

 1111 CALL ERRORS('DIHIST',ERROR)
      CALL EXITS('DIHIST')
      RETURN 1
      END


