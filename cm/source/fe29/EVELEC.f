      SUBROUTINE EVELEC(IBT,IDO,INP,ISIZE_PHIH,ISIZE_MFI,
     '  LD,LIST,NBH,NBJ,NDDATA,NEELEM,NELIST,NENP,
     '  NENQ,NHE,NHP,NHQ,NKH,NKHE,NKJE,NPF,NP_INTERFACE,NPLIST,NPLIST3,
     '  NPLIST4,NPNE,NPNODE,NPNY,NQLIST,NQNE,NQNY,NQS,NQSCNB,NQXI,NRE,
     '  NRLIST,NRLIST2,NVHE,NVJE,NVHP,NW,NXLIST,NYNE,NYNP,NYNR,NYQNR,
     '  CURVCORRECT,MFI,PHI_H,SE,WD,XA,XE,XID,XP,XQ,YP,YQ,YQS,ZA,ZD,ZE,
     &  ZP,STRING,ERROR,*)

C#### Subroutine: EVELEC
C###  Description:
C###    <HTML>
C###    EVELEC evaluates electrodes (data) from the dependent
C###    variable field or other arrays and ouputs the results
C###    to a (cmiss) signal file.
C###    <BR>
C###    <B>WARNING: </B>There are inconsistencies with data points and
C###    regions as each electrode is considered to be a data
C###    point. All points are placed in region 1 at present.
C###    </HTML>

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'

      ! Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISIZE_PHIH(2),ISIZE_MFI(3,NSSM),LD(NDM),LIST(0:NLISTM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NENQ(0:8,NQM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPLIST(0:NPM),NPLIST3(0:NPM),NPLIST4(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM,NXM),NQLIST(0:NQM),
     '  NQNE(NEQM,NQEM),NQNY(2,NYQM,0:NRCM,NXM),NQS(NEQM),NQSCNB(NQSCM),
     '  NQXI(0:NIM,NQSCM),NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVHP(NHM,NPM,NCM,0:NRM),NW(NEM,3,NXM),
     '  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),MFI(NDM,NTSM,3,NSSM),
     &  PHI_H(NY_TRANSFER_M,NTSM),
     '  SE(NSM,NBFM,NEM),WD(NJM,NDM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IY_YQSM,GP_YQSM
      PARAMETER(IY_YQSM=250)
      PARAMETER(GP_YQSM=1)

      INTEGER i,IBEG,IBEG1,IEND,IEND1,N3CO,na,nb,nd,ne,ni,NITB,niqV,
     '  NIQLIST(0:1),NIQSLIST(0:IY_YQSM),NIYLIST(0:16),nj,
     '  nodata,noelem,nolist,nonode,nonr,
     '  nonrlist,notime,np,nq,nqe,nqe1,nqe2,nqe3,nqq,nr,nss,nts,
     &  NUMTIMEDATA,NUMTIMEDATA1,nx,nxc,nytr,PHI_LIMIT,SCHEME,
     '  IY,IFROMC,IY_YQS(IY_YQSM),IY_YQST,GP_YQS(GP_YQSM),GP_YQST
      REAL*8 CALC_TIME_FROM_SAMPLE,PXI,SIGNALMAX(9),SIGNALMIN(9),
     '  RFROMC,TEND,TIME,TSTART,XI_3,YPMAX(16),YPMIN(16)
      CHARACTER ERROR1*255,FILEFORMAT*6,HISTORYFNAME*(MXCH),
     '  SIGNALFILE*(MXCH),TYPE*14,
     '  VARIABLE*3
      LOGICAL ABBREV,ALL_REGIONS,AT_DATA,CBBREV,
     '  ENDFILE,ENDFILE2,FOUND_NODE,INLIST,OPFILE,YPDATA,YQDATA,
     '  YQSDATA

      CALL ENTERS('EVELEC',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------
C#### Command: FEM evaluate electrodes<;FILENAME>
C###  Parameter:        <history FILENAME[$current]>
C###    Specify the history filename to use to evaluate the electrodes
C###    from.
C###  Parameter: <from (data/extended_nodes/grid/gvariables/nodes)[data]>
C###    Specify where the geometric locations of the electrodes are
C###    to be obtained from. If "data" is selected the electrodes
C###    correspond to the currently defined data points in the specified
C###    regions. If "nodes" is selected the electrode locations
C###    correspond to the currently defined node locations in the
C###    specified region. If "extended_nodes" is selected the electrode
C###    locations correspond to the currently defined nodal locations
C###    (as for "node") plus an extra location at the center of each
C###    element. If "grid" is selected the electrode locations
C###    correspond to the collocation grid locations in the specified
C###    regions.
C###  Parameter:        <electrodes (#s/all)[all]>
C###    Specify the electrode numbers to be evaluated.
C###  Parameter:        <(yq/yqs)[yq]>
C###    From problems evaluated from grids or grid_varables, specify
C###    whether YQ or YQS is to be used.
C###  Parameter:        <grid_points (#s/all)[1]>
C###    If you are evaluating electrodes from grid points and the YQS
C###    array, then this parameter specifies which grid points to use.
C###  Parameter:        <iy (#s/all)[1]>
C###    Specify the iy number to be used in index the YP or YQ
C###    arrays when evaluating the electrodes. If you are using YQS then
C###    multiple iy's can be specified.
C###  Parameter:        <tstart (#/beginning)[beginning]>
C###    Specify start time
C###  Parameter:        <tend (#/end)[end]>
C###    Specify end time
C###  Parameter:        <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary (.binsig) or
C###    ascii (.ipsign) file.
C###  Parameter:        <coupled>
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <intoregion (#s/all)[same]>
C###  Parameter:        <using (fit/solve)[solve]>
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Creates a signal file from a history file by evaluating the
C###    history file at specfied 'electrode' locations. The electrodes
C###    can either by specified by data points, nodes or extended_nodes
C###    (node locations plus an 'extra' node in the middle of the
C###    element).

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<history FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<from (data/extended_nodes/grid/'
     '    //'gvariables/nodes)[data]>'
        OP_STRING(4)=BLANK(1:15)//'<(yq/yqs)[yq]>'
        OP_STRING(5)=BLANK(1:15)//'<*grid_points (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<iy (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<electrodes (#s/all)[all]>'
        OP_STRING(8)=BLANK(1:15)//'<tstart (#/beginning)[beginning]>'
        OP_STRING(9)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(10)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(11)=BLANK(1:15)//'<coupled>'
        OP_STRING(12)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(13)=BLANK(1:15)//'<intoregion (#s/same)[same]>'
        OP_STRING(14)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(15)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C--------------------------------------------------------------------

C#### Command: FEM evaluate electrodes<;FILENAME>
C###  Parameter:        <from (PHI_H/MFI)[-]>
C###  Parameter:        <tstart (#/beginning)[beginning]>
C###    Specify start time
C###  Parameter:        <tend (#/end)[end]>
C###    Specify end time
C###  Parameter:        <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <intoregion (#s/same)[same]>
C###  Parameter:        <using (fit/solve)[solve]>
C###  Parameter:        <at (nodes/data)[nodes]>
C###    Specify the electrodes locations using node or data points.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <nss #[1]>
C###    Specify the signal set (for MFI) to use.
C###  Description:
C###    Creates a signal file from the inversely computed PHI_H array.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<from PHI_H>'
        OP_STRING(3)=BLANK(1:15)//'<tstart (#/beginning)[beginning]>'
        OP_STRING(4)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<intoregion (#s/same)[same]>'
        OP_STRING(8)=BLANK(1:15)//'<using (fit/solve)[solve]>'
        OP_STRING(9)=BLANK(1:15)//'<at (nodes/data)[nodes]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(11)=BLANK(1:15)//'<nss #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVELEC',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          SIGNALFILE=COQU(noco,1)
          OPFILE=.TRUE.
          IOFI=IOFILE1
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF


        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

C LKC 13-JAN-98 Changes the goto tag to *9997 for these as if nx=0
C               crashes when it tries to close the HISTORY and access
C               the array NYNR(0,0,1,0,nx)
        IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this fit class',
     '        ERROR,*9997)
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9997)
          ENDIF
        ELSE
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9997)
        ENDIF

C LKC needs to be down here so nx is allocated
        CALL ASSERT(NJT+1.LE.NJM,'>>Increase NJM',ERROR,*9999)

        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
          IF(ABBREV(CO(N3CO+1),'EXTENDED_NODES',2)) THEN
            TYPE='EXTENDED_NODES'
          ELSE IF(ABBREV(CO(N3CO+1),'GRID',2)) THEN
            TYPE='GRID'
          ELSE IF(ABBREV(CO(N3CO+1),'GVARIABLES',2)) THEN
            TYPE='GVAR'
          ELSE IF(ABBREV(CO(N3CO+1),'MFI',3)) THEN
            TYPE='MFI'            
            IF(CBBREV(CO,'NSS',3,noco+1,NTCO,N3CO)) THEN
              nss=IFROMC(CO(N3CO+1))
            ELSE
              nss=1
            ENDIF
          ELSE IF(ABBREV(CO(N3CO+1),'NODES',2)) THEN
            TYPE='NODES'
          ELSE IF(ABBREV(CO(N3CO+1),'PHI_H',2)) THEN
            TYPE='PHI_H'
C LKC adding the option to specify electro positions
            IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
              IF(ABBREV(CO(N3CO+1),'DATA',2)) THEN
                AT_DATA=.TRUE.
              ELSE
                AT_DATA=.FALSE.
              ENDIF
            ELSE
              AT_DATA=.FALSE.
            ENDIF
          ELSE
            TYPE='DATA'
          ENDIF
        ELSE
          TYPE='DATA'
        ENDIF

        
C *** DPN 15 February 2000 - need to do this after we know what array
C *** we are using
C        IF(CBBREV(CO,'IY',2,noco+1,NTCO,N3CO)) THEN
C          iy=IFROMC(CO(N3CO+1))
C        ELSE
C          iy=1
C        ENDIF

        IF(TYPE(1:4).EQ.'GRID'.OR.TYPE(1:4).EQ.'GVAR') THEN
          IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
            VARIABLE='YQS'
          ELSE
            VARIABLE='YQ'
          ENDIF
        ELSE
          VARIABLE='YP'
        ENDIF

        IF(CBBREV(CO,'ELECTRODES',2,noco+1,NTCO,N3CO)) THEN

C*** Evaluate from DATA
          IF(TYPE(1:4).EQ.'DATA') THEN
            CALL ASSERT(CALL_XI,'>> Calculate XI positions',
     '        ERROR,*9999)

            CALL PARSIL(CO(N3CO+1),NDM,LIST(0),LIST(1),ERROR,*9999)
            DO nolist=1,LIST(0)
              nd=LIST(nolist)
              NDDATA(nolist,1)=nd
              CALL ASSERT(nd.GT.0.AND.nd.LE.NDT,
     '          '>>Data point not defined',ERROR,*9999)
C LKC 21-JUN-1998 WD for electrode not set
              WD(NJT+1,nd)=1.0d0

C LKC 24-APR-1998 Check there is an element assoc. with nd
              IF(LD(nd).EQ.0) THEN
                OP_STRING(1)='WARNING: No element for data point'
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF

            ENDDO !nonlist (nd)
            NDDATA(0,1)=LIST(0)
C*** Evaluate from EXTENEDED NODES
          ELSE IF(TYPE(1:14).EQ.'EXTENDED_NODES') THEN
            CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),
     '        ERROR,*9999)
            CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '        *9999)
            nd=0
            NELIST(0)=0
            DO ne=1,NEM
              NELIST(ne)=0
            ENDDO
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              FOUND_NODE=.FALSE.
              DO nonrlist=1,NRLIST(0)
                nr=NRLIST(nonrlist)
                IF(INLIST(np,NPNODE(1,nr),NPNODE(0,nr),i)) THEN
                  FOUND_NODE=.TRUE.
                ENDIF
              ENDDO
              CALL ASSERT(FOUND_NODE,
     '          '>>Node number not defined',ERROR,*9999)
              nd=nd+1
              IF(nd.LE.NDM) THEN
                NDDATA(nd,1)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0
                CALL CALC_NP_XI(IBT,INP,NBJ,ne,NENP,np,NPNE,NRLIST,
     '            XID(1,nd),ERROR,*9999)
                LD(nd)=ne
                DO nj=1,NJT
                  ZD(nj,nd)=XP(1,1,nj,np)
                ENDDO !nj
                IF(NIM.GE.3) XI_3=XID(3,nd)
                DO nonrlist=1,NRLIST(0)
                  nr=NRLIST(nonrlist)
                  DO noelem=1,NENP(np,0,nr)
                    ne=NENP(np,noelem,nr)
                    IF((.NOT.INLIST(ne,NELIST(1),NELIST(0),i))
     '                .AND.(NRE(ne).EQ.NP_INTERFACE(np,1))) THEN
                      nd=nd+1
                      IF(nd.LE.NDM) THEN
                        NDDATA(nd,1)=nd
                        DO nj=1,NJT
                          WD(nj,nd)=1.0d0
                        ENDDO
                        WD(NJT+1,nd)=1.0d0
                        DO ni=1,NIM
                          XID(ni,nd)=0.5d0
                        ENDDO
                        IF(NIM.GE.3) XID(3,nd)=XI_3
                        LD(nd)=ne
                        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '                    NPF(1,1),NPNE(1,1,ne),
     '                    NRE(ne),NVJE(1,1,1,ne),
     '                    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                        DO nj=1,NJT
                          nb=NBJ(nj,ne)
                          ZD(nj,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,XID(1,nd),XE(1,nj))
                        ENDDO !nj
                        NELIST(0)=NELIST(0)+1
                        NELIST(NELIST(0))=ne
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)
C*** Evaluate from GRID points
          ELSE IF(TYPE(1:4).EQ.'GRID') THEN
            CALL PARSIL(CO(N3CO+1),NQT,NQLIST(0),NQLIST(1),
     '        ERROR,*9999)
            CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '        *9999)
            CALL ASSERT(NQLIST(0).LE.NDM,'>>Increase NDM',ERROR,
     '        *9999)
            CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
            CALL ASSERT(niqV.NE.0,'>>Invalid niqV',ERROR,*9999)
            nd=0
            DO nolist=1,NQLIST(0)
              nq=NQLIST(nolist)
              nd=nd+1
              IF(nd.LE.NDM) THEN
                NDDATA(nd,1)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0

                ne=NENQ(1,nq)
                CALL ASSERT(ne.GT.0,'>>No element found for grid point',
     '            ERROR,*9999)

                LD(nd)=ne
                SCHEME=NQS(ne)
                NITB=NIT(NQSCNB(SCHEME))
C*** Find XID
                IF(NITB.EQ.1) THEN
                  DO nqe1=1,NQXI(1,SCHEME)
                    nqq=NQNE(ne,nqe1)
                    IF(nqq.EQ.nq) THEN
                      XID(1,nd)=DBLE(nqe1-1)/DBLE(NQXI(1,SCHEME)-1)
                    ENDIF
                  ENDDO
                ELSEIF(NITB.EQ.2) THEN
                  DO nqe2=1,NQXI(2,SCHEME)
                    DO nqe1=1,NQXI(1,SCHEME)
                      nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))
                      nqq=NQNE(ne,nqe)
                      IF(nqq.EQ.nq) THEN
                        XID(1,nd)=DBLE(nqe1-1)/DBLE(NQXI(1,SCHEME)-1)
                        XID(2,nd)=DBLE(nqe2-1)/DBLE(NQXI(2,SCHEME)-1)
                      ENDIF
                    ENDDO
                  ENDDO
                ELSEIF(NITB.EQ.3) THEN
                  DO nqe3=1,NQXI(3,SCHEME)
                    DO nqe2=1,NQXI(2,SCHEME)
                      DO nqe1=1,NQXI(1,SCHEME)
                        nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))+((nqe3-1)*
     '                    NQXI(1,SCHEME)*NQXI(2,SCHEME))
                        nqq=NQNE(ne,nqe)
                        IF(nqq.EQ.nq) THEN
                          XID(1,nd)=DBLE(nqe1-1)/DBLE(NQXI(1,SCHEME)-1)
                          XID(2,nd)=DBLE(nqe2-1)/DBLE(NQXI(2,SCHEME)-1)
                          XID(3,nd)=DBLE(nqe3-1)/DBLE(NQXI(3,SCHEME)-1)
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ELSE
                  CALL ASSERT(.FALSE.,'>>Incorrect # of xi directions',
     '              ERROR,*9999)
                ENDIF

                DO nj=1,NJT
                  ZD(nj,nd)=XQ(nj,nq)
                ENDDO !nj
                LIST(nd)=nq
              ENDIF
            ENDDO !nonlist (nq)

            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)

C*** Evaluate from GRID VARIABLES points
          ELSE IF(TYPE(1:4).EQ.'GVAR') THEN
            CALL PARSIL(CO(N3CO+1),NDM,LIST(0),LIST(1),ERROR,*9999)
            CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '        *9999)
            IF(VARIABLE(1:3).EQ.'YQS') THEN
              CALL ASSERT(CALL_CELL,'>>Must define a cell first',
     '          ERROR,*9999)
              NDT=NIQST
              NDDATA(0,1)=NDT
              CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)
              DO nd=1,NDT
                NDDATA(nd,1)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0d0
                  ZD(nj,nd)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0
                DO ni=1,NIM
                  XID(ni,nd)=0.0d0
                ENDDO
                LD(nd)=1
              ENDDO
            ELSE
              ERROR='>>Not implemented'
              GOTO 9999
            ENDIF

C*** Evaluate from NODES
          ELSE IF(TYPE(1:5).EQ.'NODES') THEN
            CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),
     '        ERROR,*9999)
            CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '        *9999)
            CALL ASSERT(NPLIST(0).LE.NDM,'>>Increase NDM',ERROR,
     '        *9999)
            nd=0
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              FOUND_NODE=.FALSE.
              DO nonrlist=1,NRLIST(0)
                nr=NRLIST(nonrlist)
                IF(INLIST(np,NPNODE(1,nr),NPNODE(0,nr),i)) THEN
                  FOUND_NODE=.TRUE.
                ENDIF
              ENDDO
              CALL ASSERT(FOUND_NODE,
     '          '>>Node number not defined',ERROR,*9999)
              nd=nd+1
              IF(nd.LE.NDM) THEN
                NDDATA(nd,1)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0
                CALL CALC_NP_XI(IBT,INP,NBJ,ne,NENP,np,NPNE,NRLIST,
     '            XID(1,nd),ERROR,*9999)
                LD(nd)=ne
                DO nj=1,NJT
                  ZD(nj,nd)=XP(1,1,nj,np)
                ENDDO !nj
              ENDIF
            ENDDO !nonlist (np)
            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)
          ENDIF

C*** Default to use all electrodes
        ELSE
          IF(TYPE(1:4).EQ.'DATA') THEN
            NDDATA(0,1)=NDT
            DO nd=1,NDT
              NDDATA(nd,1)=nd
C LKC 21-JUN-1998 WD for electrode not set
              WD(NJT+1,nd)=1.0d0
            ENDDO !nd
          ELSE IF(TYPE(1:14).EQ.'EXTENDED_NODES') THEN
            nd=0
            NELIST(0)=0
            DO ne=1,NEM
              NELIST(ne)=0
            ENDDO
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '          *9999)
              CALL ASSERT(NPNODE(0,nr).LE.NDM,'>>Increase NDM',ERROR,
     '          *9999)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                nd=nd+1
                IF(nd.LE.NDM) THEN
                  NDDATA(nd,1)=nd
                  DO nj=1,NJT
                    WD(nj,nd)=1.0d0
                  ENDDO
                  WD(NJT+1,nd)=1.0d0
                  CALL CALC_NP_XI(IBT,INP,NBJ,ne,NENP,np,NPNE,NRLIST,
     '              XID(1,nd),ERROR,*9999)
                  LD(nd)=ne
                  DO nj=1,NJT
                    ZD(nj,nd)=XP(1,1,nj,np)
                  ENDDO !nj
                  IF(NIM.GE.3) XI_3=XID(3,nd)
                  DO noelem=1,NENP(np,0,nr)
                    ne=NENP(np,noelem,nr)
                    IF((.NOT.INLIST(ne,NELIST(1),NELIST(0),i))
     '                .AND.(NRE(ne).EQ.NP_INTERFACE(np,1))) THEN
                      nd=nd+1
                      IF(nd.LE.NDM) THEN
                        NDDATA(nd,1)=nd
                        DO nj=1,NJT
                          WD(nj,nd)=1.0d0
                        ENDDO
                        WD(NJT+1,nd)=1.0d0
                        DO ni=1,NIM
                          XID(ni,nd)=0.5d0
                        ENDDO
                        IF(NIM.GE.3) XID(3,nd)=XI_3
                        LD(nd)=ne
                        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),
     '                    NPF(1,1),NPNE(1,1,ne),
     '                    NRE(ne),NVJE(1,1,1,ne),
     '                    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                        DO nj=1,NJT
                          nb=NBJ(nj,ne)
                          ZD(nj,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,1,XID(1,nd),XE(1,nj))
                        ENDDO !nj
                        NELIST(0)=NELIST(0)+1
                        NELIST(NELIST(0))=ne
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)

          ELSE IF(TYPE(1:4).EQ.'GRID') THEN
            CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '        *9999)
            CALL ASSERT(NQT.LE.NDM,'>>Increase NDM',ERROR,
     '        *9999)
            CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)

            nd=0
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              DO nq=NQR(1,nr),NQR(2,nr)
                nd=nd+1
                IF(nd.LE.NDM) THEN
                  NDDATA(nd,1)=nd
                  DO nj=1,NJT
                    WD(nj,nd)=1.0d0
                  ENDDO
                  WD(NJT+1,nd)=1.0d0

                  ne=NENQ(1,nq)
                  CALL ASSERT(ne.GT.0,
     '              '>>No element found for grid point',ERROR,*9999)

                  LD(nd)=ne
                  SCHEME=NQS(ne)
                  NITB=NIT(NQSCNB(SCHEME))
C*** Find XID
                  IF(NITB.EQ.1) THEN
                    IF(NQXI(1,SCHEME).GT.1) THEN
                      DO nqe1=1,NQXI(1,SCHEME)
                        nqq=NQNE(ne,nqe1)
                        IF(nqq.EQ.nq) THEN
                          XID(1,nd)=DBLE(nqe1-1)/DBLE(NQXI(1,SCHEME)-1)
                        ENDIF
                      ENDDO !nqq
                    ELSE !NQXI(1,SCHEME)==0
                      XID(1,nd)=0.0d0
                    ENDIF !NQXI(1,SCHEME)
                  ELSEIF(NITB.EQ.2) THEN
                    DO nqe2=1,NQXI(2,SCHEME)
                      DO nqe1=1,NQXI(1,SCHEME)
                        nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))
                        nqq=NQNE(ne,nqe)
                        IF(nqq.EQ.nq) THEN
                          XID(1,nd)=DBLE(nqe1-1)/DBLE(NQXI(1,SCHEME)-1)
                          XID(2,nd)=DBLE(nqe2-1)/DBLE(NQXI(2,SCHEME)-1)
                        ENDIF
                      ENDDO
                    ENDDO
                  ELSEIF(NITB.EQ.3) THEN
                    DO nqe3=1,NQXI(3,SCHEME)
                      DO nqe2=1,NQXI(2,SCHEME)
                        DO nqe1=1,NQXI(1,SCHEME)
                          nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))+((nqe3-1)*
     '                      NQXI(1,SCHEME)*NQXI(2,SCHEME))
                          nqq=NQNE(ne,nqe)
                          IF(nqq.EQ.nq) THEN
                            XID(1,nd)=DBLE(nqe1-1)/
     '                        DBLE(NQXI(1,SCHEME)-1)
                            XID(2,nd)=DBLE(nqe2-1)/
     '                        DBLE(NQXI(2,SCHEME)-1)
                            XID(3,nd)=DBLE(nqe3-1)/
     '                        DBLE(NQXI(3,SCHEME)-1)
                          ENDIF
                        ENDDO
                      ENDDO
                    ENDDO
                  ELSE
                    CALL ASSERT(.FALSE.,
     '                '>>Incorrect # of xi directions',ERROR,*9999)
                  ENDIF

                  DO nj=1,NJT
                    ZD(nj,nd)=XQ(nj,nq)
                  ENDDO !nj
                  LIST(nd)=nq
                ENDIF
              ENDDO
            ENDDO !nonlist (nq)

            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)

C*** Evaluate from GRID VARIABLES points
          ELSE IF(TYPE(1:4).EQ.'GVAR') THEN
            LIST(0)=1
            LIST(1)=NQR(1,nr)
            CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '        *9999)
            IF(VARIABLE(1:3).EQ.'YQS') THEN
              CALL ASSERT(CALL_CELL,'>>Must define a cell first',
     '          ERROR,*9999)
              NDT=NIQST
              NDDATA(0,1)=NDT
              CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)
              DO nd=1,NDT
                NDDATA(nd,1)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0d0
                  ZD(nj,nd)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0
                DO ni=1,NIM
                  XID(ni,nd)=0.0d0
                ENDDO
                LD(nd)=1
              ENDDO
            ELSE
              ERROR='>>Not implemented'
              GOTO 9999
            ENDIF

          ELSE IF(TYPE(1:5).EQ.'NODES') THEN
            nd=0
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              CALL ASSERT(USE_DATA.GT.0,'>>Set USE_DATA to 1',ERROR,
     '          *9999)
              CALL ASSERT(NPNODE(0,nr).LE.NDM,'>>Increase NDM',ERROR,
     '          *9999)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                nd=nd+1
                IF(nd.LE.NDM) THEN
                  NDDATA(nd,1)=nd
                  DO nj=1,NJT
                    WD(nj,nd)=1.0d0
                  ENDDO
                  WD(NJT+1,nd)=1.0d0
                  CALL CALC_NP_XI(IBT,INP,NBJ,ne,NENP,np,NPNE,NRLIST,
     '              XID(1,nd),ERROR,*9999)
                  LD(nd)=ne
                  DO nj=1,NJT
                    ZD(nj,nd)=XP(1,1,nj,np)
                  ENDDO !nj
                ENDIF
              ENDDO !nonode
            ENDDO !nr
            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)


          ELSE IF(TYPE(1:3).EQ.'MFI') THEN

C LKC 13-APR-2006 Exporting MFI array to a signal file. Need
C             to export 3 components * #electrodes            
            CALL ASSERT(3*ISIZE_MFI(1,nss).LE.NDM,
     &        '>>Increase NDM to be 3 times #sensors',ERROR,*9999)
            nd=0
            DO nj=1,3
              DO nodata=1,ISIZE_MFI(1,nss)
                nd=nd+1
                NDDATA(nd,1)=nd
                DO ni=1,NJT
                  WD(nj,ni)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0
                LD(nd)=1
                ZD(nj,nd)=ZD(nj,nodata) !the first nelec will be duplicates                
              ENDDO !nelec
            ENDDO !nj
            NDT=nd
            NDDATA(0,1)=NDT

            
          ELSE IF(TYPE(1:5).EQ.'PHI_H') THEN
C            CALL ASSERT(EVALUATE_INVERSE.OR.EVALUATE_PHI_NEAREST
C     '        ,'>>Evaluate inverse or phi/nearest  first',ERROR,*9999)
C JMB 15/11/00 Setup up appropriate arrays for potential inverse
C            IF(EVALUATE_INVERSE) THEN
C GBS 31-Jan-2001 Added electrode coordinates using NPLIST3
C              WRITE(OP_STRING,'('' Warning: Electrode coordinates not '
C     '          //'calculated correctly yet'')')
C              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            nd=0

C LKC 6-SEP-2002 I think don't this is appropriate - say i want to
C export a signal from PHI_H without having readin the T_BH
C information yet... Should be using the ISIZE_PHIH array instead
C
CC GBS 12-Nov-2001 Added appropriate sizing of array
C              IF (PHI_H_IS_HEART) THEN
C                PHI_LIMIT=ISIZE_TBH(2) !number of heart nodes
C              ELSE
C                PHI_LIMIT=ISIZE_TBH(1) !number of torso nodes
C           ENDIF

            PHI_LIMIT=ISIZE_PHIH(1)

            DO nytr=1,PHI_LIMIT
              nd=nd+1
              IF(nd.LE.NDM) THEN
                NDDATA(nd,1)=nd
                DO nj=1,NJT
                  WD(nj,nd)=1.0d0
                ENDDO
                WD(NJT+1,nd)=1.0d0
                LD(nd)=1

C LKC 10-SEP-2002 adding the option to skip this section if we
C                 don't want to associate the electrodes with nodes

                IF(.NOT.AT_DATA) THEN
C GBS 31-Jan-2001 Added electrode coordinates from first surface nodes
                  IF (PHI_H_IS_HEART) THEN
                    np=NPLIST3(nd)
                  ELSE
                    np=NPLIST4(nd)
                  ENDIF
                  DO nj=1,NJT
                    ZD(nj,nd)=XP(1,1,nj,np)
                  ENDDO !nj
                ENDIF !AT_DATA
              ENDIF !nd
            ENDDO !nonode
            NDT=nd
            NDDATA(0,1)=NDT
            CALL ASSERT(NDT.LE.NDM,'>>Increase NDM',ERROR,*9999)
C            ENDIF
          ENDIF !type
        ENDIF !electrodes

C        IF(TYPE(1:4).EQ.'GRID') THEN
C          IF(CBBREV(CO,'YQS',3,noco+1,NTCO,N3CO)) THEN
C            VARIABLE='YQS'
C          ELSE
C            VARIABLE='YQ'
C          ENDIF
C        ENDIF

C *** DPN 18 February 2000 - Want to enable multiple variables for
C *** grid/yqs cases.
        IF(TYPE(1:4).EQ.'GRID'.AND.VARIABLE(1:3).EQ.'YQS') THEN
          IF(CBBREV(CO,'IY',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),IY_YQSM,IY_YQST,IY_YQS,ERROR,*9999)
          ELSE
            IY_YQST=1
            IY_YQS(1)=1
          ENDIF
        ELSE !not using grid and yqs
          IF(CBBREV(CO,'IY',2,noco+1,NTCO,N3CO)) THEN
            iy=IFROMC(CO(N3CO+1))
          ELSE
            iy=1
          ENDIF
        ENDIF !IF(TYPE(1:4).EQ.'GRID'.AND.VARIABLE(1:3).EQ.'YQS')

C *** DPN 18 February 2000 - If evaluating at grid points and using
C *** the YQS array, want to be able to specify grid points.
        IF(TYPE(1:4).EQ.'GRID'.AND.VARIABLE(1:3).EQ.'YQS') THEN
          IF(CBBREV(CO,'GRID_POINTS',6,noco+1,NTCO,N3CO)) THEN
C!!! LKC 26-FEB-2007 I believe this is only implemented for 1 grid_point
C!!!   need to double check            
            CALL PARSIL(CO(N3CO+1),GP_YQSM,GP_YQST,GP_YQS,ERROR,*9999)
            
          ELSE
            GP_YQST=1
            GP_YQS(1)=1
          ENDIF
        ENDIF !not using grid and yqs

        IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          HISTORYFNAME=CO(N3CO+1)(IBEG:IEND)
        ELSE
          CALL STRING_TRIM(FILE00,IBEG,IEND)
          HISTORYFNAME=FILE00(IBEG:IEND)
        ENDIF

        IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
          TSTART=RFROMC(CO(N3CO+1))
        ELSE
          TSTART=-RMAX
        ENDIF

        IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
          TEND=RFROMC(CO(N3CO+1))
        ELSE
          TEND=RMAX
        ENDIF

C*** Fileformat
        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'INTOREGION',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NRM,NRLIST2(0),NRLIST2(1),
     '      ERROR,*9999)
        ELSE
          NRLIST2(0)=NRLIST(0)
          DO nonr=1,NRLIST(0)
            NRLIST2(nonr)=NRLIST(nonr)
          ENDDO
        ENDIF


C***
C*** Start doing the actual exporting of the data here ....
C***        
        
        IF(OPFILE) THEN
          SIGNAL_HEADER(IOFILE1)=' '
          SIGNAL_HOWGENERATED(IOFILE1)=0
C***      Open up OUTPUT the signal file
C MLB time was not set
          TIME=0.0d0

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,SIGNALFILE,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C *** Set-up the electrode data
C *** DPN 19 February 2000 - Using variables/currents as data points
C *** if using grid/yqs
          IF(TYPE(1:4).EQ.'GRID'.AND.VARIABLE(1:3).EQ.'YQS'.AND.
     '      IY_YQST.GT.1) THEN
            !using grid and YQS array and want to use multiple
            !variables/currents
            SIGNAL_NUMREGIONS(IOFILE1)=1
            SIGNAL_REGNAME(1,IOFILE1)='Andre'
            SIGNAL_NUMELEC(1,IOFILE1)=IY_YQST
            SIGNAL_REGTYPE(1,IOFILE1)=0
            SIGNAL_ELEMLOC(IOFILE1)=0
            SIGNAL_NUMXI(1,IOFILE1)=0
            !need to set the number of "data points" and
            !set-up NDDATA and WD,ZD
            CALL ASSERT(NDM.GE.IY_YQST,
     '        'Need to increase NDM',ERROR,*9999)
            NDDATA(0,0)=IY_YQST
            NDDATA(0,1)=IY_YQST
            DO nd=1,IY_YQST
              NDDATA(nd,1)=nd
              DO nj=1,NJT+1
                WD(nj,nd)=1.0d0
              ENDDO
              DO nj=1,NJT
                ZD(nj,nd)=0.0d0
              ENDDO
            ENDDO
          ELSE
            SIGNAL_NUMREGIONS(IOFILE1)=1
            SIGNAL_REGNAME(1,IOFILE1)=' '
            SIGNAL_NUMELEC(1,IOFILE1)=NDT
            SIGNAL_REGTYPE(1,IOFILE1)=0 !irregular rectangular cartesian
            SIGNAL_ELEMLOC(IOFILE1)=1
            SIGNAL_NUMXI(1,IOFILE1)=0
            DO nonrlist=1,NRLIST(0)
              nr=NRLIST(nonrlist)
              IF(NIT(NBJ(1,NEELEM(1,nr))).GT.SIGNAL_NUMXI(1,IOFILE1))
     '          SIGNAL_NUMXI(1,IOFILE1)=NIT(NBJ(1,NEELEM(1,nr)))
            ENDDO
          ENDIF

C LKC 23-FEB-2000 Ensure there are some electrodes to export
C   Need to loop over regions later
          CALL ASSERT(SIGNAL_NUMELEC(1,IOFILE1).GE.1,
     '      'No electodes to export',ERROR,*9999)

C*** Write electrode data
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'WRITE',FILEFORMAT,SIGNALFILE,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        ENDIF !OPFILE


        IF(TYPE(1:5).EQ.'PHI_H') THEN !PHI_H

C LKC 6-SEP-2002 this has already been setup above
C
CC GBS  9-Aug-2000 Phi_h has heart/torso information
C          IF (PHI_H_IS_HEART) THEN
C            PHI_LIMIT=ISIZE_TBH(2) !number of heart nodes
C          ELSE
C            PHI_LIMIT=ISIZE_TBH(1) !number of torso nodes
C         ENDIF

C GBS  9-Aug-2000 Mapping stored in region 0
          IF(EVALUATE_PHI_NEAREST) THEN
            nr=0
          ELSE
            nr=1
          ENDIF

C LKC 9-SEP-2002 should be using ISIZE_PHIH not ISIZE_PHI
C
C         DO nts=1,ISIZE_PHI(2) !number of time steps
          DO nts=1,ISIZE_PHIH(2) !number of time steps
            TIME=CALC_TIME_FROM_SAMPLE(nts)
            IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
              DO nd=1,PHI_LIMIT
                ZD(NJT+1,nd)=PHI_H(nd,nts)
              ENDDO
              IF(OPFILE) THEN
                CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '            FILEFORMAT,SIGNALFILE,'SIGNAL_DATA',ENDFILE,
     '            .TRUE.,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'(/'' Electrode potentials at t='','
     '            //'D12.4)') TIME
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(1X,5D13.5,/:(1X,5D13.5)))')
     '            (ZD(NJT+1,NDDATA(nodata,nr)),nodata=1,NDDATA(0,nr))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF !time
          ENDDO !nts


          
        ELSEIF(TYPE(1:3).EQ.'MFI') THEN

          DO nts=1,ISIZE_MFI(2,nss) !number of time steps
            TIME=CALC_TIME_FROM_SAMPLE(nts)          
            IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN

              nodata=0
              DO nj=1,3 !3d magnetic field
                DO nd=1,ISIZE_MFI(1,nss)
                  nodata=nodata+1
                  ZD(NJT+1,nodata)=MFI(nd,nts,nj,nss)
                ENDDO
              ENDDO
              IF(OPFILE) THEN
                CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '            SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '            FILEFORMAT,SIGNALFILE,'SIGNAL_DATA',ENDFILE,
     '            .TRUE.,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'(/'' Electrode potentials at t='','
     '            //'D12.4)') TIME
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(1X,5D13.5,/:(1X,5D13.5)))')
     '            (ZD(NJT+1,NDDATA(nodata,nr)),nodata=1,NDDATA(0,nr))
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF !time
          ENDDO !nts
         

        ELSE !all others

          
C***      Open up the history file
          IF(VARIABLE(1:2).EQ.'YP') THEN
            NIYLIST(0)=1
            NIYLIST(1)=iy
            NIQLIST(0)=0
            NIQSLIST(0)=0
          ELSE IF(VARIABLE(1:3).EQ.'YQS') THEN
            IF(TYPE(1:4).EQ.'GVAR') THEN
              CALL ASSERT(NDT.LE.99,'>>Increase NIQSLIST dimension',
     '          ERROR,*9999)
              NIQSLIST(0)=NDT
              DO nd=1,NDT
                NIQSLIST(nd)=nd
              ENDDO
              NIYLIST(0)=0
              NIQLIST(0)=0
            ELSE
C *** DPN 18 February 2000 - Adding ability to have multiple variables
C *** from a grid point
c            NIQSLIST(0)=1
c            NIQSLIST(1)=iy
              NIQSLIST(0)=IY_YQST
              DO iy=1,IY_YQST
                NIQSLIST(iy)=IY_YQS(iy)
              ENDDO
              NIYLIST(0)=0
              NIQLIST(0)=0
            ENDIF
          ELSE IF(VARIABLE(1:2).EQ.'YQ') THEN
            NIQLIST(0)=1
            NIQLIST(1)=iy
            NIYLIST(0)=0
            NIQSLIST(0)=0
          ELSE
            ERROR='>>Invalid variable'
            GOTO 9999
          ENDIF

          na=1
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA1,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,HISTORYFNAME,
     '      'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)

          IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
            ENDFILE=.FALSE.
            NUMTIMEDATA=0
            DO WHILE(.NOT.ENDFILE)
              CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),
     '          YQS,'READ',
     '          FILEFORMAT,HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) THEN
                IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN

                  IF(TYPE(1:4).NE.'GRID') THEN
                    DO nonrlist=1,NRLIST(0)
                      nr=NRLIST(nonrlist)
                      CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '                  NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     '                  NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                    ENDDO !nonrlist

                    DO nodata=1,NDDATA(0,1)
                      nd=NDDATA(nodata,1)
                      ne=LD(nd)
                      IF(ne.NE.0) THEN
                        nb=NBH(NH_LOC(1,nx),1,ne)
                        CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),
     '                    NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '                    NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '                    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '                    ZE,ZP,ERROR,*9999)
                        ZD(NJT+1,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,1,XID(1,nd),ZE)
                      ENDIF
                    ENDDO !nolist (nd)
                  ELSE !Grid points
C *** DPN 18 February 2000 - this is bad!! Need to check is used and
C *** actually use the variable specified
c                    DO nodata=1,NDDATA(0,1)
c                      nd=NDDATA(nodata,1)
c                      nq=LIST(nd)
c                      ZD(NJT+1,nd)=YQ(nq,niqV,na,nx)
c                    ENDDO
                    IF(VARIABLE(1:3).EQ.'YQS') THEN
                      IF(IY_YQST.GT.1) THEN
                        CALL ASSERT(NDM.GE.IY_YQST,
     '                    'Need to increase NDM',ERROR,*9999)
                        !need to set the number of "data points" and
                        !set-up NDDATA and WD,ZD
                        NDDATA(0,0)=IY_YQST
                        NDDATA(0,1)=IY_YQST
                        DO nd=1,IY_YQST
                          iy=IY_YQS(nd)
                          nq=GP_YQS(1)
                          ZD(NJT+1,nd)=YQS(iy,nq)
                          NDDATA(nd,1)=nd
                          DO nj=1,NJT+1
                            WD(nj,nd)=1.0d0
                          ENDDO
                        ENDDO
                        SIGNAL_NUMREGIONS(IOFILE1)=1
                        SIGNAL_REGNAME(1,IOFILE1)='Andre'
                        SIGNAL_NUMELEC(1,IOFILE1)=IY_YQST
                        SIGNAL_REGTYPE(1,IOFILE1)=0
                        SIGNAL_ELEMLOC(IOFILE1)=0
                        SIGNAL_NUMXI(1,IOFILE1)=0
                      ELSE
                        DO nodata=1,NDDATA(0,1)
                          nd=NDDATA(nodata,1)
                          nq=LIST(nd)
                          iy=IY_YQS(1)
                          ZD(NJT+1,nd)=YQS(iy,nq)
                        ENDDO
                      ENDIF
                    ELSE
                      DO nodata=1,NDDATA(0,1)
                        nd=NDDATA(nodata,1)
                        nq=LIST(nd)
                        ZD(NJT+1,nd)=YQ(nq,iy,na,nx)
                      ENDDO !YQS
                    ENDIF
                  ENDIF

                  IF(OPFILE) THEN

                    CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '                SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '                FILEFORMAT,SIGNALFILE,'SIGNAL_DATA',ENDFILE,
     '                .TRUE.,ERROR,*9999)

                  ELSE
                    WRITE(OP_STRING,'(/'' Electrode potentials at t='','
     '                //'D12.4)') TIME
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'(1X,5D13.5,/:(1X,5D13.5)))')
     '                (ZD(NJT+1,NDDATA(nodata,1)),nodata=1,NDDATA(0,1))
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO

          ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
            DO notime=1,NUMTIMEDATA1

C***        Read the history data

              CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,
     '          NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,
     '          NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,
     '          YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)

              IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
                IF(TYPE(1:4).EQ.'GVAR') THEN
                  IF(VARIABLE(1:3).EQ.'YQS') THEN
                    nq=LIST(1)
                    DO nodata=1,NDDATA(0,1)
                      nd=NDDATA(nodata,1)
                      ZD(NJT+1,nd)=YQS(nd,nq)
                    ENDDO
                  ELSE
                    nq=LIST(1)
                    DO nodata=1,NDDATA(0,1)
                      nd=NDDATA(nodata,1)
                      ZD(NJT+1,nd)=YQ(nq,iy,na,nx)
                    ENDDO
                  ENDIF
                ELSE IF(TYPE(1:4).NE.'GRID') THEN
                  DO nonrlist=1,NRLIST(0)
                    nr=NRLIST(nonrlist)
                    CALL YPZP(iy,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '                NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,
     '                NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                  ENDDO !nonrlist

                  DO nodata=1,NDDATA(0,1)
                    nd=NDDATA(nodata,1)
                    ne=LD(nd)
                    IF(ne.NE.0) THEN
                      nb=NBH(NH_LOC(1,nx),1,ne)
                      CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
     '                  NPF(1,1),NPNE(1,1,ne),NRE(ne),NVHE(1,1,1,ne),
     '                  NW(ne,1,nx),nx,CURVCORRECT(1,1,1,ne),SE(1,1,ne),
     '                  ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
                      ZD(NJT+1,nd)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,XID(1,nd),ZE)
                    ENDIF
                  ENDDO !nolist (nd)
                ELSE !Grid points
                  IF(VARIABLE(1:3).EQ.'YQS') THEN
C *** DPN 18 February 2000 - if more than one variable is specified
C *** each variable is a data point, otherwise simply use the specified
C *** variable for all grid points.xxxxxxx
                    IF(IY_YQST.GT.1) THEN
                      DO nd=1,IY_YQST
                        iy=IY_YQS(nd)
                        nq=GP_YQS(1)
                        ZD(NJT+1,nd)=YQS(iy,nq)
                      ENDDO
                    ELSE
                      DO nodata=1,NDDATA(0,1)
                        nd=NDDATA(nodata,1)
                        nq=LIST(nd)
                        iy=IY_YQS(1)
                        ZD(NJT+1,nd)=YQS(iy,nq)
                      ENDDO
                    ENDIF
                  ELSE
                    DO nodata=1,NDDATA(0,1)
                      nd=NDDATA(nodata,1)
                      nq=LIST(nd)
                      ZD(NJT+1,nd)=YQ(nq,iy,na,nx)
                    ENDDO !YQS
                  ENDIF
                ENDIF !GRID

                IF(OPFILE) THEN

                  CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,
     '              SIGNALMAX,SIGNALMIN,TIME,WD,XID,ZD,'WRITE',
     '              FILEFORMAT,SIGNALFILE,'SIGNAL_DATA',ENDFILE,
     '              .TRUE.,ERROR,*9999)

                ELSE
                  WRITE(OP_STRING,'(/'' Electrode potentials at t='','
     '              //'D12.4)') TIME
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'(1X,5D13.5,/:(1X,5D13.5)))')
     '              (ZD(NJT+1,NDDATA(nodata,1)),nodata=1,NDDATA(0,1))
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF

            ENDDO !notime

          ENDIF

          ENDFILE2=.TRUE.
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA1,
     '      nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),
     '      YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      HISTORYFNAME,
     '      ' ',ENDFILE2,YPDATA,YQDATA,YQSDATA,ENDFILE,ERROR,*9999)

        ENDIF !TYPE

        IF(OPFILE) THEN
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGNALFILE,
     '      ' ',ENDFILE,.TRUE.,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVELEC')
      RETURN
 9999 CALL ERRORS('EVELEC',ERROR)
      CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,
     '  NIYLIST,NPNY(0,1,0,nx),
     '  NQNY(1,1,0,nx),NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),
     '  NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,
     '  'CLOSE',FILEFORMAT,HISTORYFNAME,' ',ENDFILE,.TRUE.,YPDATA,
     '  YQDATA,YQSDATA,ERROR1,*9998)
 9998 IF(OPFILE) THEN
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGNALFILE,
     '    ' ',ENDFILE,.TRUE.,ERROR1,*9997)
      ENDIF
 9997 CALL EXITS('EVELEC')
      RETURN 1
      END


