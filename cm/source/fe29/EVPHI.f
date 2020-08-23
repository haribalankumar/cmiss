      SUBROUTINE EVPHI(ISIZE_PHI,ISIZE_PHIH,LD,NBJ,NDDATA,NHP,NHQ,
     '  NKH,NPLIST3,
     '  NPLIST4,NPLIST5,NPNY,NQNY,NRLIST,NVHP,NXLIST,NYNP,NYNR,NYQNR,
     '  PHI,PHI_H,PHI_H_EXACT,SIGMA_PHI,U_PHI,VT_PHI,
     '  WD,WK1_INV,WK3_INV,XID,XP,YP,YQ,YQS,ZD,STRING,ERROR,*)

C#### Subroutine: EVPHI
C###  Description:
C###    EVPHI evaluates the PHI array and properties of it.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER ISIZE_PHI(2),ISIZE_PHIH(2),LD(NDM),NBJ(NJM,NEM),
     '  NDDATA(0:NDM,0:NRM),
     '  NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPLIST3(0:NPM),NPLIST4(0:NPM),NPLIST5(0:NPM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NQNY(2,NYQM,0:NRCM,NXM),NRLIST(0:NRM),
     '  NVHP(NHM,NPM,NCM,0:NRM),
     '  NXLIST(0:NXM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 PHI(NY_TRANSFER_M,NTSM),PHI_H(NY_TRANSFER_M,NTSM),
     '  PHI_H_EXACT(NY_TRANSFER_M,NTSM),
     '  SIGMA_PHI(NY_TRANSFER_M),U_PHI(NY_TRANSFER_M,NY_TRANSFER_M),
     '  VT_PHI(NTSM,NTSM),WD(NJM,NDM),
     '  WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M),
     '  WK3_INV(NY_TRANSFER_M,NY_TRANSFER_M),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER CALC_SAMPLE_FROM_TIME,ERR,IBEG,IBEG1,icol,iDUMMY,IEND,
     '  IEND1,IFAIL,irow,isize,
     '  LWORK,MAXROWS,na,N3CO,nd,ndd,ndnp,nh,nhx,NIQLIST(0:1),
     '  NIQSLIST(0:1),NIYLIST(0:16),nj,nk,nnp,noelec,
     '  nolist,notime,np,nr,nr1,nts,nonr,NUMTIMEDATA,NUMTIMEDATA1,nv,nx,
     '  nxc,ny,TRSF_NRLIST2(0:9)
      REAL*8 CONDITION_NUMBER,DELTATIME,DIST,FREQUENCY,MINDIST,PREVTIME,
     '  RFROMC,SIGNALMAX(9),SIGNALMIN(9),SUMDIST,TEND,TIME,
     '  TSTART,YPMAX(16),YPMIN(16)
      CHARACTER DATAFNAME*(MXCH),ERROR_DUMMY*(MXCH),
     '  FILE*(MXCH),FILEFORMAT*6,HISTORYFNAME*(MXCH),SIGNALFNAME*(MXCH)
      LOGICAL ALL_ELECTRODES,ALL_REGIONS,AT_NODES,CBBREV,ENDFILE,
     '  ENDFILE2,HISTORY,INLIST,NEAREST_MAP,OPDATA,OPFILE,PHIH,
     '  PHIHEXACT,SIGNAL,SVD,SVD_CUTOFF,WARNING,YPDATA,YQDATA,YQSDATA

      CALL ENTERS('EVPHI',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM evaluate phi<;FILENAME>
C###  Parameter:        <history FILENAME[$current]>
C###    Specify the history filename (binhis or iphist) to
C###    evaluate the PHI matrix from.
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
C###    Change the region numbers to read into.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Creates the PHI(nytr,nts) matrix from a history file.  The
C###    columns of PHI correspond to each time sample and the rows
C###    correspond to torso electrode/ny.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<history FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<tstart (#/beginning)[beginning]>'
        OP_STRING(4)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<intoregion (#s/same)[same]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate phi<;FILENAME>
C###  Parameter:        <signal FILENAME[$current]>
C###    Specify the input signal files (binsig or ipsign)
C###     from which the PHI/PHI_H_EXACT matrix is created.
C###  Parameter:        <(phi/phi_h_exact)[phi]>
C###    Specify whether phi (torso signals vs time) or phi_h_exact
C###    (exact heart signals vs time) to evaluate.
C###  Parameter:        <tstart (#/beginning)[beginning]>
C###    Specify start time
C###  Parameter:        <tend (#/end)[end]>
C###    Specify end time
C###  Parameter:        <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:        <at (nodes/electrodes)[nodes]>
C###    By default, signals are at nodes, but may be specified to be
C###    at electrode (data) locations.
C###  Parameter:        <nearest_map>
C###    Signals are not at nodes, and require mapping of nearest nodes
C###    to the electrodes.
C###  Parameter:        <all_electrodes>
C###    If calculating nearest mapping, use all electrodes (default is
C###    to use only non-rejected electrodes).
C###  Description:
C###    Creates the PHI(nytr,nts) matrix or PHI_H_EXACT(nytr,nts)
C###    matrix from a signal file.  The columns of PHI/PHI_H_EXACT
C###    correspond to each time sample and the rows correspond to
C###    torso/heart electrode/ny.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<signal FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(phi/phi_h/phi_h_exact)[phi]>'
        OP_STRING(4)=BLANK(1:15)//'<tstart (#/beginning)[beginning]>'
        OP_STRING(5)=BLANK(1:15)//'<tend (#/end)[end]>'
        OP_STRING(6)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(9)=BLANK(1:15)//'<at (nodes/electrodes)[nodes]>'
        OP_STRING(10)=BLANK(1:15)//'<nearest_map>'
        OP_STRING(11)=BLANK(1:15)//'<all_electrodes>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate phi<;FILENAME>
C###  Parameter:        <svd>
C###    Perform a Singular Value Decomposition on the PHI matrix.
C###  Parameter:        <datafile FILENAME>
C###    Outputs the singular values and normalised values
C###    to a text data file (*.dat).
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Computes the SVD of PHI

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<svd>'
        OP_STRING(3)=BLANK(1:15)//'<datafile FILENAME>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate phi<;FILENAME>
C###  Parameter:        <svd_cutoff SVD_CUTOFF_RATIO#>
C###    Specify a value below which the singular values
C###    are ignored.
C###  Parameter:        <datafile FILENAME>
C###    Outputs the singular values and normalised values
C###    to a text data file (*.dat).
C###  Parameter:        <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Defines the cutoff ratio (sigma/sigma_max) of the singular
C###    values of PHI. Sets the variable svd_cutoff and
C###    svd_cutoff_ratio.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<svd_cutoff SVD_CUTOFF_RATIO#>'
        OP_STRING(3)=BLANK(1:15)//'<datafile FILENAME>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVPHI',ERROR,*9999)
      ELSE
        SIGNAL=.FALSE.
        HISTORY=.FALSE.

        CALL ASSERT(USE_TIME.EQ.1,
     '    '>>Set USE_TIME to 1 in parameter file',ERROR,*9999)
        CALL ASSERT(USE_TRANSFER.EQ.1,
     '    '>>Set USE_TRANFER to 1 in parameter file',ERROR,*9999)

        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opphi','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        OPDATA=.FALSE.
        IF(CBBREV(CO,'DATAFILE',3,noco+1,NTCO,N3CO)) THEN
          OPDATA=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          DATAFNAME=CO(N3CO+1)(IBEG:IEND)
        ENDIF

        IF(CBBREV(CO,'HISTORY',3,noco+1,NTCO,N3CO)) THEN
          HISTORY=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          HISTORYFNAME=CO(N3CO+1)(IBEG:IEND)
        ELSEIF(CBBREV(CO,'SIGNAL',2,noco+1,NTCO,N3CO)) THEN
          SIGNAL=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          SIGNALFNAME=CO(N3CO+1)(IBEG:IEND)
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

        IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
          FILEFORMAT='BINARY'
        ELSE
          FILEFORMAT='ASCII'
        ENDIF

        IF(CBBREV(CO,'PHI_H_EXACT',6,noco+1,NTCO,N3CO)) THEN
          PHIHEXACT=.TRUE.
        ELSE
          PHIHEXACT=.FALSE.
        ENDIF

        IF(CBBREV(CO,'PHI_H',4,noco+1,NTCO,N3CO)) THEN
          PHIH=.TRUE.
        ELSE
          PHIH=.FALSE.
        ENDIF

        IF(CBBREV(CO,'AT',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'ELECTRODES',1,n3co+1,n3co+2,N3CO)) THEN
            AT_NODES=.FALSE.
          ELSE
            AT_NODES=.TRUE.
          ENDIF
        ELSE
          AT_NODES=.TRUE.
        ENDIF

        IF(CBBREV(CO,'NEAREST_MAP',4,noco+1,NTCO,N3CO)) THEN
          CALL ASSERT(KTYP95.LE.2,
     '      '>>Cannot map nodes for transmembrane to epi matrix',
     '      ERROR,*9999)
          NEAREST_MAP=.TRUE.
        ELSE
          NEAREST_MAP=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ALL_ELECTRODES',3,noco+1,NTCO,N3CO)) THEN
          ALL_ELECTRODES=.TRUE.
        ELSE
          ALL_ELECTRODES=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SVD',3,noco+1,NTCO,N3CO)) THEN
          SVD=.TRUE.
        ELSE
          SVD=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SVD_CUTOFF',4,noco+1,NTCO,N3CO)) THEN
          SVD_CUTOFF=.TRUE.
          SVD_CUTOFF_RATIO=RFROMC(CO(N3CO+1))
        ELSE
          SVD_CUTOFF=.FALSE.
          SVD_CUTOFF_RATIO=-RMAX
        ENDIF
C***    Open up the history file

        IF(HISTORY) THEN
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          CALL ASSERT(CALL_TRANSFER,
     '      '>>Define the transfer first',ERROR,*9999)

          nr1=TRSF_NR_OUTER !outer region number of body
          CALL ASSERT(NRLIST(1).EQ.nr1,
     '      '>>PHI region does not match transfer region',ERROR,*9999)
          NIYLIST(0)=1
          NIYLIST(1)=1
          NIQLIST(0)=0
          NIQSLIST(0)=0

          IF(CBBREV(CO,'INTOREGION',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),9,TRSF_NRLIST2(0),TRSF_NRLIST2(1),
     '        ERROR,*9999)
          ELSE
            TRSF_NRLIST2(0)=NRLIST(0)
            DO nonr=1,NRLIST(0)
              TRSF_NRLIST2(nonr)=NRLIST(nonr)
            ENDDO
          ENDIF

C         Check the frequency of the signal.
          PREVTIME=-RMAX
          DELTATIME=0.0d0
          WARNING=.FALSE.

          na=1
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),
     '      NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),
     '      YQS,
     '      'READ',FILEFORMAT,HISTORYFNAME,'OPEN',ENDFILE,.TRUE.,
     '      YPDATA,YQDATA,YQSDATA,ERROR,*9999)

          ERR=0
          IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
            ENDFILE=.FALSE.
            NUMTIMEDATA=0
            DO WHILE(.NOT.ENDFILE)
              CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,TRSF_NRLIST2,
     '          NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,HISTORYFNAME,'TIME_DATA',ENDFILE,.TRUE.,
     '          YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) THEN
                IF(PREVTIME.GT.-RMAX) THEN
                  IF(DELTATIME.EQ.0.0d0) THEN
                    DELTATIME=TIME-PREVTIME
                  ELSE
C LKC 20-FEB-2003 Need to change to LOOSE_TOL as errors
C                   were creaping in with ASCII files.
C                    
C                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.ZERO_TOL)
                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.LOOSE_TOL)
     '                WARNING=.TRUE.
                  ENDIF
                ENDIF
                PREVTIME=TIME
                IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
C***              Find nys of outer surface and put the corresponding
C***              time values into PHI
                  nts=CALC_SAMPLE_FROM_TIME(TIME,ERR,ERROR)
                  IF(ERR.NE.0) GOTO 9999
                  irow=0
                  DO nolist=1,NPLIST4(0) !list of body nodes
                    np=NPLIST4(nolist)
                    DO nhx=1,NHP(np,nr1,nx)
                      nh=NH_LOC(nhx,nx)
                      DO nv=1,NVHP(nh,np,1,nr1)
                        DO nk=1,
     '                    MAX(NKH(nh,np,1,nr1)-KTYP93(1,nr1),1)
                          irow=irow+1
                          CALL ASSERT(irow.LE.NY_TRANSFER_M,
     '                      '>>Increase ny_transfer_m',ERROR,*9999)
                          ny=NYNP(nk,nv,nh,np,0,1,nr1)
C                           global ny # for this variable
                          PHI(irow,nts)=YP(ny,1,nx)
                        ENDDO !nk
                      ENDDO !nv
                    ENDDO !nh
                  ENDDO !nolist (np)
                  MAXROWS=irow
                ENDIF !time
              ENDIF !end of file
            ENDDO !end while
            NTST=nts

          ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
            DO notime=1,NUMTIMEDATA1

C***          Read the history data

              CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,
     '          NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,TRSF_NRLIST2,
     '          NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,HISTORYFNAME,'TIME_DATA',ENDFILE,
     '          .TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
              IF(TIME.GE.TSTART.AND.TIME.LE.TEND) THEN
                IF(PREVTIME.GT.-RMAX) THEN
                  IF(DELTATIME.EQ.0.0d0) THEN
                    DELTATIME=TIME-PREVTIME
                  ELSE
C LKC 20-FEB-2003 Need to change to LOOSE_TOL as errors
C                   were creaping in with ASCII files.
C                    
C                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.ZERO_TOL)
                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.LOOSE_TOL)
     '                WARNING=.TRUE.
                  ENDIF
                ENDIF
                PREVTIME=TIME
C***            Find nys of outer surface and put the corresponding time
C***            values into PHI
                nts=CALC_SAMPLE_FROM_TIME(TIME,ERR,ERROR)
                IF(ERR.NE.0) GOTO 9999
                irow=0
                DO nolist=1,NPLIST4(0) !list of body nodes
                  np=NPLIST4(nolist)
                  DO nhx=1,NHP(np,nr1,nx)
                    nh=NH_LOC(nhx,nx)
                    DO nv=1,NVHP(nh,np,1,nr1)
                      DO nk=1,
     '                  MAX(NKH(nh,np,1,nr1)-KTYP93(1,nr1),1)
                        irow=irow+1
                        CALL ASSERT(irow.LE.NY_TRANSFER_M,
     '                    '>>Increase ny_transfer_m',ERROR,*9999)
                        ny=NYNP(nk,nv,nh,np,0,1,nr1)
C                           global ny # for this variable
                        PHI(irow,nts)=YP(ny,1,nx)
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !nh
                ENDDO !nolist (np)
                MAXROWS=irow
              ENDIF !time
            ENDDO !notime
            NTST=nts

          ENDIF

          ENDFILE2=.TRUE.
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      HISTORYFNAME,' ',ENDFILE2,YPDATA,YQDATA,
     '      YQSDATA,ENDFILE,ERROR,*9999)

          IF(WARNING) THEN
            WRITE(OP_STRING,'('' >>WARNING: History time steps are not '
     '        //'consistent'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(DABS(DELTATIME).GT.ZERO_TOL) THEN
            FREQUENCY=1.0d0/DELTATIME
          ELSE
            ERROR='>>Zero delta time'
            GOTO 9999
          ENDIF
C KAT Transfer frequency is in kHz, but FREQUENCY units seem to be undefined.
CC cpb 16/8/02 Transfer frequency is now in kHz, FREQUENCY is in Hz.
          IF(DABS(FREQUENCY-TRSF_FREQUENCY).GT.ZERO_TOL) THEN
            WRITE(OP_STRING,'('' >>WARNING: History frequency of '','
     '        //'D12.4,'' and transfer frequency of '',D12.4,'
     '        //''' kHz are not consistent'')') FREQUENCY,TRSF_FREQUENCY
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          ISIZE_PHI(1)=MAXROWS
          ISIZE_PHI(2)=NTST
          IF(OPFILE) THEN
            WRITE(OP_STRING,'(/'' PHI'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO irow=1,MAXROWS
              WRITE(OP_STRING,'('' Row '',I5,'' :'',5D12.4,'
     '          //'/:(12X,5D12.4))') irow,
     '          (PHI(irow,nts),nts=1,NTST)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ELSE IF(SIGNAL) THEN
          IF(AT_NODES.AND..NOT.NEAREST_MAP) THEN
            WRITE(OP_STRING,
     '        '(/'' >>WARNING: Signals are assumed to be at nodes'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' They are also assumed to be in the same order'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,
     '        '('' as the nodal list used in the detrans/deinve'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

C         Check the frequency of the signal.
          PREVTIME=-RMAX
          DELTATIME=0.0d0
          WARNING=.FALSE.

C***      Open up the signal file

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '      'OPEN',ENDFILE,.TRUE.,ERROR,*9999)

C***      Read electrode geometry etc.

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '      'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

          CALL ASSERT(SIGNAL_NUMREGIONS(IOFILE1).EQ.1,
     '      '>>More than one region in file - where are torso sigs?',
     '      ERROR,*9999)

          nr=1

          CALL ASSERT(SIGNAL_NUMELEC(1,IOFILE1).LE.NY_TRANSFER_M,
     '      '>>Increase NY_TRANSFER_M',ERROR,*9999)

C***      Find nearest node to each electrode and store into NPLIST4
          IF(NEAREST_MAP) THEN
            SUMDIST=0.0d0
            nd=0
            DO ndd=1,SIGNAL_NUMELEC(nr,IOFILE1)
C             Only create transfer matrix for non-rejected electrodes
C               unless ALL_ELECTRODES specified
C             Store electrode numbers in NDDATA(nd,0)
              IF(WD(NJT+1,ndd).GT.ZERO_TOL.OR.ALL_ELECTRODES) THEN
                MINDIST=1.0d10
                ndnp=0
                DO nnp=1,NPLIST5(0)
                  np=NPLIST5(nnp)
                  DIST=0.0d0
                  DO nj=1,NJT
                    DIST=DIST+(ZD(nj,ndd)-XP(1,1,nj,np))**2
                  ENDDO
                  IF(DIST.LT.MINDIST) THEN
                    MINDIST=DIST
                    ndnp=np
                  ENDIF
                ENDDO
                IF (INLIST(ndnp,NPLIST4(1),nd,iDUMMY)) THEN
                  WRITE(OP_STRING,'('' Warning: Node '',I7,
     '              '' already used'')') ndnp
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
                nd=nd+1
                NPLIST4(nd)=ndnp
                NDDATA(nd,0)=ndd
                SUMDIST=SUMDIST+MINDIST
                IF(OPFILE) THEN
                  WRITE(OP_STRING,'('' Electrode '',I5,'' to node '',I5,
     '              '': Dist = '',F15.5)') nd,ndnp,SQRT(MINDIST)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO
            NPLIST4(0)=nd
            NDDATA(0,0)=nd
            EVALUATE_PHI_NEAREST=.TRUE.
            IF(OPFILE) THEN
              WRITE(OP_STRING,'('' RMS distance: '',F15.5)')
     '          SQRT(SUMDIST/NPLIST4(0))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE IF (AT_NODES) THEN
            DO nd=1,NDDATA(0,nr)
              NDDATA(nd,0)=NDDATA(nd,nr)
            ENDDO
            NDDATA(0,0)=NDDATA(0,nr)
          ELSE    ! PHI is stored at electrodes
            nd=0
            DO ndd=1,SIGNAL_NUMELEC(nr,IOFILE1)
C             Only create transfer matrix for non-rejected electrodes
C               unless ALL_ELECTRODES specified
C             Store electrode numbers in NDDATA(nd,0)
              IF(WD(NJT+1,ndd).GT.ZERO_TOL.OR.ALL_ELECTRODES) THEN
                nd=nd+1
                NPLIST4(nd)=ndd
                NDDATA(nd,0)=ndd
              ENDIF
            ENDDO
            NPLIST4(0)=nd
            NDDATA(0,0)=nd
          ENDIF !nearest_map

C***      IF ASCII need to calculate NUMTIMEDATA
          IF(FILEFORMAT.EQ.'ASCII') THEN
            NUMTIMEDATA=0
            ENDFILE=.FALSE.
            DO WHILE(.NOT.ENDFILE)
              CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '          SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '          'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
              IF(.NOT.ENDFILE) THEN
                NUMTIMEDATA=NUMTIMEDATA+1
              ENDIF
            ENDDO
          ENDIF

          CALL ASSERT(NUMTIMEDATA.LE.NTSM,'>>Increase NTSM',ERROR,*9999)

C***      Reset the input file
          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '      'RESET',ENDFILE,.TRUE.,ERROR,*9999)

C***      Read in the electrode data
          DO notime=1,NUMTIMEDATA !Loop over the times
C***        Read the signal data in
            CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '        SIGNALMIN,TIME,WD,XID,ZD,'READ',FILEFORMAT,SIGNALFNAME,
     '        'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

            IF(PREVTIME.GT.-RMAX) THEN
              IF(DELTATIME.EQ.0.0d0) THEN
                DELTATIME=TIME-PREVTIME
              ELSE
C LKC 20-FEB-2003 Need to change to LOOSE_TOL as errors
C                   were creaping in with ASCII files.
C                    
C                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.ZERO_TOL)
                    IF(DABS((TIME-PREVTIME)-DELTATIME).GT.LOOSE_TOL)
     '            WARNING=.TRUE.
              ENDIF
            ENDIF
            PREVTIME=TIME

            IF(PHIHEXACT) THEN
C***          Put signals into the PHI_H_EXACT array.  The signals are
C***          assumed to be at nodal locations and in the  '
C***          same order as the list of nodes in NPLIST3
              DO noelec=1,SIGNAL_NUMELEC(1,IOFILE1)
                nd=NDDATA(noelec,1)
                PHI_H_EXACT(nd,notime)=ZD(NJT+1,nd)
              ENDDO !noeles

C***          set ISIZE_PHIH for the heart phi matrix
              ISIZE_PHIH(1)=SIGNAL_NUMELEC(1,IOFILE1)
              ISIZE_PHIH(2)=NUMTIMEDATA

            ELSEIF(PHIH) THEN
              DO noelec=1,SIGNAL_NUMELEC(1,IOFILE1)
                nd=NDDATA(noelec,1)
                PHI_H(nd,notime)=ZD(NJT+1,nd)
              ENDDO !noeles
C***          set ISIZE_PHIH for the heart phi matrix
              ISIZE_PHIH(1)=SIGNAL_NUMELEC(1,IOFILE1)
              ISIZE_PHIH(2)=NUMTIMEDATA

            ELSE
C***          Put signals into the PHI array.  The signals are assumed
C***          to be at nodal locations and in the same order as the
C***          list of nodes in NPLIST4
              DO noelec=1,NDDATA(0,0)
                nd=NDDATA(noelec,0)
                PHI(noelec,notime)=ZD(NJT+1,nd)
              ENDDO !noeles

C***          set ISIZE_PHI for the torso surface phi matrix
C              !same as SIGNAL_NUMELEC(1,IOFILE1) ?
              ISIZE_PHI(1)=NDDATA(0,0)
              ISIZE_PHI(2)=NUMTIMEDATA

            ENDIF !PHI types
          ENDDO !nts
          NTST=NUMTIMEDATA

          IF(PHIHEXACT) THEN
C*** Check num electrodes = Num first layer nodes (heart surface)
            EVALUATE_PHI_H_EXACT=.TRUE.

            CALL ASSERT(SIGNAL_NUMELEC(1,IOFILE1).EQ.NPLIST3(0),
     '        '>> Num electrodes inconsistent with transfer setup',
     '        ERROR,*9999)

          ENDIF

          CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '       SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGNALFNAME,
     '       ' ',ENDFILE,.TRUE.,ERROR,*9999)

          IF(WARNING) THEN
            WRITE(OP_STRING,'('' >>WARNING: Signal time steps are not '
     '        //'consistent'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
          IF(DABS(DELTATIME).GT.ZERO_TOL) THEN
            FREQUENCY=1.0d0/DELTATIME
          ELSE
            ERROR='>>Zero delta time'
            GOTO 9999
          ENDIF
          IF(DABS(FREQUENCY-TRSF_FREQUENCY).GT.ZERO_TOL) THEN
            WRITE(OP_STRING,'('' >>WARNING: Signal frequency of '','
     '        //'D12.4,'' kHz and transfer frequency of '',D12.4,'
     '        //''' kHz are not consistent'')') FREQUENCY,TRSF_FREQUENCY
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(SVD) THEN

          CALL ASSERT(ISIZE_PHI(1).LE.NY_TRANSFER_M,
     '      '>> Increase NY_TRANSFER_M',ERROR,*9999)
          CALL ASSERT(ISIZE_PHI(2).LE.NY_TRANSFER_M,
     '      '>> Increase NY_TRANSFER_M',ERROR,*9999)

          DO irow=1,ISIZE_PHI(1)
            DO icol=1,ISIZE_PHI(2)
              WK1_INV(irow,icol)=PHI(irow,icol)
            ENDDO
          ENDDO
          IFAIL=1

          LWORK=MAX(3*MIN(ISIZE_PHI(1),NTST)+MAX(ISIZE_PHI(1),NTST),
     '      5*MIN(ISIZE_PHI(1),NTST)-4)
          CALL ASSERT(LWORK.LE.NY_TRANSFER_M*NY_TRANSFER_M,
     '      '>>Work array too small for SVD of PHI',ERROR,*9999)
          CALL DGESVD('A','N',ISIZE_PHI(1),NTST,WK1_INV,
     '      NY_TRANSFER_M,SIGMA_PHI,U_PHI,NY_TRANSFER_M,VT_PHI,NTSM,
     '      WK3_INV,NY_TRANSFER_M*NY_TRANSFER_M,IFAIL)

          IF(IFAIL.NE.0)THEN
            WRITE(OP_STRING,*)' IFAIL=',IFAIL,' in SVD of PHI'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ERROR='>>IFAIL<>0 in DGESVD'
            GOTO 9999
          ENDIF
          ISIZE=MIN(ISIZE_PHI(1),ISIZE_PHI(2))
          IF(SIGMA_PHI(ISIZE).GT.RDELTA)THEN
            CONDITION_NUMBER=SIGMA_PHI(1)/SIGMA_PHI(ISIZE)
          ELSE
            CONDITION_NUMBER=RMAX !default value
          ENDIF
          IF(OPFILE)THEN
            WRITE(OP_STRING,'(/'' Condition number of PHI :'',D12.4)')
     '        CONDITION_NUMBER
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Singular values of PHI'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5D12.4,'
     '          //'/:(13X,5D12.4))')
     '          (SIGMA_PHI(icol),icol=1,ISIZE)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Normalised Singular values of PHI'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)


            IF(DABS(SIGMA_PHI(1)).GT.ZERO_TOL) THEN
              IF(INT(ISIZE/5).LT.25) THEN
                WRITE(OP_STRING,'(13X,5D12.4,'
     '            //'/:(13X,5D12.4))')
     '            (SIGMA_PHI(icol)/SIGMA_PHI(1),icol=1,ISIZE)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE
                WRITE(OP_STRING,'('' First 125 of '',I8,
     '            '' Singular values of PHI'')') ISIZE
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(13X,5D12.4,'
     '            //'/:(13X,5D12.4))')
     '            (SIGMA_PHI(icol)/SIGMA_PHI(1),icol=1,100)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF ! 25
            ELSE
              WRITE(OP_STRING,'('' WARNING: SIGNMA_PHI(1) is 0'')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF ! Divide by 0
          ENDIF

          WRITE(OP_STRING,'(/'' Condition number of PHI :'',D12.4)')
     '      CONDITION_NUMBER
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          IF(INT(ISIZE/5).LT.25) THEN
            WRITE(OP_STRING,'('' Singular values of PHI'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5E14.4E3,'
     '        //'/:(13X,5E14.4E3))')
     '        (SIGMA_PHI(icol),icol=1,ISIZE)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' First 125 of '',I8,
     '        '' Singular values of PHI'')') ISIZE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5E14.4E3,'
     '        //'/:(13X,5E14.4E3))')
     '        (SIGMA_PHI(icol),icol=1,100)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          IF(INT(ISIZE/5).LT.25) THEN
            WRITE(OP_STRING,'('' Normalised Singular values of PHI'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5E14.4E3,'
     '        //'/:(13X,5E14.4E3))')
     '        (SIGMA_PHI(icol)/SIGMA_PHI(1),icol=1,ISIZE)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' First 125 of '',I8,
     '        '' Normalised Singular values of PHI'')') ISIZE
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(13X,5E14.4E3,'
     '        //'/:(13X,5E14.4E3))')
     '        (SIGMA_PHI(icol)/SIGMA_PHI(1),icol=1,100)
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF

          EVALUATE_PHI_SVD=.TRUE.
        ELSEIF(SVD_CUTOFF) THEN
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF


C*** Output the singular values and the normalised singular values
C*** to a data file for use by gnuplot

        IF( OPDATA ) THEN
          ISIZE=MIN(ISIZE_PHI(1),ISIZE_PHI(2))

C*** Open file
          CALL STRING_TRIM(DATAFNAME,IBEG1,IEND1)
          CALL OPENF(IOFILE1,'DISK',DATAFNAME(IBEG1:IEND1)//'.dat',
     '      'NEW','SEQUEN','FORMATTED',132,ERROR,*9999)

C*** Output information
          WRITE(OP_STRING(1),'(''# Singular values of matrix'')')
          WRITE(OP_STRING(2),'(''#          '')')
          WRITE(OP_STRING(3),'(''# Condition number: '',E12.5)')
     '      SIGMA_PHI(1)/SIGMA_PHI(ISIZE)
          CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)


C*** Calcuate the effective rank
          IF(SVD_CUTOFF) THEN

            irow=0
            DO icol=1,ISIZE
              IF(SVD_CUTOFF_RATIO.LT.SIGMA_PHI(icol)/SIGMA_PHI(1)) THEN
                irow=irow+1
              ENDIF
            ENDDO

            WRITE(OP_STRING(1),'(''# Rank           : '',I12)')
     '        ISIZE
            WRITE(OP_STRING(2),'(''# Effective rank : '',I12)')
     '        irow
            WRITE(OP_STRING(3),'(''# Cutoff value   : '',F12.5)')
     '        SVD_CUTOFF_RATIO

          ELSE
            WRITE(OP_STRING(1),'(''# Rank           : '',I12)')
     '        ISIZE
            WRITE(OP_STRING(2),'(''# Effective rank : full'')')
            WRITE(OP_STRING(3),'(''# Cutoff value   : none'')')
          ENDIF
          CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)


          WRITE(OP_STRING(1),'(''#          '')')
          WRITE(OP_STRING(2),'(''# Rank    SV           norm SV'')')
          CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)

          DO icol=1,ISIZE
            WRITE(OP_STRING,'(I5,''  '',E12.5,''  '',E12.5)')
     '        icol,SIGMA_PHI(icol),SIGMA_PHI(icol)/SIGMA_PHI(1)
            CALL WRITES(IOFILE1,OP_STRING,ERROR,*9999)
          ENDDO

C*** Close file
          CALL CLOSEF(IOFILE1,ERROR,*9999)

        ENDIF !OPDATA

        EVALUATE_PHI=.TRUE.

      ENDIF

      CALL EXITS('EVPHI')
      RETURN
 9999 CALL ERRORS('EVPHI',ERROR)
C***  Close files correctly
      IF(SIGNAL) THEN
        CALL IOSIGN(IOFILE1,LD,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID,ZD,'CLOSE',FILEFORMAT,SIGNALFNAME,
     '    ' ',ENDFILE,.TRUE.,ERROR_DUMMY,*9998)
      ENDIF
 9998 IF(HISTORY) THEN
        ENDFILE2=.TRUE.
        CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,
     '    NIYLIST,NPNY(0,1,0,nx),
     '    NQNY(1,1,0,nx),NRLIST,TRSF_NRLIST2,NUMTIMEDATA1,nx,
     '    NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,
     '    YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,HISTORYFNAME,' ',
     '    ENDFILE2,YPDATA,YQDATA,YQSDATA,ENDFILE,ERROR_DUMMY,*9997)
      ENDIF
 9997 CALL EXITS('EVPHI')
      RETURN 1
      END

      
