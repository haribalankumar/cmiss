      SUBROUTINE APTRSF(ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '  NHP,NHQ,NKH,NPLIST3,
     '  NPLIST4,NPLIST5,NPNY,NQNY,NVHP,NXLIST,NYNP,NYNR,NYQNR,
     '  PHI_H,T_BH,YP,YQ,YQS,STRING,ERROR,*)

C#### Subroutine: APTRSF
C###  Description:
C###    APTRSF applies the transfer matrix T_BH to a set of
C###    nodal values.
C**** List of nodes to which T_bh is to be applied are stored in
C**** NPLIST1, results are for nodes in NPLIST2.  Define transfer
C**** should have been called to set these up.
C**** T_bh is applied to values stored in YP - new calculated values
C**** are also stored in YP (old values at this location are stored
C**** in YP(ny,7) - the analytic location).

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List

      INTEGER ISIZE_PHI(2),ISIZE_PHIH(2),ISIZE_TBH(2),
     '  NHP(NPM,0:NRM,NXM),NHQ(NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),
     '  NPLIST3(0:NPM),NPLIST4(0:NPM),NPLIST5(0:NPM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NQNY(2,NYQM,0:NRCM,NXM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
      REAL*8  PHI_H(NY_TRANSFER_M,NTSM),
     '  T_BH(NY_TRANSFER_M,NY_TRANSFER_M),
     '  YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM)
      CHARACTER STRING*(MXCH),ERROR*(*)

!     Local Variables
      INTEGER CALC_SAMPLE_FROM_TIME,DERIVATIVE_TYPE,ERR,i,IBEG,IBEG1,
     '  icol,IEND,IEND1,
     '  irow,ITIME_LENGTH,j,N1LIST,N3CO,na,
     '  nh1,nh2,nhx1,nhx2,NIQLIST(0:1),NIQSLIST(0:1),NIYLIST(0:16),
     '  nk1,nk2,nolist1,nolist2,nonr,notime,no_nynr,
     '  np1,np2,nr_inner,NRLIST_LOC(0:9),NRLIST2_LOC(0:9),nr_outer,
     '  nts,NUMTIMEDATA,NUMTIMEDATA1,nv1,nv2,nx,nxc,ny1,
     '  ny2,nytr,start_sample,stop_sample,TRSF_NRLIST2(0:9)
      REAL*8 CALC_FORWARD_ACTN,CALC_TIME_FROM_SAMPLE,RFROMC,
     '  SUM,TEND,TSTART,TIME,YPMAX(16),YPMIN(16)
      CHARACTER ERROR1*255,FILE*(MXCH),FILEFORMAT*6,
     '  INFILE*(MXCH),OUTFILE*(MXCH)
      LOGICAL ACTIVATION,CBBREV,ENDFILE,FINISHED,HISTORY,
     '  INLIST,INTEGRATE,OUTARRAY,YPDATA,YQDATA,YQSDATA

      CALL ENTERS('APTRSF',*9999)
      CALL STRING_TRIM(FILE00,IBEG1,IEND1)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM apply transfer
C###  Parameter:     <class #[1]>
C###  Description:
C###    Apply the transfer matrix to a set of nodal values defined by
C###    the class number (of solve type).

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM apply transfer activation
C###  Parameter:      <tstart #[0.0]>
C###    Defines the lower time limit for which the potential is
C###    calulated.
C###  Parameter:      <tend #[end]>
C###    Defines the upper time limit for which the potential is
C###    calulated.
C###  Parameter:      <outfile FILENAME[$current_new]>
C###    Specifies the output filename
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the element file is stored as a binary or
C###    ascii file.
C###  Parameter:      <(direct/integrate)[direct]>
C###    Specifies whether to calculate the body surface potentials
C###    directly or to calculate the time derivatives and then
C###    integrate.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Apply the transfer matrix (double layer) to a set of nodal
C###    activation times (stored in YP(ny,1))
C###    and write out the resulting potentials (PHI) in a history file.
C###    tstart and tend are the limits of the time interval for
C###    which the potentials are calculated.

        OP_STRING(1)=STRING(1:IEND)//' activation'
        OP_STRING(2)=BLANK(1:15)//'<tstart #[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<tend #[end]>'
        OP_STRING(4)=BLANK(1:15)//'<outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(5)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(6)=BLANK(1:15)//'<(direct/integrate)[direct]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM apply transfer activation outarray
C###  Parameter:      <tstart #[0.0]>
C###    Defines the lower time limit for which the potential is
C###    calulated.
C###  Parameter:      <tend #[end]>
C###    Defines the upper time limit for which the potential is
C###    calulated.
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the file is stored as a binary or
C###    ascii file.
C###  Parameter:      <(direct/integrate)[direct]>
C###    Specifies whether to calculate the body surface potentials
C###    directly or to calculate the time derivatives and then
C###    integrate.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Apply the transfer matrix (double layer) to a set of nodal
C###    activation times (stored in YP(ny,1))
C###    and store the resulting potentials in PHI_H.
C###    tstart and tend are the limits of the time interval for
C###    which the potentials are calculated.

        OP_STRING(1)=STRING(1:IEND)//' activation outarray'
        OP_STRING(2)=BLANK(1:15)//'<tstart #[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<tend #[end]>'
        OP_STRING(4)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(5)=BLANK(1:15)//'<(direct/integrate)[direct]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM apply transfer history
C###  Parameter:      <infile FILENAME[$current]>
C###    Specifies the input filename of the nodal values
C###  Parameter:      <outfile FILENAME[$current_new]>
C###    Specifies the output filename of the changed nodal values
C###  Parameter:      <(ascii/binary)[ascii]>
C###    Specify whether the element file is stored as a binary or
C###    ascii file.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Apply the transfer matrix T_BH to a set of nodal
C###    values stored in a history file and write out the result
C###    to another history file.

        OP_STRING(1)=STRING(1:IEND)//' history <infile FILENAME['
     '    //FILE00(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<outfile FILENAME['
     '    //FILE00(IBEG1:IEND1)//'_new]>'
        OP_STRING(3)=BLANK(1:15)//'<(ascii/binary)[ascii]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','APTRSF',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.aptrsf','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          IOFI=IOOP
        ENDIF

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        CALL ASSERT(NPLIST3(0).GT.0,'>>Call define transfer first',
     '    ERROR,*9999)
        CALL ASSERT(NPLIST4(0).GT.0,'>>No body surface nodes ??',
     '    ERROR,*9999)
        CALL ASSERT(EVALUATE_TRANSFER,'>>Evaluate transfer first',
     '    ERROR,*9999)

        IF(CBBREV(CO,'ACTIVATION',1,noco+1,NTCO,N3CO)) THEN
          ACTIVATION=.TRUE.
        ELSE
          ACTIVATION=.FALSE.
        ENDIF
        IF(ACTIVATION) THEN
          IF(CBBREV(CO,'OUTARRAY',4,noco+1,NTCO,N3CO)) THEN
            OUTARRAY=.TRUE.
          ELSE
            OUTARRAY=.FALSE.
          ENDIF
          IF(CBBREV(CO,'TSTART',2,noco+1,NTCO,N3CO)) THEN
            TSTART=RFROMC(CO(N3CO+1))
          ELSE
            TSTART=0.0d0
          ENDIF
          IF(CBBREV(CO,'TEND',2,noco+1,NTCO,N3CO)) THEN
            TEND=RFROMC(CO(N3CO+1))
          ELSE
C LKC 11-SEPT-2002 this is not the way to compute TEND
C           TEND=DBLE(ISIZE_PHI(2)) ! Number of timesteps in PHI
            TEND=CALC_TIME_FROM_SAMPLE(ISIZE_PHI(2))
          ENDIF
          IF(CBBREV(CO,'OUTFILE',4,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            OUTFILE=FILE00(IBEG1:IEND1)//'_new'
          ENDIF

          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)

          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF

          IF(CBBREV(CO,'INTEGRATE',2,noco+1,NTCO,N3CO)) THEN
            INTEGRATE=.TRUE.
          ELSE
            INTEGRATE=.FALSE.
          ENDIF

        ENDIF

        HISTORY=.FALSE.
        IF(CBBREV(CO,'HISTORY',2,noco+1,NTCO,N3CO)) THEN
          HISTORY=.TRUE.
          IF(CBBREV(CO,'INFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            INFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            INFILE=FILE00(IBEG1:IEND1)
          ENDIF

          IF(CBBREV(CO,'OUTFILE',2,noco+1,NTCO,N3CO)) THEN
            CALL STRING_TRIM(CO(N3CO+1),IBEG1,IEND1)
            OUTFILE=CO(N3CO+1)(IBEG1:IEND1)
          ELSE
            OUTFILE=FILE00(IBEG1:IEND1)//'_new'
          ENDIF

          CALL STRING_TRIM(INFILE,IBEG,IEND)
          CALL STRING_TRIM(OUTFILE,IBEG1,IEND1)
          IF(IBEG.EQ.IBEG1.AND.IEND.EQ.IEND1.AND.
     '      INFILE(IBEG:IEND).EQ.OUTFILE(IBEG1:IEND1)) THEN
            ERROR='>>Input and output file names must be different'
            GOTO 9999
          ENDIF
          IF(CBBREV(CO,'BINARY',2,noco+1,NTCO,N3CO)) THEN
            FILEFORMAT='BINARY'
          ELSE
            FILEFORMAT='ASCII'
          ENDIF
        ENDIF

        nr_inner=TRSF_NR_FIRST
        nr_outer=TRSF_NR_SECOND

        CALL ASSERT(NIYM.GE.7,'>>NIYM must be at least 7',ERROR,
     '    *9999)

        IF(ACTIVATION.OR.HISTORY) THEN
C***      Move existing solution (if any) to the analytic solution
C***      location
          DO no_nynr=1,NYNR(0,0,1,nr_outer,nx)
            ny1=NYNR(no_nynr,0,1,nr_outer,nx)
            IF(NPNY(0,ny1,0,nx).EQ.1) THEN
              np1=NPNY(4,ny1,0,nx)
              IF(INLIST(np1,NPLIST5(1),NPLIST5(0),N1LIST)) THEN
                !np1 is on the outer surface
                YP(ny1,7,nx)=YP(ny1,1,nx)
                YP(ny1,1,nx)=0.0d0
              ENDIF !inlist
            ELSE
              ERROR='>>No np1 found for this ny1'
              GOTO 9999
            ENDIF
          ENDDO !no_nynr
        ENDIF !activation or history

        IF(ACTIVATION) THEN
C***      Will construct time-varying potentials on the second
C***      surface using the standard UDL transfer matrix eqtn:
C***          phi(y,t) = int(surface)[T_BH*phi_m] dSx
C***      The efficient way of doing this is to differentiate the above
C***      (e.g. if phi_m is the heaviside step function, H, then it
C***      differentiates to a dirac delta) then construct a temporal
C***      matrix of the derivatives (stored in phi) then integrate this
C***      and write to a history file.

          CALL ASSERT(USE_TIME.EQ.1,'>>USE_TIME must be set to 1',
     '      ERROR,*9999)

C***      Check that T_BH is not all zero (e.g. has not been read in
C***      from an incorrect file)

          i=1
          j=1
          FINISHED=.FALSE.
          DO WHILE ((.NOT.TBH_OK).AND.(.NOT.FINISHED))
            TBH_OK=(DABS(T_BH(i,j)).GT.ZERO_TOL)
            i=i+1
            IF(i.GT.ISIZE_TBH(1)) THEN
              i=1
              j=j+1
            ENDIF
            IF(j.GT.ISIZE_TBH(2)) THEN
              FINISHED=.TRUE.
            ENDIF
          ENDDO
          CALL ASSERT(TBH_OK,'>>T_BH matrix is zero?',ERROR,*9999)

          CALL ASSERT(TEND-TSTART.GT.RDELTA,
     '      '>>Specified activation period is 0 - check TEND,TSTART',
     '      ERROR,*9999)

          CALL ASSERT(TSTART.LT.TEND,'>>Start time after end time?',
     '      ERROR,*9999)

          ERR=0
          start_sample=CALC_SAMPLE_FROM_TIME(TSTART,ERR,ERROR)
          IF(ERR.NE.0) GOTO 9999
          stop_sample=CALC_SAMPLE_FROM_TIME(TEND,ERR,ERROR)
          IF(ERR.NE.0) GOTO 9999

          ITIME_LENGTH=stop_sample-start_sample+1
          CALL ASSERT(ITIME_LENGTH.LE.NTSM,'>>NTSM too small',
     '      ERROR,*9999)

C CPB 17/9/01 It is just PHI_H that is initialised here. Why the IF????

          IF(OUTARRAY) THEN
C***        PHI_H is used to temporarily store PHI values (since this
C***        may have been set by e.g. the optimisation
            DO nytr=1,ISIZE_TBH(1)
              DO nts=1,ITIME_LENGTH
                PHI_H(nytr,nts)=0.0d0
              ENDDO !nts
            ENDDO !nytr
          ELSE
C***        Initialise PHI

            DO nytr=1,ISIZE_TBH(1)
              DO nts=1,NTSM        ! Make sure whole array is zeroed
                PHI_H(nytr,nts)=0.0d0
              ENDDO !nts
            ENDDO !nytr

          ENDIF !OUTARRAY
          PHI_H_IS_HEART=.FALSE.

C LKC 9-SEP-2002 need to setup ISIZE_PHIH
          ISIZE_PHIH(1)=ISIZE_TBH(1) !electrodes
          ISIZE_PHIH(2)=ITIME_LENGTH !time steps



C***      Loop over inner surface nodes and find activation time for
C***      that node.

          IF(INTEGRATE) THEN
            DERIVATIVE_TYPE=1
          ELSE
            DERIVATIVE_TYPE=0
          ENDIF

          ERR=0
          DO nts=start_sample,stop_sample
            TIME=CALC_TIME_FROM_SAMPLE(nts)
            DO nytr=1,ISIZE_TBH(1)
              PHI_H(nytr,nts)=CALC_FORWARD_ACTN(DERIVATIVE_TYPE,0,0,0,
     '          nytr,NHP(1,0,nx),NKH,NPLIST3,NVHP,nx,NYNP,
     '          0.5d0/TRSF_FREQUENCY,T_BH,TIME,YP,
     '          TRSF_ACTN_WAVE_INTERPOLATE,ERR,ERROR)
              IF(ERR.NE.0) GOTO 9999
            ENDDO !nytr
          ENDDO !nts

C***      If integrate is .TRUE. has been used then d(phi)/dt is
C***      stored in PHI, otherwise PHI contains the surface potentials.

          NIYLIST(0)=1
          NIYLIST(1)=1
          NIQLIST(0)=0
          NIQSLIST(0)=0
          na=1

          IF(.NOT.OUTARRAY) THEN

            YPDATA=.TRUE.
            YQDATA=.FALSE.
            YQSDATA=.FALSE.

C***        Open the output history file
            CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST,
     '        NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '        FILEFORMAT,OUTFILE,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)

          ENDIF
          ENDFILE=.FALSE.

          DO nts=start_sample,stop_sample

            TIME=CALC_TIME_FROM_SAMPLE(nts)
            IF(INTEGRATE) THEN
C***          Integrate using simple trapezoidal scheme
              IF(nts.EQ.start_sample) THEN
                DO nytr=1,ISIZE_TBH(1)
                  PHI_H(nytr,nts)=0.5d0*(PHI_H(nytr,nts+1)+
     '              PHI_H(nytr,nts))
                ENDDO !nytr
              ELSE IF(nts.LE.stop_sample-1) THEN
                DO nytr=1,ISIZE_TBH(1)
                  PHI_H(nytr,nts)=0.5d0*(PHI_H(nytr,nts+1)+
     '              PHI_H(nytr,nts))+PHI_H(nytr,nts-1)
                ENDDO !nytr
              ELSE
                DO nytr=1,ISIZE_TBH(1)
                  PHI_H(nytr,nts)=PHI_H(nytr,nts-1)
                ENDDO !nytr
              ENDIF !nts
            ENDIF !integrate
            IF(.NOT.OUTARRAY) THEN
C***          Evaluate YP array for second surface.
              irow=0
              DO nolist2=1,NPLIST4(0) !list of second surface nodes
                np2=NPLIST4(nolist2)
                DO nhx2=1,NHP(np2,nr_outer,nx)
                  nh2=NH_LOC(nhx2,nx)
                  DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                    DO nk2=1,MAX(NKH(nh2,np2,1,nr_outer)-
     '                KTYP93(1,nr_outer),1)  !body column number of GK
                      irow=irow+1
                      ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
                      YP(ny2,1,nx)=PHI_H(irow,nts)
                    ENDDO !nk2
                  ENDDO !nv2
                ENDDO !nh2
              ENDDO !nolist2

C***          Write out new solution
              YPDATA=.TRUE.
              YQDATA=.FALSE.

              CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST,
     '          NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '          FILEFORMAT,OUTFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '          YQDATA,YQSDATA,ERROR,*9999)
            ENDIF !outarray
          ENDDO !End of times

          IF(.NOT.OUTARRAY) THEN
C***        Close history file

            YPDATA=.TRUE.
            YQDATA=.FALSE.
            CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '        NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST,
     '        NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '        TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',
     '        FILEFORMAT,OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,
     '        YQSDATA,ERROR,*9999)
          ENDIF

        ELSE IF(HISTORY) THEN
C***      Applying a transfer matrix to a history file.

          TRSF_NRLIST2(0)=TRSF_NRLIST(0)
          DO nonr=1,TRSF_NRLIST(0)
            TRSF_NRLIST2(nonr)=TRSF_NRLIST(nonr)
          ENDDO

          NIYLIST(0)=1
          NIYLIST(1)=1
          NIQLIST(0)=0
          NIQSLIST(0)=0

C***      Open the input history file
          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',FILEFORMAT,
     '      INFILE,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

C***      Check on the number of regions in the history file and
C***      decide what to do if different.
          IF(NRLIST_HIST(0,IOFILE1).LT.TRSF_NRLIST(0)) THEN
C***        Less regions than current setup in the history file.  This
C***        can be obtained by e.g. fitting measured signals on the
C***        outside of a homogeneous torso and then applying a
C***        multi-region inverse matrix to it.
            IF(NRLIST_HIST(0,IOFILE1).EQ.1) THEN
C***          In this case assume that the data should be read into the
C***          last location
              NRLIST_LOC(0)=1
              NRLIST_LOC(1)=NRLIST_HIST(1,IOFILE1)
            ELSE
              WRITE(OP_STRING,
     '          '('' >> WARNING: Do not know where to store data'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ELSEIF(NRLIST_HIST(0,IOFILE1).GT.TRSF_NRLIST(0)) THEN
C***        More regions were used to generate the data than the
C***        current setup has.
            IF(NRLIST_HIST(0,IOFILE1).EQ.TRSF_NRLIST(0)+1) THEN
C***          Will occur when e.g. the data has been generated using
C***          a moving dipole inside the heart - the interior of the
C***          heart being the extra region not needed in the transfer
              CALL ASSERT(NRLIST_HIST(0,IOFILE1).LE.9,
     '          '>>NRLIST_LOC not big enough',ERROR,*9999)
              NRLIST_LOC(0)=TRSF_NRLIST(0)
              DO nolist1=1,NRLIST_LOC(0)
                NRLIST_LOC(nolist1)=TRSF_NRLIST(nolist1)
              ENDDO
              DO nolist1=NRLIST_LOC(0)+1,NRLIST_HIST(0,IOFILE1)
C               !initialising rest of NRLIST_LOC - needed in IOHIST.
                NRLIST_LOC(nolist1)=0
              ENDDO
            ELSE

C LKC 1-MAR-2000 Need to do something when the transfer and history
C  is mismatched by more that one.
C
C***  This may occur when applying a epi2torso transfer matrix to a
C***  history file which has been generated by a
C***  transmembrane2epicardial transfer matrix - this will
C***  contain all (?) the regions of the model.


              WRITE(OP_STRING,
     '          '('' >> WARNING: Do not know where to store data'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING(1),
     '          '('' >>  '',I2,'' Regions contained in transfer'')')
     '          TRSF_NRLIST(0)
              WRITE(OP_STRING(2),
     '          '('' >>  '',I2,'' Regions contained in history'')')
     '          NRLIST_HIST(0,IOFILE1)
              WRITE(OP_STRING(3),
     '          '('' >>   Using the '',I2,'' transfer regions'')')
     '          TRSF_NRLIST(0)
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)


              CALL ASSERT(NRLIST_HIST(0,IOFILE1).LE.9,
     '          '>>NRLIST_LOC not big enough',ERROR,*9999)
              NRLIST_LOC(0)=TRSF_NRLIST(0)
              DO nolist1=1,NRLIST_LOC(0)
                NRLIST_LOC(nolist1)=TRSF_NRLIST(nolist1)
              ENDDO
              DO nolist1=NRLIST_LOC(0)+1,NRLIST_HIST(0,IOFILE1)
C               !initialising rest of NRLIST_LOC - needed in IOHIST.
                NRLIST_LOC(nolist1)=0
              ENDDO

            ENDIF
          ELSE
C***        Problem setup and history file have the same number of
C***        regions.  Could check that the number of nys etc are the
C***        same, but for now just assume it is okay.
            NRLIST_LOC(0)=TRSF_NRLIST(0)
            DO nolist1=1,NRLIST_LOC(0)
              NRLIST_LOC(nolist1)=TRSF_NRLIST(nolist1)
            ENDDO
          ENDIF
          na=1

          NRLIST2_LOC(0)=NRLIST_LOC(0)
          DO nonr=1,NRLIST_LOC(0)
            NRLIST2_LOC(nonr)=NRLIST_LOC(nonr)
          ENDDO

C***      Open the output history file
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,
     '      OUTFILE,'OPEN',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

          IF(FILEFORMAT(1:5).EQ.'ASCII') THEN
            ENDFILE=.FALSE.
            NUMTIMEDATA=0
C***        Read YP for each time series
            DO WHILE(.NOT.ENDFILE)

              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,INFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '          YQDATA,YQSDATA,ERROR,*9999)
              IF(.NOT.ENDFILE) THEN
C***            Apply T_bh to first surface nodes.
                irow=0
                DO nolist2=1,NPLIST4(0) !list of torso nodes
                  np2=NPLIST4(nolist2)
                  DO nhx2=1,NHP(np2,nr_outer,nx)
                    nh2=NH_LOC(nhx2,nx)
                    DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                      DO nk2=1,MAX(NKH(nh2,np2,1,nr_outer)-
     '                  KTYP93(1,nr_outer),1)
                        ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
                        !body column number of GK
                        irow=irow+1
                        icol=0
                        SUM=0.0d0
                        DO nolist1=1,NPLIST3(0) !list of heart nodes
                          np1=NPLIST3(nolist1)
                          DO nhx1=1,NHP(np1,nr_inner,nx)
                            nh1=NH_LOC(nhx1,nx)
                            DO nv1=1,NVHP(nh1,np1,1,nr_inner)
                              DO nk1=1,MAX(NKH(nh1,np1,1,nr_inner)-
     '                          KTYP93(1,nr_inner),1)
                                ny1=NYNP(nk1,nv1,nh1,np1,0,1,nr_inner)
                                !heart column number of GK
                                icol=icol+1
                                SUM=SUM+T_BH(irow,icol)*YP(ny1,1,nx)
                              ENDDO !nk1
                            ENDDO !nv1
                          ENDDO !nh1
                        ENDDO !np1
                        YP(ny2,1,nx)=SUM
                      ENDDO !nk2
                    ENDDO !nv2
                  ENDDO !nh2
                ENDDO !np2

C***            Write out new solution

                YPDATA=.TRUE.
                YQDATA=.FALSE.
                CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,
     '            NIYLIST,NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,
     '            TRSF_NRLIST2,NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),
     '            NYQNR(0,0,1,0,nx),TIME,YP(1,1,nx),YPMAX,YPMIN,
     '            YQ(1,1,1,nx),YQS,'WRITE',FILEFORMAT,OUTFILE,
     '            'TIME_DATA',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '            ERROR,*9999)

              ENDIF !not end of file
            ENDDO !Hit end of file

          ELSE IF(FILEFORMAT(1:6).EQ.'BINARY') THEN
            DO notime=1,NUMTIMEDATA
C***          Read the history data

              CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '          NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'READ',
     '          FILEFORMAT,INFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '          YQDATA,YQSDATA,ERROR,*9999)

C***          Apply T_bh to first surface nodes.
              irow=0
              DO nolist2=1,NPLIST4(0) !list of torso nodes
                np2=NPLIST4(nolist2)
C                CALL NODE_CHANGE(np2,.FALSE.,ERROR,*9999)
                DO nhx2=1,NHP(np2,nr_outer,nx)
                  nh2=NH_LOC(nhx2,nx)
                  DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                    DO nk2=1,MAX(NKH(nh2,np2,1,nr_outer)-
     '                KTYP93(1,nr_outer),1)
                      ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
                      !body column number of GK
                      irow=irow+1
                      icol=0
                      SUM=0.0d0
                      DO nolist1=1,NPLIST3(0) !list of heart nodes
                        np1=NPLIST3(nolist1)
                        DO nhx1=1,NHP(np1,nr_inner,nx)
                          nh1=NH_LOC(nhx1,nx)
                          DO nv1=1,NVHP(nh1,np1,1,nr_inner)
                            DO nk1=1,MAX(NKH(nh1,np1,1,nr_inner)-
     '                        KTYP93(1,nr_inner),1)
                              ny1=NYNP(nk1,nv1,nh1,np1,0,1,nr_inner)
                              !heart column number of GK
                              icol=icol+1
                              SUM=SUM+T_BH(irow,icol)*YP(ny1,1,nx)
                            ENDDO !nk1
                          ENDDO !nv1
                        ENDDO !nh1
                      ENDDO !np1
                      YP(ny2,1,nx)=SUM
                    ENDDO !nk2
                  ENDDO !nv2
                ENDDO !nh2
              ENDDO !np2

C***          Write out new solution

              YPDATA=.TRUE.
              YQDATA=.FALSE.
              CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '          NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '          NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),
     '          TIME,YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'WRITE',
     '          FILEFORMAT,OUTFILE,'TIME_DATA',ENDFILE,.TRUE.,YPDATA,
     '          YQDATA,YQSDATA,ERROR,*9999)

            ENDDO !notime
          ENDIF !fileformat

C***      Close both history files

          CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '      NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      INFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR,*9999)
          CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '      NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '      NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '      YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '      OUTFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)


C*** not activation or history file
        ELSE !not activation or history file

C***      Apply T_bh to first surface nodes.

          irow=0
          DO nolist2=1,NPLIST4(0) !list of torso nodes
            np2=NPLIST4(nolist2)
C            CALL NODE_CHANGE(np2,.FALSE.,ERROR,*9999)
            DO nhx2=1,NHP(np2,nr_outer,nx)
              nh2=NH_LOC(nhx2,nx)
              DO nv2=1,NVHP(nh2,np2,1,nr_outer)
                DO nk2=1,
     '            MAX(NKH(nh2,np2,1,nr_outer)-KTYP93(1,nr_outer),1)
! ny2=NYNP(nk2,nv2,nh2,np2,1,1,nr_outer) !body row number of GK

C LKC not required?
                  ny2=NYNP(nk2,nv2,nh2,np2,0,1,nr_outer)
                  !body column number of GK
                  irow=irow+1
                  icol=0
                  SUM=0.0d0
                  DO nolist1=1,NPLIST3(0) !list of heart nodes
                    np1=NPLIST3(nolist1)
                    DO nhx1=1,NHP(np1,nr_inner,nx)
                      nh1=NH_LOC(nhx1,nx)
                      DO nv1=1,NVHP(nh1,np1,1,nr_inner)
                        DO nk1=1,MAX(NKH(nh1,np1,1,nr_inner)-
     '                    KTYP93(1,nr_inner),1)
                          ny1=NYNP(nk1,nv1,nh1,np1,0,1,nr_inner)
                          !heart column number of GK
                          icol=icol+1
                          SUM=SUM+T_BH(irow,icol)*YP(ny1,1,nx)
                        ENDDO !nk1
                      ENDDO !nv1
                    ENDDO !nh1
                  ENDDO !np1
                  YP(ny2,1,nx)=SUM
                ENDDO !nk2
              ENDDO !nv2
            ENDDO !nh2
          ENDDO !np2
        ENDIF
        APPLY_TRANSFER=.TRUE.
      ENDIF


      CALL EXITS('APTRSF')
      RETURN
 9999 CALL ERRORS('APTRSF',ERROR)
C*** Close history files correctly
      CALL IOHIST(IOFILE1,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '  NPNY(0,1,0,nx),NQNY(1,1,0,nx),NRLIST_LOC,NRLIST2_LOC,
     '  NUMTIMEDATA,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '  YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '  INFILE,' ',ENDFILE,.TRUE.,YPDATA,YQDATA,YQSDATA,ERROR1,*9998)
 9998 YQDATA=.TRUE.
      YQSDATA=.FALSE.
      YPDATA=.FALSE.
      CALL IOHIST(IOFILE2,na,NHQ(1,nx),NIQLIST,NIQSLIST,NIYLIST,
     '  NPNY(0,1,0,nx),NQNY(1,1,0,nx),TRSF_NRLIST,TRSF_NRLIST2,
     '  NUMTIMEDATA1,nx,NYNR(0,0,1,0,nx),NYQNR(0,0,1,0,nx),TIME,
     '  YP(1,1,nx),YPMAX,YPMIN,YQ(1,1,1,nx),YQS,'CLOSE',FILEFORMAT,
     '  OUTFILE,' ',ENDFILE,.TRUE.,YQDATA,YQSDATA,YPDATA,ERROR1,*9997)
 9997 CALL EXITS('APTRSF')
      RETURN 1
      END


