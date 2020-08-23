      SUBROUTINE FITSIG_SPLINE(IBT,IDO,INP,LD,LN,NBH,NBJ,NDDATA,NEELEM,
     '  nelec,NHQ,NKJE,NPF,NPNE,NPNODE,NPNY,NQNY,nslice,nr,NRE,NRLIST,
     '  NRLIST2,NVJE,nx,NYNR,NYQNR,SE,START_NODE,SNODE_NUM,TEND,TSTART,
     '  WD,XID,XE,XP,YP,YQ,YQS,ZD,FILEFORMAT,ERROR,*)

C#### Subroutine: FITSIG_SPLINE
C###  Description:
C###    FITSIG fits field variables to signal data using cubic splines

C###  Comment: Cubic Spline Signal Fitting
C###  Description:
C###    <HTML>
C###    Note : In this routine the terminology
C###    refering to nodes and data points have been reverse. It must be
C###    remembered that the nodes which we are trying to fit to have
C###    been defined as data points and their xi locations are the
C###    locations of the fitted nodes (nr) in relation to the
C###    spline mesh (nr_spline). The data points from the
C###    sampled signal information has been defined into a bi-linear
C###    finite element mesh.
C###    <BR>
C###    <H4>Some variable specific to this routine</H4>
C###
C###    <UL>
C###      <LI> P(NKNOTSM,NKNOTSM) - stores 2nd deriv information
C###              for each horizontal spline
C###      <LI> P2(NKNOTSM) - stores tempoary 2nd deriv info
C###              for vert splines
C###      <LI> UPPERDIAG(NKNOTSM) - upper diagonal of matrix
C###      <LI> DIAG(KNOTSM)       - diagonal of matrix
C###    </UL>

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'sign00.cmn'
      INCLUDE 'time02.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LD(NDM),LN(0:NEM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NEELEM(0:NE_R_M,0:NRM),NHQ(NRM),
     '  NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM),NQNY(2,NYQM,0:NRCM),nr,
     '  NRE(NEM),NRLIST(0:NRM),NRLIST2(0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYQNR(0:NYQM,0:NRCM,NCM,0:NRM),
     '  SNODE_NUM,START_NODE
      REAL*8 SE(NSM,NBFM,NEM),TEND,TSTART,WD(NJM,NDM),
     '  XA(NAM,NJM,NEM),XID(NIM,NDM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),YQ(NYQM,NIQM,NAM),
     '  YQS(NIQSM,NQM),ZD(NJM,NDM)
      CHARACTER ERROR*(*),FILEFORMAT*(*)

!     Local Parameters
      INTEGER NKNOTSM
      PARAMETER(NKNOTSM=40)

!     Local Variables
      INTEGER layer,LD_SIGN(NDM),na,nb,nd,nelec,nh,nhj,nhx,ni,nj,
     '  nn,nnfit,nn_bottom,nn_spline,nr_spline,nslice,NIQLIST(0:1),
     '  NIQSLIST(0:1),NIYLIST(0:1),notime,NUMTIMEDATA,NUMTIMEDATA1,row
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      REAL*8 DIAG(NKNOTSM),PXI,P(NKNOTSM,NKNOTSM),P2(NKNOTSM),
     '  UPPERDIAG(NKNOTSM),SIGNALMAX(NSIGNALREGIONSMX),
     '  SIGNALMIN(NSIGNALREGIONSMX),TIME,XI(NJM),XID_SIGN(NIM,NDM),
     '  YPMAX(1),YPMIN(1)
      LOGICAL ENDFILE,NEXTTIME,YPDATA,YQDATA,YQSDATA
      CHARACTER CERROR*50

! Functions
      REAL*8 SPLINE_INTERP

      CALL ENTERS('FITSIG_SPLINE',*9999)
      TIME=0.D0


      IF(nr.EQ.DATA_REGION) THEN
        WRITE(OP_STRING,'(''WARNING: Fit region == Data Region'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      nr_spline=DATA_REGION   !region where fitting mesh is defined
      CALL ASSERT(CALL_DATA,'>>Define data first',ERROR,*9999)
      CALL ASSERT(NJM.GE.NJT+1,'>>Increase NJM',ERROR,*9999)
      CALL ASSERT(NDM.GE.NDT+NSLICE,
     '  '>>Increase NDM to NDT+NSLICE',ERROR,*9999)
C     NOTE : allows for temp storage of vertical splines

      CALL ASSERT(NELEC*NSLICE.EQ.NPNODE(0,nr_spline),
     '  '>>NELEC*NSLICE.NE.NPNODE(0,nr_spline)',ERROR,*9999)



C*** Need to linearfield defined correctly for element field and
C    electrode region

      DO nhj=1,NUM_FIT(1)
        nhx=NLH_FIT(nhj,3,1)
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh,1,LN(1))
        DO ni=1,NIT(nb)
          CALL ASSERT(IBT(1,ni,nb).EQ.1
     '      .OR.IBT(1,ni,nb).EQ.5
     '      .OR.IBT(1,ni,nb).EQ.6,
     '      '>>Use Lagrange for ELFD ',ERROR,*9999)
          CALL ASSERT(IBT(2,ni,nb).EQ.1,
     '      '>>Use Linear bfn for ELFD ',ERROR,*9999)
        ENDDO !ni
      ENDDO !nhj

      nb=NBJ(1,NEELEM(1,nr_spline) )
      DO ni=1,NIT(nb)
        CALL ASSERT(IBT(1,ni,nb).EQ.1,
     '    '>>Use Lagrange for ELFD ',ERROR,*9999)
        CALL ASSERT(IBT(2,ni,nb).EQ.1,
     '    '>>Use Linear bfn for ELFD ',ERROR,*9999)
      ENDDO


      IF(NJ_LOC(NJL_GEOM,0,nr).GT.2) THEN !assuming a 3d torso with tops closed
        IF(NDT.GT.NPNODE(0,nr)-START_NODE+1) THEN
          WRITE(OP_STRING,
     '      '(''WARNING: More Data points than required '')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(NDT.LT.NPNODE(0,nr)-START_NODE+1) THEN
          WRITE(OP_STRING,
     '      '(''WARNING: Insufficient Data points '')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE !assuming a 2d slice
        IF(NDT.GT.NPNODE(0,nr)-START_NODE+1) THEN
          WRITE(OP_STRING,
     '      '(''WARNING: More Data points than required '')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        IF(NDT.LT.NPNODE(0,nr)-START_NODE+1) THEN
          WRITE(OP_STRING,
     '      '(''WARNING: Insufficient Data points '')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      DO nd=1,NDT
        CALL ASSERT(LD(nd).GT.0,'>> No element for nd',ERROR,*9999)
        CALL ASSERT(LD(nd).LE.NEELEM(NEELEM(0,nr_spline),nr_spline),
     '    '>> Invalid element projection',ERROR,*9999)
      ENDDO

C***  Open up signal file & read electrode data
      CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID_SIGN,ZD,'READ',FILEFORMAT,INFILENAME,
     '  'OPEN',ENDFILE,.TRUE.,ERROR,*9999)
      CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID_SIGN,ZD,'READ',FILEFORMAT,INFILENAME,
     '  'ELECTRODE_DATA',ENDFILE,.TRUE.,ERROR,*9999)

      CALL ASSERT(SIGNAL_NUMELEC(nr,IOFILE1).EQ.NPNODE(0,nr_spline)
     '  ,'>>#Signal NE #Nodes',
     '  ERROR,*9999)

C***  Calculate the number of time data points
      IF(FILEFORMAT.EQ.'ASCII') THEN
        NUMTIMEDATA=0
        ENDFILE=.FALSE.
        DO WHILE(.NOT.ENDFILE) !finds the number of time samples
          CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,
     '      SIGNALMAX,
     '      SIGNALMIN,TIME,WD,XID_SIGN,ZD,'READ',FILEFORMAT,INFILENAME,
     '      'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)
          IF(.NOT.ENDFILE) NUMTIMEDATA=NUMTIMEDATA+1
        ENDDO
      ENDIF

      CALL ASSERT(NUMTIMEDATA.GT.0,'>>No signal data',ERROR,*9999)
      CALC_XI=SIGNAL_ELEMLOC(IOFILE1).EQ.1
      CALL ASSERT(CALC_XI,'>>Define xi positions first',ERROR,*9999)

C***  Open up the history file
      NRLIST(0)=1
      NRLIST(1)=nr
      NRLIST2(0)=1
      NRLIST2(1)=nr
      NIYLIST(0)=1
      NIYLIST(1)=1
      NIQLIST(0)=0
      NIQSLIST(0)=0
      NEXTTIME=.TRUE.
      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      na=1
      CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,YPMAX,YPMIN,
     '  YQ,YQS,'WRITE',FILEFORMAT,OUTFILENAME,'OPEN',ENDFILE,NEXTTIME,
     '  YPDATA,YQDATA,YQSDATA,ERROR,*9999)


C***    Reset both files
      CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID_SIGN,ZD,'READ',FILEFORMAT,INFILENAME,
     '  'RESET',ENDFILE,.TRUE.,ERROR,*9999)
      CALL ASSERT(NUMTIMEDATA.GT.0,
     '  '>>No signal data, NUMTIMEDATA = 0',ERROR,*9999)

      IF(NIT(NBJ(1,NEELEM(1,nr_spline))).EQ.2) THEN
C       !2d Surface so assume body with closed top and bottom
        NNFIT=NPNODE(0,nr)-1
        YP(NPNODE(0,nr),1)=1.D0
      ELSE
        NNFIT=NPNODE(0,nr)
      ENDIF

      CALL CPU_TIMER(CPU_USER,TIME_START)

C**** Create Horizontal Splines
      DO layer=1,nslice
        CALL MAKESPLINE('HORIZONTAL',nelec,layer,ZD,
     '    DIAG,UPPERDIAG,P(1,layer),
     '    ERROR,*9999)
        IF(SIGNAL_OUTPUT.GE.1) THEN
          WRITE(*,'(''Creating Horizont Spline '',I2,'' '')')layer
          WRITE(*,'(''  U '',12(E11.3))')
     '      (UPPERDIAG(nn),nn=2,nelec-1)
          WRITE(*,'(''  P '',12E11.3)')(P(nn,layer),nn=2,nelec)
          WRITE(*,*)
        ENDIF
      ENDDO

      IF(SIGNAL_OUTPUT.GE.4) THEN
        WRITE(*,*)
     '    'nn Layer        X positions           Potentials'
        WRITE(*,*) 'Making Vertical Splines'
      ENDIF


C***  Loop over the number of signals

      DO notime=1,NUMTIMEDATA
C***    Read in the signal data
        CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '    SIGNALMIN,TIME,WD,XID_SIGN,ZD,'READ',FILEFORMAT,INFILENAME,
     '    'SIGNAL_DATA',ENDFILE,.TRUE.,ERROR,*9999)

        IF(MOD(notime,SIGNAL_SKIPSAMPLE).EQ.0.AND.
     '    TIME.GE.TSTART.AND.TIME.LE.TEND) THEN

          IF(SIGNAL_OUTPUT.GE.4) THEN
            WRITE(*,*)
            WRITE(*,*) ' ------ at notime # ------',notime
          ENDIF

C*** Setup for Vertical Splines
          DO nn=1,NNFIT

            IF(nn.LT.START_NODE) THEN !remove unwanted nodes from fit
              YP(nn,1)=1.D0
            ELSE

C*** Calculate the bottom node of vertical spline
              nn_bottom=MOD(nn-START_NODE+1,nelec) !for slice
              IF(nn_bottom.EQ.0) THEN !an end node
                row=INT((nn-START_NODE+1)/nelec) !finds which row
                nn_bottom=INT((nn-START_NODE+1)/row)

C 8-JUL-1998 wrong missing a bit
C                nn_bottom=INT(nn/nn_bottom)
C 22-APR new
C                nn_bottom=nn_bottom-START_NODE+1
C              ELSE
C                nn_bottom=nn_bottom+START_NODE-1
C                nn_bottom=nn_bottom
              ENDIF
              CALL ASSERT(nn_bottom.GT.0,'>>Invalid Node',ERROR,*9999)

C*** Work way up the vertical spline
              DO layer=1,nslice !each horiz spline
                nn_spline=(nn_bottom+(layer-1)*nelec)

                CALL ASSERT(LD(nn_spline).GT.0,
     '            '>> Error going up spline',ERROR,*9999)

C*** Interpolate to obtain the global coordinates from data
                CALL XPXE(NBJ(1,LD(nn_spline)),
     '            NKJE(1,1,1,LD(nn_spline)),NPF(1,1),
     '            NPNE(1,1,LD(nn_spline)),NRE(LD(nn_spline)),
     '            NVJE(1,1,1,LD(nn_spline)),SE(1,1,LD(nn_spline)),
     '            XA(1,1,LD(nn_spline)),XE,XP,ERROR,*9999)

                DO nj=1,NIT(NBJ(1,1)) !all NBJ's should be the same
                  XI(nj)=XID(nj,nn_spline)
                ENDDO

                DO nj=1,NJ_LOC(NJL_GEOM,0,nr_spline)
                  nb=NBJ(nj,LD(nn_spline))
                  ZD(nj,NDT+layer)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,1,XI,XE(1,nj))
                ENDDO !nj

C*** Interp to get new potentials
                ZD(NJT+1,NDT+layer)=
     '            SPLINE_INTERP('HORIZ',START_NODE,NBJ,nn,
     '            nelec,layer,
     '            P(1,layer),ZD,LD(nn-START_NODE+1),
     '            XID,1,NPNE,SNODE_NUM,NEELEM(1,nr_spline))

                IF(SIGNAL_OUTPUT.GE.4) THEN
                  WRITE(*,'(I4,I4,F12.4,F12.4,'' '',E12.5)')
     '              nn,layer,ZD(1,NDT+layer),ZD(1,NDT+layer)
     '              ,ZD(NJT+1,NDT+layer)
                ENDIF
              ENDDO !layer

C*** Create vertical Splines
              IF(nslice.EQ.1) THEN !only one slice ie no vert splines
                YP(nn,1)=ZD(NJT+1,NDT+nslice)
              ELSE
                CALL MAKESPLINE('VERTICAL',nslice,nn,ZD,
     '            DIAG,UPPERDIAG,P2(layer),ERROR,*9999)

                YP(nn,1)=SPLINE_INTERP('VERTI',START_NODE,NBJ,nn,
     '            nelec,layer,
     '            P2,ZD,LD(nn-START_NODE+1),XID,2,
     '            NPNE,SNODE_NUM,NEELEM(1,nr_spline))

                IF(SIGNAL_OUTPUT.GE.4) THEN
                  write (*,'(''Soln'',E12.5)') YP(nn,1)
                  write (*,*)
                ENDIF
              ENDIF !nslice
            ENDIF !in the fit
          ENDDO !nn

          NEXTTIME=.TRUE.
          YPDATA=.TRUE.
          YQDATA=.FALSE.
          YQSDATA=.FALSE.
          na=1
          CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,
     '      NQNY,NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,
     '      YPMAX,YPMIN,YQ,YQS,'WRITE',FILEFORMAT,OUTFILENAME,
     '      'TIME_DATA',ENDFILE,NEXTTIME,YPDATA,YQDATA,YQSDATA,
     '      ERROR,*9999)

          IF(SIGNAL_OUTPUT.GE.3) THEN
            WRITE(OP_STRING,'('' Fitting signal '',I5,'//
     '        ' '' at time ='',D12.4)') notime,TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !if fitting this time sample
      ENDDO !notime

      CALL CPU_TIMER(CPU_USER,TIME_STOP)
      ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)
      IF(SIGNAL_OUTPUT.GE.2) THEN
        WRITE(OP_STRING,'(/'' Total CPU time for signal fitting '
     '    //'problem: '',D11.4,'' s'')') ELAPSED_TIME
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Average CPU time per signal: '',D11.4,'
     '    //''' s'')') ELAPSED_TIME/NUMTIMEDATA
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  '')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C***  Close the history and signal file
      CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID_SIGN,ZD,'CLOSE',FILEFORMAT,INFILENAME,
     '  ' ',ENDFILE,.TRUE.,ERROR,*9998)
      NEXTTIME=.TRUE.
      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      na=1
      CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,YPMAX,YPMIN,
     '  YQ,YQS,'CLOSE',FILEFORMAT,OUTFILENAME,' ',ENDFILE,NEXTTIME,
     '  YPDATA,YQDATA,YQSDATA,ERROR,*9999)

      CALL EXITS('FITSIG_SPLINE')
      RETURN
 9999 CALL ERRORS('FITSIG_SPLINE',ERROR)

C***  Close the history and signal file
      CALL IOSIGN(IOFILE1,LD_SIGN,NBJ,NDDATA,NUMTIMEDATA,0,SIGNALMAX,
     '  SIGNALMIN,TIME,WD,XID_SIGN,ZD,'CLOSE',FILEFORMAT,INFILENAME,
     '  ' ',ENDFILE,.TRUE.,CERROR,*9998)
 9998 NEXTTIME=.TRUE.
      YPDATA=.TRUE.
      YQDATA=.FALSE.
      YQSDATA=.FALSE.
      na=1
      CALL IOHIST(IOFILE2,na,NHQ,NIQLIST,NIQSLIST,NIYLIST,NPNY,NQNY,
     '  NRLIST,NRLIST2,NUMTIMEDATA1,nx,NYNR,NYQNR,TIME,YP,YPMAX,YPMIN,
     '  YQ,YQS,'CLOSE',FILEFORMAT,OUTFILENAME,' ',ENDFILE,NEXTTIME,
     '  YPDATA,YQDATA,YQSDATA,CERROR,*9997)
 9997 CALL EXITS('FITSIG_SPLINE')
      RETURN 1
      END


