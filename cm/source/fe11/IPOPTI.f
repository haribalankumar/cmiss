      SUBROUTINE IPOPTI(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '  IDO,INP,ISIZE_MFI,ISIZE_PHI,
     '  LD,LDR,NBH,NBJ,NDIPOLES,NEELEM,NENP,NHP,
     '  NKB,NKHE,NKH,NKJ,NLNO,NMNO,NNB,NONL,NONM,NONY,
     '  NP1OPT,NP2OPT,NP3OPT,NPL,NPLIST4,
     '  NPNE,NPNODE,NPNY,nr,NVHE,NVHP,NVJP,
     '  nx_opt,nx_sol,NXI,NYNE,NYNO,NYNP,NYNR,NYNY,PAOPTY,
     '  CONY,CYNO,CYNY,DIPOLE_CEN,DIPOLE_DIR,
     '  DL,PAOPTI,PBOPTI,PMAX,PMIN,TIME,XP,YP,FIX,ERROR,*)

C#### Subroutine: IPOPTI
C###  Description:
C###    IPOPTI defines optimisation parameters for region nr.
C###    If KTYP1B=2 (derivative as angle) then there are n-1
C###    optimising angle variables instead of n components.
C###    If there are 3 components, angles are taken as
C###    Euler azimuth and angle, if there are only two
C###    components then the angle is taken between the
C###    two components.
C###    If KTYP1B=3 (non-normalised derivative) then there are no
C###    constraints on the derivative magnitude.  Components are
C###    constrained individually.

C#### Variable: NT_RES
C###  Type: INTEGER
C###  Set_up: IPOPTI,GLOBALO,NAGMINA
C###  Description:
C###    NT_RES is the number of residuals.

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp100.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER
     '  DIPOLE_CEN_NTIME(NDIPOLEM,NRM),DIPOLE_DIR_NTIME(NDIPOLEM,NRM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),ISIZE_MFI(3,NSSM),
     '  ISIZE_PHI(2),LD(NDM),LDR(0:NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDIPOLES(NRM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NHP(NPM,0:NRM,NXM),
     '  NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),
     '  NKH(NHM,NPM,NCM),NKJ(NJM,NPM),NLNO(NOPM,NXM),
     '  NMNO(1:2,0:NOPM,NXM),NNB(4,4,4,NBFM),
     '  NONL(NLM,NXM),NONM(NMM,NPM,NXM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NP1OPT(NOPM),NP2OPT(NOPM),NP3OPT(NOPM),NPL(5,0:3,NLM),
     '  NPLIST4(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),nr,NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  nx_opt,nx_sol,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYNY(0:NYYM,NYM,NRM,NXM),
     '  PAOPTY(NOPM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),CYNY(0:NYYM,NYM,NRM,NXM),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DL(3,NLM),PAOPTI(*),PBOPTI(*),PMAX(*),PMIN(*),TIME,
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)

!     Local Variables
      INTEGER i,IBEG,IBEG1,IBEG2,IBEG3,IBEG4,IBEG5,
     '  ICHAR,IEND,IEND1,IEND2,IEND3,IEND4,
     '  IEND5,IIY,INFO,IOSTAT,IREC,IY,
     '  K,NDATA,n,N1,nb,ncritical,ndat,ndipole,ne,nh,nhx,nj,njj,
     '  nj2,nj3,NJH(3),nk,nk2,nk3,NL,nm,nmlist,nn,no,no_aero,
     '  noelem,nojoin,
     '  nonode,nonode1,NONODE2,noopti,NOQUES,nores,no_wake,
     '  no_nynr,np,NP1,NP2,np3,npr,nps,nr1,nr2,nrc,
     '  nv,nv2,nv3,ny,ny2,ny3,ny_first,PRTLIM,var
      REAL*8 RAD,DATAN_MOD
      CHARACTER CHAR*2,CHAR1*30,CHAR2*1,CHAR3*2,
     '  CHAR4*4,CHAR5*6,CHAR7*6,CHAR8*6,CHAR9*6,CHAR10*2,CHAR11*3
      LOGICAL FILEIP,INLIST,ISFIXED,OPENED,SHARE,CROSSDERIV,DERIV

      CALL ENTERS('IPOPTI',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

C LKC 24-MAR-1999 Initialise MAX_MAJOR_ITER
C   Currently only used for inverse proceedures
      MAX_MAJOR_ITER=-1

C LKC 16-NOV-1999 Initialise MAX_MINOR_ITER
C   Currently only used for inverse proceedures
      MAX_MINOR_ITER=-1


C CPB 1/6/94 Adding MINOS
      FORMAT='('' Specify optimisation package [1]:'''//
     '  '/''   (1) NPSOL '''//
     '  '/''   (2) MINOS'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP29
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP29=IDATA(1)
      IF(KTYP29.EQ.1) THEN
        CALL ASSERT(USE_NPSOL.EQ.1,'>>Set USE_NPSOL to 1 to use NPSOL',
     '    ERROR,*9999)
      ELSE IF(KTYP29.EQ.2) THEN
        CALL ASSERT(USE_MINOS.EQ.1,'>>Set USE_MINOS to 1 to use MINOS',
     '    ERROR,*9999)
        FORMAT='($,'' Is the problem sparse [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(SPARSEJAC) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(2)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          SPARSEJAC=(ADATA(1).EQ.'Y')
        ENDIF
        IF(.NOT.SPARSEJAC) THEN
          CALL ASSERT(NZ_MINOSM.GE.(NOPM*(NCOM+NLCM)),
     '      '>>Increase NZ_MINOS, s.b. > NOPM*(NCOM+NLCM)',
     '      ERROR,*9999)
        ENDIF
      ENDIF

      FORMAT='('' Specify whether optimising [1]:'''//
     '  '/''   (1) Material parameters '''//
     '  '/''   (2) Geometric parameters'''//
     '  '/''   (3) Micro-structure parameters'''//
     '  '/''   (4) Holmes constitutive law constants'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP26
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP26=IDATA(1)

      IF(KTYP26.EQ.1) THEN      !optimising material parameters
        CALL ASSERT(CALL_MATE,'>>Define materials first',ERROR,*9999)
        FORMAT='('' Specify type of objective function [2]: '''//
     '    '/''   (1) Square of max principle stress difference'''//
     '    '/''   (2) Sum of squared reaction differences'''//
     '    '/''   (3) Zero flux differences'''//
     '    '/''   (4) Hydrostatic pressure condition'''//
     '    '/''   (5) Geometric least squares'''//
     '    '/''   (6) Reaction differences - unknown geometry'''//
     '    '/''   (7) Unused'''//
     '    '/''   (8) Unused'''//
     '    '/''   (9) Unused'''//
     '      '/$,''    '',I1)'
      ELSE IF(KTYP26.EQ.2) THEN !optimising geometric parameters
        FORMAT='('' Specify type of objective function [2]: '''//
     '    '/''   (1) Squared area of trapezoids'''//
     '    '/''   (2) Curvature of stripes'''//
     '    '/''   (3) Zero flux differences'''//
     '    '/''   (4) Hydrostatic pressure condition'''//
     '    '/''   (5) Data Fitting'''//
     '    '/''   (6) Fluid interface condition'''//
     '    '/''   (7) Aerofoil wake (& sail stress)'''//
     '    '/''   (8) Aerofoil lift and wake'''//
     '    '/''   (9) Boundary layer thickness condition'''//
     '    '/''  (10) Customisation parameters'''//
     '    '/''  (11) Difference in second moments'''//
     '    '/''  (12) Activation times'''//
     '    '/''  (13) Dipole source'''//
     '      '/$,''   '',I2)'

      ELSE IF(KTYP26.EQ.3) THEN !optimising micro-structparameters
        FORMAT='('' Specify type of objective function [1]: '''//
     '    '/''   (1) Minimise fibre stress, t11'''//
     '    '/''   (2) Minimise fibre stress and '
     '    //'transmural fibre stress gradient'''//
     '    '/''   (3) Minimise transmural fibre '
     '    //'stress gradient'''//
     '    '/''   (4) Maximise inter-sheet shear strain, e23'''//
     '    '/''   (5) Minimise sheet stress, t22'''//
     '      '/$,''   '',I2)'
      ELSE IF(KTYP26.EQ.4) THEN !optimising Holmes constants
        FORMAT='('' Specify parameters to optimise [1]: '''//
     '    '/''   (1) C1 '''//
     '    '/''   (2) C2 '''//
     '    '/''   (3) C1 and C2 '''//
     '    '/''   (4) C3 '''//
     '    '/''   (5) C1, C2 and C3'''//
     '      '/$,''   '',I2)'
      ENDIF
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP27
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,13,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP27=IDATA(1)

C LKC 20-DEC-1999 adding output timing
      FORMAT='('' Specify output from optimisation [1]: '''//
     '  '/''   (1) None'''//
     '  '/''   (2) Total Solution Time'''//
     '  '/''  *(3) Unused'''//
     '  '/$,''   '',I2)'
      IDEFLT(1)=1
C XSL 23Aug2010 IWRIT5 is defined in iwrit00.cmn as IWRIT5(nr,nx)
C But for some cases here nr, nx_opt were not initialised
C Hence rename IWRIT5 to IWRIT5_OPTI
      IF(IOTYPE.EQ.3) IDATA(1)=IWRIT5_OPTI
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) IWRIT5_OPTI=IDATA(1)


      IF(KTYP26.EQ.1) THEN !optimising material parameters
        IF(IOTYPE.EQ.3) IDATA(1)=NMNO(1,0,nx_opt)
        FORMAT=
     '    '($,'' Enter number of material params in fit [1]: '',I2)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NOPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NMNO(1,0,nx_opt)=IDATA(1)
C AJP Fixing material fitting. 17/3/96
        CHAR1=' '
        DO nm=1,NMNO(1,0,nx_opt)
          IDEFLT(nm)=nm
          IBEG=2*nm-1
          IEND=2*nm
          WRITE(CHAR1(IBEG:IEND),'(I2)') IDEFLT(nm)
        ENDDO
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
cC PJH 14Nov95 If index string>30 use 'default'
c        IF(IEND-IBEG.GT.30) THEN
c          CHAR1='default'
c          IBEG=1
c          IEND=7
c        ENDIF

        FORMAT='($,'' Enter list of material params in fit ['
     '    //CHAR1(IBEG:IEND)//']: '',20I3)'
        IF(IOTYPE.EQ.3) THEN
          DO nmlist=1,NMNO(1,0,nx_opt)
            IDATA(nmlist)=NMNO(1,nmlist,nx_opt) !needs checking AJP
          ENDDO
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    NMNO(1,0,nx_opt),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,1,NMM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nmlist=1,NMNO(1,0,nx_opt)
            NMNO(1,nmlist,nx_opt)=IDATA(nmlist) !needs checking AJP
          ENDDO
        ENDIF

C Set up NMNO and NONM in GLOBALO. AJP 17/3/96
cC PJH 9MAR96 allow for spatially varying material params
c          NMNO(1,0)=0
c          DO n=1,NTOPTI !list of material params in fit
c            il=IDATA(n) !is material parameter
c            IF(ILP(il,1,nr,nx_sol).EQ.1) THEN !material param constant
c              NMNO(1,NMNO(0)+n)=IDATA(n)
c              NMNO(1,0)=NMNO(1,0)+1
c            ELSE IF(ILP(il,1,nr,nx_sol).EQ.2) THEN !piecewise constant
c              DO noelem=1,NEELEM(0,nr)
c                NMNO(1,NMNO(0)+noelem)=IDATA(n)
c              ENDDO
c              NMNO(1,0)=NMNO(1,0)+NEELEM(0,nr)
c            ELSE IF(ILP(il,1,nr,nx_sol).EQ.3) THEN !piecewise linear
c              DO nonode=1,NPNODE(0,nr)
c                NMNO(1,NMNO(0)+nonode)=IDATA(n)
c              ENDDO
c              NMNO(1,0)=NMNO(1,0)+NPNODE(0,nr)
c            ELSE IF(ILP(il,1,nr,nx_sol).EQ.4) THEN !Gauss points
cC             needs completing
c            ENDIF !ilp
c          ENDDO !n
c        ENDIF !iotype

        IF((KTYP27.NE.5).AND.(KTYP27.NE.6)) THEN
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP28
          FORMAT='('' Enter type of optimisation [1]:'''
     '      //'/''   (0) Calculate residuals during fit  '''
     '      //'/''   (1) Use n existing sets of residuals'''
     '      //'/$,''    '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN 
            KTYP28=IDATA(1)
          ELSE
            KTYP28=0
          ENDIF
        ENDIF

        RDEFLT(1)= 1.0D0
        RDEFLT(2)= 0.0D0
        RDEFLT(3)=10.0D0
        DO nmlist=1,NMNO(1,0,nx_opt)
          WRITE(CHAR,'(I2)') NMNO(1,nmlist,nx_opt)
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=PAOPTI(nmlist)
            RDATA(2)=PMIN(nmlist)
            RDATA(3)=PMAX(nmlist)
          ENDIF
          FORMAT='($,'' Enter init,min,max vals for param #'
     '     //CHAR(1:2)//' [1 0 10]:'',3E13.6)'
          !The question is asked only once for each nm.  If nodally
          !based fitting is used the defaults will be applied to
          !each nm at each node.
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            PAOPTI(nmlist)=RDATA(1)
            PMIN(nmlist)  =RDATA(2)
            PMAX(nmlist)  =RDATA(3)
          ENDIF
        ENDDO
        
        IF((KTYP27.EQ.5).OR.(KTYP27.EQ.6)) THEN
          KTYP29B=1
          FORMAT='($,'' Enter the command file name [OPTIMISE]: '',A)'
          CDEFLT(1)='OPTIMISE'
          IF(IOTYPE.EQ.3) CDATA(1)=COM_FILE
          CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) COM_FILE=CDATA(1)(1:255)
        ELSE
          DO K=2,KTYP28  !to read sets of residuals beyond the first set
            IY=3+2*K
            IIY=0
            IF(IOTYPE.EQ.1.OR.IOTYPE.EQ.3) THEN
              FORMAT='('' Experiment number '',I2)'
              INQUIRE(UNIT=IFILE,OPENED=OPENED,NEXTREC=IREC)
              WRITE(UNIT=IFILE,FMT=FORMAT,REC=IREC,IOSTAT=IOSTAT) IY
              FORMAT='(1X,10E13.5)'
              INQUIRE(UNIT=IFILE,OPENED=OPENED,NEXTREC=IREC)
              WRITE(UNIT=IFILE,FMT=FORMAT,REC=IREC,IOSTAT=IOSTAT)
     '          (YP(ny,IY,nx_opt),NY=1,NYT(2,1,nx_opt))
              INQUIRE(UNIT=IFILE,OPENED=OPENED,NEXTREC=IREC)
              WRITE(UNIT=IFILE,FMT=FORMAT,REC=IREC,IOSTAT=IOSTAT)
     '          (YP(ny,IY+1,nx_opt),NY=1,NYT(2,1,nx_opt))
            ELSE IF(IOTYPE.EQ.2) THEN
              FORMAT='('' Experiment number '',I2)'
              INQUIRE(UNIT=IFILE,OPENED=OPENED,NEXTREC=IREC)
              READ(UNIT=IFILE,FMT=FORMAT,REC=IREC,IOSTAT=IOSTAT) IIY
              FORMAT='(1X,10E13.5)'
              INQUIRE(UNIT=IFILE,OPENED=OPENED,NEXTREC=IREC)
              READ(UNIT=IFILE,FMT=FORMAT,REC=IREC,IOSTAT=IOSTAT)
     '          (YP(ny,IY,nx_opt),NY=1,NYT(2,1,nx_opt))
              INQUIRE(UNIT=IFILE,OPENED=OPENED,NEXTREC=IREC)
              READ(UNIT=IFILE,FMT=FORMAT,REC=IREC,IOSTAT=IOSTAT)
     '          (YP(ny,IY+1,nx_opt),NY=1,NYT(2,1,nx_opt))
            ENDIF
          ENDDO !k
        ENDIF
C ***   Set up mapping betw optimising params and physical params

        CALL GLOBALO(IDO,INP,ISIZE_MFI,ISIZE_PHI,
     '    LD,LDR,NBH,NBJ,NEELEM,NENP,NKB,NKHE,NKH,
     '    NLNO(1,nx_opt),NMNO(1,0,nx_opt),NNB,
     '    NONL(1,nx_opt),NONM(1,1,nx_opt),NONY(0,1,1,nr,nx_opt),
     '    NPL,NPLIST4,NPNE,NPNODE,NPNY(0,1,0,nx_opt),nr,NVHE,NVHP,
     '    nx_opt,nx_sol,NXI,NYNE,NYNO(0,1,1,nr,nx_opt),NYNP,NYNR,
     '    NYNY(0,1,nr,nx_opt),PAOPTY,CONY(0,1,1,nr,nx_opt),
     '    CYNO(0,1,1,nr,nx_opt),CYNY(0,1,nr,nx_opt),PAOPTI,PMAX,PMIN,XP,
     '    FIX,ERROR,*9999)


        IF (KTYP27.EQ.5) THEN
          CALL ASSERT(nx_sol.GT.0,'>>no nx defined for solving',
     '      ERROR,*9999)
          NT_RES=NDT*3
C Should really be NJ_LOC(NJL_GEOM,0,nr) not sure how to get it in here though.
        ELSE IF (KTYP27.NE.6) THEN
          CALL ASSERT(nx_sol.GT.0,'>>no nx defined for solving',
     '      ERROR,*9999)
          NT_RES=NOT(1,1,nr,nx_sol)
        ELSE
          IDEFLT(1)=1
          FORMAT='($,'' Enter the number of residuals [1]: '',I4)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=NT_RES
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NREM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NT_RES=IDATA(1)
        ENDIF

      ELSE IF(KTYP26.EQ.2) THEN !optimising geometric parameters

C AJP 13-2 98 What is this??? nr has been passed through
C        nr=1 !needs fixing

        IF(KTYP27.EQ.1) THEN !Squared area of trapezoids
          IF(IOTYPE.EQ.3) IDATA(1)=NMNO(1,0,nx_opt)
          FORMAT='($,'' Enter number of geometric parameters in fit '
     '      //'[1]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,20,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NMNO(1,0,nx_opt)=IDATA(1)
C AJP 17/3/96  NTOPTI=NMNO(1,0,nx_opt)

          DO nm=1,NMNO(1,0,nx_opt)
            NMNO(1,nm,nx_opt)=nm
            WRITE(CHAR,'(I2)') NMNO(1,nm,nx_opt)
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(nm)
              RDATA(2)=PMIN(nm)
              RDATA(3)=PMAX(nm)
            ENDIF
            FORMAT=
     '        '($,'' Enter initial, minimum & maximum values for '
     '        //'parameter number '//CHAR(1:2)//': '',3E11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RZERO(1),-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(nm)=RDATA(1)
              PMIN(nm)  =RDATA(2)
              PMAX(nm)  =RDATA(3)
            ENDIF
          ENDDO

          CALL GLOBALO(IDO,INP,ISIZE_MFI,ISIZE_PHI,
     '      LD,LDR,NBH,NBJ,NEELEM,NENP,NKB,NKHE,NKH,
     '      NLNO(1,nx_opt),NMNO(1,0,nx_opt),NNB,
     '      NONL(1,nx_opt),NONM(1,1,nx_opt),NONY(0,1,1,nr,nx_opt),
     '      NPL,NPLIST4,NPNE,NPNODE,NPNY(0,1,0,nx_opt),nr,NVHE,NVHP,
     '      nx_opt,nx_sol,NXI,NYNE,NYNO(0,1,1,nr,nx_opt),NYNP,NYNR,
     '      NYNY(0,1,nr,nx_opt),PAOPTY,CONY(0,1,1,nr,nx_opt),
     '      CYNO(0,1,1,nr,nx_opt),CYNY(0,1,nr,nx_opt),PAOPTI,PMAX,PMIN,
     '      XP,FIX,ERROR,*9999)

        ELSE IF(KTYP27.EQ.2) THEN !Curvature of stripes

        ELSE IF(KTYP27.EQ.3) THEN !Zero flux differences
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP100
          FORMAT='('' Enter type of problem [1]:'''
     '      //'/''   (1) Flow from cavities   '''
     '      //'/''   (2) Flow around cavities '''
     '      //'/$,''    '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP100=IDATA(1)
          IF(KTYP100.EQ.2) THEN
            IF(IOTYPE.EQ.3) RDATA(1)=THETA_SAT
            RDEFLT(1)=1.0d0
            FORMAT='($,'' Enter value of theta_sat: '',E11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THETA_SAT=RDATA(1)
          ENDIF
C ***     Assume that there are two regions specified.  Locate all
C ***     nodes on the interface of the regions and store in NPJOIN.
C ***     Thus NPJOIN(I,1),i=1,..,NPJOIN(0,1) are the nodes on the
C ***     interface
C ***     (EXCLUDING those nodes which also lie on the fixed surface).
C ***     NPJOIN is a two dimensional array to allow for the extension
C ***     to more than one interface.
C ***     Locate all nodes that belong to both regions.
          IF(IOTYPE.EQ.3) IDATA(1)=NPSHARE(0)
          FORMAT='($,'' Enter the number of nodes lying on both the '
     '      //'fixed and free surfaces [0]:'',I1)'
          IDEFLT(1)=0
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NPSHARE(0)=IDATA(1)
          IF(NPSHARE(0).GT.0) THEN
            IF(IOTYPE.EQ.3) THEN
              DO nps=1,NPSHARE(0)
                IDATA(nps)=NPSHARE(nps)
              ENDDO
            ENDIF
            WRITE(CHAR1,'(I1)') NPSHARE(0)
            FORMAT='(/$,'' Enter the '//CHAR1(1:1)//' node numbers: '','
     '        //'10I3)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NPSHARE(0),
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,0,NPM,
     '        LDATA,LDEFLT,RDATA,RZERO(1),RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nps=1,NPSHARE(0)
                NPSHARE(nps)=IDATA(nps)
              ENDDO
            ENDIF
          ENDIF
          nr1=1 !Assumed FE region
          nr2=2 !Assumed BE region
          NPJOIN(0,1)=0
          DO nonode1=1,NPNODE(0,nr1) !First region
            NP1=NPNODE(nonode1,nr1)
            DO NONODE2=1,NPNODE(0,nr2)
              NP2=NPNODE(NONODE2,nr2)
              IF(NP1.EQ.NP2) THEN
                SHARE=.FALSE.
                DO nps=1,NPSHARE(0)
                  IF(NP1.EQ.NPSHARE(nps))SHARE=.TRUE.
                ENDDO
                IF(.NOT.SHARE) THEN
                  NPJOIN(0,1)=NPJOIN(0,1)+1
                  NPJOIN(NPJOIN(0,1),1)=NP1
                  GOTO 1234
                ENDIF
              ENDIF
            ENDDO
1234      ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,*)' Number of nodes on the interface=',
     '        NPJOIN(0,1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)' Node numbers :'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nonode1=1,NPJOIN(0,1)
              WRITE(OP_STRING,*)NPJOIN(nonode1,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

          NTOPTI=NPJOIN(0,1)
          DO noopti=1,NTOPTI
            PAOPTI(noopti)=1.0d0  !is scaling on coordinates
            PBOPTI(noopti)=1.0d0  !is scaling on coordinates (saved)
          ENDDO
C ***     Need to locate all nodes not on a fixed surface which are
C ***     on the same radial line as each free surface node.
C ***     This is stored in NPRAD:
C ***     NPRAD(I,noopti,nr), I=1,..,NPRAD(0,noopti,nr) are the nodes
C ***     on the (approximate) same radial line as the current (noopti)
C ***     free surface node which are in the FE region and not on any
C ***     fixed surface. The node on the fixed surface is identified by
C ***     NPRAD(-1,noopti,nr).  This node number is needed in the mesh
C ***     updating.
C ***     Due to the presence of fixed surfaces the user is prompted
C ***     for a list of nodes on the same radial line
          DO nonode=1,NPJOIN(0,1)
            WRITE(CHAR4,'(I4)') NPJOIN(nonode,1)
            FORMAT='('' For node number '//CHAR4(1:4)//'''/,'' Enter '
     '        //'the number of nodes on the same radial line''/$,'
     '        //''' (exclude nodes on fixed surfaces and interfaces:) '
     '        //''',I4)'
            IF(IOTYPE.EQ.3)IDATA(1)=NPRAD(0,nonode,1)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     &        NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3)NPRAD(0,nonode,1)=IDATA(1)
            CALL ASSERT(IDATA(1).LE.10,'>>NPRAD array not big enough',
     '        ERROR,*9999)
            FORMAT='($,''  Enter FIRST the node number on the fixed '
     '        //'surface''/,''    and then the other node numbers on '
     '        //'the same radial line '',10I4)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=NPRAD(-1,nonode,nr1)
              DO npr=1,NPRAD(0,nonode,nr1)
                IDATA(npr+1)=NPRAD(npr,nonode,nr1)
              ENDDO
            ENDIF
            NDATA=NPRAD(0,nonode,nr1)+1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NDATA,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     &        1,NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              NPRAD(-1,nonode,nr1)=IDATA(1)
              DO npr=1,NPRAD(0,nonode,nr1)
                NPRAD(npr,nonode,nr1)=IDATA(npr+1)
              ENDDO
            ENDIF
          ENDDO !End of loop over interface nodes
          IF(DOP) THEN
            DO nonode1=1,NPJOIN(0,1)
              WRITE(OP_STRING,*)' Interface node =',NPJOIN(nonode1,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,*)' Corresp. fixed surface node =',
     '          NPRAD(-1,nonode1,nr1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,*)'   Node nums on same radial line :'
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO npr=1,NPRAD(0,nonode1,nr1)
                WRITE(OP_STRING,*)NPRAD(npr,nonode1,nr1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDDO
          ENDIF

C          NTOPTI=0
C 100      IF(IOTYPE.EQ.3) THEN
C            IDATA(1)=NP1OPT(NTOPTI+1)
C            IDATA(2)=NP2OPT(NTOPTI+1)
C          ENDIF
C          FORMAT='($,'' Enter fem node # & corresponding bem node #'//
C     '      ' [exit]: '',2I4)'
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NPT(1),
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IDATA(1).GT.0) THEN
C            NTOPTI=NTOPTI+1
C            IF(IOTYPE.NE.3) THEN
C              NP1=IDATA(1)
C              NP2=IDATA(2)
C              NP1OPT(NTOPTI)=NP1    !is finite element node number
C              NP2OPT(NTOPTI)=NP2    !is boundary element node number
C              PAOPTI(NTOPTI)=1.0D0  !is scaling on coordinates
C              PBOPTI(NTOPTI)=1.0D0  !is scaling on coordinates (saved)
C            ENDIF
C            GO TO 100
C          ENDIF
C          IF(IOTYPE.EQ.3) THEN
C            IDATA(1)=NRAELE
C            IDATA(2)=NTHELE
C          ENDIF
C          FORMAT='($,'' Enter number of elements in radial and theta'//
C     '      ' directions: '',2I4)'
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NPT(1),
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            NRAELE=IDATA(1)
C            NTHELE=IDATA(2)
C          ENDIF
C
        ELSE IF(KTYP27.EQ.4) THEN !Hydrostatic pressure condition
          NTOPTI=0
 200      IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NP1OPT(NTOPTI+1)
          ENDIF
          FORMAT='($,'' Enter fem node number on free boundary'//
     '      ' [exit]: '',I4)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NPT(1),
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).GT.0) THEN
            NTOPTI=NTOPTI+1
            IF(IOTYPE.NE.3) THEN
              NP1=IDATA(1)
              NP1OPT(NTOPTI)=NP1  !is f.e. node number of free bdry
            ENDIF
            WRITE(CHAR,'(I2)') NTOPTI
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(NTOPTI)
              RDATA(2)=PMIN(NTOPTI)
              RDATA(3)=PMAX(NTOPTI)
            ENDIF
            FORMAT=
     '     '($,'' Enter initial, minimum & maximum values for '
     '        //'parameter number '//CHAR(1:2)//': '',3E11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(NTOPTI)=RDATA(1)  !is scaling on vertical coord
              PBOPTI(NTOPTI)=1.0D0     !is scaling on vert. coords (saved)
              PMIN(NTOPTI)  =RDATA(2)
              PMAX(NTOPTI)  =RDATA(3)
            ENDIF
            GO TO 200
          ENDIF

        ELSE IF(KTYP27.EQ.5) THEN !Data fitting
          nr=1 !This may need fixing
          CALL ASSERT(CALL_FIT,'>>Define fit first',ERROR,*9999)
          CALL ASSERT(USE_NONLIN.EQ.1,'>>Define fit first',ERROR,*9999)

          IF(KTYP29.EQ.1) THEN
            FORMAT='('' Specify residual(s) in optimisation [1]:'''//
     '        '/''   (1) Components of data projections'''//
     '        '/''   (2) Euclidean norm of data projection magnitudes'''
     '        //'/''   (3) Squares of data projection magnitudes'''
     '        //'/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP29B
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP29B=IDATA(1)
          ENDIF
          IF(KTYP29.EQ.2.OR.KTYP29B.EQ.3) THEN
            FORMAT=
     '        '($,'' Are lines to be included in the fit [N]?:'',A)'
            IF(IOTYPE.EQ.3) THEN
              IF(G1SCALING) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            ADEFLT(1)='N'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
                G1SCALING=.TRUE.
              ELSE
                G1SCALING=.FALSE.
              ENDIF
            ENDIF
          ELSE
            G1SCALING=.FALSE.
          ENDIF

          DO nj=1,NJT
            RDEFLT(1)=25.0d0
            RDEFLT(2)=25.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=NPMIN(nj)
              RDATA(2)=NPMAX(nj)
            ENDIF
            WRITE(CHAR1,'(I1)') nj
            FORMAT='($,'' Enter min & max change limit for nodes for '
     '        //'direction '//CHAR1(1:1)
     '        //' [25.0,25.0]: '',2(1X,D11.3))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              NPMIN(nj)=RDATA(1)
              NPMAX(nj)=RDATA(2)
            ENDIF
          ENDDO
C GMH 18/3/96 Check for derivative as angle
          DERIV=.FALSE.
          DO nb=1,NBFT
            DO nn=1,NNT(nb)
              IF(NKT(nn,nb).GE.2) DERIV=.TRUE.
            ENDDO
          ENDDO
          IF(DERIV) THEN
            FORMAT='('' Treat derivatives as [1]: '''//
     '        '/''   (1) Unit length, optimise components'''//
     '        '/''   (2) Unit length, optimise angles'''//
     '        '/''   (3) Non-unit length, optimise components'''//
     '        '/$,''   '',I2)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP1B
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP1B=IDATA(1)
            IF(KTYP1B.EQ.3) THEN !non-unit length first derivs
              RDEFLT(1)=-1.0d0
              RDEFLT(2)=1.0d0
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=FDMIN
                RDATA(2)=FDMAX
              ENDIF
              FORMAT='($,'' Enter min & max limits for first'
     '          //' derivatives [-1.0,1.0]: '',2D11.3)'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                FDMIN=RDATA(1)
                FDMAX=RDATA(2)
              ENDIF
            ELSE
              FDMIN=-1.0d0
              FDMAX=1.0d0
            ENDIF !
          ELSE !No deriv
            KTYP1B=1
          ENDIF
          CROSSDERIV=.FALSE.
          DO nb=1,NBFT
            DO nn=1,NNT(nb)
              IF(NKT(nn,nb).GE.4) CROSSDERIV=.TRUE.
            ENDDO
          ENDDO
          IF(CROSSDERIV) THEN
            RDEFLT(1)=-1.0d-1
            RDEFLT(2)=1.0d-1
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=CDMIN
              RDATA(2)=CDMAX
            ENDIF
            FORMAT='($,'' Enter min & max limits for cross'
     '         //' derivatives [-0.1,0.1]: '',2D11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              CDMIN=RDATA(1)
              CDMAX=RDATA(2)
            ENDIF
          ENDIF
          IF(G1SCALING) THEN
            RDEFLT(1)=0.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=LINEMAX
            ELSE
              DO nl=1,NLT
                IF(DL(3,nl).GT.RDEFLT(1)) RDEFLT(1)=DL(3,nl)
              ENDDO
              RDEFLT(1)=3.0d0*RDEFLT(1)
            ENDIF
            WRITE(CHAR5,'(F6.1)') RDEFLT(1)
            FORMAT='($,'' Enter maximum limit for lines ['//CHAR5//
     '        ']: '',D11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              LINEMAX=RDATA(1)
            ENDIF
          ENDIF

C ***     Set up mapping betw optimising params and physical params
          CALL GLOBALO(IDO,INP,ISIZE_MFI,ISIZE_PHI,
     '      LD,LDR,NBH,NBJ,NEELEM,NENP,NKB,NKHE,NKH,
     '      NLNO(1,nx_opt),NMNO(1,0,nx_opt),NNB,
     '      NONL(1,nx_opt),NONM(1,1,nx_opt),NONY(0,1,1,nr,nx_opt),
     '      NPL,NPLIST4,NPNE,NPNODE,NPNY(0,1,0,nx_opt),nr,NVHE,NVHP,
     '      nx_opt,nx_sol,NXI,NYNE,NYNO(0,1,1,nr,nx_opt),NYNP,NYNR,
     '      NYNY(0,1,nr,nx_opt),PAOPTY,CONY(0,1,1,nr,nx_opt),
     '      CYNO(0,1,1,nr,nx_opt),CYNY(0,1,nr,nx_opt),PAOPTI,PMAX,PMIN,
     '      XP,FIX,ERROR,*9999)

C ***     Set the initial parameter values and their bounds
          DO noopti=1,NTOPTI
            IF(PAOPTY(noopti).EQ.1) THEN
              ny=NYNO(1,noopti,2,nr,nx_opt)
              nk=NPNY(1,ny,0,nx_opt)
              nv=NPNY(2,ny,0,nx_opt)
              nj=NPNY(3,ny,0,nx_opt)
              np=NPNY(4,ny,0,nx_opt)
              IF(nk.EQ.1) THEN
                PAOPTI(noopti)=XP(nk,nv,nj,np)
                PMIN(noopti)=PAOPTI(noopti)-NPMIN(nj)
                PMAX(noopti)=PAOPTI(noopti)+NPMAX(nj)
              ELSE IF(nk.EQ.2.OR.nk.EQ.3) THEN
                IF(KTYP1B.EQ.2) THEN !derivative as angle
                  IF(NYNO(0,noopti,2,nr,nx_opt).EQ.3) THEN
                    ny2=NYNO(2,noopti,2,nr,nx_opt)
                    ny3=NYNO(3,noopti,2,nr,nx_opt)
                    nk2=NPNY(1,ny2,0,nx_opt)
                    nk3=NPNY(1,ny3,0,nx_opt)
                    nv2=NPNY(2,ny2,0,nx_opt)
                    nv3=NPNY(2,ny3,0,nx_opt)
                    nj2=NPNY(3,ny2,0,nx_opt)
                    nj3=NPNY(3,ny3,0,nx_opt)
                    np2=NPNY(4,ny2,0,nx_opt)
                    np3=NPNY(4,ny3,0,nx_opt)
C                   Which angle are we
                    IF(NONY(1,ny,2,nr,nx_opt).EQ.noopti) THEN !theta
                      PAOPTI(noopti)=DATAN_MOD(XP(nk,nv,nj,np),
     '                  XP(nk2,nv2,nj2,np2))
                      PMIN(noopti)=PAOPTI(noopti)-0.5d0*PI
                      PMAX(noopti)=PAOPTI(noopti)+0.5d0*PI
                    ELSE IF(NONY(2,ny,2,nr,nx_opt).EQ.noopti) THEN !phi
                      RAD=DSQRT(XP(nk,nv,nj,np)*XP(nk,nv,nj,np)+
     '                  XP(nk2,nv2,nj2,np2)*XP(nk2,nv2,nj2,np2))
                      PAOPTI(noopti)=DATAN_MOD(XP(nk3,nv3,nj3,np3),
     '                  RAD)
                      PMIN(noopti)=0.0d0
                      PMAX(noopti)=PI
                    ELSE
                      ERROR='>>Invalid angle'
                      GOTO 9999
                    ENDIF
                  ELSEIF(NYNO(0,noopti,2,nr,nx_opt).EQ.2) THEN
                    ny2=NYNO(2,noopti,2,nr,nx_opt)
                    nk2=NPNY(1,ny2,0,nx_opt)
                    nv2=NPNY(2,ny2,0,nx_opt)
                    nj2=NPNY(3,ny2,0,nx_opt)
                    np2=NPNY(4,ny2,0,nx_opt)
C                   We must be theta
                    PAOPTI(noopti)=DATAN_MOD(XP(nk,nv,nj,np),
     '                XP(nk2,nv2,nj2,np2))
                    PMIN(noopti)=PAOPTI(noopti)-0.5d0*PI
                    PMAX(noopti)=PAOPTI(noopti)+0.5d0*PI
                  ELSE
                    ERROR='>>Angle must correspond '
     '                //'to at least 2 components'
                    GOTO 9999
                  ENDIF
C                 Calculate the angle
                ELSE
                  PAOPTI(noopti)=XP(nk,nv,nj,np)
                  PMIN(noopti)=FDMIN
                  PMAX(noopti)=FDMAX
                ENDIF
              ELSE IF(nk.EQ.4) THEN
                PAOPTI(noopti)=XP(nk,nv,nj,np)
                PMIN(noopti)=CDMIN
                PMAX(noopti)=CDMAX
              ENDIF
            ELSE IF(PAOPTY(noopti).EQ.2) THEN
              PAOPTI(noopti)=DL(1,NLNO(noopti,nx_opt))
              PMIN(noopti)=0.0d0
              PMAX(noopti)=LINEMAX
            ELSE
              ERROR='>> Invalid PAOPTY for data fitting'
              GOTO 9999
            ENDIF
          ENDDO !noopti

        ELSE IF(KTYP27.EQ.6) THEN !Fluid interface condition
          IF(IOTYPE.EQ.3) IDATA(1)=NPJOIN(0,1)
          FORMAT='($,'' Enter number of interface nodes in fit [1]: '','
     '      //'I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,50,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NPJOIN(0,1)=IDATA(1)

          IF(IOTYPE.EQ.3) THEN
            DO nojoin=1,NPJOIN(0,1)
              IDATA(nojoin)=NPJOIN(nojoin,1)
            ENDDO
          ENDIF
          FORMAT='($,'' Enter list of interface nodes in fit [1]: '','
     '      //'50I3)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      NPJOIN(0,1),
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,50,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nojoin=1,NPJOIN(0,1)
              NPJOIN(nojoin,1)=IDATA(nojoin)
            ENDDO
          ENDIF

          NTOPTI=2*NPJOIN(0,1)
          NMNO(1,0,nx_opt)=NTOPTI
C!!!!Should be done in GLOBALO
          DO noopti=1,NMNO(1,0,nx_opt)
            NMNO(1,noopti,nx_opt)=noopti
          ENDDO
          DO noopti=1,NTOPTI
            WRITE(CHAR,'(I2)') noopti
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(noopti)
              RDATA(2)=PMIN(noopti)
              RDATA(3)=PMAX(noopti)
            ENDIF
            FORMAT=
     '     '($,'' Enter initial, minimum & maximum values for '
     '        //'parameter number '//CHAR(1:2)//': '',3E11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RZERO(1),-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(noopti)=RDATA(1)
              PBOPTI(noopti)=RDATA(1)
              PMIN(noopti)  =RDATA(2)
              PMAX(noopti)  =RDATA(3)
            ENDIF
          ENDDO

        ELSE IF(KTYP27.EQ.7) THEN !Aerofoil wake dPHI residuals
!         & (if sail stress eqtns defined) sail stress residuals
          IDEFLT(1)=NL_WAKE(0,1)
          WRITE(CHAR,'(I2)') IDEFLT(1)
          FORMAT='($,'' Enter #wake parameters ['//CHAR(1:2)//']: '','
     '      //'I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=N_OPTI(1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,50,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) N_OPTI(1)=IDATA(1)

! wake node positions
          IF(N_OPTI(1).GE.NL_WAKE(0,1)) THEN
!           Store node#s of wake nodes excluding trailing edge
            DO no_wake=1,NL_WAKE(0,1)
              NP1OPT(no_wake)=NPL(3,1,NL_WAKE(no_wake,1)) !upper surface
              NP2OPT(no_wake)=NPL(3,1,NL_WAKE(no_wake,2)) !lower surface
            ENDDO
            NTOPTI=NL_WAKE(0,1)
          ELSE
            NTOPTI=0 !no wake parameters
          ENDIF

! trailing edge velocity constraint (nolonger needed for 2D)
c         IF(N_OPTI(1).EQ.NL_WAKE(0,1)+1) THEN !constraint applied
c           FORMAT='('' Enter constraint type [0]: '''//
c    '        '/''   (0) none'''//
c    '        '/''   (1) equal trailing edge velocities'''//
c    '        '/''   (2) '''//
c    '        '/$,''    '',I1)'
c           IF(IOTYPE.EQ.3) IDATA(1)=CONSTRAINT_TYPE
c           CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
c    '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,1,
c    '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c           IF(IOTYPE.NE.3) THEN
c             CONSTRAINT_TYPE=IDATA(1)
c             PMIN(NTOPTI+1)=0.0d0
c             PMAX(NTOPTI+1)=0.0d0
c           ENDIF
c         ENDIF
          CONSTRAINT_TYPE=0

          IF(N_OPTI(1).GE.NL_WAKE(0,1)) THEN
            IF(IOTYPE.EQ.1) THEN !prompted input
              WRITE(OP_STRING,'(/'' Wake node pairs are:'','
     '          //'10(I5,1X,I5,5X))')
     '          (NP1OPT(no_wake),NP2OPT(no_wake),no_wake=1,NL_WAKE(0,1))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF

          DO noopti=1,NTOPTI
            np1=NP1OPT(noopti)
            PBOPTI(noopti)=1.0d0 !not used
            WRITE(CHAR,'(I2)') noopti
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(noopti)
              RDATA(2)=PMIN(noopti)
              RDATA(3)=PMAX(noopti)
            ENDIF
c           IF(noopti.LT.NTOPTI.OR.N_OPTI(1).EQ.NL_WAKE(0,1)) THEN
              RDEFLT(1)=XP(1,1,2,np1)
              RDEFLT(2)=0.0d0
              RDEFLT(3)=2.0d0*XP(1,1,2,np1)
              FORMAT='($,'' Enter initial, min & max values'
     '          //' for node pair #'//CHAR(1:2)//' [xp,0,2*xp]: '','
     '          //'3E11.3)'
c           ELSE IF(N_OPTI(1).GE.NL_WAKE(0,1)+1) THEN
c             RDEFLT(1)= 0.0d0
c             RDEFLT(2)=-1.0d0
c             RDEFLT(3)= 1.0d0
c             FORMAT='($,'' Enter initial, min & max values'
c    '          //' for deriv pair #'//CHAR(1:2)//' [0,-1,1]: '','
c    '          //'3E11.3)'
c           ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &          *9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(noopti)=RDATA(1)
              PMIN(noopti)  =RDATA(2)
              PMAX(noopti)  =RDATA(3)
            ENDIF
          ENDDO

! sail stress parameters
          nr=2 !is assumed region# for sail stress
          IF(NPNODE(0,nr).GT.0) THEN !sail stress nodes defined
            IDEFLT(1)=0
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np,nr,nx_opt)
                nh=NH_LOC(nhx,nx_opt)
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,NKH(nh,NP,1)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    IF(.NOT.FIX(ny,3,nx_opt)) THEN !Include unfixed vars
                      IDEFLT(1)=IDEFLT(1)+1
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !nonode(np)
            WRITE(CHAR,'(I2)') IDEFLT(1)
            FORMAT='($,'' Enter #sail stress parameters ['//CHAR(1:2)
     '        //']: '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=N_OPTI(2)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        50,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) N_OPTI(2)=IDATA(1)
          ELSE
            N_OPTI(2)=0
          ENDIF !npnode(0,nr)>0

          IF(N_OPTI(2).GT.0) THEN !sail stress parameters included
            IF(IOTYPE.EQ.1) THEN !prompted input
              WRITE(OP_STRING,'(/'' Sail degrees of freedom:'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              WRITE(CHAR1,'(I5)') np
              CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
              DO nhx=1,NHP(np,nr,nx_opt)
                nh=NH_LOC(nhx,nx_opt)
                WRITE(CHAR2,'(I1)') nh
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,NKH(nh,NP,1)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' np='',I5,'' nh='',I1,'
     '                  //''' nk='',I1,'' ny='',I5)') np,nh,nk,NY
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(.NOT.FIX(ny,3,nx_opt)) THEN !Include unfixd vars in fit
                      NTOPTI=NTOPTI+1
                      NP1OPT(NTOPTI)=NP
                      NP2OPT(NTOPTI)=NH
                      NP3OPT(NTOPTI)=NK
                      WRITE(CHAR3,'(I1)') nk
                      IF(IOTYPE.EQ.3) THEN
                        RDATA(1)=PAOPTI(NTOPTI)
                        RDATA(2)=PMIN(NTOPTI)
                        RDATA(3)=PMAX(NTOPTI)
                      ENDIF
                      RDEFLT(1)=       XP(nk,nv,nh,np)
                      RDEFLT(2)=-2.0D0*XP(nk,nv,nh,np)
                      RDEFLT(3)= 2.0D0*XP(nk,nv,nh,np)
                      FORMAT='($,'' Enter initial, min & max values for'
     '                   //' node '//CHAR1(IBEG1:IEND1)//' variable '
     '                   //CHAR2(1:1)//' deriv '//CHAR3(1:1)
     '                   //' [x,-2x,2x]: '',3E11.3)'
                      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                  FILEIP,FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                  IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '                  -RMAX,RMAX,INFO,ERROR,*9999)
                      IF(IOTYPE.NE.3) THEN
                        PAOPTI(NTOPTI)=RDATA(1)
                        PMIN(NTOPTI)  =RDATA(2)
                        PMAX(NTOPTI)  =RDATA(3)
                      ENDIF
                    ENDIF
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !nonode
          ENDIF !n_opti(2)>0

          NT_RES=NTOPTI !#residuals
          IF(CONSTRAINT_TYPE.GT.0) THEN
            NTCNTR=1 !#constraints
          ELSE
            NTCNTR=0
          ENDIF

          NMNO(1,0,nx_opt)=NTOPTI
          DO noopti=1,NMNO(1,0,nx_opt)
            NMNO(1,noopti,nx_opt)=noopti
          ENDDO

        ELSE IF(KTYP27.EQ.8) THEN !Aerofoil lift & wake press diff.
          IDEFLT(1)=NL_WAKE(0,1)
          WRITE(CHAR,'(I2)') IDEFLT(1)
          FORMAT='($,'' Enter #wake parameters ['//CHAR(1:2)//']: '','
     '      //'I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=N_OPTI(1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,50,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) N_OPTI(1)=IDATA(1)

          IF(N_OPTI(1).GE.NL_WAKE(0,1)) THEN
!           Store node#s of wake nodes excluding trailing edge
            DO no_wake=1,NL_WAKE(0,1)
              NP1OPT(no_wake)=NPL(3,1,NL_WAKE(no_wake,1)) !upper surface
              NP2OPT(no_wake)=NPL(3,1,NL_WAKE(no_wake,2)) !lower surface
            ENDDO
            IF(IOTYPE.EQ.1) THEN !prompted input
              WRITE(OP_STRING,'(/'' Wake node pairs are: '','
     '          //'10(2I5,3X))')
     '          (NP1OPT(no_wake),NP2OPT(no_wake),no_wake=1,NL_WAKE(0,1))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            DO noopti=1,NL_WAKE(0,1)
              np1=NP1OPT(noopti)
              PBOPTI(noopti)=1.d0 !not used
              WRITE(CHAR,'(I2)') noopti
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=PAOPTI(noopti)
                RDATA(2)=PMIN(noopti)
                RDATA(3)=PMAX(noopti)
              ENDIF
              RDEFLT(1)=XP(1,1,2,np1)
              RDEFLT(2)=-2.0D0*XP(1,1,2,np1)
              RDEFLT(3)= 2.0D0*XP(1,1,2,np1)
              FORMAT='($,'' Enter initial, min & max values'
     '          //' for node pair #'//CHAR(1:2)//' [xp,-2*xp,2*xp]: '','
     '          //'3E11.3)'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,3,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                PAOPTI(noopti)=RDATA(1)
                PMIN(noopti)  =RDATA(2)
                PMAX(noopti)  =RDATA(3)
              ENDIF
            ENDDO
            NTOPTI=NL_WAKE(0,1)
          ELSE
            NTOPTI=0 !no wake parameters
          ENDIF

!         Store node#s of aerofoil nodes excluding ends
          DO no_aero=1,NL_AERO(0,1)-1
            NP1OPT(N_OPTI(1)+no_aero)=NPL(3,1,NL_AERO(no_aero,1)) !up
            NP2OPT(N_OPTI(1)+no_aero)=NPL(3,1,NL_AERO(no_aero,2)) !lo
          ENDDO
          IF(IOTYPE.EQ.1) THEN !prompted input
            WRITE(OP_STRING,'( '' Aero node pairs are: '',10(2I5,3X))')
     '        (NP1OPT(N_OPTI(1)+no_aero),
     '         NP2OPT(N_OPTI(1)+no_aero),no_aero=1,N_OPTI(1)-1)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          NTOPTI=NTOPTI+NL_AERO(0,1)-1
          DO noopti=N_OPTI(1)+1,NTOPTI
            np1=NP1OPT(noopti)
            PBOPTI(noopti)=1.0d0 !not used
            WRITE(CHAR,'(I2)') noopti
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(noopti)
              RDATA(2)=PMIN(noopti)
              RDATA(3)=PMAX(noopti)
            ENDIF
            RDEFLT(1)=XP(1,1,2,np1)
            RDEFLT(2)=0.0d0
            RDEFLT(3)=2.0d0*XP(1,1,2,np1)
            FORMAT='($,'' Enter initial, min & max values'
     '        //' for node pair #'//CHAR(1:2)//' [xp,0,2*xp]: '','
     '        //'3E11.3)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(noopti)=RDATA(1)
              PMIN(noopti)  =RDATA(2)
              PMAX(noopti)  =RDATA(3)
            ENDIF
          ENDDO

          RDEFLT(1)=4.0d0
          WRITE(CHAR1,'(F5.2)') RDEFLT(1)
          IF(IOTYPE.EQ.3) RDATA(1)=UPPER_BOUND
          FORMAT='($,'' Enter upper bound for lift ['//CHAR1(1:5)
     '      //']: '',E12.5)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            UPPER_BOUND=RDATA(1)
          ENDIF

          FORMAT='('' Enter constraint type [0]: '''//
     '      '/''   (0) none'''//
     '      '/''   (1) wetted area'''//
     '      '/''   (2) crosssection area'''//
     '      '/''   (3) '''//
     '      '/''   (4) '''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=CONSTRAINT_TYPE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,1,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) CONSTRAINT_TYPE=IDATA(1)

          IF(CONSTRAINT_TYPE.GT.0) THEN
            RDEFLT(1)=AERO_PERIM
            WRITE(CHAR1,'(E12.5)') RDEFLT(1)
            IF(IOTYPE.EQ.3) RDATA(1)=CONSTRAINT_VALUE
            FORMAT='($,'' Enter constraint value ['//CHAR1(1:12)
     '        //']: '',E12.5)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              CONSTRAINT_VALUE=RDATA(1)
              PMIN(NTOPTI+1)=RDATA(1)
              PMAX(NTOPTI+1)=RDATA(1)
            ENDIF
          ENDIF

          NT_RES=N_OPTI(1)+1 !#residuals
          IF(CONSTRAINT_TYPE.GT.0) THEN
            NTCNTR=1 !#constraints
          ELSE
            NTCNTR=0
          ENDIF

        ELSE IF(KTYP27.EQ.10) THEN !Customisation parameters
           FORMAT='('' Enter type of customisation: '''//
     '      '/''   (1) Use Circumferential Measurements '''//
     '      '/''   (2) Use Volumes'''//
     '      '/''   (3) Use Elliptic Approximation '''//
     '      '/''   (4) Use 2-D Measurements'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=CUSTOMISATION_TYPE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) CUSTOMISATION_TYPE=IDATA(1)

          FORMAT='($,'' Enter the number of polynomial coeffs [4]: '''
     '      //',I2)'
          IDEFLT(1)=4
          IF(IOTYPE.EQ.3) IDATA(1)=NPC
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NPC=IDATA(1)

          DO noopti=1,NPC
            WRITE(CHAR1,'(I2)') noopti
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(noopti)
            ENDIF
            FORMAT='($,'' Enter the initial value for optimisation '
     '        //'variable '//CHAR1(1:2)//' [1.0]: '',D12.4))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(noopti)=RDATA(1)
            ENDIF
            RDEFLT(1)=-25.0d0
            RDEFLT(2)=25.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PMIN(noopti)
              RDATA(2)=PMAX(noopti)
            ENDIF
            FORMAT='($,'' Enter min & max bounds for optimisation '
     '        //'variable '//CHAR1(1:2)
     '        //' [-25.0,25.0]: '',2(1X,D11.3))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PMIN(noopti)=RDATA(1)
              PMAX(noopti)=RDATA(2)
            ENDIF
          ENDDO !noopti

          FORMAT='($,'' Enter the number of cosine coeffs [4]: '''
     '      //',I2)'
          IDEFLT(1)=4
          IF(IOTYPE.EQ.3) IDATA(1)=NCC
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NCC=IDATA(1)

          DO noopti=NPC+1,NPC+NCC
            WRITE(CHAR1,'(I2)') noopti
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(noopti)
            ENDIF
            FORMAT='($,'' Enter the initial value for optimisation '
     '        //'variable '//CHAR1(1:2)//' [1.0]: '',D12.4))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PAOPTI(noopti)=RDATA(1)
            ENDIF
            RDEFLT(1)=-25.0d0
            RDEFLT(2)=25.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PMIN(noopti)
              RDATA(2)=PMAX(noopti)
            ENDIF
            FORMAT='($,'' Enter min & max bounds for optimisation '
     '        //'variable '//CHAR1(1:2)
     '        //' [-25.0,25.0]: '',2(1X,D11.3))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PMIN(noopti)=RDATA(1)
              PMAX(noopti)=RDATA(2)
            ENDIF
          ENDDO !noopti
          NTOPTI=NPC+NCC

          FORMAT='('' Enter type of length scaling: '''//
     '           '/''   (1) Simple scaling '''//
     '           '/''   (2) Variable length scaling'''//
     '           '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=SCLTYPE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) SCLTYPE=IDATA(1)

          IF(SCLTYPE.EQ.1) THEN
            FORMAT='($,'' Enter the total length '
     '        //'[1.0]: '',D12.3)'
            RDEFLT(1)=1.0d0
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) RDATA(1)=TORSO_LENGTHS(0,1)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) TORSO_LENGTHS(0,1)=RDATA(1)

          ELSEIF(SCLTYPE.EQ.2) THEN
            FORMAT='($,'' Enter the number of length measurements '
     '        //'[2]: '',I2)'
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) IDATA(1)=NTL
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        99,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NTL=IDATA(1)
            DO n=1,NTL
              FORMAT='($,'' Enter height on generic model and '
     '       //'corresponding actual height [1.0,1.0]: '',D12.4,D12.4)'
              IF(IOTYPE.EQ.3) THEN
                DO ndat=1,2
                  RDATA(ndat)=TORSO_LENGTHS(ndat,n)
                ENDDO !ndat
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NPM,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                DO ndat=1,2
                 TORSO_LENGTHS(ndat,n)=RDATA(ndat)
                ENDDO !ndat
              ENDIF
            ENDDO !n
          ENDIF
          IF(CUSTOMISATION_TYPE.EQ.1) THEN
            FORMAT='($,'' Enter the number of circumferential '
     '        //'measurements [3]: '',I2)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=NT_RES
            ENDIF
            IDEFLT(1)=3
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NMEASUREMX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              NT_RES=IDATA(1)
            ENDIF
            DO nores=1,NT_RES
              WRITE(CHAR1,'(I2)') nores
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=CIRMEASURE(1,nores)
              ENDIF
              FORMAT='($,'' Enter the circumference for measurement '
     '          //CHAR1(1:2)//' [1.0]: '',D12.4))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                CIRMEASURE(1,nores)=RDATA(1)
              ENDIF
              RDEFLT(1)=DBLE(nores)/DBLE(NT_RES)
              WRITE(CHAR4,'(F4.1)') RDEFLT(1)
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=CIRMEASURE(2,nores)
              ENDIF
              FORMAT='($,'' Enter the height of measurement '
     '          //CHAR1(1:2)//' ['//CHAR4//']: '',D12.4))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                CIRMEASURE(2,nores)=RDATA(1)
              ENDIF
            ENDDO !nores
          ELSEIF(CUSTOMISATION_TYPE.EQ.2) THEN
            NT_RES=2
            FORMAT='($,'' Enter the region number for outer volume '
     '       //'[3]: '',I2)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=vreg1
            ENDIF
            IDEFLT(1)=1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NMEASUREMX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              vreg1=IDATA(1)
            ENDIF
            FORMAT='($,'' Enter the region number for inner volume '
     '       //'[3]: '',I2)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=vreg2
            ENDIF
            IDEFLT(1)=2
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NMEASUREMX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              vreg2=IDATA(1)
            ENDIF
              RDEFLT(1)=20.0d0
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=PERCENT_FAT
              ENDIF
            FORMAT='($,'' Enter the percentage volume fat for torso '
     '        //' [20.0]: '',D12.4))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              PERCENT_FAT=RDATA(1)
            ENDIF
          ELSEIF(CUSTOMISATION_TYPE.EQ.3.OR.CUSTOMISATION_TYPE.EQ.4)
     '      THEN
            FORMAT='($,'' Enter the number of 2-d measurements '
     '       //'[3]: '',I2)'
            IF(IOTYPE.EQ.3) THEN
              IDATA(1)=NT_RES
            ENDIF
            IDEFLT(1)=3
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NMEASUREMX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              NT_RES=IDATA(1)*2
            ENDIF
            DO nores=1,NT_RES/2
              WRITE(CHAR1,'(I2)') nores
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=WDMEASURE(1,nores)
              ENDIF
              FORMAT='($,'' Enter the 2d measurement in x direction '
     '          //CHAR1(1:2)//' [1.0]: '',D12.4))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                WDMEASURE(1,nores)=RDATA(1)
              ENDIF
              RDEFLT(1)=1.0d0
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=WDMEASURE(2,nores)
              ENDIF
              FORMAT='($,'' Enter the 2d measurement in y direction '
     '          //CHAR1(1:2)//' [1.0]: '',D12.4))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                WDMEASURE(2,nores)=RDATA(1)
              ENDIF
              RDEFLT(1)=DBLE(nores)/DBLE(NT_RES)
              WRITE(CHAR4,'(F4.1)') RDEFLT(1)
              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=WDMEASURE(3,nores)
              ENDIF
              FORMAT='($,'' Enter the height of the measurements '
     '          //CHAR1(1:2)//' ['//CHAR4//']: '',D12.4))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                WDMEASURE(3,nores)=RDATA(1)
              ENDIF
            ENDDO
            IF(CUSTOMISATION_TYPE.EQ.3) THEN
              NT_RES=NT_RES/2
              DO nores=1,NT_RES
                CIRMEASURE(1,nores)=2.d0*PI*SQRT(0.5d0*(
     '             (0.5d0*WDMEASURE(1,nores))**2.d0+
     '             (0.5d0*WDMEASURE(2,nores))**2d0))
                CIRMEASURE(2,nores)=WDMEASURE(3,nores)
              ENDDO
            ENDIF
          ENDIF

        ELSE IF(KTYP27.EQ.11) THEN   !Second moments

          NTOPTI=3
          DO noopti=1,NTOPTI
            WRITE(CHAR1,'(I2)') noopti
            RDEFLT(1)=1.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PAOPTI(noopti)
            ENDIF
            FORMAT='($,'' Enter initial angle of rotation about '
     '        //'axis '//CHAR1(1:2)//' [1.0]: '',D12.4))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
               PAOPTI(noopti)=RDATA(1)
            ENDIF
            RDEFLT(1)=-180.0d0
            RDEFLT(2)=180.0d0
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=PMIN(noopti)
              RDATA(2)=PMAX(noopti)
            ENDIF
            FORMAT='($,'' Enter min & max bounds for optimisation '
     '        //'variable '//CHAR1(1:2)//' [-180,180]: '',2(1X,D11.3))'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              PMIN(noopti)=RDATA(1)
              PMAX(noopti)=RDATA(2)
            ENDIF
          ENDDO !noopti

          NT_RES=NJT+2+NJT-3


        ELSE IF(KTYP27.EQ.12) THEN !activation times

          FORMAT='($,'' Enter the number of critical points [1]: '',I2)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=NPJOIN(0,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NPJOIN(0,1)=IDATA(1)
          DO ncritical=1,NPJOIN(0,1)
            WRITE(CHAR1(1:2),'(I2)') ncritical
            FORMAT='($,'' Enter the node number of critical point '
     '        //CHAR1(1:2)//': '',I5)'
            IF(IOTYPE.EQ.3) IDATA(1)=NPJOIN(ncritical,1)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NPM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NPJOIN(ncritical,1)=IDATA(1)
            FORMAT='($,'' Enter the critical time of this point'
     '        //' [ms]: '',D10.3)'
            IF(IOTYPE.EQ.3) RDATA(1)=RVALUE(ncritical)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,DBLE(NTSM),INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) RVALUE(ncritical)=RDATA(1)
            !Find ny for this np and fix it.
            ny=NYNP(1,1,1,NPJOIN(ncritical,1),0,1,TRSF_NR_FIRST)
            YP(ny,1,nx_sol)=RVALUE(ncritical)
            FIX(ny,1,nx_sol)=.TRUE.
          ENDDO !# of critical points
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=ACTN_MIN(1)
            RDATA(2)=ACTN_MAX(1)
          ENDIF
          FORMAT='($,'' Enter min & max change limit for activation '
     '      //'times (ms) [10.0,10.0]: '',2(1X,D11.3))'
          RDEFLT(1)=10.0d0
          RDEFLT(2)=10.0d0
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            ACTN_MIN(1)=RDATA(1)
            ACTN_MAX(1)=RDATA(2)
          ENDIF
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=ACTN_MIN(2)
            RDATA(2)=ACTN_MAX(2)
          ENDIF
          RDEFLT(1)=0.0d0
          RDEFLT(2)=100.0d0
          FORMAT='($,'' Enter the min & max bounds on the activation '
     '      //'interval (ms) [0.0,100.0]: '',2(1X,D11.3))'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            ACTN_MIN(2)=RDATA(1)
            ACTN_MAX(2)=RDATA(2)
          ENDIF
          FORMAT='($,'' Do you want to optimise any wavefront '
     '      //'parameters [N]? '',A)'
          ADEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF(ACTN_OPTI_WAVE_PARAM) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ACTN_OPTI_WAVE_PARAM=ADATA(1).EQ.'Y'

          ACTN_MIN(3)=0.0d0
          ACTN_MAX(3)=0.0d0
          ACTN_MIN(4)=0.0d0
          ACTN_MAX(4)=0.0d0
          IF(ACTN_OPTI_WAVE_PARAM) THEN

            IF(TRSF_ACTN_WAVE_TYPE.EQ.1) THEN !Heaviside step function

              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=ACTN_MIN(3)
                RDATA(2)=ACTN_MAX(3)
              ENDIF
              RDEFLT(1)=5.0d0
              RDEFLT(2)=5.0d0
              FORMAT='($,'' Enter the min & max bounds on the '
     '          //'transmembrane jump (mV) [5.0,5.0]: '',2(1X,D11.3))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                ACTN_MIN(3)=RDATA(1)
                ACTN_MAX(3)=RDATA(2)
              ENDIF

            ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.2) THEN !Sigmoidal function

              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=ACTN_MIN(3)
                RDATA(2)=ACTN_MAX(3)
              ENDIF
              RDEFLT(1)=5.0d0
              RDEFLT(2)=5.0d0
              FORMAT='($,'' Enter the min & max bounds on the '
     '          //'transmembrane jump (mV) [5.0,5.0]: '',2(1X,D11.3))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                ACTN_MIN(3)=RDATA(1)
                ACTN_MAX(3)=RDATA(2)
              ENDIF

              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=ACTN_MIN(4)
                RDATA(2)=ACTN_MAX(4)
              ENDIF
              RDEFLT(1)=1.0d0
              RDEFLT(2)=1.0d0
              FORMAT='($,'' Enter the min & max bounds on the '
     '          //'function window(ms) [1.0,1.0]: '',2(1X,D11.3))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                ACTN_MIN(4)=RDATA(1)
                ACTN_MAX(4)=RDATA(2)
              ENDIF

            ELSE IF(TRSF_ACTN_WAVE_TYPE.EQ.3) THEN !Arctan function

              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=ACTN_MIN(3)
                RDATA(2)=ACTN_MAX(3)
              ENDIF
              RDEFLT(1)=5.0d0
              RDEFLT(2)=5.0d0
              FORMAT='($,'' Enter the min & max bounds on the '
     '          //'transmembrane jump (mV) [5.0,5.0]: '',2(1X,D11.3))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                ACTN_MIN(3)=RDATA(1)
                ACTN_MAX(3)=RDATA(2)
              ENDIF

              IF(IOTYPE.EQ.3) THEN
                RDATA(1)=ACTN_MIN(4)
                RDATA(2)=ACTN_MAX(4)
              ENDIF
              RDEFLT(1)=1.0d0
              RDEFLT(2)=1.0d0
              FORMAT='($,'' Enter the min & max bounds on the '
     '          //'function window(ms) [1.0,1.0]: '',2(1X,D11.3))'
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                ACTN_MIN(4)=RDATA(1)
                ACTN_MAX(4)=RDATA(2)
              ENDIF

            ELSE
              ERROR='>>Invalid activation wave type'
              GOTO 9999
            ENDIF
          ENDIF

          IF(IOTYPE.NE.3) THEN
            SS_OBJ_WEIGHT=1.0d0
            CC_OBJ_WEIGHT=0.0d0
          ENDIF
          FORMAT='($,'' Enter the weighting for sum-of-squares '
     '      //'objective [1.0]: '',D11.3)'
          RDEFLT(1)=1.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=SS_OBJ_WEIGHT
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) SS_OBJ_WEIGHT=RDATA(1)
          FORMAT='($,'' Enter the weighting for correlation '
     '      //'coefficient objective [0.0]: '',D11.3)'
          RDEFLT(1)=0.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=CC_OBJ_WEIGHT
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) CC_OBJ_WEIGHT=RDATA(1)

C*** Extra 'constraints' on the objective function
          IDEFLT(1)=1
          FORMAT='(/'' Specify any additional objective function '
     '      //'components [1]:'''//
     '      '/''   (1) None'''//
     '      '/''   (2) Surface Laplacian'''//
     '      '/''   (3) Unused'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=ACTN_IREGULARISE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ACTN_IREGULARISE=IDATA(1)

          IF(ACTN_IREGULARISE.EQ.2) THEN
            RDEFLT(1)=5.0d-1
            FORMAT='($,'' Enter the regularisation parameter'
     '        //' [0.5] : '',D11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=ACTN_REG_PARAM_LAPLACE
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-1.D-9,RMAX,INFO,
     &        ERROR,*9999)
            IF(IOTYPE.NE.3) ACTN_REG_PARAM_LAPLACE=RDATA(1)
          ENDIF !ACTN_IREGULARISE = 2

          IDEFLT(1)=1
          FORMAT='(/'' Specify the residual function calculation '
     '      //'[1]:'''//
     '      '/''   (1) Normal'''//
     '      '/''   (2) Optimised'''//
     '      '/''   (3) Unused'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=ACTN_IRESFUN
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ACTN_IRESFUN=IDATA(1)

C***      Set up mapping betw optimising params and physical params

          CALL GLOBALO(IDO,INP,ISIZE_MFI,ISIZE_PHI,
     '      LD,LDR,NBH,NBJ,NEELEM,NENP,NKB,NKHE,NKH,
     '      NLNO(1,nx_opt),NMNO(1,0,nx_opt),NNB,
     '      NONL(1,nx_opt),NONM(1,1,nx_opt),NONY(0,1,1,nr,nx_opt),
     '      NPL,NPLIST4,NPNE,NPNODE,NPNY(0,1,0,nx_opt),nr,NVHE,NVHP,
     '      nx_opt,nx_sol,NXI,NYNE,NYNO(0,1,1,nr,nx_opt),NYNP,NYNR,
     '      NYNY(0,1,nr,nx_opt),PAOPTY,CONY(0,1,1,nr,nx_opt),
     '      CYNO(0,1,1,nr,nx_opt),CYNY(0,1,nr,nx_opt),PAOPTI,PMAX,PMIN,
     '      XP,FIX,ERROR,*9999)

C***      Set the initial parameter values and their bounds
          DO noopti=1,NTOPTI
            IF(PAOPTY(noopti).EQ.1) THEN
              ny=NYNO(1,noopti,2,nr,nx_opt)
              nk=NPNY(1,ny,0,nx_sol)
              nv=NPNY(2,ny,0,nx_sol)
              nj=NPNY(3,ny,0,nx_sol)
              np=NPNY(4,ny,0,nx_sol)
              IF(nk.EQ.1) THEN
                PAOPTI(noopti)=YP(ny,1,nx_sol)
                PMIN(noopti)=YP(ny,1,nx_sol)-ACTN_MIN(1)
                IF(PMIN(noopti).LE.ACTN_MIN(2)) PMIN(noopti)=ACTN_MIN(2)
                PMAX(noopti)=YP(ny,1,nx_sol)+ACTN_MAX(1)
                IF(PMAX(noopti).GE.ACTN_MAX(2)) PMAX(noopti)=ACTN_MAX(2)
                !activation times must lie between ACTN_MIN(2) and
                !ACTN_MAX(2)
              ELSE
                ERROR='>>Can''t handle nk>1 yet'
                GOTO 9999
              ENDIF
            ELSE IF(PAOPTY(noopti).EQ.2) THEN !Transmembrane jump
              PAOPTI(noopti)=TRSF_ACTN_WAVE_JUMP
              PMIN(noopti)=TRSF_ACTN_WAVE_JUMP-ACTN_MIN(3)
              PMAX(noopti)=TRSF_ACTN_WAVE_JUMP+ACTN_MAX(3)
            ELSE IF(PAOPTY(noopti).EQ.3) THEN !Wave width
              PAOPTI(noopti)=TRSF_ACTN_WAVE_WIDTH
              PMIN(noopti)=MAX(TRSF_ACTN_WAVE_WIDTH-ACTN_MIN(4),0.0d0)
              PMAX(noopti)=TRSF_ACTN_WAVE_WIDTH+ACTN_MAX(4)
            ELSE
              ERROR='>>Invalid PAOPTY for activation fitting'
              GOTO 9999
            ENDIF
          ENDDO !noopti

C LKC 20-MAR-1999 Adding maximum major iteration question
          FORMAT='($,'' Do you wish to specify a maximum number'//
     '      ' of major iterations [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'N') THEN
            MAX_MAJOR_ITER=-1
          ELSE
            FORMAT=
     '        '($,'' Enter the maximum major iterations [50]: '',I3)'
            IDEFLT(1)=50
            IF(IOTYPE.EQ.3) IDATA(1)=MAX_MAJOR_ITER
c cpb 9/8/2000 Increasing max number of iterations
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        9999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '        *9999)
            IF(IOTYPE.NE.3) MAX_MAJOR_ITER=IDATA(1)
          ENDIF

C LKC 16-NOV-1999 Adding maximum minor iteration question
          FORMAT='($,'' Do you wish to specify a maximum number'//
     '      ' of minor iterations [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'N') THEN
            MAX_MINOR_ITER=-1
          ELSE
            FORMAT=
     '        '($,'' Enter the maximum minor iterations [50]: '',I3)'
            IDEFLT(1)=50
            IF(IOTYPE.EQ.3) IDATA(1)=MAX_MAJOR_ITER
C LKC 6-AUG-2000 increasing the max number of minor iterations
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
C     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
C     '        999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
C     '        *9999)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        3000,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '        *9999)
            IF(IOTYPE.NE.3) MAX_MINOR_ITER=IDATA(1)
          ENDIF




        ELSE IF(KTYP27.EQ.13) THEN !dipole

          IDEFLT(1)=1
          FORMAT='('' Specify the residual field(s) [1]:'''//
     '      '/''   (1) Magnetic'''//
     '      '/''   (2) Potential'''//
     '      '/''   (3) Both Potential and Magnetic'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=KTYP27B
          ENDIF

          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP27B=IDATA(1)

          ndipole=1
          CALL ASSERT(NDIPOLES(nr).EQ.1,
     '      '>> Only implemented for single dipole',
     '      ERROR,*9999)

          IF(KTYP27B.EQ.3) THEN
            FORMAT=
     '        '($,'' Enter the data file for potential electrodes: '',)'
            CDEFLT(1)='electrodes'
            IF(IOTYPE.EQ.3) CDATA(1)=DATASET1
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) DATASET1=CDATA(1)(1:100)

            FORMAT=
     '        '($,'' Enter the data file for magnetic sensors: '',A)'
            CDEFLT(1)='electrodes'
            IF(IOTYPE.EQ.3) CDATA(1)=DATASET2
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) DATASET2=CDATA(1)(1:100)
          ENDIF !KTYP27B


          CALL GLOBALO(IDO,INP,ISIZE_MFI,ISIZE_PHI,
     '      LD,LDR,NBH,NBJ,NEELEM,NENP,NKB,NKHE,NKH,
     '      NLNO(1,nx_opt),NMNO(1,0,nx_opt),NNB,
     '      NONL(1,nx_opt),NONM(1,1,nx_opt),NONY(0,1,1,nr,nx_opt),
     '      NPL,NPLIST4,NPNE,NPNODE,NPNY(0,1,0,nx_opt),nr,NVHE,NVHP,
     '      nx_opt,nx_sol,NXI,NYNE,NYNO(0,1,1,nr,nx_opt),NYNP,NYNR,
     '      NYNY(0,1,nr,nx_opt),PAOPTY,CONY(0,1,1,nr,nx_opt),
     '      CYNO(0,1,1,nr,nx_opt),CYNY(0,1,nr,nx_opt),PAOPTI,PMAX,PMIN,
     '      XP,FIX,ERROR,*9999)


C LKC 17-MAY-2002 replaced with GETDIPOLE
C          IF(DIPOLE_CEN_NTIME(ndipole,nr).EQ.0) THEN
C            DO nj=1,NJT
C              PAOPTI(nj)=DIPOLE_CEN(nj,0,ndipole,nr)
C           ENDDO
C          ELSE
C            DO nj=1,NJT
C              PAOPTI(nj)=DIPOLE_CEN(nj,STIME,ndipole,nr)
C            ENDDO
C            WRITE(OP_STRING(1),
C     '        '('' Time varying direction: using nt=1 as estimate'')')
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          ENDIF


          CALL GETDIPOLE(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,ndipole,
     '      nr,PAOPTI(1),DIPOLE_CEN,DIPOLE_DIR,PAOPTI(NJT+1),TIME,ERROR,
     '      *9999)

C LKC Bounds must be positive (the minimum bound is subtracted from the
C initial guess and the max bound is added to the initial guess).
          FORMAT='($,'' Enter min & max change limit for dipole '
     '      //'center [10.0,10.0]: '',2(1X,D11.3))'
          RDEFLT(1)=10.0d0
          RDEFLT(2)=10.0d0
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nj=1,NJT
              PMIN(nj)=PAOPTI(nj)-RDATA(1)
              PMAX(nj)=PAOPTI(nj)+RDATA(2)
            ENDDO
          ENDIF

          IF(DOP) THEN
            WRITE(OP_STRING(1),'('' Dipole Center  '',3F12.3)')
     '        (PAOPTI(nj),nj=1,NJT)
            WRITE(OP_STRING(2),'(''   Upper Bound  '',3F12.3)')
     '        (PMAX(nj),nj=1,NJT)
            WRITE(OP_STRING(3),'(''   Lower Bound  '',3F12.3)')
     '        (PMIN(nj),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C LKC 17-MAY-2002 replaced with GETDIPOLE
C          IF(DIPOLE_DIR_NTIME(ndipole,nr).EQ.0) THEN
C            DO nj=1,NJT
C              PAOPTI(nj+NJT)=DIPOLE_DIR(nj,0,ndipole,nr)
C            ENDDO
C          ELSE
C            DO nj=1,NJT
C              PAOPTI(nj+NJT)=DIPOLE_DIR(nj,1,ndipole,nr)
C            ENDDO
C            WRITE(OP_STRING(1),
C     '        '('' Time varying direction: using nt=1 as estimate'')')
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C          ENDIF

          FORMAT='($,'' Enter min & max change limit for dipole '
     '      //'components [1.0,1.0]: '',2(1X,D11.3))'
          RDEFLT(1)=1.0d0
          RDEFLT(2)=1.0d0
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nj=1,NJT
              PMIN(nj+NJT)=PAOPTI(nj+NJT)-RDATA(1)
              PMAX(nj+NJT)=PAOPTI(nj+NJT)+RDATA(2)
            ENDDO
          ENDIF


          IF(DOP) THEN
            WRITE(OP_STRING(1),'('' Dipole Orient. '',3F12.3)')
     '        (PAOPTI(nj+NJT),nj=1,NJT)
            WRITE(OP_STRING(2),'(''   Upper Bound  '',3F12.3)')
     '        (PMAX(nj+NJT),nj=1,NJT)
            WRITE(OP_STRING(3),'(''   Lower Bound  '',3F12.3)')
     '        (PMIN(nj+NJT),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF


          FORMAT='($,'' Do you wish to specify a maximum number'//
     '      ' of major iterations [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'N') THEN
            MAX_MAJOR_ITER=-1
          ELSE
            FORMAT=
     '        '($,'' Enter the maximum major iterations [50]: '',I3)'
            IDEFLT(1)=50
            IF(IOTYPE.EQ.3) IDATA(1)=MAX_MAJOR_ITER
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '        *9999)
            IF(IOTYPE.NE.3) MAX_MAJOR_ITER=IDATA(1)
          ENDIF

          FORMAT='($,'' Do you wish to specify a maximum number'//
     '      ' of minor iterations [N]? '',A)'
          IF(IOTYPE.EQ.3) ADATA(1)='N'
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'N') THEN
            MAX_MINOR_ITER=-1
          ELSE
            FORMAT=
     '        '($,'' Enter the maximum minor iterations [50]: '',I3)'
            IDEFLT(1)=50
            IF(IOTYPE.EQ.3) IDATA(1)=MAX_MAJOR_ITER
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        999,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '        *9999)
            IF(IOTYPE.NE.3) MAX_MINOR_ITER=IDATA(1)
          ENDIF

        ELSE
          ERROR='Unknown KTYP27'
          GOTO 9999
        ENDIF !KTYP27


      ELSE IF(KTYP26.EQ.3) THEN !optimising micro-struct parameters
        IF(IOTYPE.EQ.3) THEN
          IDEFLT(1)=NH_LOC(0,nx_opt)
        ENDIF
        FORMAT='($,'' Enter the number of '
     '    //'optimisation variables [1]: '',I2)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NJM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) nhx=IDATA(1)

        DO i=1,nhx
          WRITE(CHAR7,'(I2)') i
          CALL STRING_TRIM(CHAR7,IBEG1,IEND1)
          IDEFLT(1)=4
          IF(IOTYPE.EQ.3) THEN
C LKC 7-DEC-2000 NJH(i) has not been set -- may as well set IDELT to i
C            IDEFLT(1)=NJH(i)
            IDEFLT(1)=i
          ENDIF
          FORMAT='($,'' Enter optimisation variable number '
     '      //CHAR7(IBEG1:IEND1)//' [4]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NJM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            nj=IDATA(1)
            NJH(i)=nj
          ENDIF
        ENDDO !i

        NH_LOC(0,nx_opt)=nhx
        DO i = 1, NH_LOC(0,nx_opt)
          NH_LOC(0,0)=NH_LOC(0,0)+1
          NH_LOC(i,nx_opt)=NH_LOC(0,0)
        ENDDO

        CALL CALC_NY_MAPS_IND_OPT(NJH,NKJ,NPNODE,NPNY,nr,NVJP,
     '    nx_opt,NYNP,NYNR,ERROR,*9999)

C***    Initialise arrays
        IF(IOTYPE.NE.3) THEN
C***      Set all ny's in the region to being fixed (ie out of the opt)
          DO no_nynr=1,NYNR(0,0,1,nr,nx_opt) !loop over global variables
            ny=NYNR(no_nynr,0,1,nr,nx_opt) !is global variable number
            DO iy=1,NIYFIXM
              FIX(ny,iy,nx_opt)=.TRUE.
            ENDDO !iy
          ENDDO !no_nynr
C***     Set all ny's in the region and in the nj_opt to being free (ie
C***     in the opt)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nhx=1,NH_LOC(0,nx_opt)
              nh=NH_LOC(nhx,nx_opt)
              nj=NJH(nhx)
              DO nv=1,NVJP(nj,np)
                DO nk=1,NKJ(nj,np)
                  ny=NYNP(nk,nv,nh,np,0,1,nr)
                  DO iy=1,NIYFIXM
                    FIX(ny,iy,nx_opt)=.FALSE.
                  ENDDO !iy
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhx
          ENDDO !nonode (np)
        ENDIF ! IOTYPE.NE.3

        nonode=0
 6100   FORMAT='(/$,'' Enter node #s/name to fix [EXIT]: '',I4)'
        IF(IOTYPE.EQ.3) THEN
          nonode=nonode+1
          IF(nonode.LE.NPNODE(0,nr)) THEN
            np=NPNODE(nonode,nr)
            IDATA(1)=np
          ELSE
            IDATA(1)=0
          ENDIF
        ENDIF
 6500   CDATA(1)='NODES' !for use with group input
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IZERO,0,NPT(nr),
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).NE.0) THEN !not default exit
          NPLIST4(0)=IDATA(0)
          DO n=1,IDATA(0)
            NPLIST4(n)=IDATA(n)
            np=IDATA(n)
            IF(.NOT.INLIST(np,NPNODE(1,nr),
     '        NPNODE(0,nr),N1)) THEN
              WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '          //'belong to the current region'')') np
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              GOTO 6500
            ENDIF
          ENDDO !n

C         Define optimise condition for first node in group
          np=NPLIST4(1) !rest of group filled at end of njj loop
          DO njj=1,NH_LOC(0,nx_opt)
            WRITE(CHAR1,'(I10)') njj
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            DO nhx=1,NH_LOC(0,nx_opt)
              nh=NH_LOC(nhx,nx_opt)
              nj=NJH(nhx)
              WRITE(CHAR2,'(I1)') nhx
              ADEFLT(1)='N'
              IF(NVJP(nj,np).EQ.1.AND.NKJ(nj,np).EQ.1) THEN
                FORMAT='($,'' Is variable '//CHAR2(1:1)
     '            //' of optimisation '
     '            //'variable '//CHAR1(IBEG1:IEND1)//' fixed '
     '            //'[N]?: '',A)'
              ELSE
                FORMAT='($,'' Are any variables for variable '
     '            //CHAR2(1:1)//' of optimisation variable '
     '            //CHAR1(IBEG1:IEND1)//' fixed [N]?: '',A)'
              ENDIF
              IF(IOTYPE.EQ.3) THEN
                ISFIXED=.FALSE.
                DO nv=1,NVJP(nj,np)
                  DO nk=1,NKJ(nj,np)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    IF(FIX(ny,1,nx_opt)) ISFIXED=.TRUE.
                  ENDDO !nk
                ENDDO !nv
                IF(ISFIXED) THEN
                  ADATA(1)='Y'
                ELSE
                  ADATA(1)='N'
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
                  ISFIXED=.TRUE.
                ELSE
                  ISFIXED=.FALSE.
                ENDIF
              ENDIF
              IF(ISFIXED) THEN
                IF(NVJP(nj,np).EQ.1.AND.NKJ(nj,np).EQ.1)
     '            THEN
C*** Don't need to ask a question as it has already been answered
                  ny=NYNP(1,1,nh,np,0,1,nr)
                  FIX(ny,1,nx_opt)=.TRUE.
                ELSE
                  DO nv=1,NVJP(nj,np) !loop over versions
                    IF(NVJP(nj,np).GT.1) THEN
                      WRITE(CHAR3,'(I2)') nv
                      FORMAT='('' For version '
     '                  //'number '//CHAR3(1:2)
     '                  //':'')'
                      CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     '                  FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '                  ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '                  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                    ENDIF
                    DO nk=1,NKJ(nj,np)
                      ny=NYNP(nk,nv,nh,np,0,1,nr)
                      ADEFLT(1)='N'
                      IF(nk.EQ.1) THEN
                        FORMAT='($,'' Is variable '//CHAR2(1:1)
     '                    //' of optimisation variable '
     '                    //CHAR1(IBEG1:IEND1)
     '                    //' fixed [N]?: '',A)'
                      ELSE IF(nk.GT.1) THEN
                        WRITE(CHAR3,'(I1)') nk
                        FORMAT='($,'' Is variable '//CHAR2(1:1)
     '                    //' of optimisation variable '
     '                    //CHAR1(IBEG1:IEND1)//' derivative '
     '                    //CHAR3(1:1)//' fixed [N]?: '',A)'
                      ENDIF
                      IF(IOTYPE.EQ.3) THEN
                        IF(FIX(ny,1,nx_opt)) THEN
                          ADATA(1)='Y'
                        ELSE
                          ADATA(1)='N'
                        ENDIF
                      ENDIF
                      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,
     '                  FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '                  ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '                  RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                      IF(IOTYPE.NE.3) THEN
                        IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
                          FIX(ny,1,nx_opt)=.TRUE.
                        ELSE
                          FIX(ny,1,nx_opt)=.FALSE.
                        ENDIF
                      ENDIF
                    ENDDO !nk
                  ENDDO !nv
                ENDIF
              ENDIF !isfixed
            ENDDO !nhx
          ENDDO !njj

C         Apply opt conditions to rest of nodes group
          DO n=2,NPLIST4(0)
            np2=NPLIST4(n)
            DO nhx=1,NH_LOC(0,nx_opt)
              nh=NH_LOC(nhx,nx_opt)
              nj=NJH(nhx)
              DO nv=1,NVJP(nj,np2)
                DO nk=1,NKJ(nj,np2)
                  ny=NYNP(nk,nv,nh,np2,0,1,nr)
                  IF(NVJP(nj,np2).NE.NVJP(nj,np)) THEN
                    ! this will do for now
                    WRITE(OP_STRING,
     '                '('' >>WARNING: fixing higher versions too'')')
                      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    ny_first=NYNP(nk,1,nh,np,0,1,nr)
                  ELSE
                    ny_first=NYNP(nk,1,nh,np,0,1,nr)
                  ENDIF
                  IF(ny.NE.0) FIX(ny,1,nx_opt)=FIX(ny_first,1,nx_opt)
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhx
          ENDDO !n

          GO TO 6100 !for more nodes
        ENDIF !default exit

        no=0
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NH_LOC(0,nx_opt)
            WRITE(CHAR2,'(I1)') nhx
            nh=NH_LOC(nhx,nx_opt)
            nj=NJH(nhx)
            var=0
            DO nv=1,NVJP(nj,np)
              DO nk=1,NKJ(nj,np)
                var=var+1
                ny=NYNP(nk,nv,nh,np,0,1,nr)
                IF(.NOT.FIX(ny,1,nx_opt)) THEN
                  no=no+1
                  IF(IOTYPE.EQ.3) THEN
                    RDEFLT(1)=PAOPTI(no)*180.0d0/PI
                    RDEFLT(2)=PMIN(no)*180.0d0/PI
                    RDEFLT(3)=PMAX(no)*180.0d0/PI
                  ELSE
                    IF (nk.EQ.1) THEN
                      RDEFLT(1)=XP(nk,nv,nj,np)*180.0d0/PI
                      RDEFLT(2)=XP(nk,nv,nj,np)*180.0d0/PI-15.0d0
                      RDEFLT(3)=XP(nk,nv,nj,np)*180.0d0/PI+15.0d0
                    ELSE
                      RDEFLT(1)=XP(nk,nv,nj,np)
                      RDEFLT(2)=-1.0d0
                      RDEFLT(3)=1.0d0
                    ENDIF
                  ENDIF
                  RDATA(1)=RDEFLT(1)
                  RDATA(2)=RDEFLT(2)
                  RDATA(3)=RDEFLT(3)
                  WRITE(CHAR7,'(F6.1)') RDATA(1)
                  WRITE(CHAR8,'(F6.1)') RDATA(2)
                  WRITE(CHAR9,'(F6.1)') RDATA(3)
                  WRITE(CHAR10,'(I2)') var
                  WRITE(CHAR11,'(I3)') np
                  CALL STRING_TRIM(CHAR7,IBEG1,IEND1)
                  CALL STRING_TRIM(CHAR8,IBEG2,IEND2)
                  CALL STRING_TRIM(CHAR9,IBEG3,IEND3)
                  CALL STRING_TRIM(CHAR10,IBEG4,IEND4)
                  CALL STRING_TRIM(CHAR11,IBEG5,IEND5)
                  FORMAT=
     '              '($,'' Enter init min & max values for '
     '              //'var # '//CHAR10(IBEG4:IEND4)
     '              //' of opt var  '//CHAR2(1:1)
     '              //' at node '//CHAR11(IBEG5:IEND5)
     '              //' ['//CHAR7(IBEG1:IEND1)
     '              //','//CHAR8(IBEG2:IEND2)
     '              //','//CHAR9(IBEG3:IEND3)//']: '',3E11.3)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,
     &              FILEIP,FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) THEN
                    IF (RDATA(2).NE.RDATA(3)) THEN
                      PAOPTI(no)=RDATA(1)*PI/180.0d0
                      PMIN(no)  =RDATA(2)*PI/180.0d0
                      PMAX(no)  =RDATA(3)*PI/180.0d0
                    ENDIF
                  ENDIF
                  IF (RDATA(2).NE.RDATA(3)) THEN
                    DO nrc=1,2
                      NONY(0,ny,nrc,nr,nx_opt)=1
                      NONY(1,ny,nrc,nr,nx_opt)=no
                      NYNO(0,no,nrc,nr,nx_opt)=1
                      NYNO(1,no,nrc,nr,nx_opt)=ny
                    ENDDO
                  ELSE
                    no=no-1
                    FIX(ny,1,nx_opt)=.TRUE.
                  ENDIF
                ENDIF
              ENDDO !nk
            ENDDO !nv
          ENDDO !nhx
        ENDDO !nonode (np)
        NTOPTI=no

        FORMAT='('' Specify the form of the residuals [1]: '''//
     '    '/''   (1) Nearest Gauss points '
     '    //'in connected elements'''//
     '    '/''   (2) All Gauss points in connected '
     '    //'elements'''//
     '    '/''   (3) All Gauss points in connected '
     '    //'elements weighted'''//
     '    '/$,''   '',I2)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP29C
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP29C=IDATA(1)

        IF(KTYP29C.EQ.1) THEN
          NT_RES=0
          DO no=1,NTOPTI
            ny=NYNO(1,no,1,nr,nx_opt)
            np=NPNY(4,ny,0,nx_opt)
            nk=NPNY(1,ny,0,nx_opt)
            IF(nk.EQ.1) THEN
              DO noelem=1,NENP(np,0,nr)
                ne=NENP(np,noelem,nr)
                nb=NBJ(1,ne)
                DO nn=1,NNT(nb)
                  np2=NPNE(nn,nb,ne)
                  IF(np2.EQ.np) THEN
                    NT_RES=NT_RES+1
                  ENDIF
                ENDDO ! nn
              ENDDO ! noelem
            ENDIF
          ENDDO ! no
        ELSE IF((KTYP29C.EQ.2).OR.(KTYP29C.EQ.3)) THEN
          NT_RES=0
          DO no=1,NTOPTI
            ny=NYNO(1,no,1,nr,nx_opt)
            np=NPNY(4,ny,0,nx_opt)
            nk=NPNY(1,ny,0,nx_opt)
            IF(nk.EQ.1) THEN
              DO noelem=1,NENP(np,0,nr)
                ne=NENP(np,noelem,nr)
                NT_RES=NT_RES+NGT(NBJ(1,ne))
              ENDDO ! noelem
            ENDIF
          ENDDO ! no
        ENDIF

        FORMAT='($,'' Enter the command file name [OPTIMISE]: '',A)'
        CDEFLT(1)='OPTIMISE'
        IF(IOTYPE.EQ.3) CDATA(1)=COM_FILE
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) COM_FILE=CDATA(1)(1:255)


      ELSE IF(KTYP26.EQ.4) THEN !optimising holmes constants
        IF(IOTYPE.EQ.3) THEN
          IDEFLT(1)=NH_LOC(0,nx_opt)
        ENDIF

        IF(KTYP27.EQ.1) THEN
          nhx=1
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C1 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(1)=RDATA(1)
              PMIN(1)  =RDATA(2)
              PMAX(1)  =RDATA(3)
            ENDIF
          ENDIF

        ELSE IF(KTYP27.EQ.2) THEN
          nhx=1
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C2 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(1)=RDATA(1)
              PMIN(1)  =RDATA(2)
              PMAX(1)  =RDATA(3)
            ENDIF
          ENDIF

        ELSE IF(KTYP27.EQ.3) THEN
          nhx=2
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C1 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(1)=RDATA(1)
              PMIN(1)  =RDATA(2)
              PMAX(1)  =RDATA(3)
            ENDIF
          ENDIF
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C2 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(2)=RDATA(1)
              PMIN(2)  =RDATA(2)
              PMAX(2)  =RDATA(3)
            ENDIF
          ENDIF

        ELSE IF(KTYP27.EQ.4) THEN
          nhx=1
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C3 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(1)=RDATA(1)
              PMIN(1)  =RDATA(2)
              PMAX(1)  =RDATA(3)
            ENDIF
          ENDIF

        ELSE IF(KTYP27.EQ.5) THEN
          nhx=3
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C1 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(1)=RDATA(1)
              PMIN(1)  =RDATA(2)
              PMAX(1)  =RDATA(3)
            ENDIF
          ENDIF
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C2 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(2)=RDATA(1)
              PMIN(2)  =RDATA(2)
              PMAX(2)  =RDATA(3)
            ENDIF
          ENDIF
          RDEFLT(1)=1
          RDEFLT(2)=1
          RDEFLT(3)=1
          RDATA(1)=RDEFLT(1)
          RDATA(2)=RDEFLT(2)
          RDATA(3)=RDEFLT(3)
          FORMAT=
     '     '($,'' Enter init min & max values for C3 [1,1,1]:'',3E11.3)'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '    FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '    IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '    RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF (RDATA(2).NE.RDATA(3)) THEN
              PAOPTI(3)=RDATA(1)
              PMIN(3)  =RDATA(2)
              PMAX(3)  =RDATA(3)
            ENDIF
          ENDIF
        ENDIF

        NH_LOC(0,nx_opt)=nhx
        DO i = 1, NH_LOC(0,nx_opt)
          NH_LOC(0,0)=NH_LOC(0,0)+1
          NH_LOC(i,nx_opt)=NH_LOC(0,0)
        ENDDO
        NTOPTI=nhx

        IDEFLT(1)=1
        FORMAT='($,'' Enter the number of residuals [1]: '',I2)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NT_RES
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NREM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NT_RES=IDATA(1)

        FORMAT='($,'' Enter the command file name [OPTIMISE]: '',A)'
        CDEFLT(1)='OPTIMISE'
        IF(IOTYPE.EQ.3) CDATA(1)=COM_FILE
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) COM_FILE=CDATA(1)(1:255)
      ENDIF

      IF(   (KTYP26.EQ.1.AND.KTYP27.EQ.2)       !mat pars sqd reac diffs
     '  .OR.(KTYP26.EQ.1.AND.KTYP27.EQ.5)       !mat pars sqd data diffs
     '  .OR.(KTYP26.EQ.1.AND.KTYP27.EQ.6)       !mat pars reaction
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.5)       !Data Fitting 
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.7)       !aerofoil wake & stress
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.8)       !aerofoil lift
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.10)      !torso measurements
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.11)      !second moments
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.12)      !activation times
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.13)      !dipole inverses
     '  .OR.(KTYP26.EQ.3)                       !micro-struct fibre stress
     '  .OR.(KTYP26.EQ.4)) THEN                 !holmes constants
C***    optimisations that call E04UPF or MINOS

        IDEFLT(1)=10
        FORMAT='($,'' Enter print level required [10]: '',I5)'
        IF(IOTYPE.EQ.3) IDATA(1)=IPPLEV
        IF(KTYP29.EQ.1) THEN !NPSOL
          PRTLIM=30
        ELSE IF(KTYP29.EQ.2) THEN !MINOS
          PRTLIM=11111
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,PRTLIM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) IPPLEV=IDATA(1)

        IF(KTYP29.EQ.2) THEN !MINOS
          FORMAT='($,'' Do you want a summary file [N]? '',A)'
          CDEFLT(1)='N'
          IF(IOTYPE.EQ.3) THEN
            IF(SUMMFILE) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            SUMMFILE=(ADATA(1).EQ.'Y')
          ENDIF
        ENDIF

        RDEFLT(1)=1.0D-7
        FORMAT='($,'' Enter the optimality tolerance [1.0d-7]: '','
     '    //'E11.3)'
        IF(IOTYPE.EQ.3) RDATA(1)=OPTTOL
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) OPTTOL=RDATA(1)

        IF((KTYP26.EQ.1.AND.KTYP27.EQ.2).OR.  !mat pars, sqd reac diffs
     '    (KTYP26.EQ.2.AND.KTYP27.EQ.5).OR.   !data fitting
     '    (KTYP26.EQ.3)) THEN !Micro Structure, fibre
          RDEFLT(1)=0.9D0
        ELSE
          RDEFLT(1)=0.1D0
        ENDIF
        WRITE(CHAR1,'(F3.1)') RDEFLT(1)
        FORMAT='($,'' Enter the line search tolerance ['//CHAR1(1:3)
     '    //']: '',F7.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=LNSTOL
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,1.0D0,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) LNSTOL=RDATA(1)

        RDEFLT(1)=2.0D0
        WRITE(CHAR1,'(F3.1)') RDEFLT(1)
        FORMAT='($,'' Enter the step limit ['//CHAR1(1:3)
     '    //']: '',F7.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=STEPLIM
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) STEPLIM=RDATA(1)

        IF(IOTYPE.NE.3) NLFTOL=0.0D0
        IF(NTCNTR.GT.0) THEN
          RDEFLT(1)=1.0D-5
          FORMAT='($,'' Enter the nonlinear feasibility '
     '      //'tolerance [1.0d-5]: '',D11.3)'
          IF(IOTYPE.EQ.3) RDATA(1)=NLFTOL
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NLFTOL=RDATA(1)
        ENDIF

        IF(KTYP27.EQ.5) THEN !Data fitting
          IDEFLT(1)=3
        ELSEIF(KTYP27.EQ.12) THEN !inverse activation setup
          IDEFLT(1)=3
        ELSE
          IDEFLT(1)=0
        ENDIF
        WRITE(CHAR,'(I1)') IDEFLT(1)
        FORMAT='($,'' Enter derivative level ['//CHAR(1:1)//']: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=IPDLEV
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,-1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) IPDLEV=IDATA(1)

        IF(IOTYPE.NE.3) DIFFER=0.0D0
c        IF((IPDLEV.LT.3).OR.
c     '    (IPDLEV.LE.2.AND.(KTYP26.NE.1.OR.KTYP27.NE.2))) THEN
        IF((IPDLEV.EQ.0).OR.
     '    (IPDLEV.LE.2.AND.(KTYP26.NE.1.OR.KTYP27.NE.2))) THEN
C ***     For material params + squared reaction diffs only use
C ***     difference interval if derivative level=0
          RDEFLT(1)=0.0D0
          FORMAT='($,'' Enter the difference interval [computed]: '','
     '      //'D10.3)'
          IF(IOTYPE.EQ.3) RDATA(1)=DIFFER
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) DIFFER=RDATA(1)
        ENDIF

        IF(KTYP26.EQ.3) THEN ! Micro-structure
          IF(IOTYPE.NE.3) CEN_DIFFER=0.0D0
          IF(IPDLEV.EQ.0) THEN
            RDEFLT(1)=0.0D0
            FORMAT='($,'' Enter the central '
     '        //'difference interval [computed]: '','
     '        //'D10.3)'
            IF(IOTYPE.EQ.3) RDATA(1)=CEN_DIFFER
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) CEN_DIFFER=RDATA(1)
          ENDIF
        ENDIF

        IDEFLT(1)=-1
        FORMAT='($,'' Enter verify level [-1]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=IPVLEV
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,-1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) IPVLEV=IDATA(1)

        IF(IOTYPE.NE.3) THEN
          ISROCV=0
          ISPOCV=0
        ENDIF
        IF(IPVLEV.EQ.1.OR.IPVLEV.EQ.3) THEN
          IDEFLT(1)=1
          FORMAT='($,'' Enter the starting objective check variable '
     '      //'[  1]: '',I3)'
          IF(IOTYPE.EQ.3) IDATA(1)=ISROCV
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NTOPTI,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ISROCV=IDATA(1)
          IDEFLT(1)=NTOPTI
          WRITE(CHAR1,'(I3)') NTOPTI
          FORMAT='($,'' Enter the stopping objective check variable ['
     '      //CHAR1(1:3)//']: '',I3)'
          IF(IOTYPE.EQ.3) IDATA(1)=ISPOCV
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,ISROCV,
     '      NTOPTI,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '      *9999)
          IF(IOTYPE.NE.3) ISPOCV=IDATA(1)
        ENDIF

        IF(IOTYPE.NE.3) THEN
          ISRCCV=0
          ISPCCV=0
        ENDIF
        IF(NTCNTR.GT.0) THEN
          IF(IPVLEV.EQ.2.OR.IPVLEV.EQ.3) THEN
            IDEFLT(1)=1
            FORMAT='($,'' Enter the starting constraint check '
     '        //'variable [  1]: '',I3)'
            IF(IOTYPE.EQ.3) IDATA(1)=ISRCCV
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        NTOPTI,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) ISRCCV=IDATA(1)
            IDEFLT(1)=NTOPTI
            WRITE(CHAR1,'(I3)') NTOPTI
            FORMAT='($,'' Enter the stopping constraint check '
     '        //'variable ['//CHAR1(1:3)//']: '',I3)'
            IF(IOTYPE.EQ.3) IDATA(1)=ISPCCV
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        ISRCCV,NTOPTI,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     &        ERROR,*9999)
            IF(IOTYPE.NE.3) ISPCCV=IDATA(1)
          ENDIF
        ENDIF
        IF(KTYP29.EQ.2) THEN !MINOS
          IDEFLT(1)=0
          FORMAT='($,'' Enter the debug level [0]: '',I3)'
          IF(IOTYPE.EQ.3) IDATA(1)=DBGLEV
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,100,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) DBGLEV=IDATA(1)
          IDEFLT(1)=3*NTCNTR+10*NTOPTI
          WRITE(CHAR1,'(I5)') IDEFLT(1)
          FORMAT='($,'' Enter the iterations limit '','
     '      //'''['//CHAR1(1:5)//']: '',I5)'
          IF(IOTYPE.EQ.3) IDATA(1)=ITERLM
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) ITERLM=IDATA(1)
        ENDIF
        IF(KTYP27.EQ.12) THEN !inverse activation setup
          RDEFLT(1)=1.0d-3
          FORMAT='($,'' Enter function precision [1.0d-3]: '',E11.3)'
          IF(IOTYPE.EQ.3) RDATA(1)=FUNPREC
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RDELTA,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) FUNPREC=RDATA(1)
        ENDIF
      ENDIF


C!!! LKC 18-MAr-2002 Don't think this is for
C       "activation inverse modelling" -- is for cell modelling
C         ITYP2=2cellular based modelling
C         ITYP5=2 == Time integration
      IF(nx_sol.GT.0) THEN !check for activation modelling
        IF(ITYP5(nr,nx_sol).EQ.2.AND.ITYP2(nr,nx_sol).EQ.9) THEN
          NT_RES=NDT
        ENDIF
      ENDIF !nx_sol>0

      CALL_OPTI=.TRUE.

      CALL EXITS('IPOPTI')
      RETURN
 9999 CALL ERRORS('IPOPTI',ERROR)
      CALL EXITS('IPOPTI')
      RETURN 1
      END


