      SUBROUTINE EVSOLU(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  ISIZE_MFI,NBH,NBJ,NDDATA,NDET,NDIPOLES,NEELEM,NENP,NGAP,NHE,
     '  NHP,NKH,NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,NPNODE,NQET,NQNE,
     '  NQS,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,CE,
     '  CURVCORRECT,DET,DIPOLE_CEN,DIPOLE_DIR,DL,DRDN,PG,MFI,
     '  RAD,RD,RG,SE,WD,WG,XA,XE,XG1,XIG,XN,XP,XQ,XR,YD,YP,YQ,
     '  ZA,ZD,ZE,ZF,ZP,STRING,RET_ERROR,*)

C#### Subroutine: EVSOLU
C###  Description:
C###    EVSOLU calculates the solution at mulitple domain points
C###    for boundary element solutions.
C**** Created by Martin Buist, June 2001

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'time02.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),ISIZE_MFI(3,NSSM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDDATA(0:NDM,0:NRM),
     '  NDET(NBFM,0:NNM),NDIPOLES(NRM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),NGAP(NIM,NBM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),
     '  NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),NQNE(NEQM,NQEM),NQS(NEQM),
     &  NRE(NEM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     &  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     &  NXLIST(0:NXM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      
      REAL*8 CE(NMM,NEM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  DET(NBFM,0:NNM,NGM,6),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),DL(3,NLM),DRDN(NGM),
     '  MFI(NDM,NTSM,3,NSSM),
     '  PG(NSM,NUM,NGM,NBM),RAD(NGM),RD(NGM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WD(NJM,NDM),WG(NGM,NBM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG1(NJM,NUM,NGM),
     '  XIG(NIM,NGM,NBM),XN(NJM,NGM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  XR(NJM,NGM),YD(NHM),YP(NYM,NIYM,NXM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZE(NSM,NHM),ZF(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER RET_ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER ERR,IBEG,IEND,IFROMC,iy,N3CO,nd,ne,ng,nh,nhx,nj,nk,
     '  noelem,nonode,np,nq,nqe,nr,NRLIST_LOCAL(0:99),nrr,ns,nss,nt,
     '  nu,nv,nx,nxc,nxgrid,ny,SQUID_CONFIG
      REAL*8 BASELINE,CPCUTOFF,LOC_TIME,ZD_TMP(3)
      REAL*8 SQUID_PT1(3),RAD_COIL,CURRENT,HEIGHT_COIL,HEIGHT_COIL2,
     '  BFIELD(3),BFIELD_2(3)
      REAL*8 DET_ADAPT(99,0:64,99),DRDN_(99),GRADPHI(3),GRADPHI2(3),
     '  PG_J(64,11,99),PG_Q(64,11,99),PG_U(64,11,99),
     '  RAD_(99),RD_(99),RG_(99),XE_(64,9),XG1_(9,11,99),
     '  XIG_J(3,99),XIG_Q(3,99),XIG_U(3,99),
     '  XN_(9,99),XR_(9,99),YD_LOCAL(9,9),ZE_(64,9),ZF_(64,9)
      REAL ELAPSED_TIME,TIME_START(1),TIME_STOP(1)
      CHARACTER ERROR*(ERRSTRLEN),EVALTYPE*8
      LOGICAL ALL_REGIONS,ANALYTIC,CBBREV,ELECTRIC,ERROR_FLAG,GRADIENT,
     '  MAGDIPOLE,MAGVOLUME,OPTIML,POTENTIAL,
     '  VIEW,VPOTENTIAL,COIL

C*** List of channels which record more than one vector
C*** They are numbered using cmiss numbers, 1..19)
      INTEGER MULTICHANNEL(5) 
      DATA MULTICHANNEL / 1,8,11,14,17 / 
      
!     Functions
      INTEGER CALC_SAMPLE_FROM_TIME
      REAL*8 RFROMC
      LOGICAL INLIST
      
      CALL ENTERS('EVSOLU',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate solution grid electric
C###  Parameter:        <optimise>
C###    Use optimised code
C###  Parameter:        <iy #[1]>
C###    Grid solution index to write answers into
C###  Parameter:        <class #[1]>
C###    Boundary element solution class
C###  Parameter:        <grid_class #[2]>
C###    Boundary element solution class
C###  Parameter:        <region #[1]>
C###    Boundary element solution region
C###  Parameter:        <element #[1]>
C###    Element number of the single grid element
C###  Parameter:        <cpcutoff #[0.25]>
C###    All coefficients below this value are set to zero
C###  Parameter:        <potential>
C###    Calculate the domain potential field
C###  Parameter:        <gradient>
C###    Calculate the potential gradient field
C###  Description:
C###    Evaluate boundary element domain solution at grid
C###    points.

        OP_STRING(1)=STRING(1:IEND)//' grid electric'
        OP_STRING(2)=BLANK(1:15)//'<optimise>'
        OP_STRING(3)=BLANK(1:15)//'<iy #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<grid_class #[2]>'
        OP_STRING(6)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<element #[1]>'
        OP_STRING(8)=BLANK(1:15)//'<cpcutoff #[0.25]>'
        OP_STRING(9)=BLANK(1:15)//'<potential>'
        OP_STRING(10)=BLANK(1:15)//'<gradient>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM evaluate solution magnetic
C###  Parameter:        <optimise>
C###    Use optimised code
C###  Parameter:        <from (grid/data/MFI)[grid]>
C###    Evaluate the solution at <grid> points or <data> points
C###    or the magnetic field intensity matrix <MFI>.
C###  Parameter:        <iy #[1]>
C###    Grid solution index to write answers into
C###  Parameter:        <nss #[1]>
C###    MFI index to write answers into
C###  Parameter:        <class #[1]>
C###    Boundary element solution class
C###  Parameter:        <grid_class #[2]>
C###    Boundary element solution class
C###  Parameter:        <region #[1]>
C###    Boundary element solution region
C###  Parameter:        <element #[1]>
C###    Element number of the single grid element
C###  Parameter:        <analytic>
C###    Calculate the analytic solution field for a dipole
C###    in the x direction on the z axis
C###  Parameter:        <time #[0.0]>
C###    The time at which the solution is being evaluated
C###  Parameter:        <frequency #[-]>
C###    The frequency to map between time and time-step
C###    (indices in the storage arrays)
C###  Parameter:        <dipole/volume/both [both]>
C###    Specify if calculating the magnetic field from any
C###    dipole sources, just the volume conductor, or both.
C###  Parameter:        <vector_potential>
C###    This is used as an alternative to using dipole or
C###    volume and cannot be used with those qualifiers. It
C###    calculates the vector potential of the solution field.
C###    Note this does not include the vector potential of any
C###    free-space sources.
C###  Parameter:        <noview>
C###    Do not diplay information to the screen
C###  Parameter:        <squid_config #[0]>
C###    Specify the SQUID configuration.
C###    If configuration 0, then compute all components for
C###        all electrodes.
C###    If configuration 1 (Vanderbilt SQUID)
C###    Only compute 3 components for channels 1,8,11,14,17.
C###    All other channels (total of 19) compute only the z component.
C###    If configuration 2 (y component only)
C###    Only compute only the y component for all channels. The y
C###    component corresponds to the anterior/posterior direction
C###  Parameter:        <baseline #[0.0]>
C###    Calculate gradients of the MFI about a specified baseline.
C###    If the default baseline of 0 is used, then absolute fields
C###    are calculated.
C###  Description:
C###    Evaluate boundary element domain solution at grid
C###    points.
C###  Parameter:        <coil>
C###    Specify if calculating the magnetic field due to a circular
C###    current filament (a one turn coil)
C###  Parameter:        <rad_coil #[0.0]>
C###    Specify the radius of the one turn coil
C###  Parameter:        <current #[0.0]>
C###    Specify the current running through the one turn coil
C###  Parameter:        <height_coil #[0.0]>
C###    Specify the height of the 1 turn coil above the surface of
C###    body

        OP_STRING(1)=STRING(1:IEND)//' magnetic'
        OP_STRING(2)=BLANK(1:15)//'<optimise>'
        OP_STRING(3)=BLANK(1:15)//'<from (grid/data/MFI)[grid]>'
        OP_STRING(4)=BLANK(1:15)//'<nss #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<iy #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(7)=BLANK(1:15)//'<grid_class #[2]>'
        OP_STRING(8)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(9)=BLANK(1:15)//'<element #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<analytic>'
        OP_STRING(11)=BLANK(1:15)//'<time #[0.0]>'
        OP_STRING(12)=BLANK(1:15)//'<frequency #[-]>'
        OP_STRING(13)=BLANK(1:15)//'<dipole/volume/both [both]>'
        OP_STRING(14)=BLANK(1:15)//'<vector_potential>'
        OP_STRING(15)=BLANK(1:15)//'<noview>'
        OP_STRING(16)=BLANK(1:15)//'<squid_config #[0]>'
        OP_STRING(17)=BLANK(1:15)//'<baseline #[0.0]>'
        OP_STRING(18)=BLANK(1:15)//'<coil>'
        OP_STRING(19)=BLANK(1:15)//'<rad_coil #[0.0]>'
        OP_STRING(20)=BLANK(1:15)//'<current #[0.0]>'
        OP_STRING(21)=BLANK(1:15)//'<height_coil #[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVSOLU',ERROR,*9999)
      ELSE

C LKC 7-MAR-2002 No longer need to have a grid for magnetic fields
C
C        TYPE=' '
C        IF(CBBREV(CO,'GRID',2,noco+1,NTCO,N3CO)) THEN
C          TYPE='GRID'
C        ENDIF
C        IF(TYPE(1:1).EQ.' ') THEN
C          CALL STRING_TRIM(STRING,IBEG,IEND)
C          CALL STRING_TRIM(CO(noco),IBEG1,IEND1)
C          CTEMP=CO(noco)(IBEG1:IEND1)
C          CALL STRING_TRIM(CTEMP,IBEG1,IEND1)
C          STRING=STRING(1:IEND)//' '//CTEMP(IBEG1:IEND1)
C          CO(noco+1)='?'
C          GO TO 1
C        ENDIF

        IF(CBBREV(CO,'ELECTRIC',3,noco+1,NTCO,N3CO)) THEN
          ELECTRIC=.TRUE.
        ELSEIF(CBBREV(CO,'MAGNETIC',3,noco+1,NTCO,N3CO)) THEN
          ELECTRIC=.FALSE.
        ELSE
          ELECTRIC=.TRUE.
        ENDIF
      
C    HB 07-OCT-2004 Define whether we are solving for a coil
        IF(CBBREV(CO,'COIL',3,noco+1,NTCO,N3CO)) THEN
          COIL=.TRUE.
        ELSE
          COIL=.FALSE.
        ENDIF
        
        IF(.NOT.COIL) THEN
           CALL ASSERT(CALL_SOLV,'>>Must have a solution',ERROR,*9999)
        ENDIF

        IF(ELECTRIC) THEN
          CALL ASSERT(CALL_GRID,
     '      '>>Must have a grid defined for electric fields',
     '      ERROR,*9999)

          EVALTYPE='GRID' ! potential only implemented for grid

        ELSE

          IF(CBBREV(CO,'GRID',3,noco+1,NTCO,N3CO)) THEN
            EVALTYPE='GRID'
          ELSEIF(CBBREV(CO,'DATA',3,noco+1,NTCO,N3CO)) THEN
            EVALTYPE='DATA'
          ELSEIF(CBBREV(CO,'MFI',3,noco+1,NTCO,N3CO)) THEN
            EVALTYPE='MFI'
            CALL ASSERT(USE_MAGNETIC.EQ.1,'>> Set USE_MAGNETIC to 1',
     '        ERROR,*9999)

          ELSE
            EVALTYPE='UNKNOWN'
          ENDIF

        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
           nr=NRLIST(1)

        IF(.NOT.COIL) THEN
           CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
           nxc=NXLIST(1)
           CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
           CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'IY',2,noco+1,NTCO,N3CO)) THEN
          iy=IFROMC(CO(N3CO+1))
        ELSE
          iy=1
        ENDIF

        IF(CBBREV(CO,'NSS',3,noco+1,NTCO,N3CO)) THEN
          nss=IFROMC(CO(N3CO+1))
        ELSE
          nss=1
        ENDIF

        IF(EVALTYPE(1:4).EQ.'GRID') THEN
          IF(CBBREV(CO,'ELEMENT',4,noco+1,NTCO,N3CO)) THEN
            ne=IFROMC(CO(N3CO+1))
          ELSE
            ne=1
          ENDIF

C LKC 19-JULY-2006 Check that grid scheme has been defined in
C  that element
          IF(CALL_GRID) THEN
            CALL ASSERT(NQS(ne).GT.0,
     &        '>>ERROR: No grid scheme defined for element',
     &        ERROR,*9999)
          ELSE
            ERROR='Grids have not been defined yet'
            GOTO 9999
          ENDIF
        ENDIF
        
        IF(CBBREV(CO,'OPTIMISE',2,noco+1,NTCO,N3CO)) THEN
          OPTIML=.TRUE.
        ELSE
          OPTIML=.FALSE.
        ENDIF

        IF(CBBREV(CO,'CPCUTOFF',2,noco+1,NTCO,N3CO)) THEN
          CPCUTOFF=RFROMC(CO(N3CO+1))
        ELSE
          CPCUTOFF=0.25d0
        ENDIF

        IF(CBBREV(CO,'POTENTIAL',3,noco+1,NTCO,N3CO)) THEN
          POTENTIAL=.TRUE.
        ELSE
          POTENTIAL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'GRADIENT',2,noco+1,NTCO,N3CO)) THEN
          GRADIENT=.TRUE.
          CALL ASSERT(NIQM.GE.(iy+3),'>>ERROR: Increase NIQM',
     &      ERROR,*9999)
        ELSE
          GRADIENT=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ANALYTIC',3,noco+1,NTCO,N3CO)) THEN
          ANALYTIC=.TRUE.
        ELSE
          ANALYTIC=.FALSE.
        ENDIF

        IF(CBBREV(CO,'DIPOLE',3,noco+1,NTCO,N3CO)) THEN
          MAGDIPOLE=.TRUE.
          MAGVOLUME=.FALSE.
        ELSEIF(CBBREV(CO,'VOLUME',3,noco+1,NTCO,N3CO)) THEN
          MAGDIPOLE=.FALSE.
          MAGVOLUME=.TRUE.
        ELSE
          MAGDIPOLE=.TRUE.
          MAGVOLUME=.TRUE.
        ENDIF

        IF(CBBREV(CO,'VECTOR_POTENTIAL',4,noco+1,NTCO,N3CO)) THEN
          VPOTENTIAL=.TRUE.
          MAGDIPOLE=.FALSE.
          MAGVOLUME=.FALSE.
        ELSE
          VPOTENTIAL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'TIME',3,noco+1,NTCO,N3CO)) THEN
          LOC_TIME=RFROMC(CO(N3CO+1))
        ELSE
          LOC_TIME=0.0d0
        ENDIF

C Over-ride the frequency possibly set in the iptrsf file
        IF(CBBREV(CO,'FREQUENCY',3,noco+1,NTCO,N3CO)) THEN
          TRSF_FREQUENCY=RFROMC(CO(N3CO+1))
        ENDIF

        IF(CBBREV(CO,'GRID_CLASS',6,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=2
        ENDIF

        IF(ELECTRIC.OR.EVALTYPE(1:4).EQ.'GRID') THEN
          CALL NX_LOC(NX_INQUIRE,nxc,nxgrid,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nxgrid.NE.0,
     &      '>>No nx defined for this grid class',ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'NOVIEW',3,noco+1,NTCO,N3CO)) THEN
          VIEW=.FALSE.
        ELSE
          VIEW=.TRUE.
        ENDIF

        IF(CBBREV(CO,'SQUID_CONFIG',3,noco+1,NTCO,N3CO)) THEN
          SQUID_CONFIG=IFROMC(CO(N3CO+1))
        ELSE
          SQUID_CONFIG=0
        ENDIF

        IF(CBBREV(CO,'BASELINE',3,noco+1,NTCO,N3CO)) THEN
          BASELINE=RFROMC(CO(N3CO+1))

          IF(DABS(BASELINE).LT.ZERO_TOL) THEN !no baseline
            BASELINE=0.d0
          ELSE
            CALL ASSERT(BASELINE.GT.0.d0,
     &        '>>Baseline must be greater than 0.0',ERROR,*9999)
            GRADIENT=.TRUE.
          ENDIF !no baseline
          
        ELSE
          BASELINE=0.d0
C LKC 27-MAR-2003 Shouldn't need to set this because gradient will
C         already be set by the "Potential" Gradient flags.
C          !hopefully this doesn't conflict with "potential" gradients
C          GRADIENT=.FALSE. 
        ENDIF
  
C    HB 07-OCT-04 Define the size of the coil radius in mm 
        IF(CBBREV(CO,'RAD_COIL',3,noco+1,NTCO,N3CO)) THEN
            RAD_COIL=RFROMC(CO(N3CO+1))
        ENDIF

C     Define the current going into the coil (can be -ve or +ve)
          IF(CBBREV(CO,'CURRENT',3,noco+1,NTCO,N3CO)) THEN
             CURRENT=RFROMC(CO(N3CO+1))
          ENDIF

C     Define the height of the coil above the surface of the body
          IF(CBBREV(CO,'HEIGHT_COIL',3,noco+1,NTCO,N3CO)) THEN
             HEIGHT_COIL=RFROMC(CO(N3CO+1))
          ENDIF

C**** Start calculations ....

        CALL CPU_TIMER(CPU_USER,TIME_START)

        CALL ASSERT(NBFM.LE.99,'>>Increase NBFM hard coded '
     &    //'dimension in EVSOLU',ERROR,*9999)
        CALL ASSERT(NHM.LE.9,'>>Increase NHM hard coded dimension '
     &    //'in EVSOLU',ERROR,*9999)
        CALL ASSERT(NIM.LE.3,'>>Increase NIM hard coded dimension '
     &    //'in EVSOLU',ERROR,*9999)
        CALL ASSERT(NJM.LE.9,'>>Increase NJM hard coded dimension '
     &    //'in EVSOLU',ERROR,*9999)
        CALL ASSERT(NGM.LE.99,'>>Increase NGM hard coded dimension '
     &    //'in EVSOLU',ERROR,*9999)
        CALL ASSERT(NKM.LE.8,'>>Increase NKM hard coded dimension '
     &    //'in EVSOLU',ERROR,*9999)
        CALL ASSERT(NNM.LE.64,'>>Increase NNM hard coded dimension '
     &    //'in EVSOLUEVSOLU',ERROR,*9999)
        CALL ASSERT(NSM.LE.64,'>>Increase NSM hard coded dimension '
     &    //'in ',ERROR,*9999)
        CALL ASSERT(NUM.LE.11,'>>Increase NUM hard coded dimension '
     &    //'in EVSOLU',ERROR,*9999)


C!!! LKC 19-MAR-2002 can probably reorder these loops for more speed
        DO nj=1,9
          DO ns=1,64
            XE_(ns,nj)=0.0d0
          ENDDO
        ENDDO
        DO ng=1,99
          DO nu=1,11
            DO ns=1,64
              PG_J(ns,nu,ng)=0.0d0
              PG_Q(ns,nu,ng)=0.0d0
              PG_U(ns,nu,ng)=0.0d0
            ENDDO
          ENDDO
        ENDDO

        IF(ELECTRIC) THEN
          IF(OPTIML) THEN
            IF(NJT.EQ.2) THEN
              IF(POTENTIAL) THEN
                DO nqe=1,NQET(1)
                  nq=NQNE(ne,nqe)
                  CALL DOMSOLOPTI2D(IBT,IDO,INP,NBH,NBJ,NEELEM,NGAP,
     '              NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),NKHE,NKJE,NLL,
     '              NPF,NP_INTERFACE,NPNE,NPNODE,nr,NVHE,NVHP(1,1,1,nr),
     '              NVJE,NW(1,1,nx),nx,NYNP,CE(1,1,nx),CPCUTOFF,
     '              CURVCORRECT,DET,DET_ADAPT,DL,DRDN,PG,PG_J,PG_Q,PG_U,
     '              RAD,RD,RG,SE,WG,XE,XG1,XIG,XIG_J,XIG_Q,XIG_U,XN,XP,
     '              XQ(1,nq),XR,YQ(nq,iy,1,nxgrid),YP(1,1,nx),ZA,ZE,ZF,
     '              ZP,ERROR,*9999)
                ENDDO !nqe
              ENDIF
              IF(GRADIENT) THEN
                DO nqe=1,NQET(1)
                  nq=NQNE(ne,nqe)
                  CALL GRADDOMSOLOPTI2D(NBH,NBJ,NEELEM,NHE(1,nx),
     '              NHP(1,nr,nx),NKH(1,1,1,nr),NKHE,NKJE,NLL,NPF,
     '              NP_INTERFACE,NPNE,NPNODE,nr,NVHE,NVHP(1,1,1,nr),
     '              NVJE,NW(1,1,nx),nx,NYNE,NYNP,CE(1,1,nx),CURVCORRECT,
     '              DET,DL,DRDN,GRADPHI,PG,RAD,RG,SE,WG,XE,XG1,XN,XP,
     '              XQ(1,nq),XR,YP(1,1,nx),ZA,ZE,ZF,ZP,ERROR,*9999)
                    YQ(nq,iy+1,1,nxgrid)=GRADPHI(1)
                    YQ(nq,iy+2,1,nxgrid)=GRADPHI(2)
                ENDDO !nqe
              ENDIF
            ELSE
              IF(POTENTIAL) THEN
                DO noelem=1,NEELEM(0,nr)
                  IF(IGREN(nr).EQ.2) CE(1,NEELEM(noelem,nr),nx)=1.0d0
                ENDDO
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,
     '            NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                ERROR_FLAG=.FALSE.

C$OMP           PARALLEL DO
C$OMP&          PRIVATE(DET_ADAPT,DRDN_,nq,nqe,PG_J,PG_Q,PG_U,RAD_,RD_,
C$OMP&            RG_,XE_,XG1_,XIG_J, XIG_Q,XIG_U,XN_,XR_,ZE_,ZF_,
C$OMP&            ERROR),
C$OMP&          SHARED(IBT,IDO,INP,NBH,NBJ,NDET,NEELEM,NGAP,NHE,NHP,NKH,
C$OMP&            NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,NPNODE,NQNE,nr,
C$OMP&            NVHE,NVHP,NVJE,NW,nx,NYNP,CE,CPCUTOFF,CURVCORRECT,DET,
C$OMP&            DL,PG,SE,WG,XIG,XP,XQ,YQ,YP,ZA,ZP,ERROR_FLAG)
                DO nqe=1,NQET(1)
                  IF(.NOT.ERROR_FLAG) THEN
                    nq=NQNE(ne,nqe)
                    CALL DOMSOLOPTI3D(IBT,IDO,INP,NBH,NBJ,NDET,NEELEM,
     '                NGAP,NHE(1,nx),NKH(1,1,1,nr),NKHE,NKJE,NLL,NPF,
     '                NP_INTERFACE,NPNE,nr,NVHE,NVJE,NW(1,1,nx),nx,
     '                CE(1,1,nx),CPCUTOFF,CURVCORRECT,DET,DET_ADAPT,DL,
     '                DRDN_,PG,PG_J,PG_Q,PG_U,RAD_,RD_,RG_,SE,WG,XE_,
     '                XG1_,XIG,XIG_J,XIG_Q,XIG_U,XN_,XP,XQ(1,nq),XR_,
     '                YQ(nq,iy,1,nxgrid),ZA,ZE_,ZF_,ZP,
     '                ERROR,*100)
                    GO TO 102
C                   This statement is designed to be skipped if no
C                   errors occur. However if a error occurs within
C                   a subroutine the alternate return points to line
C                   100 and sets the flag
 100                CONTINUE
C$OMP               CRITICAL(EVSOLU_1)
                    ERROR_FLAG=.TRUE.
                    WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                    CALL WRITES(IOER,OP_STRING,ERROR,*101)
                    WRITE(OP_STRING,'(/'' >>An error occurred - '
     '                //'results may be unreliable!'')')
                    CALL WRITES(IOER,OP_STRING,ERROR,*101)
 101                CONTINUE
C$OMP               END CRITICAL(EVSOLU_1)
 102                CONTINUE
                  ENDIF
                ENDDO !nqe
C$OMP           END PARALLEL DO
              ENDIF
              IF(GRADIENT) THEN
                DO noelem=1,NEELEM(0,nr)
                  IF(IGREN(nr).EQ.2) CE(1,NEELEM(noelem,nr),nx)=1.0d0
                ENDDO
                CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,
     '            NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
                ERROR_FLAG=.FALSE.

C$OMP           PARALLEL DO
C$OMP&          PRIVATE(DRDN_,GRADPHI,RAD_,RG_,XE_,XG1_,XN_,XR_,ZE_,ZF_,
C$OMP&            nqe,nq,ERROR),
C$OMP&          SHARED(NBH,NBJ,NEELEM,NHE,NKH,NKHE,NKJE,NPF,
C$OMP&            NP_INTERFACE,NPNE,nr,NVHE,NVJE,NW,nx,CE,CURVCORRECT,
C$OMP&            DET,PG,SE,WG,XP,XQ,ZA,ZP,iy,ne,NQNE,nxgrid,
C$OMP&            ERROR_FLAG)
                DO nqe=1,NQET(1)
                  IF(.NOT.ERROR_FLAG) THEN
                    nq=NQNE(ne,nqe)
                    CALL GRADDOMSOLOPTI3D(NBH,NBJ,NEELEM,NHE(1,nx),
     '                NKH(1,1,1,nr),NKHE,NKJE,NPF,NP_INTERFACE,NPNE,nr,
     '                NVHE,NVJE,NW(1,1,nx),nx,CE(1,1,nx),CURVCORRECT,
     '                DET,DRDN_,GRADPHI,PG,RAD_,RG_,SE,WG,XE_,XG1_,XN_,
     '                XP,XQ(1,nq),XR_,ZA,ZE_,ZF_,ZP,ERROR,*200)
                    YQ(nq,iy+1,1,nxgrid)=GRADPHI(1)
                    YQ(nq,iy+2,1,nxgrid)=GRADPHI(2)
                    YQ(nq,iy+3,1,nxgrid)=GRADPHI(3)

                    GO TO 202
C                   This statement is designed to be skipped if no
C                   errors occur. However if a error occurs within
C                   a subroutine the alternate return points to line
C                   200 and sets the flag
 200                CONTINUE
C$OMP               CRITICAL(EVSOLU_2)
                    ERROR_FLAG=.TRUE.
                    WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                    CALL WRITES(IOER,OP_STRING,ERROR,*201)
                    WRITE(OP_STRING,'(/'' >>An error occurred - '
     '                //'results may be unreliable!'')')
                    CALL WRITES(IOER,OP_STRING,ERROR,*201)
 201                CONTINUE
C$OMP               END CRITICAL(EVSOLU_2)
 202                CONTINUE
                  ENDIF
                ENDDO !nqe
C$OMP           END PARALLEL DO
              ENDIF
            ENDIF
          ELSE
            IF(POTENTIAL) THEN
              DO nqe=1,NQET(1)
                nq=NQNE(ne,nqe)
                CALL DOMSOL2(NBH,NBJ,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,
     '            NPNODE,nr,NRE,NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,
     '            NYNE,NYNP,CE(1,1,nx),CURVCORRECT,DET,DL,DRDN,PG,RAD,
     '            RD,RG,SE,WG,XA,XE,XG1,XN,XP,XQ(1,nq),XR,YD,YP(1,1,nx),
     '            ZA,ZE,ZF,ZP,ERROR,*9999)
                YQ(nq,iy,1,nxgrid)=YD(1)
              ENDDO !nqe
            ENDIF
            IF(GRADIENT) THEN
              DO nqe=1,NQET(1)
                nq=NQNE(ne,nqe)
                CALL GRADDOMSOL(NBH,NBJ,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '            NKH(1,1,1,nr),NKHE,NKJE,NLL,NPF,NP_INTERFACE,NPNE,
     '            NPNODE,nr,NRE,NVHE,NVHP(1,1,1,nr),NVJE,NW(1,1,nx),nx,
     '            NYNE,NYNP,CE(1,1,nx),CURVCORRECT,DET,DL,DRDN,GRADPHI,
     '            PG,RAD,RG,SE,WG,XA,XE,XG1,XN,XP,XQ(1,nq),XR,
     '            YP(1,1,nx),ZA,ZE,ZF,ZP,ERROR,*9999)
                DO nj=1,NJT
                  YQ(nq,iy+nj,1,nxgrid)=GRADPHI(nj)
                ENDDO !nj
              ENDDO !nqe
            ENDIF
          ENDIF



C*** Magnetic field solution
        ELSE !magnetic based solution
          IF(EVALTYPE(1:4).EQ.'GRID') THEN
            CALL ASSERT(NIQM.GE.(iy+3),'>>ERROR: Increase NIQM',
     '        ERROR,*9999)
          ELSEIF(EVALTYPE(1:4).EQ.'DATA') THEN
            CALL ASSERT(NJM.GE.(NJT+3),'>>ERROR: Increase NJM >= NJT+3',
     '        ERROR,*9999)
            CALL ASSERT(USE_DATA.EQ.1,'>>Set USE_DATA to 1',ERROR,*9999)

            CALL ASSERT(.NOT.OPTIML,
     '        '>>Optimal code not implemented for DATA',
     '        ERROR,*9999)

          ENDIF


          CALL ASSERT(NJT.EQ.3,'>>ERROR: Magnetic fields are 3D only',
     '      ERROR,*9999)

          IF(ANALYTIC) THEN
            DO nqe=1,NQET(1)
              nq=NQNE(ne,nqe)
              CALL MAGSOLANALYTIC(NPNODE(0,nr),DIPOLE_CEN(1,0,1,nr,nx),
     '          DIPOLE_DIR(1,0,1,nr,nx),GRADPHI,XP,XQ(1,nq),MAGDIPOLE,
     '          MAGVOLUME,VPOTENTIAL,ERROR,*9999)
              DO nj=1,NJT
                YQ(nq,iy+nj,1,nxgrid)=GRADPHI(nj)
              ENDDO !nj
            ENDDO !nqe

          ELSE IF(OPTIML.AND.(.NOT.VPOTENTIAL)) THEN
C***        Reordered and a bit more unrolled to help speed and
C           parallelism.
            nhx=1
            nh=NH_LOC(nhx,nx)

C***        Check to see if the BE problem is a coupled problem and
C           sort out which region boundaries need to be dealt with
            IF(IS_COUPLED(nx)) THEN
              DO nrr=1,COUP_NRLIST(0,nx)
                NRLIST_LOCAL(nrr)=COUP_NRLIST(nrr,nx)
              ENDDO !nrr
              NRLIST_LOCAL(0)=COUP_NRLIST(0,nx)
            ELSE
              NRLIST_LOCAL(0)=1
              NRLIST_LOCAL(1)=NRLIST(1)
            ENDIF

            DO nqe=1,NQET(1)
C LKC 10-APR-2003 Missing the nq initialisation              
              nq=NQNE(ne,nqe)
              
              YQ(nq,iy+1,1,nxgrid)=0.0d0
              YQ(nq,iy+2,1,nxgrid)=0.0d0
              YQ(nq,iy+3,1,nxgrid)=0.0d0
            ENDDO

            DO nrr=1,NRLIST_LOCAL(0)
              nr=NRLIST_LOCAL(nrr)
              ERROR_FLAG=.FALSE.

C***          Transfer dependent variable to ZP (YPZP)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,MAX(NKH(nh,np,1,nr)-KTYP93(1,nr),1)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    ZP(nk,nv,nh,np,1)=YP(ny,1,nx)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nonode (np)

C$OMP         PARALLEL DO
C$OMP&        PRIVATE(DET_ADAPT,GRADPHI,nq,nqe,PG_J,PG_Q,PG_U,RAD_,RG_,
C$OMP&          XE_,XG1_,XIG_J,XIG_Q,XIG_U,XN_,XR_,YD_LOCAL,ZE_,ERROR),
C$OMP&        SHARED(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,NBH,
C$OMP&          NBJ,NDET,NDIPOLES,ne,NEELEM,NENP,NGAP,NKH,NKHE,NKJE,
C$OMP&          NP_INTERFACE,NLL,NPNE,NQNE,nr,NVHE,NVJE,NW,nx,CE,
C$OMP&          CURVCORRECT,DET,DIPOLE_CEN,DIPOLE_DIR,DL,LOC_TIME,PG,SE,
C$OMP&          WG,XIG,XP,XQ,ZP,MAGDIPOLE,MAGVOLUME,ERROR_FLAG)
              DO nqe=1,NQET(1)
                IF(.NOT.ERROR_FLAG) THEN
                  nq=NQNE(ne,nqe)
                  CALL MAGSOLOPTI(DIPOLE_CEN_NTIME(1,1,nx),
     '              DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBH,NBJ,NDET,
     '              NDIPOLES(1,nx),NEELEM,NENP,NGAP,NKH,NKHE,NKJE,
     '              NP_INTERFACE,NLL,NPNE,nr,NRLIST_LOCAL,NVHE,NVJE,
     '              NW(1,1,nx),nx,CE(1,1,nx),CURVCORRECT,DET,DET_ADAPT,
     '              DIPOLE_CEN(1,0,1,1,nx),DIPOLE_DIR(1,0,1,1,nx),DL,
     '              GRADPHI,LOC_TIME,PG,PG_J,PG_Q,PG_U,RAD_,RG_,SE,WG,
     '              XE_,XG1_,XIG,XIG_J,XIG_Q,XIG_U,XN_,XP,XQ(1,nq),XR_,
     '              YD_LOCAL,ZE_,ZP,MAGDIPOLE,MAGVOLUME,ERROR,*300)
                  DO nj=1,NJT
                    YQ(nq,iy+nj,1,nxgrid)=YQ(nq,iy+nj,1,nxgrid)+
     '                GRADPHI(nj)
                  ENDDO !nj

                  GO TO 302
C                 This statement is designed to be skipped if no
C                 errors occur. However if a error occurs within
C                 a subroutine the alternate return points to line
C                 300 and sets the flag
 300              CONTINUE
C$OMP             CRITICAL(EVSOLU_3)
                  ERROR_FLAG=.TRUE.
                  WRITE(OP_STRING,'(/'' >>ERROR: '',A)') ERROR
                  CALL WRITES(IOER,OP_STRING,ERROR,*301)
                  WRITE(OP_STRING,'(/'' >>An error occurred - '
     '              //'results may be unreliable!'')')
                  CALL WRITES(IOER,OP_STRING,ERROR,*301)
 301              CONTINUE
C$OMP             END CRITICAL(EVSOLU_3)
 302              CONTINUE
                ENDIF !error flag
              ENDDO !nqe
C$OMP         END PARALLEL DO
            ENDDO !nrr

C            WRITE(*,*) YQ(1,iy+1,1,nxgrid),YQ(1,iy+2,1,nxgrid),
C     '        YQ(1,iy+3,1,nxgrid)


          ELSE ! electric or magnetic
            
            IF(EVALTYPE(1:4).EQ.'GRID') THEN
              DO nqe=1,NQET(1)
                nq=NQNE(ne,nqe)
                CALL MAGSOL(DIPOLE_CEN_NTIME(1,1,nx),
     '            DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBH,NBJ,NDET,
     '            NDIPOLES(1,nx),NEELEM,NENP,NGAP,NHE(1,nx),NHP(1,0,nx),
     '            NKH,NKHE,NKJE,NP_INTERFACE,NLL,NPF,NPNE,NPNODE,NRLIST,
     '            NRE,NVHE,NVHP,NVJE,NW(1,1,nx),nx,NYNE,NYNP,CE(1,1,nx),
     '            CURVCORRECT,DET,DET_ADAPT,DIPOLE_CEN(1,0,1,1,nx),
     '            DIPOLE_DIR(1,0,1,1,nx),DL,GRADPHI,LOC_TIME,PG,
     '            PG_J,PG_Q,PG_U,RAD,RG,SE,WG,XA,XE,XG1,XIG,XIG_J,
     '            XIG_Q,XIG_U,XN,XP,XQ(1,nq),XR,YD_LOCAL,YP(1,1,nx),
     '            ZA,ZE,ZP,MAGDIPOLE,MAGVOLUME,VPOTENTIAL,ERROR,*9999)
                DO nj=1,NJT
                  YQ(nq,iy+nj,1,nxgrid)=GRADPHI(nj)
                ENDDO !nj

C                WRITE(*,*) ne,YQ(1,iy+1,1,nxgrid)
                
              ENDDO !nqe


            ELSEIF(EVALTYPE(1:3).EQ.'MFI') THEN
  
C LKC 16-APR-2002 Put in correct mapping between time step and real time
C              nt=INT(LOC_TIME)+1

              nt=CALC_SAMPLE_FROM_TIME(LOC_TIME,ERR,ERROR)
              CALL ASSERT(nt.LE.NTSM,'>>Increase NTSM',ERROR,
     &          *9999)
              
              IF(NDDATA(0,nr).EQ.0) THEN
                WRITE(OP_STRING,'(/'' >>WARNING: '',A)')
     '            'No data points to evaluate solution at'
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF

              ISIZE_MFI(1,nss)=NDDATA(0,nr) !electrodes
              IF(ISIZE_MFI(2,nss).LT.nt) THEN
                ISIZE_MFI(2,nss)=nt !time
              ENDIF
              ISIZE_MFI(3,nss)=NJT !components

C!!! debug debug

C              WRITE(OP_STRING,'('' Evaluating solution at index '',
C     '          I3,'' and time '',F8.2)') nt,LOC_TIME
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C              IF(DIPOLE_CEN_NTIME(1,1,nx).EQ.0) THEN
C                WRITE(*,*) 'EVSOLU: Static dipole'
C                WRITE(OP_STRING,'(''   Dipole Center: '',3D18.10)')
C     '            (DIPOLE_CEN(nj,0,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                WRITE(OP_STRING,'(''   Dipole vector: '',3D18.10)')
C     '            (DIPOLE_DIR(nj,0,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              ELSE
C                WRITE(*,*) 'EVSOLU: Time dependent dipole'
C                WRITE(OP_STRING,'(''   Dipole Center: '',3D18.10)')
C     '            (DIPOLE_CEN(nj,nt-1,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C                WRITE(OP_STRING,'(''   Dipole vector: '',3D18.10)')
C     '            (DIPOLE_DIR(nj,nt-1,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              ENDIF


C              WRITE(OP_STRING,              
C     &          '('' EVSOLU: at index '',I3,'' and time '',F8.2)') nt,
C     &          LOC_TIME
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C              IF(DIPOLE_CEN_NTIME(1,1,nx).EQ.0) THEN
C                WRITE(*,*) 'EVSOLU: Static dipole'
C                WRITE(OP_STRING,'(''   Dipole Center: '',3F15.9)')
C     '            (DIPOLE_CEN(nj,0,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                WRITE(OP_STRING,'(''   Dipole vector: '',3F15.9)')
C     '            (DIPOLE_DIR(nj,0,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              ELSE
C                WRITE(*,*) 'EVSOLU: Time dependent dipole'
C                WRITE(OP_STRING,'(''   Dipole Center: '',3F15.9)')
C     '            (DIPOLE_CEN(nj,nt-1,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C
C                WRITE(OP_STRING,'(''   Dipole vector: '',3F15.9)')
C     '            (DIPOLE_DIR(nj,nt-1,1,nr,nx),nj=1,NJT)
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C              ENDIF

              

C*** HB 11-NOV-2004 Calculate the coil solution if coil is stated

C*** If the parameter 'coil' is specified, it will calculate the magnetic 
C*** field or the magnetic field gradient (if 'baseline' is specified) due 
C*** to a circular current flowing in a coil placed on the surface of the 
C*** body. This is used to validate SQUID measurements by comparing them 
C*** with those generated from a non-magnetic field.
C*** As the magnetic field is symmetric about the coil axes, by assuming 
C*** that the coils are all at the centre of the SQUID we calculate the 
C*** magnetic field at a distance specified by the distance (radius) between 
C*** the centre of the SQUID and the coil.
C*** Note that the coil positions are required to be entered into an .ipdata 
C*** file and that the the coil solution uses the full squid configuration
C*** of 19 squid data points

              IF(COIL) THEN

C***  Read in the dewar centre data
                DO nj=1,3
                  SQUID_PT1(nj)=ZD(nj,1)
                ENDDO

                IF (BASELINE.GT.ZERO_TOL) THEN
C***  Calculate the coil height at the offset height
                  HEIGHT_COIL2 = HEIGHT_COIL + BASELINE
                ENDIF

C***  For each coil sensor
                DO nd=1,NDDATA(0,nr)
       
C***  Read in the coil co-ordinates
                  DO nj=1,3
                    ZD_TMP(nj)=ZD(nj,nd)                
                  ENDDO

C***  Calculate the magnetic field at the sensor position
                  CALL MAGSOLCOIL(BFIELD,SQUID_PT1,ZD_TMP,
     '              RAD_COIL,CURRENT,HEIGHT_COIL,ERROR,*9999)

C***  If a baseline is specified, calculate the magnetic gradient between
C***  the field specified at the coil height and the field removed by a 
C***  distance BASELINE from that height.

                  IF (BASELINE.GT.ZERO_TOL) THEN

                    CALL MAGSOLCOIL(BFIELD_2,SQUID_PT1,
     '                ZD_TMP,RAD_COIL,CURRENT,HEIGHT_COIL2,
     '                ERROR,*9999)

C***  Write out the magnetic field gradient solution to MFI

                    DO nj=1,3
                      MFI(nd,nt,nj,nss)=BFIELD(nj)-BFIELD_2(nj)
                    ENDDO !nj
  
                  ELSE

C***  Write out the magnetic field solution to MFI
                    DO nj=1,3
                      MFI(nd,nt,nj,nss)=BFIELD(nj)
                    ENDDO !nj

                  ENDIF            
                ENDDO !nd
                  
              ELSE ! .NOT.COILS

                DO nd=1,NDDATA(0,nr)

                  IF(.NOT.GRADIENT) THEN
                    CALL MAGSOL(DIPOLE_CEN_NTIME(1,1,nx),
     '                DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBH,NBJ,NDET,
     '                NDIPOLES(1,nx),NEELEM,NENP,NGAP,NHE(1,nx),NHP(1,0,
     '              nx),NKH,NKHE,NKJE,NP_INTERFACE,NLL,NPF,NPNE,NPNODE,
     '              NRLIST,NRE,NVHE,NVHP,NVJE,NW(1,1,nx),nx,NYNE,NYNP,
     '              CE(1,1,nx),CURVCORRECT,DET,DET_ADAPT,DIPOLE_CEN(1,0,
     '              1,1,nx),DIPOLE_DIR(1,0,1,1,nx),DL,GRADPHI,LOC_TIME,
     '              PG,PG_J,PG_Q,PG_U,RAD,RG,SE,WG,XA,XE,XG1,XIG,XIG_J,
     '              XIG_Q,XIG_U,XN,XP,ZD(1,nd),XR,YD_LOCAL,YP(1,1,nx),
     '              ZA,ZE,ZP,MAGDIPOLE,MAGVOLUME,VPOTENTIAL,ERROR,*9999)
 
                  ELSE

C*** Absolute MFI calculates in the first direction
                    DO nj=1,3
                      ZD_TMP(nj)=ZD(nj,nd)
                    ENDDO
                    ZD_TMP(2)=ZD_TMP(2)
                    
                    CALL MAGSOL(DIPOLE_CEN_NTIME(1,1,nx),
     '                DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBH,NBJ,NDET,
     '                NDIPOLES(1,nx),NEELEM,NENP,NGAP,NHE(1,nx),NHP(1,0,
     '                nx),NKH,NKHE,NKJE,NP_INTERFACE,NLL,NPF,NPNE,
     &                NPNODE,NRLIST,NRE,NVHE,NVHP,NVJE,NW(1,1,nx),nx,
     &                NYNE,NYNP,CE(1,1,nx),CURVCORRECT,DET,DET_ADAPT,
     &                DIPOLE_CEN(1,0,1,1,nx),DIPOLE_DIR(1,0,1,1,nx),DL,
     &                GRADPHI,LOC_TIME,PG,PG_J,PG_Q,PG_U,RAD,RG,SE,WG,
     &                XA,XE,XG1,XIG,XIG_J,XIG_Q,XIG_U,XN,XP,ZD_TMP,XR,
     &                YD_LOCAL,YP(1,1,nx),ZA,ZE,ZP,MAGDIPOLE,MAGVOLUME,
     &                VPOTENTIAL,ERROR,*9999)

                  
C*** Absolute MFI calculates in the second direction.
                    ZD_TMP(2)=ZD_TMP(2)-BASELINE
                    
                    CALL MAGSOL(DIPOLE_CEN_NTIME(1,1,nx),
     '                DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBH,NBJ,NDET,
     '                NDIPOLES(1,nx),NEELEM,NENP,NGAP,NHE(1,nx),NHP(1,0,
     '                nx),NKH,NKHE,NKJE,NP_INTERFACE,NLL,NPF,NPNE,
     &                NPNODE,NRLIST,NRE,NVHE,NVHP,NVJE,NW(1,1,nx),nx,
     &                NYNE,NYNP,CE(1,1,nx),CURVCORRECT,DET,DET_ADAPT,
     &                DIPOLE_CEN(1,0,1,1,nx),DIPOLE_DIR(1,0,1,1,nx),DL,
     &                GRADPHI2,LOC_TIME,PG,PG_J,PG_Q,PG_U,RAD,RG,SE,WG,
     &                XA,XE,XG1,XIG,XIG_J,XIG_Q,XIG_U,XN,XP,ZD_TMP,XR,
     &                YD_LOCAL,YP(1,1,nx),ZA,ZE,ZP,MAGDIPOLE,MAGVOLUME,
     &                VPOTENTIAL,ERROR,*9999)

C*** Simple difference over the baseline. I think it should be larger
C*** value less the smaller value (closer field - the further field)
                    DO nj=1,NJT
                      GRADPHI(nj)=(GRADPHI(nj)-GRADPHI2(nj))
                    ENDDO
                    
                  ENDIF !Gradient
       
C*** Assign the appropriate computed magnetic fields (absolute or
C*** gradients) according to the SQUID_CONFIG specifications
               
                  IF(SQUID_CONFIG.EQ.0) THEN
                    
                    DO nj=1,NJT
                      MFI(nd,nt,nj,nss)=GRADPHI(nj)
                    ENDDO !nj
                    
                  ELSEIF(SQUID_CONFIG.EQ.1) THEN ! Vanderbilt SQUID
                  
C!!! LKC 17-FEB-2003
C NOTE: we are using the "cmiss" electrode numbers here -
C There should be a a total of 19 electrodes with the central electrode
C number 1 and numbered going clockwise and outwards
C                  
C This could be improved by renaming the five 3 component sensors 
C to be at the start or end of the data list -- worry 
C about this when changing the code to use the optimised code.

                  CALL ASSERT(NDDATA(0,nr).EQ.19,'>> Should only be '
     '              //'19 electrodes when using SQUID option',
     '              ERROR,*9999)

C LKC 25-MAR-2002 Tidyup to use INLIST
C                   IF(nd.EQ.8.OR.nd.EQ.11.OR.nd.EQ.14.OR.nd.EQ.17) THEN
                  IF(INLIST(nd,MULTICHANNEL,5,ng)) THEN
                    DO nj=1,NJT
                      MFI(nd,nt,nj,nss)=GRADPHI(nj)
                    ENDDO !nj
                  ELSE

C Main "z" component channels
C Probably don't need to reinitialise 1 and 3 components ... 
C NOTE: cmiss y == vanderbilt z component
                    MFI(nd,nt,1,nss)=0.d0
                    MFI(nd,nt,2,nss)=GRADPHI(2)
                    MFI(nd,nt,3,nss)=0.d0
                  ENDIF !INLIST

                ELSEIF(SQUID_CONFIG.EQ.2) THEN ! "y" components for all channels
                  
                  MFI(nd,nt,1,nss)=0.d0
                  MFI(nd,nt,2,nss)=GRADPHI(2)
                  MFI(nd,nt,3,nss)=0.d0

                ELSE
                  ERROR='Unknown SQUID Configuration'
                  GOTO 9999
                ENDIF !SQUID_CONFIG

              ENDDO !nd
           ENDIF  !Coil

C !!! debug debug
C              nd=1
C              WRITE(OP_STRING,'('' EVSOLU  nd: '',I2,
C     '          '' MFI : '',3E20.8)') nd,(MFI(nd,nt,nj,nss),nj=1,3)
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              
            ELSEIF(EVALTYPE(1:4).EQ.'DATA') THEN

              IF(GRADIENT) THEN
                ERROR='Gradient calculations not calculated for DATA'
                GOTO 9999
              ENDIF
              
              IF(NDDATA(0,nr).EQ.0) THEN
                WRITE(OP_STRING,'(/'' >>WARNING: '',A)')
     '            'No data points to evaluate solution at'
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
              DO nd=1,NDDATA(0,nr)
                CALL MAGSOL(DIPOLE_CEN_NTIME(1,1,nx),
     '            DIPOLE_DIR_NTIME(1,1,nx),IBT,IDO,INP,NBH,NBJ,NDET,
     '            NDIPOLES(1,nx),NEELEM,NENP,NGAP,NHE(1,nx),NHP(1,0,nx),
     '            NKH,NKHE,NKJE,NP_INTERFACE,NLL,NPF,NPNE,NPNODE,NRLIST,
     '            NRE,NVHE,NVHP,NVJE,NW(1,1,nx),nx,NYNE,NYNP,CE(1,1,nx),
     '            CURVCORRECT,DET,DET_ADAPT,DIPOLE_CEN(1,0,1,1,nx),
     '            DIPOLE_DIR(1,0,1,1,nx),DL,GRADPHI,LOC_TIME,PG,
     '            PG_J,PG_Q,PG_U,RAD,RG,SE,WG,XA,XE,XG1,XIG,XIG_J,
     '            XIG_Q,XIG_U,XN,XP,ZD(1,nd),XR,YD_LOCAL,YP(1,1,nx),
     '            ZA,ZE,ZP,MAGDIPOLE,MAGVOLUME,VPOTENTIAL,ERROR,*9999)
                
                DO nj=1,NJT
                  ZD(NJT+nj,nd)=GRADPHI(nj)
                  WD(NJT+nj,nd)=1.d0
                ENDDO !nj
              ENDDO

            ENDIF !EVALTYPE

          ENDIF

          CALL CPU_TIMER(CPU_USER,TIME_STOP)
          ELAPSED_TIME=TIME_STOP(1)-TIME_START(1)

          IF(VIEW) THEN
            WRITE(OP_STRING,'('' Time '',E12.4)')  ELAPSED_TIME
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF


          IF(DOP) THEN
            IF(EVALTYPE(1:4).EQ.'GRID'.OR.ELECTRIC) THEN
              IF(NJT.EQ.2) THEN
                DO nqe=1,NQET(1)
                  nq=NQNE(ne,nqe)
                  WRITE(OP_STRING,
     '              '('' Grid '',I8,'' potential'',4E12.4)')
     '              nq,YQ(nq,iy,1,nxgrid),XQ(1,nq),XQ(2,nq),
     '              (XQ(1,nq)**2-XQ(2,nq)**2)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO !nqe
              ELSE
                DO nqe=1,NQET(1)
                  nq=NQNE(ne,nqe)
                  WRITE(OP_STRING,
     '              '('' Grid '',I8,'' potential'',4E12.4)')
     '              nq,YQ(nq,iy,1,nxgrid),XQ(1,nq),XQ(2,nq),
     '              (XQ(1,nq)**2+XQ(2,nq)**2-2.0d0*XQ(3,nq)**2)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDDO !nqe
              ENDIF !NJT
            ENDIF
          ENDIF !DOP

        ENDIF !electric or magnetic based solution

      ENDIF

      CALL EXITS('EVSOLU')
      RETURN
 9999 CALL ERRORS('EVSOLU',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('EVSOLU')
      RETURN 1
      END


