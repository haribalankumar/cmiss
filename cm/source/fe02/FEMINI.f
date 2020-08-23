      SUBROUTINE FEMINI(ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,LD,NBJ,
     '  NDIPOLES,NEELEM,NFFACE,NHP,NHQ,NKH,NLLINE,NORD,NPF,NPNE,NPNODE,
     &  NQET,NTIME_NR,NVHP,NVJP,NWP,NYNO,NYNR,CE,ERROR,*)

C#### Subroutine: FEMINI
C###  Description:
C###    FEMINI initialises essential logicals and variables.

C**** AJP 7/7/95  Removed most array initalisation.  These should be
C**** done when the arrays are set up.

      IMPLICIT NONE
      INCLUDE 'alig00.cmn'
      INCLUDE 'anal00.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'fbgr00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'hist00.cmn'
      INCLUDE 'host00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'inver00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'linc00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'maqloc00.cmn'
      INCLUDE 'nqloc00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'scal00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'sour00.cmn'
      INCLUDE 'tol00.cmn'
      INCLUDE 'trsf00.cmn'

!     Parameter List
      INTEGER ISIZE_MFI(3,NSSM),ISIZE_PHI(2),ISIZE_PHIH(2),
     '  ISIZE_TBH(2),LD(NDM),NBJ(NJM,NEM),NDIPOLES(NRM,NXM),
     '  NEELEM(0:NE_R_M,0:NRM),NFFACE(0:NF_R_M,NRM),NHP(NPM,0:NRM,NXM),
     '  NHQ(NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NLLINE(0:NL_R_M,0:NRM),
     &  NORD(5,NE_R_M),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     &  NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),NTIME_NR(0:NTIMEVARSM,NRM),
     &  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),NYNO(0:NYOM,NOOPM,NRCM,
     &  0:NRM,NXM),NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NWP(NPM,2)

      REAL*8 CE(NMM,NEM,NXM) !new AJP 20/3/97
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NR_MX
      PARAMETER (NR_MX=9)
      INTEGER i,j,nb,nc,nd,ne,nf,nj,njj1,njj2,nh,nhx,nhost,
     '  niq,nm,nn,no,np,nq,nr,nrc,nx
      REAL*8 DLAMCH

      CALL ENTERS('FEMINI',*9999)

      ZERO_TOL=DLAMCH('EPS')*5.0d0
      CONVERG_TOL=DLAMCH('EPS')*5.0d0
      LOOSE_TOL=DSQRT(DLAMCH('EPS'))

C *** Check array dimensions are large enough in geom00.cmn & ityp00.cmn
      CALL ASSERT(NBM.LE.99,'Inc size of NBM dep arrays in geom00.cmn',
     '  ERROR,*9999)
      CALL ASSERT(NRM.LE.99,'Inc size of NRM dep arrays in geom00.cmn',
     '  ERROR,*9999)
      CALL ASSERT(NRM.LE.NR_MX,'Inc size of NRM dep arrays in '
     '  //'ityp00.cmn',ERROR,*9999)
      CALL ASSERT(NXM.LE.9,'Inc size of NXM dep arrays in ityp00.cmn',
     '  ERROR,*9999)
      CALL ASSERT(NCM.LE.11,'Increase first dim of NYT in geom00.cmn '
     '  //'to NCM',ERROR,*9999)
      CALL ASSERT(NCM.LE.4,'Increase dim of KTYP93 in ktyp90.cmn '
     '  //'to NCM',ERROR,*9999)
      CALL ASSERT(NMAQM.LE.99,'Inc dimension of arrays in maqloc00.cmn',
     '  ERROR,*9999)

      ALIGNMENT_ON=.FALSE.
      CALL_ACTI=.FALSE.
      CALL_AERO=.FALSE.
      CALL_ANAL=.FALSE.
      CALL_ANAL_TRSF=.FALSE.
      CALL_BASE=.FALSE.
      CALL_CALC_LAPL=.FALSE.
      CALL_CELL=.FALSE.
C LKC 16-JUN-2000 was uninitialised
      CALL_CELL_MATE=.FALSE.
C*** 08/10/08 JHC was uninitialised
      CALL_CONT=.FALSE.
      CALL_COUP=.FALSE.
      CALL_DATA=.FALSE.
      CALL_ELEM=.FALSE.
      CALL_ELFB=.FALSE.
      CALL_ELFD=.FALSE.
      CALL_EQUA=.FALSE.
      CALL_EXPO=.FALSE.
      CALL_FACE=.FALSE.
      CALL_FIBR=.FALSE.
      CALL_FIEL=.FALSE.
      CALL_FIT =.FALSE.
      CALL_GROW=.FALSE.
C LKC 16-JUN-2000 was uninitialised
      CALL_GMAP=.FALSE.
      CALL_INIT=.FALSE.
      CALL_INVE=.FALSE.
      CALL_IMPO=.FALSE.
C LKC 16-JUN-2000 was uninitialised
      CALL_LEAD=.FALSE.
      CALL_LINE=.FALSE.
      CALL_MATE=.FALSE.
      CALL_MAPPING=.FALSE.
      CALL_MESH=.FALSE.
      CALL_MOTI=.FALSE.
      CALL_NODE=.FALSE.
      CALL_NOIS=.FALSE.
      CALL_OBJE=.FALSE.
      CALL_OPTI=.FALSE.
      CALL_SOLV=.FALSE.
      CALL_TIME=.FALSE.
      CALL_TREE=.FALSE.
      CALL_UPGAUS_EIK=.FALSE.
      CALL_UPNODE_MAT=.FALSE.
      CALC_SHEET=.FALSE.
      CALC_XI=.FALSE.

      CELL_SPATIALLY_VARYING=.FALSE.
      FBGRAF=.FALSE.
      DO nx=1,NXM
        IS_COUPLED(nx)=.FALSE.
        SOLVELU_NZA(nx)=0
      ENDDO
      LFTRAN=.FALSE.
      RHTRAN=.FALSE.
      IF(USE_BEM.NE.0) THEN
        ADAPINT=.FALSE.
        COMPLEX=.FALSE.
        HERMITE=.FALSE.
        HYP=.FALSE.
      ENDIF
      IF(USE_DATA.NE.0) THEN
        CALL_DATA_FIBRE=.FALSE.
        CALL_DATA_FIELD=.FALSE.
        CALL_DATA_SHEET=.FALSE.
C*** 08/10/08 JHC Added new variable when readin in contact pressure from a data file
        CALL_DATA_CONT=.FALSE.
      ENDIF
      IF(USE_GRID.NE.0) THEN
        CALL_GRID=.FALSE.
        CALL_GMAP=.FALSE.
      ENDIF
      IF(USE_TRANSFER.NE.0) THEN
        APPLY_INVERSE=.FALSE.
        APPLY_TRANSFER=.FALSE.
        EVALUATE_INVERSE=.FALSE.
        EVALUATE_TRANSFER=.FALSE.
      ENDIF
      CALL_DEREFI=.FALSE.
      MAKE_DATA=.FALSE.

      IF(USE_TIME.NE.0) THEN
        EVALUATE_PHI=.FALSE.
      ENDIF
C     Initialize parameters for scale factor iterative calculation.
      LINC_ITMAX=50
      LINC_TOL=CONVERG_TOL
      
      DO nb=1,NBM
        NAT(nb)=0
      ENDDO
      NBT=0
      DO nr=1,NRM
        NTIME_NR(0,nr)=0 !initialise time-dep boundary condition variables per region
        DO nx=1,NXM
          NCT(nr,nx)=0
C LKC 20-JUN-2001 Also initialise NYNR so that lippara will work
C         without problems being setup
          DO nc=1,NCM
            NYNR(0,0,nc,nr,nx)=0
            NYNR(0,0,nc,0,nx)=0
          ENDDO
        ENDDO
      ENDDO
      NDT=0
      NET(0)=0
      NFT=0
      DO nb=1,NBM
        NGT(nb)=0
        NIT(nb)=0
        NKT(0,nb)=0
      ENDDO
      NLT=0
      DO nb=1,NBM
        NNT(nb)=0
      ENDDO
      DO nx=1,NXM
        DO nr=0,NRT
          DO nrc=1,NRCM
            DO nc=1,4
              NOT(nrc,nc,nr,nx)=0
              NOQT(nrc,nc,nr,nx)=0
            ENDDO
          ENDDO !nrc
        ENDDO !nr
      ENDDO !nx
      NPT(0)=0
      NQT=0
      NRT=1
      DO nb=1,NBM
        NST(nb)=0
        NUT(nb)=0
      ENDDO
      DO nx=1,NXM
        DO nrc=1,NRCM
          DO nc=1,NCM
            NYT(nrc,nc,nx)=0
            NYQT(nrc,nc,nx)=0
          ENDDO
        ENDDO
      ENDDO
      NJT=2
      NBFT=0
      NQSCT=0
      DO nq=1,NQSCM
        NQET(nq)=0
      ENDDO !nq
! also need to initialise NYNO for ippara
      DO no=1,NOM
        DO nr=0,NRT
          DO nrc=1,NRCM
            DO nx=1,NXM
              NYNO(0,no,nrc,nr,nx)=0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO no=1,99
        TOTALNQ(no)=0
        TOTALNY(no)=0
      ENDDO
      NTIMEPOINTST=0
      NTIMEVARST=0

C cpb 19/4/93 Initialise NJ_LOC
      DO nr=0,9
        DO njj1=0,3
          DO njj2=0,NJ_LOC_MX
            NJ_LOC(njj1,njj2,nr)=0
          ENDDO !njj2
        ENDDO !njj1
      ENDDO
      DO njj1=1,3*NJ_LOC_MX
        DO njj2=1,2
          NJ_TYPE(njj1,njj2)=0
        ENDDO !njj2
      ENDDO !njj1
C CPB 8/4/94 Adding NJ_LOC
      DO nj=1,NJT
        NJ_LOC(NJL_GEOM,nj,1)=nj
      ENDDO
      NJ_LOC(NJL_GEOM,0,1)=NJT
      NJ_LOC(NJL_GEOM,0,0)=NJT
      NJ_LOC(0,0,1)=NJT
      DO nj=1,NJT
        NJ_TYPE(nj,1)=NJL_GEOM
        NJ_TYPE(nj,2)=nj
      ENDDO
      NJ_LOC(0,0,0)=NJT

      DO nr=1,NRM
        DO nx=1,NXM
          ITYP1(nr,nx)=1
          ITYP2(nr,nx)=0
          ITYP3(nr,nx)=0
          ITYP4(nr,nx)=0
          ITYP5(nr,nx)=0
          ITYP6(nr,nx)=0
          ITYP7(nr,nx)=0
          ITYP8(nr,nx)=0
          ITYP9(nr,nx)=0
          ITYP21(nr)=1
          IF(USE_GRID.GT.0) NHQ(nr,nx)=0
        ENDDO !nx
        ITYP10(nr)=1
        ITYP11(nr)=1
      ENDDO !nr

      JTYP1=1
      JTYP2A=0
      JTYP2B=0
      JTYP2C=0
      JTYP3=0
      JTYP4=1
      JTYP5=1
      JTYP6=1
      JTYP7=0
      JTYP8=0
      JTYP9=0
      JTYP10=1
      JTYP11=0
      JTYP12=1
      FOCUS=1.0d0
      ZMINI=0.0d0
      ZMAXI=1.0d0
      ZDIFF=1.0d0
      KTYP4=0
      DO nr=1,NRM
        KTYP50(nr)=0
        KTYP51(nr)=0
        KTYP52(nr)=0
        KTYP53(nr)=0
        KTYP54(nr)=0
        KTYP55(nr)=0
        KTYP56(nr)=0
        KTYP57(nr)=0
        KTYP58(nr)=0
        KTYP59(nr)=0
        KTYP5A(nr)=0
        KTYP5B(nr)=0
        KTYP5C(nr)=0
        KTYP5D(nr)=0
        KTYP5E(nr)=0
        KTYP5F(nr)=0
        KTYP5G(nr)=0
        KTYP5H(nr)=0
        KTYP5I(nr)=0
      ENDDO
      KTYP90=0
      KTYP94=0
      DO nr=1,NRM
        DO nc=1,NCM
          KTYP93(nc,nr)=0
        ENDDO !nc
      ENDDO !nr
      KTYP11=0

c cpb 17/3/95 Adding NH_LOC and NX_LIST
      NX_LIST(0)=0
      NH_LOC(0,0)=0
      DO nx=1,9
        NX_CLASS(nx)=0
        NX_LIST(nx)=0
        NX_LOCKS(nx)=NX_NOLOCK
        NX_TYPE(nx)=0
        DO nhx=0,NH_LOC_MX
          NH_LOC(nhx,nx)=0
        ENDDO !nhx
      ENDDO !nx
      DO i=1,2
        DO j=1,19
          NH_TYPE(j,i)=0
        ENDDO !j
      ENDDO !i

C MLB 15/7/97 Adding NIQ_ arrays
      NIQ_LIST(0)=0
      MAQ_LIST(0)=0
      DO niq=1,99
        NIQ_LIST(niq)=0
        NIQ_LOCKS(niq)=0
        NIQ_SUBTYPE(niq)=0
        NIQ_TYPE(niq)=0
      ENDDO !niq

      NUM_DIP_TIMES=0
      DO nr=1,NRM
        DO nx=1,NXM
          NDIPOLES(nr,nx)=0
        ENDDO
      ENDDO

C MLB 24/5/99 initialise Salu in common block
      DO nx=1,NXM
        SALU_CONSISTENCY(nx)=.FALSE.
      ENDDO

C!!!  Temporary MPN 14-Feb-96: until a better place is found for this
C     Initialise arrays for zero'th region
      DO np=1,NPM
        DO nj=1,NJM
          NVJP(nj,np)=0
        ENDDO
        DO nc=1,NCM
          DO nh=1,NHM
            NKH(nh,np,nc,0)=0

C!!! LKC 20-AUG-1999 Also initialise the rest of NVHP
C!!!  this is required for fitting problems. NVHP is
C!!!  setup for the fitted regions yet the other
C!!!  regions are referenced in CALC_NY_MAPS_IND
C!!!  when fitting the region
C            NVHP(nh,np,nc,0)=0
            DO nr=0,NRM
              NVHP(nh,np,nc,nr)=0
            ENDDO

          ENDDO !nh
        ENDDO !nc
        DO nx=1,NXM
          NHP(np,0,nx)=0
        ENDDO !nx
      ENDDO !np
      DO nr=1,NRM
        NPNODE(0,nr)=0
        NEELEM(0,nr)=0
        NFFACE(0,nr)=0
        NLLINE(0,nr)=0
      ENDDO
      NPNODE(0,0)=0
      NEELEM(0,0)=0
      NLLINE(0,0)=0
C!!!  Temporary initialise NBJ because not region dependent
      DO ne=1,NEM
        DO nj=1,NJM
          NBJ(nj,ne)=0
        ENDDO
      ENDDO

!news AJP 20-Mar-97
      DO nm=1,NMM
        DO ne=1,NEM
          DO nx=1,NXM
            CE(nm,ne,nx)=0.0d0
          ENDDO
        ENDDO
      ENDDO
!newe AJP 20-Mar-97

!news MPN 13-Sep-94:
C     Initialise host info for parallel processing (in host00.cmn)
      DO nhost=1,N_HOSTS_MX
C       If the socket is open, close it.
        IF(SOCKET_OPEN(nhost)) THEN
          IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,ICONNID(nhost)).EQ.-1)
     '      GOTO 9999 !Signal to slave processes to stop
C IRET never used
C         IRET=FSKCLOSE(ICONNID(nhost)) !close socket
          SOCKET_OPEN(nhost)=.FALSE.
        ENDIF
      ENDDO
!newe

C     Free any memory used by the linear solvers
      DO NX=1,9
        CALL FREE_SOLVER(PRECON_CODE(NX),SOLVEROPTION(NX),SPARSEGKK(NX),
     '    NX,.TRUE.,ERROR,*9999)
      ENDDO

      DO i=1,99
        IWKS(i)=0
      ENDDO
      IWKDEF(0)=1
      IWKDEF(1)=1
      XMIN=-1.0
      XMAX= 1.0
      YMIN=-1.0
      YMAX= 1.0
      ZMIN=-1.0
      ZMAX= 1.0
      DIAG=SQRT(12.0)


C LKC 11-MAR-2002
C*** Initialise ISIZE_MFI, initialisation of ISIZE_PHI will come into
C*** this loop when the NSSM dependence is included for the PHI array.
      DO i=1,3
        DO j=1,NSSM
          ISIZE_MFI(i,j)=0
        ENDDO
      ENDDO

      DO i=1,2
        ISIZE_PHI(i)=0
        ISIZE_PHIH(i)=0
        ISIZE_TBH(i)=0
      ENDDO
      DO nd=1,NDM*USE_DATA
        LD(nd)=0
      ENDDO

C MLB 4 June 1997
C Added NPNE,NPF initialisations here because they can be set up
C in ipelem,ipface,ipmesh etc, and each can be called mulitple times.
      DO nf=1,NFM
        DO i=1,9
          NPF(i,nf)=0
        ENDDO
      ENDDO
      DO ne=1,NEM
        DO nb=1,NBFM
          DO nn=1,NNM
            NPNE(nn,nb,ne)=0
          ENDDO
        ENDDO
      ENDDO
C End 4 June 1997


C KAT 6May99: groups
      NTGRDA=0
      NTGREL=0
      NTGRFA=0
      NTGRGR=0
      NTGRLI=0
      NTGRNO=0


C LKC 1-DEC-1999 Initialise IREGULARISE in case deinve
C   has not been called for activation problems
      ACTN_IREGULARISE=1 ! No extra constraint.
      ACTN_REG_PARAM_LAPLACE=0.0d0
      IREGULARISE=1 ! No extra constraint.
C GBS 25-JUL-2000 Set default subroutine to use for RESFUN calculation
      ACTN_IRESFUN=1
C GBS 26 July 2000 True if transfer matrix T_BH is valid (APTRSF)
      TBH_OK=.FALSE.
C GBS 9-Aug-2000 By default, PHI_H refers to heart nodes
      PHI_H_IS_HEART=.TRUE.
C GBS 9-Aug-2000 True if electrodes are mapped to nodes (filling NDDATA)
      EVALUATE_PHI_NEAREST=.FALSE.

C LKC 18-APR-2000 adding new parameters for transfer jump
      TRSF_ACTN_WAVE_JUMP=1.0d0
      TRSF_ACTN_WAVE_REST=0.0d0
      TRSF_ACTN_WAVE_WIDTH=1.0d0
      TRSF_ACTN_WAVE_INTERPOLATE=.FALSE.
      TRSF_FREQUENCY=1.d0

C LKC 16-APR-2002 time for solvers
      T0=0.d0
      T1=0.d0


C CS 2/8/2001 initialse NWP
      DO np=1,NPM
        NWP(np,1)=0
      ENDDO !np
      IF(JTYP2C.EQ.1) THEN !hanging nodes
        DO np=1,NPNODE(0,1)
          NWP(np,1)=-1
        ENDDO !np
      ENDIF

      BEMCURVATURECORRECTION=.FALSE.

C KSB 1/07/2004 Initialise NORD
      IF(USE_LUNG.EQ.1) THEN
        DO j=1,NE_R_M
          DO i=1,5
            NORD(i,j)=0
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('FEMINI')
      RETURN
 9999 CALL ERRORS('FEMINI',ERROR)
      CALL EXITS('FEMINI')
      RETURN 1
      END


