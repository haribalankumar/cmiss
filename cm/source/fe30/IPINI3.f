      SUBROUTINE IPINI3(IBT,IDO,INP,ITHRES,NBH,NBHF,NBJ,NBJF,NEELEM,
     '  NEL,NENP,NENQ,NFF,NHE,NHP,NHQ,NKEF,NKH,NKHE,NKJE,NNF,NPF,
     &  NPLIST,NLL,NP_INTERFACE,NNL,NPL,NPNE,NPNF,NPNODE,NQET,NQLIST,
     &  NQNE,NQS,nr,NVHE,NVHF,NVHP,NVJE,NVJF,NW,NWQ,nx,NXLIST,NXI,
     '  NYNE,NYNP,NYNQ,NYNR,AQ,CE,CQ,CURVCORRECT,DL,PG,RCQS,RG,
     '  SE,SF,THRES,WG,XA,XE,XG,XIG,XP,XQ,YG,YP,YQ,YQS,ZA,ZE,ZP,
     '  TIME_VARIABLE_NAMES,ALL_REGIONS,FIX,FIXQ,FIX_ZERO,GENER,
     '  ERROR,*)

C#### Subroutine: IPINI3
C###  Description:
C###    IPINI3 inputs initial conditions and boundary conditions for
C###    FE30 problems.

C**** YP(ny,1) contains essential and flux  b.c.s,defined by FIX(ny,1).
C**** YP(ny,3)    "     initial solution.
C KAT 27Jan99: FIX(ny,5) not used
CC**** FIX(ny,5) is .true. for input from a file (FILE07)
C**** 23-9-92 AJP
C**** If point values are to be used for the flux or mixed bdry cond.
C**** flux then these are not stored in YP(ny,2) and the effective
C**** integrated values for these are calculated and stored in YP(ny,1).
C#### KTYPMBC=1 if the flux or mixed bc flux is an integrated value
C###         =2 if the flux or mixed bc flux is a point value
C**** GENER is .true. if initial conditions are generated automatically.
C**** Only implemented for coupled fe/be probs (assumes that FE domain
C**** is not the first region).
C**** 18-1-95 This routine needs checking to ensure it used the correct
C**** storage location in the FIX and YP arrays.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
C 20-MAR-1998     INCLUDE 'cmiss$reference:binf00.cmn'
      INCLUDE 'anal00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp90.cmn'
C 22-MAR-1998      INCLUDE 'cmiss$reference:mach00.cmn'
C 20-MAR-1998     INCLUDE 'cmiss$reference:mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lungex00.cmn'
      INCLUDE 'mesh03.cmn'
      INCLUDE 'mxbc00.cmn'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'pulm00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ITHRES(3,NGM,NEM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),NFF(6,NEM),
     '  NHE(NEM),NHP(NPM),NHQ,NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPLIST(0:NP_R_M),
     '  NLL(12,NEM),NP_INTERFACE(0:NPM,0:3),NNL(0:4,12,NBFM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),NQLIST(0:NQM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),nr,NVHE(NNM,NBFM,NHM,NEM),
     '  NVHF(NNM,NBFM,NHM),NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NW(NEM,3),NWQ(8,0:NQM,NAM),nx,NXLIST(0:NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNQ(NHM,NQM,0:NRCM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 AQ(NMAQM,NQM),CE(NMM,NEM),CQ(NMM,NQM),
     '  CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),
     '  RCQS(NQRM),RG(NGM),SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  THRES(3,NGM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM),XQ(NJM,NQM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TIME_VARIABLE_NAMES(NTIMEVARSM)*(*)
      LOGICAL ALL_REGIONS,FIX(NYM,NIYFIXM),FIXQ(NYQM,NIYFIXM,NXM),
     '  FIX_ZERO,GENER
!     Local Variables
      INTEGER ATIME,DPOT,i,IBEG,IBEG1,ICHAR,IEND,IEND1,INFO,
     '  ipulse,iy,loop,LOOPT,maqp1t0,maqp1t1,maqp1i,maqp2t0,maqp2t1,
     '  maqp2i,n,N1,na,nb,nbf,nc,ne,nef,NE_START,nf,ng,nh,nhx,niq,
     '  niq_bc,niq_old,niqV,niy,nj,nk,NKHF(NKM,NNM,NHM),
     '  NKJF(NKM,NNM,NJM),NMAX,nn,NN_TOT,noelem,nonode,no_nynq,no_nynr,
     '  NOQUES,NO_DATA,np,NPMC,nq,ns,NUMBC,NUMPULSE,NUMTIMES,nv,
     '  NVARINDEX,nxc,nx_ext,ny,ny_P,nyq,ny1,ny2,ny3,ny4,ny_first,
     '  ny1_first,nv_last
      REAL*8 beta,bc,COM1,COM2,dummy,DXIX(3,3),fo,GL(3,3),
     '  Go,G_TERM,GU(3,3),HEIGHT(NJT),INITIAL_HD,PALV_DEFLT,
     '  perimeter,R,Ro,PXI,RART,RCAP,RVIE,SUM
C      REAL*8 RELDIST,HS
      CHARACTER CHAR*6,CHAR1*20,CHAR2*15,CHAR3*18,CHAR4*1
      LOGICAL FACE_FIXED,FILEIP,FLUXBC,INIT_RESET,INLIST,MANUAL,
     '  NEXTNODE
C DPN 08 February 1999
      INTEGER kl,nl,nt
      REAL*8 COEFF(4)
      LOGICAL LINE_FIXED
!     External functions
      INTEGER IDIGITS

      CALL ENTERS('IPINI3',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.8) THEN
C ***   threshold model
        FORMAT='('' Enter initial activation points:'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        NMAX=MAX(NET(1),NGM)

        IF(iotype.ne.3) THEN
          DO ne=1,NEM
            DO ng=1,NGM
              YG(1,ng,ne)=999.999d0
              ITHRES(1,ng,ne)=0
              ITHRES(2,ng,ne)=0
              THRES(1,ng,ne)=0.0d0
              THRES(2,ng,ne)=0.0d0
              THRES(3,ng,ne)=0.0d0
            ENDDO
          ENDDO
 700      FORMAT=
     '      '($,'' Element and Gauss point numbers [exit]: '',2I4)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).ne.0) THEN
            ne=IDATA(1)
            ng=IDATA(2)
            FORMAT='($,'' Enter time of initial activation [0]: '','
     '        //'E12.4)'
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,
     '        LDATA,LDEFLT,RDATA,RZERO,RMIN,RMAX,INFO,ERROR,*9999)
            YG(1,ng,ne)=RDATA(1)
            GO TO 700
          ENDIF

        ELSE IF(IOTYPE.EQ.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ng=1,NGT(1)
              FORMAT=
     '          '($,'' Element and Gauss point numbers [exit]: '','
     '          //'2I4)'
              IDATA(1)=ne
              IDATA(2)=ng+18
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,2,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              FORMAT='($,'' Enter time of initial activation [0]: '','
     '          //'E12.4)'
              RDATA(1)=YG(1,ng,ne)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,
     '          LDATA,LDEFLT,RDATA,RZERO,RMIN,RMAX,INFO,ERROR,*9999)
            ENDDO
          ENDDO
          FORMAT='($,'' Element and Gauss point numbers [exit]: '','
     '      //'2I4)'
          IDATA(1)=0
          IDATA(2)=0
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ENDIF

        FORMAT='(/'' Enter Purkinje tissue Gauss points:'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        NMAX=MAX(NET(1),NGM)

        IF(iotype.ne.3) THEN
          DO ne=1,NEM
            DO ng=1,NGM
              ITHRES(2,ng,ne)=0
            ENDDO
          ENDDO
 800      FORMAT=
     '      '($,'' Element and Gauss point numbers [exit]: '',2I4)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).ne.0) THEN
            ne=IDATA(1)
            ng=IDATA(2)
            ITHRES(2,ng,ne)=1
            GO TO 800
          ENDIF

        ELSE IF(IOTYPE.EQ.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ng=1,NGT(1)
              FORMAT=
     '          '($,'' Element and Gauss point numbers [exit]: '','
     '          //'2I4)'
              IDATA(1)=ne
              IDATA(2)=ng+18
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,2,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            ENDDO
          ENDDO
          FORMAT='($,'' Element and Gauss point numbers [exit]: '','
     '      //'2I4)'
          IDATA(1)=0
          IDATA(2)=0
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ENDIF

      ELSE IF(ITYP2(nr,nx).EQ.9) THEN !cellular based modelling

C KAT 2001-04-18: leave initialization till we know which niq's require
C       initialization.
CC *** DPN 12 March 2001
CC ***   If solving any cell based model then need to initialise the YQ
CC       array
CC$OMP PARALLEL DO
CC$&     PRIVATE(nq,na,niq),
CC$&     SHARED(YQ)
C        DO nq=1,NQT
C          DO na=1,NAM
C            DO niq=1,NIQM
C              YQ(nq,niq,na,nx)=0.0d0
C            ENDDO !niq
C          ENDDO !na
C        ENDDO !nq
CC$OMP END PARALLEL DO

        IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
!electrical or coupled modelling

C Set up current injection conditions

          CALL ASSERT(CALL_CELL,'>>Must define cell first',ERROR,*9999)
          CALL ASSERT(NMAQM.GE.8,'>>Increase NMAQM, must be >= 8',
     '      ERROR,*9999)

C       Call maq_loc to get the three indicies for each pulse
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t0,MAQ_START,
     '      ERROR,*9999)
          IF(maqp1t0.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_I_PULSE1,maqp1t0,
     '        MAQ_START,ERROR,*9999)
          ENDIF

          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1t1,MAQ_STOP,
     '      ERROR,*9999)
          IF(maqp1t1.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_I_PULSE1,maqp1t1,
     '        MAQ_STOP,ERROR,*9999)
          ENDIF

          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE1,maqp1i,MAQ_CURRENT,
     '      ERROR,*9999)
          IF(maqp1i.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_I_PULSE1,maqp1i,
     '        MAQ_CURRENT,ERROR,*9999)
          ENDIF

          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t0,MAQ_START,
     '      ERROR,*9999)
          IF(maqp2t0.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_I_PULSE2,maqp2t0,
     '        MAQ_START,ERROR,*9999)
          ENDIF

          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2t1,MAQ_STOP,
     '      ERROR,*9999)
          IF(maqp2t1.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_I_PULSE2,maqp2t1,
     '        MAQ_STOP,ERROR,*9999)
          ENDIF

          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_I_PULSE2,maqp2i,MAQ_CURRENT,
     '      ERROR,*9999)
          IF(maqp2i.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_I_PULSE2,maqp2i,
     '        MAQ_CURRENT,ERROR,*9999)
          ENDIF

          IF(DOP) THEN
            WRITE(OP_STRING,'(''Allocated maqs  '',6I2)') maqp1t0,
     '        maqp1t1,maqp1i,maqp2t0,maqp2t1,maqp2i
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C       Initialise the AQ array for the activation parameters
          IF(IOTYPE.NE.3) THEN
            DO nq=1,NQT
              AQ(maqp1t0,nq)=0.0d0
              AQ(maqp1t1,nq)=0.0d0
              AQ(maqp1i,nq)=0.0d0
              AQ(maqp2t0,nq)=0.0d0
              AQ(maqp2t1,nq)=0.0d0
              AQ(maqp2i,nq)=0.0d0
            ENDDO
          ENDIF

C       Allocate an AQ index to store the grid point activation time
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,ATIME,MAQ_ACTIV_TIME,
     '      ERROR,*9999)
          IF(ATIME.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_TIME,ATIME,
     '        MAQ_ACTIV_TIME,ERROR,*9999)
          ENDIF

          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,DPOT,MAQ_M_DPOT,
     '      ERROR,*9999)
          IF(DPOT.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_TIME,DPOT,
     '        MAQ_M_DPOT,ERROR,*9999)
          ENDIF

C       Find extracellular nx
          IF(KTYP32.EQ.2) THEN !bidomain
            IF(NXLIST(0).GE.2) THEN
              nxc=NXLIST(2)
              CALL NX_LOC(NX_INQUIRE,nxc,nx_ext,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx_ext.GT.0,'>>No nx defined for solve class',
     '          ERROR,*9999)
            ELSE
              CALL ASSERT(.FALSE.,'>>You must specify 2 classes for'
     '          //' bidomain',ERROR,*9999)
            ENDIF
          ELSE
            nx_ext=0
          ENDIF
        ENDIF !ITYP19(nr,nx).EQ.1

        NMAX=NQM

C ***   If solving any real cellular model stimulus protcol is set
C ***   in the ipcell file. Otherwise, prompt for it here.
        IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).LT.5) THEN
          FORMAT='('' Enter initial activation points:'')'
          CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

          IF(IOTYPE.NE.3) THEN

 802        FORMAT=
     '        '($,'' Enter collocation point #s/name [EXIT]: '',I5)'
            CDATA(1)='GRIDS' !for use with group input
C KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
C               passed to GINOUT instead of IDATA.
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C
C          IF(IDATA(1).NE.0) THEN !not default exit
C            NQLIST(0)=IDATA(0)
C            DO n=1,IDATA(0)
C              NQLIST(n)=IDATA(n)
C            ENDDO !n
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,0,
     &        NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

            IF(NQLIST(1).NE.0) THEN !not default exit

C           Define b.c. or i.c. for first grid pt in group
              nq=NQLIST(1) !rest of group is filled in afterwards
              WRITE(CHAR,'(I6)') nq
              IDEFLT(1)=1
              FORMAT='($,'' Enter number of pulses (1-2) [1]: '',I4)'
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              NUMTIMES=IDATA(1)

              DO ipulse=1,NUMTIMES
                WRITE(CHAR1,'(I1)') ipulse
                FORMAT='($,'' Enter activation/deactivation times'//
     '            ' for pulse '//CHAR1(1:1)//' (ms) [0,0]: '',2E12.4)'
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,
     '            LDATA,LDEFLT,RDATA,RZERO,0.0d0,1.d9,INFO,ERROR,*9999)

                IF(ipulse.EQ.1) THEN
                  AQ(maqp1t0,nq)=RDATA(1)
                  AQ(maqp1t1,nq)=RDATA(2)
                ELSE
                  AQ(maqp2t0,nq)=RDATA(1)
                  AQ(maqp2t1,nq)=RDATA(2)
                ENDIF

                FORMAT='($,'' Enter current for pulse '//
     '            CHAR1(1:1)//' (uA/mm^3) [0.1D+03]: '',E12.4)'
                RDEFLT(1)=0.1d3
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,
     '            LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

                IF(ipulse.EQ.1) THEN
                  AQ(maqp1i,nq)=RDATA(1)
                ELSE
                  AQ(maqp2i,nq)=RDATA(1)
                ENDIF

              ENDDO

              DO n=2,NQLIST(0) !fill in rest of group
                nq=NQLIST(n)
                AQ(maqp1t0,nq)=AQ(maqp1t0,NQLIST(1))
                AQ(maqp1t1,nq)=AQ(maqp1t1,NQLIST(1))
                AQ(maqp1i,nq)=AQ(maqp1i,NQLIST(1))
                AQ(maqp2t0,nq)=AQ(maqp2t0,NQLIST(1))
                AQ(maqp2t1,nq)=AQ(maqp2t1,NQLIST(1))
                AQ(maqp2i,nq)=AQ(maqp2i,NQLIST(1))
              ENDDO
              GO TO 802

            ENDIF

          ELSE IF(IOTYPE.EQ.3) THEN
            DO nq=1,NQT
              IF(AQ(maqp1i,nq).GT.ZERO_TOL) THEN !At least 1 i.c specified
                FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I6)'
                IDATA(1)=nq
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,
     '            ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NMAX,
     '            LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

                FORMAT='($,'' Enter number of pulses (1-2) [1]: '',I4)'
                IF(DABS(AQ(maqp2i,nq)).LT.ZERO_TOL) THEN
                  NUMPULSE=1
                ELSE
                  NUMPULSE=2
                ENDIF
                IDATA(1)=NUMPULSE
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,
     '            ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,NMAX,
     '            LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

                DO ipulse=1,NUMPULSE
                  WRITE(CHAR1,'(I1)') ipulse
                  IF(ipulse.EQ.1) THEN
                    RDATA(1)=AQ(maqp1t0,nq)
                    RDATA(2)=AQ(maqp1t1,nq)
                  ELSE
                    RDATA(1)=AQ(maqp2t0,nq)
                    RDATA(2)=AQ(maqp2t1,nq)
                  ENDIF

                  FORMAT='($,'' Enter activation/deactivation times'//
     '              ' for pulse '//CHAR1(1:1)//' (ms) [0,0]: '',2E12.4)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,
     '              LDATA,LDEFLT,RDATA,RZERO,RMIN,RMAX,INFO,ERROR,*9999)

                  IF(ipulse.EQ.1) THEN
                    RDATA(1)=AQ(maqp1i,nq)
                  ELSE
                    RDATA(1)=AQ(maqp2i,nq)
                  ENDIF

                  FORMAT='($,'' Enter current for pulse '//
     '              CHAR1(1:1)//' (uA/mm^3) [0.1D+03]: '',E12.4)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,
     '              LDATA,LDEFLT,RDATA,RZERO,0.0d0,1.d9,INFO,ERROR,
     '              *9999)
                ENDDO

              ENDIF
            ENDDO

            FORMAT='($,'' Enter collocation point '
     '        //'#s/name [EXIT]: '',I6)'
            IDATA(1)=0
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     &        NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          ENDIF
        ENDIF !Not using a real electrical cellular model. ITYP19 and ITYP3

C Set up initial conditions

        IF(IOTYPE.NE.3) THEN
          IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
            DO nq=1,NQT
              AQ(ATIME,nq)=0.0d0 !initialise activation time
              AQ(DPOT,nq)=0.0d0
            ENDDO
          ENDIF

          CALL ASSERT(USE_CELL.EQ.1,' >>Set USE_CELL to 1 in ippara',
     '      ERROR,*9999)
          CALL ASSERT(USE_GRID.EQ.1,' >>Set USE_GRID to 1 in ippara',
     '      ERROR,*9999)

          IF(ITYP19(nr,nx).EQ.1) THEN !electrical
            IF(ITYP3(nr,nx).EQ.1) THEN !Cubic
              CALL CUBIC_INIT_GRID(NQLIST,CQ,RCQS,YQ(1,1,1,nx),
     '          YQS,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN !FHN
              CALL FHN_INIT_GRID(NQLIST,CQ,RCQS,YQ(1,1,1,nx),
     '          YQS,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.3) THEN !VCD
              CALL VCD_INIT_GRID(NQLIST,CQ,RCQS,YQ(1,1,1,nx),
     '          YQS,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.4) THEN !BR
              CALL BR_INIT_GRID(NQLIST,CQ,RCQS,YQ(1,1,1,nx),
     '          YQS,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.5) THEN !JRW
              !do nothing
              NQLIST(1)=1 !for bidomain
            ELSE IF(ITYP3(nr,nx).EQ.6) THEN !LR2
              !do nothing
              NQLIST(1)=1 !for bidomain
            ELSE IF(ITYP3(nr,nx).EQ.7) THEN !DFN
              CALL ASSERT(.FALSE.,' >>Not implemented',ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.8) THEN !N98
              !do nothing
              NQLIST(1)=1 !for bidomain
            ELSE IF(ITYP3(nr,nx).EQ.9) THEN !HH
              !do nothing
              NQLIST(1)=1 !for bidomain
            ELSE IF(ITYP3(nr,nx).EQ.10) THEN !User defined
              NQLIST(1)=1 !for bidomain
              !do nothing
            ENDIF
          ELSEIF(ITYP19(nr,nx).EQ.7) THEN
            NQLIST(1)=1 !for bidomain
          ENDIF !ITYP19(nr,nx).EQ.1

          IF((ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7)
     '         .AND.KTYP32.EQ.2) THEN !bidomain
             DO nq=1,NQT
                YQ(nq,NQLIST(1),1,nx_ext)=0.0d0
                                !Extracellular resting potential
                IF(CALL_CELL.AND.(ITYP3(nr,nx).GT.4))
     '            CQ(9,nq)=YQS(1,nq) !Resting potential
             ENDDO
          ENDIF

!Old way
          IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
!electrical or coupled model
            IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.1) THEN
            ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.2) THEN
            ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.3) THEN
            ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.4) THEN
C MLB please leave
C          ELSE IF(ITYP3(nr,nx).EQ.6) THEN
C           Luo-Rudy
C            CALL ASSERT(NIQM.GE.10,'>>Increase NIQM (10)',ERROR,*9999)
C            CALL ASSERT(NQM.GE.10,'>>Increase NQM (10)',ERROR,*9999)
C            NQLIST(0)=10
C
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(1),
C     '        NIQ_V,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(2),
C     '        NIQ_BNDRY,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(3),
C     '        NIQ_X,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(4),
C     '        NIQ_CAI,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(5),
C     '        NIQ_M,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(6),
C     '        NIQ_H,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(7),
C     '        NIQ_J,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(8),
C     '        NIQ_D,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(9),
C     '        NIQ_F,ERROR,*9999)
C            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(10),
C     '        NIQ_OLDSOLN,ERROR,*9999)
C
C            VALUES(2)=0.0d0
C
C            DO nq=1,NQT
C              VALUES(1)=CQ(9,nq) !Transmembrane resting potential
C              VALUES(10)=VALUES(1)
C
C              YQ(nq,NQLIST(1),1,nx)=VALUES(1)
C              YQ(nq,NQLIST(10),1,nx)=VALUES(10)
C
C              V=CQ(9,nq)
C              Am=(0.32d0*(V+47.13d0))/(1.d0-DEXP(-V-47.13d0))
C              Bm=0.08d0*DEXP(-V/11.d0)
C
C              IF(V.GE.-40.d0) THEN
C                Ah=0.d0
C                Aj=0.d0
C                Bh=1.d0/(0.13d0*(DEXP((V+10.66d0)/(-11.1d0))+1.d0))
C                Bj=0.3d0*DEXP(-2.535d-7*V)/
C     '            (1.d0+DEXP(-0.1d0*(V+32.d0)))
C              ELSE
C                Ah=0.135d0*DEXP((-80.d0-V)/6.8d0)
C                Aj=(-1.2714d5*DEXP(0.2444d0*V)
C     '            -3.474d-5*DEXP(-0.04391d0*V))*(V+37.78d0)/
C     '            (1.d0+DEXP(0.311d0*(V+79.23d0)))
C                Bh=3.56d0*DEXP(0.079d0*V)+3.1d5*DEXP(0.35d0*V)
C                Bj=0.1212d0*DEXP(-0.01052d0*V)/
C     '            (1.d0+DEXP(-0.1378d0*(V+40.14d0)))
C              ENDIF
C
C              Ad=0.095d0*DEXP(-0.01d0*(V-5.d0))/
C     '          (1.d0+DEXP(-0.072d0*(V-5.d0)))
C              Bd=0.07d0*DEXP(-(V+44.d0)/59.d0)/
C     '          (1.d0+DEXP(0.05d0*(V+44.d0)))
C              Af=0.012d0*DEXP(-0.008d0*(V+28.d0))/
C     '          (1.d0+DEXP(0.15d0*(V+28.d0)))
C              Bf=0.0065d0*DEXP(-0.02d0*(V+30.d0))/
C     '          (1.d0+DEXP(-0.2d0*(V+30.d0)))
C              Ax=5.d-4*DEXP(-(V+50.d0)/12.1d0)/
C     '          (1.d0+DEXP((V+50.d0)/17.5d0))
C              Bx=0.0013d0*DEXP(-0.06d0*(V+20.d0))/
C     '          (1.d0+DEXP(-0.04d0*(V+20.d0)))
C
C              !Steady-state values
C              YQ(nq,NQLIST(3),1,nx)=Ax/(Ax+Bx) !x1
C              YQ(nq,NQLIST(4),1,nx)=1.d-7      !Cai
C              YQ(nq,NQLIST(5),1,nx)=Am/(Am+Bm) !m
C              YQ(nq,NQLIST(6),1,nx)=Ah/(Ah+Bh) !h
C              YQ(nq,NQLIST(7),1,nx)=Aj/(Aj+Bj) !j
C              YQ(nq,NQLIST(8),1,nx)=Ad/(Ad+Bd) !d
C              YQ(nq,NQLIST(9),1,nx)=Af/(Af+Bf) !f
C            ENDDO
C
C            IF(KTYP32.EQ.2) THEN
C              DO nq=1,NQT
C                YQ(nq,NQLIST(1),1,nx_ext)=0.0d0
C                !Extracellular resting potential
C              ENDDO
C            ENDIF
C
            ELSE IF((ITYP19(nr,nx).EQ.1.AND.(ITYP3(nr,nx).EQ.5.OR.
     '          ITYP3(nr,nx).EQ.6.OR.
     '          ITYP3(nr,nx).EQ.8.OR.ITYP3(nr,nx).EQ.9.OR.
     '          ITYP3(nr,nx).EQ.10)).OR.
     '          (ITYP19(nr,nx).EQ.7.AND.ITYP3(nr,nx).LE.3)) THEN
C           JRW, Noble 98, Luo-Rudy, Hodgkin-Huxley, or User defined
C           electrical models or a coupled model
              CALL ASSERT(NIQM.GE.3,'>>Increase NIQM (>=3)',ERROR,*9999)
              CALL ASSERT(CALL_CELL,'>>You must define cell first',
     '          ERROR,*9999)
              CALL ASSERT(CALL_CELL_MATE,
     '          '>>You must define cell materials first',ERROR,*9999)

              CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(1),NIQ_V,
     '          ERROR,*9999)
              IF(NQLIST(1).EQ.0) THEN
                CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(1),
     '            NIQ_V,ERROR,*9999)
              ENDIF

              CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(2),NIQ_BNDRY,
     '          ERROR,*9999)
              IF(NQLIST(2).EQ.0) THEN
                CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(2),
     '            NIQ_BNDRY,ERROR,*9999)
              ENDIF

C KAT 2001-04-12: moved to after user specified initial conditions
C              CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,NQLIST(3),NIQ_OLDSOLN,
C     '          ERROR,*9999)
C              IF(NQLIST(3).EQ.0) THEN
C                CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,NQLIST(3),
C     '            NIQ_OLDSOLN,ERROR,*9999)
C              ENDIF

              DO nq=NQR(1,nr),NQR(2,nr)
                YQ(nq,NQLIST(1),1,nx)=YQS(Vm,nq)
C??? KAT 2001-04-12: Do we need to initialize boundary condition values
C???            or are they only used when FIX is set?
                YQ(nq,NQLIST(2),1,nx)=0.0d0
C                YQ(nq,NQLIST(3),1,nx)=YQ(nq,NQLIST(1),1,nx)
              ENDDO

              IF(KTYP32.EQ.2) THEN
                DO nq=1,NQT
                  YQ(nq,NQLIST(1),1,nx_ext)=0.0d0
!Extracellular resting potential
                ENDDO
              ENDIF

            ELSE IF(ITYP19(nr,nx).EQ.1.AND.ITYP3(nr,nx).EQ.7) THEN
C           Difrancesco-Noble
              CALL ASSERT(.FALSE.,' >>Not implemented',ERROR,*9999)
            ENDIF !Model type
          ENDIF !ITYP19(nr,nx).EQ.1
        ENDIF !Initialising

C Set up boundary conditions

        IF(IOTYPE.NE.3) THEN
C         Initialise the fix array
C KAT: Is this necessary: why not just initialize what will be used?
          DO nyq=1,NYQM
            DO niy=1,NIYFIXM
              FIXQ(nyq,niy,nx)=.FALSE.
            ENDDO
          ENDDO
          IF((ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7).AND.
     '         KTYP32.EQ.2) THEN !bidomain
            DO nyq=1,NYQM
              DO niy=1,NIYFIXM
                FIXQ(nyq,niy,nx_ext)=.FALSE.
              ENDDO
            ENDDO
          ENDIF
C         Initialise other parameters
C          T_POTE=.FALSE.
C          T_FLUX=.TRUE.
          E_POTE=.FALSE.
          E_FLUX=.FALSE.
          E_ANAL=.TRUE.
          E_INJ=.FALSE.
        ENDIF

        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq_bc,NIQ_BNDRY,ERROR,*9999)

C        IF((ITYP16(nr,nx).GE.2).OR.(KTYP32.EQ.2)) THEN
C          !implicit finite differences or bidomain
        CALL ASSERT(NIYFIXM.GE.3,'>>Increase NIYFIXM (>=3)',
     '    ERROR,*9999)

        IF(ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7) THEN
!electrical or coupled model
          DO nq=1,NQT
C!!! The first dimension of FIXQ is NYQM: Is nyq the same as nq?
            IF(NWQ(1,nq,1).GT.0) FIXQ(nq,2,nx)=.TRUE.
          ENDDO

          FORMAT=
     '      '(/'' Specify the type of intracellular b.c. [1]: '''//
     '      '/''   (1) Zero transmembrane flux'''//
     '      '/''   (2) Zero intracellular flux'''//
     '      '/$,''    '',I1)'
          !Note that for monodomain problems only the first of the
          !two bc types is ever used.
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP38
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP38=IDATA(1)

C          FORMAT='($,'' Do you want to fix any boundary transmembrane'
C     '      //' potentials? [N]? '',A)'
C          ADEFLT(1)='N'
C          IF(IOTYPE.EQ.3) THEN
C            IF(T_POTE) THEN
C              ADATA(1)='Y'
C            ELSE
C              ADATA(1)='N'
C            ENDIF
C          ENDIF
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C          FORMAT='($,'' Do you want to fix any boundary transmembrane'
C     '      //' potentials? [N]? '',A)'
C          ADEFLT(1)='N'
C          IF(IOTYPE.EQ.3) THEN
C            IF(T_POTE) THEN
C          IF(IOTYPE.NE.3) THEN
C            IF(ADATA(1).EQ.'Y') THEN
C              T_POTE=.TRUE.
C            ELSE IF(ADATA(1).EQ.'N') THEN
C              T_POTE=.FALSE.
C            ENDIF
C          ENDIF
C
C          IF(T_POTE) THEN
C            IF(IOTYPE.NE.3) THEN
C 803          FORMAT='($,'' Enter collocation point '
C     '          //'#s/name [EXIT]: '',I5)'
C              CDATA(1)='GRIDS' !for use with group input
CC KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
CC               passed to GINOUT instead of IDATA.
CC            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
CC     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NMAX,
CC     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
CC
CC            IF(IDATA(1).NE.0) THEN !not default exit
CC              NQLIST(0)=IDATA(0)
CC              DO n=1,IDATA(0)
CC                NQLIST(n)=IDATA(n)
CC              ENDDO !n
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,0,NMAX,
C     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C
C              IF(NQLIST(1).NE.0) THEN !not default exit
C
C                nq=NQLIST(1) !rest of group is filled in afterwards
C                FORMAT='($,'' Enter the value of the boundary condition'
C     '            //' [0.0]: '',E12.4)'
C                RDEFLT(1)=0.0d0
C                CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
C     '            RMAX,INFO,ERROR,*9999)
C                YQ(nq,niq_bc,1,nx)=RDATA(1)
C                FIXQ(nq,1,nx)=.TRUE.
C
C                DO n=2,NQLIST(0) !fill in rest of group
C                  nq=NQLIST(n)
C                  YQ(nq,niq_bc,1,nx)=YQ(NQLIST(1),niq_bc,1,nx)
C                  FIXQ(nq,1,nx)=.TRUE.
C                ENDDO
C                GO TO 803
C              ENDIF
C
C            ELSE IF(IOTYPE.EQ.3) THEN
C              DO nq=1,NQT
C                IF(FIXQ(nq,1,nx)) THEN
C                  FORMAT='($,'' Enter collocation point '
C     '              //'#s/name [EXIT]: '',I6)'
C                  IDATA(1)=nq
C                  CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IZERO,0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
C     '              INFO,ERROR,*9999)
C
C                  RDATA(1)=YQ(nq,niq_bc,1,nx)
C                  FORMAT='($,'' Enter the value of the boundary'
C     '              //' condition [0.0]: '',E12.4)'
C                  CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,
C     '              RMAX,INFO,ERROR,*9999)
C                ENDIF
C              ENDDO
C
C              FORMAT='($,'' Enter collocation point '
C     '          //'#s/name [EXIT]: '',I6)'
C              IDATA(1)=0
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NMAX,
C     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            ENDIF
C          ENDIF
C
C          FORMAT='($,'' Do you want to fix any boundary transmembrane'
C     '      //' fluxes? [Y]? '',A)'
C          ADEFLT(1)='Y'
C          IF(IOTYPE.EQ.3) THEN
C            IF(T_FLUX) THEN
C              ADATA(1)='Y'
C            ELSE
C              ADATA(1)='N'
C            ENDIF
C          ENDIF
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            IF(ADATA(1).EQ.'Y') THEN
C              T_FLUX=.TRUE.
C            ELSE IF(ADATA(1).EQ.'N') THEN
C              T_FLUX=.FALSE.
C            ENDIF
C          ENDIF
C
C          IF(T_FLUX) THEN
C            IF(IOTYPE.NE.3) THEN
C 804          FORMAT='($,'' Enter collocation point '
C     '          //'#s/name [EXIT]: '',I5)'
C              CDATA(1)='GRIDS' !for use with group input
CC KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
CC               passed to GINOUT instead of IDATA.
CC            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
CC     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NMAX,
CC     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
CC
CC            IF(IDATA(1).NE.0) THEN !not default exit
CC              NQLIST(0)=IDATA(0)
CC              DO n=1,IDATA(0)
CC                NQLIST(n)=IDATA(n)
CC              ENDDO !n
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,0,NMAX,
C     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C
C              IF(NQLIST(1).NE.0) THEN !not default exit
C
C                nq=NQLIST(1) !rest of group is filled in afterwards
C                FORMAT='($,'' Enter the value of the boundary condition'
C     '            //' [0.0]: '',E12.4)'
C                RDEFLT(1)=0.0d0
C                CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '            IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
C     '            RMAX,INFO,ERROR,*9999)
C                YQ(nq,niq_bc,1,nx)=RDATA(1)
C                FIXQ(nq,2,nx)=.TRUE.
C
C                DO n=2,NQLIST(0) !fill in rest of group
C                  nq=NQLIST(n)
C                  YQ(nq,niq_bc,1,nx)=YQ(NQLIST(1),niq_bc,1,nx)
C                  FIXQ(nq,2,nx)=.TRUE.
C                ENDDO
C                GO TO 804
C              ENDIF
C
C            ELSE IF(IOTYPE.EQ.3) THEN
C              DO nq=1,NQT
C                IF(FIXQ(nq,2,nx)) THEN
C                  FORMAT='($,'' Enter collocation point '
C     '              //'#s/name [EXIT]: '',I6)'
C                  IDATA(1)=nq
C                  CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IZERO,0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
C     '              INFO,ERROR,*9999)
C
C                  RDATA(1)=YQ(nq,niq_bc,1,nx)
C                  FORMAT='($,'' Enter the value of the boundary'
C     '              //' condition [0.0]: '',E12.4)'
C                  CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,-RMAX,
C     '              RMAX,INFO,ERROR,*9999)
C                ENDIF
C              ENDDO
C
C              FORMAT='($,'' Enter collocation point '
C     '          //'#s/name [EXIT]: '',I6)'
C              IDATA(1)=0
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NMAX,
C     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            ENDIF
C          ENDIF

          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq_old,NIQ_OLDSOLN,
     '      ERROR,*9999)
          IF(niq_old.EQ.0) THEN
            CALL NIQ_LOC(NIQ_ALLOCATE_AND_LOCK,NIQ_ION,niq_old,
     '        NIQ_OLDSOLN,ERROR,*9999)
          ENDIF
          CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
          DO nq=NQR(1,nr),NQR(2,nr)
            YQ(nq,niq_old,1,nx)=YQ(nq,niqV,1,nx)
          ENDDO

          IF(KTYP32.EQ.2) THEN !bidomain
            FORMAT='($,'' Do you want to fix any boundary extracellular'
     '        //' potentials? [N]? '',A)'
            ADEFLT(1)='N'
            IF(IOTYPE.EQ.3) THEN
              IF(E_POTE) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                E_POTE=.TRUE.
              ELSE IF(ADATA(1).EQ.'N') THEN
                E_POTE=.FALSE.
              ENDIF
            ENDIF

            IF(E_POTE) THEN
              IF(IOTYPE.NE.3) THEN
 805            FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I5)'
                CDATA(1)='GRIDS' !for use with group input
C KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
C               passed to GINOUT instead of IDATA.
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
C     '          0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '          ERROR,*9999)
C
C              IF(IDATA(1).NE.0) THEN !not default exit
C                NQLIST(0)=IDATA(0)
C                DO n=1,IDATA(0)
C                  NQLIST(n)=IDATA(n)
C                ENDDO !n
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)

                IF(NQLIST(1).NE.0) THEN !not default exit

                  nq=NQLIST(1) !rest of group is filled in afterwards
                  FORMAT='($,'' Enter the value of the boundary'
     '              //' condition [0.0]: '',E12.4)'
                  RDEFLT(1)=0.0d0
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,INFO,ERROR,*9999)
C KAT 2001-04-12: shouldn't niqV also be set here?
C                 do we even need to set niq_bc?
                  YQ(nq,niq_bc,1,nx_ext)=RDATA(1)
                  FIXQ(nq,1,nx_ext)=.TRUE.
                  DO n=2,NQLIST(0) !fill in rest of group
                    nq=NQLIST(n)
                    YQ(nq,niq_bc,1,nx_ext)=YQ(NQLIST(1),niq_bc,1,nx_ext)
                    FIXQ(nq,1,nx_ext)=.TRUE.
                  ENDDO
                  GO TO 805
                ENDIF

              ELSE IF(IOTYPE.EQ.3) THEN
                DO nq=1,NQT
                  IF(FIXQ(nq,1,nx_ext)) THEN
                    FORMAT='($,'' Enter collocation point #s/name'
     '                //' [EXIT]: '',I6)'
                    IDATA(1)=nq
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IZERO,0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)

                    RDATA(1)=YQ(nq,niq_bc,1,nx_ext)
                    FORMAT='($,'' Enter the value of the boundary'
     '                //' condition [0.0]: '',E12.4)'
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     &                -RMAX,RMAX,INFO,ERROR,*9999)
                  ENDIF
                ENDDO

                FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I6)'
                IDATA(1)=0
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
              ENDIF
            ENDIF

            FORMAT='($,'' Do you want to fix any boundary extracellular'
     '        //' fluxes? [N]? '',A)'
            ADEFLT(1)='N'
            IF(IOTYPE.EQ.3) THEN
              IF(E_FLUX) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                E_FLUX=.TRUE.
              ELSE IF(ADATA(1).EQ.'N') THEN
                E_FLUX=.FALSE.
              ENDIF
            ENDIF

            IF(E_FLUX) THEN
              WRITE(OP_STRING,
     '          '('' >>WARNING: at least one extracellular'
     '          //' potential'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''          must be specified'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)

              IF(IOTYPE.NE.3) THEN
 806            FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I5)'
                CDATA(1)='GRIDS' !for use with group input
C KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
C               passed to GINOUT instead of IDATA.
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
C     '          0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '          ERROR,*9999)
C
C              IF(IDATA(1).NE.0) THEN !not default exit
C                NQLIST(0)=IDATA(0)
C                DO n=1,IDATA(0)
C                  NQLIST(n)=IDATA(n)
C                ENDDO !n
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)

                IF(NQLIST(1).NE.0) THEN !not default exit

                  nq=NQLIST(1) !rest of group is filled in afterwards
                  FORMAT='($,'' Enter the value of the boundary'
     '              //' condition [0.0]: '',E12.4)'
                  RDEFLT(1)=0.0d0
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,INFO,ERROR,*9999)
                  YQ(nq,niq_bc,1,nx_ext)=RDATA(1)
                  FIXQ(nq,2,nx_ext)=.TRUE.

                  DO n=2,NQLIST(0) !fill in rest of group
                    nq=NQLIST(n)
                    YQ(nq,niq_bc,1,nx_ext)=YQ(NQLIST(1),niq_bc,1,nx_ext)
                    FIXQ(nq,2,nx_ext)=.TRUE.
                  ENDDO
                  GO TO 806
                ENDIF

              ELSE IF(IOTYPE.EQ.3) THEN
                DO nq=1,NQT
                  IF(FIXQ(nq,2,nx_ext)) THEN
                    FORMAT='($,'' Enter collocation point #s/name'
     '                //' [EXIT]: '',I6)'
                    IDATA(1)=nq
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IZERO,0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)

                    RDATA(1)=YQ(nq,niq_bc,1,nx_ext)
                    FORMAT='($,'' Enter the value of the boundary'
     '                //' condition [0.0]: '',E12.4)'
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RZERO,
     &                -RMAX,RMAX,INFO,ERROR,*9999)
                  ENDIF
                ENDDO

                FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I6)'
                IDATA(1)=0
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
              ENDIF
            ENDIF

            FORMAT='($,'' Do you want to fix any analytic boundary'
     '        //' ext. potentials? [Y]? '',A)'
            ADEFLT(1)='Y'
            IF(IOTYPE.EQ.3) THEN
              IF(E_ANAL) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                E_ANAL=.TRUE.
              ELSE IF(ADATA(1).EQ.'N') THEN
                E_ANAL=.FALSE.
              ENDIF
            ENDIF

            IF(E_ANAL) THEN
              WRITE(OP_STRING,'('' >>WARNING: this only works for'
     '          //' equal anisotropy'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)

              IF(IOTYPE.NE.3) THEN
 807            FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I5)'
                CDATA(1)='GRIDS' !for use with group input
C KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
C               passed to GINOUT instead of IDATA.
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
C     '          0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '          ERROR,*9999)
C
C              IF(IDATA(1).NE.0) THEN !not default exit
C                NQLIST(0)=IDATA(0)
C                DO n=1,IDATA(0)
C                  NQLIST(n)=IDATA(n)
C                ENDDO !n
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)

                IF(NQLIST(1).NE.0) THEN !not default exit

                  DO n=1,NQLIST(0)
                    nq=NQLIST(n)
                    FIXQ(nq,3,nx_ext)=.TRUE.
                  ENDDO
                  GO TO 807
                ENDIF

              ELSE IF(IOTYPE.EQ.3) THEN
                DO nq=1,NQT
                  IF(FIXQ(nq,3,nx_ext)) THEN
                    FORMAT='($,'' Enter collocation point #s/name'
     '                //' [EXIT]: '',I6)'
                    IDATA(1)=nq
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IZERO,0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                  ENDIF
                ENDDO

                FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I6)'
                IDATA(1)=0
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
              ENDIF
            ENDIF

            FORMAT='($,'' Do you want any time dependent '
     '        //'extracellular bcs [N]? '',A)'
            ADEFLT(1)='N'
            IF(IOTYPE.EQ.3) THEN
              IF(E_INJ) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                E_INJ=.TRUE.
              ELSE IF(ADATA(1).EQ.'N') THEN
                E_INJ=.FALSE.
              ENDIF
            ENDIF

            !initialise NWQ for no time variable
            DO nq=1,NQT
              NWQ(8,nq,1)=0
            ENDDO

            IF(E_INJ) THEN
              CALL ASSERT(CALL_TIME,'>>Must define time first',
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
 808            FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I5)'
                CDATA(1)='GRIDS' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)

                IF(NQLIST(1).NE.0) THEN !not default exit
                  FORMAT='($,'' Enter the name of the time '
     '              //'variable []: '',A32)'
C LKC 6-DEC-2000 Zero length string not allowe
C                  CDEFLT(1)=''
                  CDEFLT(1)='-'
                  CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IZERO,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)

                  CALL STRING_TRIM(CDATA(1),IBEG1,IEND1)
                  NVARINDEX=0
                  DO i=1,NTIMEVARST
                    CALL STRING_TRIM(TIME_VARIABLE_NAMES(i),IBEG,IEND)
                    IF((IEND1-IBEG1).EQ.(IEND-IBEG)) THEN
                      IF(CDATA(1)(IBEG1:IEND1).EQ.
     '                  TIME_VARIABLE_NAMES(i)(IBEG:IEND)) NVARINDEX=i
                    ENDIF
                  ENDDO !i
                  CALL ASSERT(NVARINDEX.NE.0,
     '              '>>Time variable not found',ERROR,*9999)

                  DO n=1,NQLIST(0)
                    nq=NQLIST(n)
                    NWQ(8,nq,1)=NVARINDEX
                  ENDDO
                  GOTO 808
                ENDIF

              ELSE IF(IOTYPE.EQ.3) THEN
                DO nq=1,NQT
                  IF(NWQ(8,nq,1).GT.0) THEN
                    FORMAT='($,'' Enter collocation point #s/name'
     '                //' [EXIT]: '',I6)'
                    IDATA(1)=nq
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IZERO,0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)

                    FORMAT='($,'' Enter the name of the time '
     '                //'variable []: '',A32)'
                    CDATA(1)=TIME_VARIABLE_NAMES(NWQ(8,nq,1))
                    CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IZERO,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &                RMAX,INFO,ERROR,*9999)
                  ENDIF
                ENDDO

                FORMAT='($,'' Enter collocation point '
     '            //'#s/name [EXIT]: '',I6)'
                IDATA(1)=0
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,NMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
              ENDIF
            ENDIF

C KAT 31Mar01:    using old solution in forming initial guess for
C                 iterative linear solvers.  There is no reason why this
C                 can't be used for collocation also.  See MARCH8.
          DO nq=NQR(1,nr),NQR(2,nr)
            YQ(nq,niq_old,1,nx_ext)=YQ(nq,niqV,1,nx_ext)
          ENDDO

          ENDIF
        ENDIF !ITYP19(nr,nx).EQ.1
C         Check the necessary boundary conditions have been applied
C        IF(ITYP16(nr,nx).GE.2) THEN
C          IF((.NOT.T_POTE).AND.(.NOT.T_FLUX))
C     '      CALL ASSERT(.FALSE.,'>>No boundary conditions applied',
C     '      ERROR,*9999)
C        ENDIF
        IF((ITYP19(nr,nx).EQ.1.OR.ITYP19(nr,nx).EQ.7).AND.
     '    KTYP32.EQ.2) THEN
          IF((.NOT.E_POTE).AND.(.NOT.E_ANAL))
     '      CALL ASSERT(.FALSE.,'>>No boundary grid points fixed',
     '      ERROR,*9999)
          IF((.NOT.E_POTE).AND.(.NOT.E_ANAL).AND.(.NOT.E_FLUX))
     '      CALL ASSERT(.FALSE.,'>>No boundary conditions applied',
     '      ERROR,*9999)
        ENDIF

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volumes also
      ELSE IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '        ITYP4(nr,nx).EQ.7) THEN !collocation or grid-based FE or grid FV
C ***   collocation 
        IF(iotype.ne.3) THEN
          DO nq=1,NQT
            DO na=1,NAM
              DO niq=1,NIQM
                YQ(nq,niq,na,nx)=0.0d0
              ENDDO !niq
            ENDDO !na
            IF(NWQ(1,nq,1).GT.0) THEN !boundary point
              DO na=1,NAM
                NWQ(5,nq,na)=2 !initialize to Neumann bdry condition
              ENDDO !na
            ELSE !interior point
              DO na=1,NAM
                NWQ(5,nq,na)=0
              ENDDO !na
            ENDIF !NWQ
          ENDDO !nq
        ENDIF !iotype

        DO loop=1,2
          IF(loop.EQ.1) THEN
            FORMAT='('' Bdry or initial conditions'//
     '        ' at collocation points:'')'
          ELSE IF(loop.EQ.2) THEN
            FORMAT='(/'' Flux conditions at collocation points:'')'
          ENDIF
          CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,-IMAX,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          RDATA(1)=0.d0 !initialize default
 500      FORMAT='($,'' Enter collocation point #s/name [EXIT]: '',I5)'
          CDATA(1)='GRIDS' !for use with group input
C KAT 12May99:  To avoid large unadjustable static arrays, NQLIST is
C               passed to GINOUT instead of IDATA.
C          CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NQM,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IDATA(1).NE.0) THEN !not default exit
C            NQLIST(0)=IDATA(0)
C            DO n=1,IDATA(0)
C              NQLIST(n)=IDATA(n)
C            ENDDO !n
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NQLIST,IZERO,0,NQM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(NQLIST(1).NE.0) THEN !not default exit
C***        Define b.c. or i.c. for first grid pt in group
            nq=NQLIST(1) !rest of group is filled in afterwards
            WRITE(CHAR,'(I6)') nq
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=YQ(nq,1,1,nx)
            ENDIF
            RDEFLT(1)=RDATA(1)
            WRITE(CHAR1,'(E12.5)') RDEFLT(1)
            IF(NWQ(1,nq,1).EQ.0) THEN !Interior point
              FORMAT='($,'' Enter i.c. at int. pt #'//CHAR(IBEG:IEND)
     '          //' ['//CHAR1(1:12)//']: '',E12.5)'
            ELSE IF(NWQ(1,nq,1).GT.0) THEN !Boundary point
              FORMAT='($,'' Enter b.c. at bdry pt #'//CHAR(IBEG:IEND)
     '          //' ['//CHAR1(1:12)//']: '',E12.5)'
            ENDIF !NWQ>0
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
C***          Apply bdry conditions to whole group
              DO n=1,NQLIST(0)
                nq=NQLIST(n)
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPINI3_1)
                  WRITE(OP_STRING,'('' Group data: n='',I3,'' nq='''
     '              //',I5)') n,nq
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(IPINI3_1)
                ENDIF
                IF(NWQ(1,nq,1).GT.0) THEN !Boundary point
C                  DO na=1,NMGT
                  DO na=1,NAM
                    NWQ(5,nq,na)=loop
                  ENDDO
                ENDIF
                IF(loop.EQ.1) THEN !store Dirichlet b.c.
C                  DO na=1,NMGT
                  DO na=1,NAM
                    YQ(nq,1,na,nx)=RDATA(1)
                  ENDDO
                ELSE IF(loop.EQ.2) THEN !store Neumann b.c.
                ENDIF !loop
              ENDDO !n
            ENDIF !iotype.ne.3

            GOTO 500
          ENDIF !idata(1).ne.0
        ENDDO !loop


      ELSE IF((ITYP4(nr,nx).EQ.3).AND.(ITYP3(nr,nx).EQ.1)) THEN
                                              ! FD & flow in elastic tubes

        FORMAT='($,'' Do you want start from a'//
     '    ' previous solution [N]? '',A)'
        ADEFLT(1)='N'
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '    IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '    ERROR,*9999)
        IF(ADATA(1).EQ.'Y') THEN
          INIT_RESET=.FALSE.
        ELSE
          INIT_RESET=.TRUE.
        ENDIF

C If a previous solution is started from the flow variables in YQ are
C left unchanged otherwise YQ is initialised so all flow and pressure
C gradients are zero

C PM 26-JUL-01 : unnecessary
C        Po=10.6D0

C PM 28-NOV-01 : This is not required as expoint is now used.
C        FORMAT='($,'' Do you want to out put data files'//
C     '    ' to write over old files [Y]? '',A)'
C        ADEFLT(1)='Y'
C        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '    IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
C     '    ERROR,*9999)

C        IF(ADATA(1).EQ.'Y') THEN
C          DATA_FILE_COUNT=0
C        ELSE
C          FORMAT='($,''Do you want to set data file counter [N]? '',A)'

C          ADEFLT(1)='N'

C          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '      IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
C     '      ERROR,*9999)

C          IF(ADATA(1).EQ.'Y') THEN
C            FORMAT='($,''Enter the initial file number[1]: '',I6)'
C            IDEFLT(1)=1
C            IF(IOTYPE.EQ.3) IDATA(1)=DATA_FILE_COUNT
C            CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,10000,
C     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C            IF(IOTYPE.NE.3) DATA_FILE_COUNT=IDATA(1)
c          ENDIF
C        ENDIF

C determines the number and naming convention of the output .exnode
C and gnuplot files

        FORMAT='($,'' Enter value initial pressure '//
     '    'Po [10.6 kPa]: '',F8.4)'
        RDEFLT(1)=10.60d0
        IF(IOTYPE.EQ.3) RDATA(1)=Po
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-50.0d0,50.0d0,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Po=RDATA(1)

C determines the inital pressure which is set to be constant throughout
C the system

C PM 26-JUL-01 : Question re-phrased
C        FORMAT='('' Enter the black box constant [0.0001]'''//
C     '    '/$,''   '',F8.4)'
C        RDEFLT(1)=0.0001d0
C        IF(IOTYPE.EQ.3) RDATA(1)=FLOW_CONST
C        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,5.0d0,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) FLOW_CONST=RDATA(1)

C sets the black box constant

        FORMAT='($,'' Do you want to solve the '//
     '    'venous network [Y]? '',A)'
        ADEFLT(1)='Y'
        IF(IOTYPE.EQ.3) ADATA(1)=VENOUS_NETWORK
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '    IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '    ERROR,*9999)
        IF(IOTYPE.NE.3) VENOUS_NETWORK=ADATA(1)


C checks whether computaion on the venous side is necessary

        IF((VENOUS_NETWORK.EQ.'Y').OR.
     '    (VENOUS_NETWORK.EQ.'y')) THEN
          FORMAT='($,'' Is venous network coupled via micro-'//
     '      'circulation network [Y]? '',A)'
          ADEFLT(1)='Y'
          IF(IOTYPE.EQ.3) ADATA(1)=MICRO_NETWORK
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '      ERROR,*9999)
          IF(IOTYPE.NE.3) MICRO_NETWORK=ADATA(1)
        ENDIF

C PM 29-NOV-01 : New questions added for solving non-identical networks
        IF((VENOUS_NETWORK.EQ.'Y').OR.
     '    (VENOUS_NETWORK.EQ.'y')) THEN
          FORMAT='(/'' Geomtery for venous network is [1] : '''//
     '      '/''   (1) Indentical to the arterial network   '''//
     '      '/''   (2) Different from the arterial network  '''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=N_VENOUS_GEOM
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '      1,2,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '      ERROR,*9999)
          IF(IOTYPE.NE.3) N_VENOUS_GEOM=IDATA(1)
        ELSE
          N_VENOUS_GEOM=0
        ENDIF

        IF((VENOUS_NETWORK.EQ.'Y').OR.
     '    (VENOUS_NETWORK.EQ.'y')) THEN
          IF(N_VENOUS_GEOM.EQ.1) THEN
            FORMAT='($,'' Enter the vien/artery radius '//
     '        'ratio [1.22]: '',F8.4)'
            RDEFLT(1)=1.2247449d0
            IF(IOTYPE.EQ.3) RDATA(1)=VIEN_RATIO
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,5.0d0,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) VIEN_RATIO=RDATA(1)
C           ELSE
C             VIEN_RATIO=1.00d0
C           ENDIF

C sets the ratio of the venous radius to the arterial radius for the
C pair of veins which parallel each artery

          ELSEIF(N_VENOUS_GEOM.EQ.2) THEN
            FORMAT='($,'' No.of terminal points in '//
     '        'each network [1]: '',I5)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=N_TERM_P(0)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) N_TERM_P(0)=IDATA(1)

            DO nn=1,N_TERM_P(0)
              WRITE(CHAR1,'(I5)') nn
              IDEFLT(1)=1
              FORMAT='($,'' Node no. for terminal '
     '          //CHAR1(1:5)//' [1]: '',I5)'
              IF(IOTYPE.EQ.3) IDATA(1)=N_TERM_P(nn)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '          IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) N_TERM_P(nn)=IDATA(1)
            ENDDO
          ENDIF
        ENDIF

C PM 26-JUL-01 :These can now be input through iptime
C          FORMAT='('' Enter initial entry point pressure '//
C     '      '[0.0 kPa]''/$,''   '',F8.4)'
C          RDEFLT(1)=0.000d0
C          IF(IOTYPE.EQ.3) RDATA(1)=ENTRY_PRESSURE0
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-50.0d0,50.0d0,INFO,ERROR,*9999)
c          IF(IOTYPE.NE.3) ENTRY_PRESSURE0=RDATA(1)

C          FORMAT='('' Enter the final entry point pressure '//
C     '      '[5.0 kPa]''/$,''   '',F8.4)'
C          RDEFLT(1)=5.0d0
c          IF(IOTYPE.EQ.3) RDATA(1)=ENTRY_PRESSURE1
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-50.0d0,50.0d0,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) ENTRY_PRESSURE1=RDATA(1)

C          FORMAT='('' Enter the initial exit point pressure '//
C     '      '[0.0 kPa]''/$,''   '',F8.4)'
C          RDEFLT(1)=0.000d0
C          IF(IOTYPE.EQ.3) RDATA(1)=ENTRY_PRESSURE0
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-50.0d0,50.0d0,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) EXIT_PRESSURE0=RDATA(1)

C          FORMAT='('' Enter the final exit point pressure '//
C     '      '[0.0 kPa]''/$,''   '',F8.4)'
C          RDEFLT(1)=0.0d0
C          IF(IOTYPE.EQ.3) RDATA(1)=ENTRY_PRESSURE1
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,-50.0d0,50.0d0,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) EXIT_PRESSURE1=RDATA(1)

C sets the increments of the inflow and outflow pressures on
C top of the initial pressure for the begining and end of the
C time specified by the ipsolve file Note: these pressures are
C linearly interpolated at the intermediate time steps

        FORMAT='($,'' Enter the value of gravity [9.81d0] '',D12.4)'
        RDEFLT(1)=9.81d0
        IF(IOTYPE.EQ.3) RDATA(1)=GRAVITY
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GRAVITY=RDATA(1)

        IF(GRAVITY.NE.0.d0) THEN
          DO nj=1,NJT
            WRITE(CHAR,'(I1)') nj
            IF(nj.LT.3) THEN
              RDEFLT(1)=0.00d0
            ELSE
              RDEFLT(1)=1.00d0
            ENDIF
            WRITE(CHAR1,'(F8.6)') RDEFLT(1)
            FORMAT='($,'' Enter the Xj'//CHAR(1:1)//
     &        ' vector direction of gravity ['//CHAR1(1:6)//
     &        '] '',D12.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=G_VECTOR(nj)
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) G_VECTOR(nj)=RDATA(1)
          ENDDO
          CALL NORMALISE(NJT,G_VECTOR,ERROR,*9999) !unit vector - direction of gravity
        ELSE
          DO nj=1,NJT
            G_VECTOR(nj)=0.d0
          ENDDO
        ENDIF
        IF(ITYP12(nr,nx).EQ.2) THEN !pulmonary Navier-Stokes solutions
          FORMAT=
     &      '($,'' Enter the pleural pressure (kPa) [-0.49d0] '',D12.4)'
          RDEFLT(1)=-0.49d0
          IF(IOTYPE.EQ.3) RDATA(1)=PLEURAL_P
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) PLEURAL_P=RDATA(1)
        ENDIF
                
        FORMAT='($,'' Enter the initial flow entry '//
     '    'element[1]: '',I6)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NENQ(1,NQ_START(nr))
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NEM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NE_START=IDATA(1)

        IF(iotype.ne.3) THEN          
          IF(NXI(-1,0,NE_START).EQ.0) THEN
            CALL ASSERT(NXI(1,0,NE_START).NE.0,'>>Element is not'
     '        //' connected',ERROR,*9999)
            NQ_START(nr)=NQNE(NE_START,1)
          ELSE

            CALL ASSERT(NXI(1,0,NE_START).EQ.0,'>>Element is'
     '        //'not at start',ERROR,*9999)
            NQ_START(nr)=NQNE(NE_START,NQET(NQS(NE_START)))
          ENDIF
          
C sets and checks the validality of the inflow outflow
C element at the top of the network
          
          IF(INIT_RESET) THEN
            DO nj=1,NJT
              XQ_IN(nj)=XQ(nj,NQ_START(nr)) !position of nq_start, pressure relative to this
            ENDDO
            
            DO nq=NQR(1,nr),NQR(2,nr)
              
C initalising black box parameters
              
              CALL USER3_CORONARY1(DENOM_BOX(1,1),dummy,
     '          0.0d0,NUM_BOX(1,1),RART,Po,0.0d0,ERROR,*9999)
              
              CALL USER3_CORONARY1(DENOM_BOX(1,3),dummy,
     '          0.0d0,NUM_BOX(1,3),RVIE,Po,0.0d0,ERROR,*9999)
              
              CALL USER3_CORONARY2(DENOM_BOX(1,2),
     '          NUM_BOX(1,2),dummy,dummy,
     '          RCAP,Po,Po,ERROR,*9999)
              
              CALL USER3_CORONARY1(DENOM_BOX(1,4),dummy,
     '          0.0d0,NUM_BOX(1,4),COM1,Po,Po,ERROR,*9999)
              
              CALL USER3_CORONARY1(DENOM_BOX(1,5),dummy,
     '          0.0d0,NUM_BOX(1,5),COM2,Po,Po,ERROR,*9999)

              CALL ASSERT(NIQM.GE.6,'>>Increase NIQM, must be >= 6',
     '          ERROR,*9999)
              YQ(NYNQ(1,nq,0),4,1,nx)=COM1*1000d0 ! C1 arterial capacitance
              YQ(NYNQ(2,nq,0),4,1,nx)=RART/1000d0 !arterial resistance
              YQ(NYNQ(3,nq,0),4,1,nx)=0.0d0    !FC capillary flow
              YQ(NYNQ(4,nq,0),4,1,nx)=COM2*1000d0 !c2 veinous capacitance
              YQ(NYNQ(5,nq,0),4,1,nx)=RVIE/1000d0 !vienous resistance

C sets the inital value for the microcirculation black box which are
C needed to calculate the finite difference time derivatives in
C the subroutine BRANCH2
              DO no_nynq=1,NHQ
                ny=NYNQ(no_nynq,nq,0)
                IF((no_nynq.EQ.1).OR.(no_nynq.EQ.4)) THEN
C set up initial values in YQ
C PM 02-OCT-01 : Initialize all relevant YQs
                  DO niq=1,NIQM
                    YQ(ny,niq,1,nx)=0.0d0
                  ENDDO

                  G_TERM=0.d0
                  DO nj=1,NJT
                    HEIGHT(nj)=XQ_IN(nj)-XQ(nj,nq)
                    G_TERM=G_TERM+HEIGHT(nj)*CQ(1,nq)*G_VECTOR(nj)
     &                *GRAVITY*1000.d0 !GRAVITY TERM (mm/s**2)
                  ENDDO
                  YQ(ny,3,1,nx)=Po+G_TERM !initial pressure
                  YQ(ny,1,1,nx)=Po+G_TERM !initial pressure
                  YQ(ny,6,1,nx)=Po+G_TERM !initial pressure
                  
                ELSEIF ((no_nynq.EQ.2).OR.(no_nynq.EQ.5)) THEN
C initial Radius and the trace component values
C PM 26-JUL-01 : This computation has been moved to upmate
C                  ne=NENQ(1,nq)
C                  nb=NBJ(1,ne)
c                  nj=NJ_LOC(NJL_FIEL,1,nr)
C                  nj2=NJ_LOC(NJL_FIEL,2,nr)
C                  IF((NXI(1,0,ne).LE.1).AND.(NXI(-1,0,ne).LE.1))THEN
C not a bifurcation thus radius is interpolated

C                    CALL XPXE(NBJ(1,ne),
C     '                NKE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
C     '                nr,NVJE(1,1,1,ne),SE(1,1,ne),XE,XP,ERROR,*9999)
C                    XI(1)=0.0d0
C                    nqele=NQNE(ne,1)
C                    XI_COUNT=1
C                    DO WHILE((nqele.NE.nq).AND.
C     '                (nqele.NE.NQNE(ne,NQET(NQS(ne)))))
C                      XI_COUNT=XI_COUNT+1
C                      nqele=NQNE(ne,XI_COUNT)
C                    ENDDO
C                    XI(1)=(1.0d0/(NQET(NQS(ne))-1.0d0))*(XI_COUNT-1.0d0)

C determines the xi position of the grid point nq in the elements ne

C                    nb=NBJ(nj,ne)
C                    Ro=PXI(IBT(1,1,nb),
C     '                IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1),XE(1,nj))

C now interpolate the trace component using the same xi location

C                    nb=NBJ(nj2,ne)
C                    TRACE_COMP=1.0d0-PXI(IBT(1,1,nb),
C     '                IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,XE(1,nj2))

C                  ELSE
C the element is adjacent to a bifurcation thus a discontinuity in radius and
C the trace exists between elements.
C                    R1=XP(1,1,nj,NPNE(1,nb,ne))
C                    R2=XP(1,1,nj,NPNE(NNT(nb),nb,ne))
C                    Ro=MIN(R1,R2)

C                    TRACE_COMP1=1.0d0-XP(1,1,nj2,NPNE(1,nb,ne))
C                    TRACE_COMP2=1.0d0-XP(1,1,nj2,NPNE(NNT(nb),nb,ne))
C                    TRACE_COMP=MAX(TRACE_COMP1,TRACE_COMP2)
C                  ENDIF
C                  CQ(4,nq)=Ro !unstressed radius
                  Ro=CQ(4,nq)

                  IF(no_nynq.EQ.5) THEN !vien
                    Ro=Ro*VIEN_RATIO
                    Go=CQ(11,nq)
                    beta=DMAX1((CQ(9,nq)+CQ(10,nq)),1.0d0)
                  ELSE !artery
                    Go=CQ(5,nq)
                    beta=DMAX1((CQ(6,nq)+CQ(12,nq)),1.0d0)
                  ENDIF

C                    IF(Po.GT.0.0D0) THEN
C PM 02-OCT-01 :
C                    YQ(ny,3,1,nx)=Ro*(((Po/Go)+1.0d0)**
C     '                (1.0d0/beta)) !Radius
C                    YQ(ny,1,1,nx)=Ro*(((Po/Go)+1.0d0)**
C     '                (1.0d0/beta))

C                      DO niq=1,NIQM
C                        YQ(ny,niq,1,nx)=Ro*(((Po/Go)+1.0d0)**
C     '                    (1.0d0/beta))
C                      ENDDO
C                      YQ(ny,3,1,nx)=Ro*(((Po/Go)+1.0d0)**
C     '                  (1.0d0/beta))
C                      YQ(ny,1,1,nx)=Ro*(((Po/Go)+1.0d0)**
C     '                  (1.0d0/beta))
C                      
C                    ELSE
C                      fo=beta*Go/CQ(7,nq)
C PM 02-OCT-01 :
C                    YQ(ny,3,1,nx)=Ro*(1.0d0-(Po/fo))**(-1.0d0/CQ(7,nq))
C                    YQ(ny,1,1,nx)=YQ(ny,3,1,nx)

C                      DO niq=1,NIQM
C                        YQ(ny,niq,1,nx)=Ro*(1.0d0-(Po/fo))**
C     '                    (-1.0d0/CQ(7,nq))
C                      ENDDO
C                      YQ(ny,3,1,nx)=Ro*(1.0d0-(Po/fo))**(-1.0d0/CQ(7,
C     '                  nq))
C                      YQ(ny,1,1,nx)=YQ(ny,3,1,nx)
C                 ENDIF
C.. KSB 2004: this modified to take gravity into account in initial conditions
                    ny_P=NYNQ(1,nq,0) !Pressure value for nq
                    IF(YQ(ny_P,3,1,nx).GT.0.0d0) THEN
                      DO niq=1,NIQM
                        YQ(ny,niq,1,nx)=Ro*(((YQ(ny_P,3,1,nx)/Go)+1.0d0)
     '                    **(1.0d0/beta))
                      ENDDO
                      YQ(ny,3,1,nx)=Ro*(((YQ(ny_P,3,1,nx)/Go)+1.0d0)**
     '                  (1.0d0/beta))
                      YQ(ny,1,1,nx)=Ro*(((YQ(ny_P,3,1,nx)/Go)+1.0d0)**
     '                  (1.0d0/beta))                      
                    ELSE
                      fo=beta*Go/CQ(7,nq)
                      DO niq=1,NIQM
                        YQ(ny,niq,1,nx)=Ro*(1.0d0-(YQ(ny_P,3,1,nx)/fo))
     '                    **(-1.0d0/CQ(7,nq))
                      ENDDO
                      YQ(ny,3,1,nx)=Ro*(1.0d0-(YQ(ny_P,3,1,nx)/fo))**
     '                  (-1.0d0/CQ(7,nq))
                      YQ(ny,1,1,nx)=YQ(ny,3,1,nx)
                      
                    ENDIF
                  
C                 YQ(ny,2,1,nx)=TRACE_COMP
                  
                ELSE IF((no_nynq.EQ.3).OR.(no_nynq.EQ.6)) THEN
C sets the inital flow velocity values and the vessel stretch values to be zero
C PM 02-OCT-01 :
C                  YQ(ny,3,1,nx)=0.0d0 !Velocity
C                  YQ(ny,2,1,nx)=1.0d0 !initial lambda
C                  YQ(ny,9,1,nx)=1.0d0 !initial lambda
C                  YQ(ny,1,1,nx)=0.0d0 !initial velocity

                  DO niq=1,NIQM
                    YQ(ny,niq,1,nx)=1.0d0
                  ENDDO
                  YQ(ny,3,1,nx)=0.0d0  !initial velocity
                  YQ(ny,1,1,nx)=0.0d0  !initial velocity

                ENDIF !pressure or radius or velocity
              ENDDO  !no_nynq
            ENDDO   !nq
          ENDIF    !INIT_RESET ie all flows are set to zero
        ENDIF     ! (iotype.ne.3)

      ELSE
C ***   All other models

!       Check whether any bdry elements

C BDRY_ELEMENT not used
C        BDRY_ELEMENT=.FALSE.
C        DO nr_loop=1,NRT
C          IF(ITYP4(nr_loop,nx).EQ.2.OR.ITYP4(nr_loop,nx).EQ.3) THEN
C            BDRY_ELEMENT=.TRUE.
C          ENDIF
C        ENDDO

        IF(IOTYPE.ne.3) THEN
C***      Need to only initialise those belonging to the
C***      current region.
          IF(NIYFIXM.GT.NIYM) THEN
            WRITE(CHAR4,'(I1)') IDIGITS(NIYFIXM)
            WRITE(ERROR,'(''>>Increase NIYM to '',I'//CHAR4//')')
     '        NIYFIXM
            GO TO 9999
          ENDIF
          DO nc=1,NCT(nr,nx)  !NCM changed NCM->NCT=actual # variables
            DO no_nynr=1,NYNR(0,0,nc) !loop over global variables
              ny=NYNR(no_nynr,0,nc) !global variable number
              DO iy=1,NIYFIXM
                FIX(ny,iy)=.FALSE.
                YP(ny,iy)=0.0d0
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C Advection-diffusion analytical solution
C DMAL 22-MAY-2002 :
        IF(GENER.AND.ITYP5(nr,nx).EQ.2) THEN ! time integration

          CALL ASSERT(ANAL_CHOICE(nr).NE.0,
     '      '>>Analytic formula not set for this region',ERROR,*9999)

        IF(ITYP2(nr,nx).EQ.3) THEN ! Advection-diffusion
          IF(ANAL_CHOICE(nr).EQ.1) THEN
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              ny1=NYNP(1,1,1,np,0,1,nr)
              ny2=NYNP(2,1,1,np,0,1,nr)
              ny3=NYNP(3,1,1,np,0,1,nr)
              ny4=NYNP(4,1,1,np,0,1,nr)
              YP(ny1,1)=0.0d0
              YP(ny2,1)=0.0d0
              YP(ny3,1)=0.0d0
              YP(ny4,1)=0.0d0
              CALL ANS001(YP(ny4,1),YP(ny2,1),YP(ny3,1),
     '          YP(ny1,1),ANAL_TIME,XP(1,1,1,np),XP(1,1,2,np),ERROR,
     '          *9999) ! calculates analytic solution
              YP(ny1,3)=YP(ny1,1) ! initial conditions
              YP(ny2,3)=YP(ny2,1) ! initial conditions
              YP(ny3,3)=YP(ny3,1) ! initial conditions
              YP(ny4,3)=YP(ny4,1) ! initial conditions
            ENDDO !nonode

            ! determine boundary node and fix them as b.c's
            NPLIST(0)=0
            DO nl=1,NLT
              IF(NEL(0,nl).EQ.1) THEN !line nl is on boundary
                IF(NPL(1,1,nl).EQ.1) THEN  !linear Lagrange basis
                  NN_TOT=2
                ELSE IF(NPL(1,1,nl).EQ.2) THEN !quadratic Lagrange
                  NN_TOT=3
                ELSE IF(NPL(1,1,nl).EQ.3) THEN !cubic Lagrange
                  NN_TOT=4
                ELSE IF(NPL(1,1,nl).EQ.4) THEN !cubic Hermite
                  NN_TOT=2
                ENDIF
                DO nn=1,NN_TOT
                  IF(DOP) THEN
                    WRITE(*,'('' nl='',I4,'' nn='',I2,'' np='',I5)')
     '                nl,nn,NPL(1+nn,1,nl)
                  ENDIF
                  IF(.NOT.INLIST(NPL(1+nn,1,nl),
     '              NPLIST(1),MIN(NPLIST(0),NPM),n1)) THEN
                    NPLIST(0)=NPLIST(0)+1
                    IF(NPLIST(0).LE.NPM) NPLIST(NPLIST(0))=
     '                NPL(1+nn,1,nl)
                    DO iy=1,NIYFIXM
                      ny1=NYNP(1,1,1,NPL(1+nn,1,nl),0,1,nr)
                      ny2=NYNP(2,1,1,NPL(1+nn,1,nl),0,1,nr)
                      ny3=NYNP(3,1,1,NPL(1+nn,1,nl),0,1,nr)
                      FIX(ny1,iy)=.TRUE. ! fix bc value
                      IF(NPL(1,0,nl).EQ.1) THEN
                        FIX(ny2,iy)=.TRUE. ! fix bc derivative
                      ELSEIF(NPL(1,0,nl).EQ.2) THEN
                        FIX(ny3,iy)=.TRUE. ! fix bc derivative
                      ENDIF
                    ENDDO
                  ELSE
                    DO iy=1,NIYFIXM
                      ny2=NYNP(2,1,1,NPLIST(n1),0,1,nr)
                      ny3=NYNP(3,1,1,NPLIST(n1),0,1,nr)
                      IF(NPL(1,0,nl).EQ.1) THEN
                        FIX(ny2,iy)=.TRUE. ! fix bc derivative
                      ELSEIF(NPL(1,0,nl).EQ.2) THEN
                        FIX(ny3,iy)=.TRUE. ! fix bc derivative
                      ENDIF
                    ENDDO !iy
                  ENDIF
                ENDDO !nn
              ENDIF !NEL
            ENDDO !nl

            CALL ASSERT(NPLIST(0).LE.NPM,'>>Increase NPM',
     '        ERROR,*9999)

          ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              ny1=NYNP(1,1,1,np,0,1,nr)
              ny2=NYNP(2,1,1,np,0,1,nr)
              YP(ny1,1)=0.0d0
              YP(ny2,1)=0.0d0
              CALL ANS002(YP(ny2,1),YP(ny1,1),0.0d0,
     '          XP(1,1,1,np),ERROR,*9999) ! calculates analytic solution
              YP(ny1,3)=YP(ny1,1) ! initial conditions
              YP(ny2,3)=YP(ny2,1) ! initial conditions
            ENDDO !nonode

            ! determine boundary node and fix them as b.c's
            NPLIST(0)=0
            DO nonode=1,NPNODE(0,nr)
              DO iy=1,NIYFIXM
                np=NPNODE(nonode,nr)
                ny1=NYNP(1,1,1,np,0,1,nr)
                IF(NENP(np,0,nr).EQ.1) THEN
                  FIX(ny1,iy)=.TRUE. ! fix bc value
                ENDIF
              ENDDO !iy
            ENDDO !nonode
            CALL ASSERT(NPLIST(0).LE.NPM,'>>Increase NPM',
     '        ERROR,*9999)
          ELSE
            ERROR='>>Not implemented for this analytic problem'
            GOTO 9999
          ENDIF
        ELSE
          ERROR='>>Not implemented for this problem'
          GOTO 9999
        ENDIF


      ELSEIF(GENER) THEN

C***      If we are dealing with more than one region then
C***      only set up the initial conditions for the nodes
C***      in the first and last regions that don't belong to
C***      an interface

         CALL ASSERT(ANAL_CHOICE(nr).NE.0,
     '      '>>Analytic formula not set for this region',ERROR,*9999)

          DO nc=1,NCT(nr,nx)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np)
                nh=nh_loc(nhx,nx)
                DO nv=1,NVHP(nh,np,nc)
                  DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                    IF(nc.EQ.1) THEN
C cpb 22/1/97 Changing test
C                      IF(.NOT.ALL_REGIONS.OR.(nr.EQ.1.AND.
C     '                  NP_INTERFACE(np,0).EQ.1).OR.(nr.EQ.NRT.AND.
C     '                  NP_INTERFACE(np,0).EQ.1)) THEN
CC***                    Only dealing with one region, or dealing with
CC***                    them all and np is in first reigon and no
CC***                    interface or in the last region and no interface
                      IF((.NOT.ALL_REGIONS).OR.
     '                  (NP_INTERFACE(np,0).EQ.1)) THEN
C                       If we are on a 'free' surface or we want to
C                       specify a particular region(s)
                        IF(nk.EQ.1) THEN
                          ny1=NYNP(1,nv,nh,np,0,1,nr)
C CPB 23/1/96 Changing dipole cases to actually use a dipole
C                          IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).GE.5.AND.
C     '                      ANAL_CHOICE(nr).LE.8) THEN
C                            R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
C     '                        XP(1,nv,3,np)**2)
C                            IF(DABS(R-MESH3_RAD(NSPHERES))/
C     '                        MESH3_RAD(NSPHERES).GT.SPHERE_RAD_TOL)
C     '                        THEN
CC***                          Not on outer sphere
C                              YP(ny1,1)=YP(ny1,7)
C                              FIX(ny1,1)=.TRUE.
C                            ENDIF
                          IF(NJT.EQ.2) THEN
                            IF(ANAL_CHOICE(nr).EQ.5) THEN !single circle
                              IF(np.EQ.ANAL_FIXEDNODE) THEN
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !mul circ
                              R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
                              IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                          MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL)
     '                          THEN
                                IF(np.EQ.ANAL_FIXEDNODE) THEN
C                                 First node on outer circle
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                              R=XP(1,nv,1,np)
                              IF(DABS(R-MESH3_RAD(1))/
     '                          MESH3_RAD(1).LE.SPHERE_RAD_TOL) THEN
C                               Inner circle
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.10.OR.
     '                          ANAL_CHOICE(nr).EQ.11) THEN
                              IF(DABS(XP(1,1,1,np)).LE.1.0d-1.OR.
     '                          DABS(XP(1,1,2,np)).LE.1.0d-1.OR.
     '                          DABS(XP(1,1,2,np)-ANISO_H).LE.1.0d-1.OR.
     '                          DABS(XP(1,1,1,np)-ANISO_L).LE.1.0d-1)
     '                          THEN
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ELSE IF(NJT.EQ.3) THEN
                            IF(ANAL_CHOICE(nr).EQ.12.OR.
     '                        ANAL_CHOICE(nr).EQ.14) THEN !single sphere
                              IF(np.EQ.ANAL_FIXEDNODE) THEN
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.13.OR.
     '                      ANAL_CHOICE(nr).EQ.15) THEN !mul sph
                              R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '                          XP(1,nv,3,np)**2)
                              IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                          MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL)
     '                          THEN
                                IF(np.EQ.ANAL_FIXEDNODE) THEN
C                                 First node on outer sphere
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ELSE
                          ny1=NYNP(nk,nv,nh,np,0,1,nr)
C                          IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).GE.5.AND.
C     '                      ANAL_CHOICE(nr).LE.8) THEN
C                            R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
C     '                        XP(1,nv,3,np)**2)
C                            IF(DABS(R-MESH3_RAD(NSPHERES))/
C     '                        MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL)
C     '                        THEN
C                              YP(ny1,1)=YP(ny1,7)
C                              FIX(ny1,1)=.TRUE.
C                            ENDIF
                          IF(NJT.EQ.2) THEN
                            IF(ANAL_CHOICE(nr).GE.5.AND.
     '                        ANAL_CHOICE(nr).LE.6) THEN
C                             Do Nothing
                            ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                              R=XP(1,nv,1,np)
                              IF(DABS(R-MESH3_RAD(1))/
     '                          MESH3_RAD(1).LE.SPHERE_RAD_TOL) THEN
C                               Inner circle
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE IF(ANAL_CHOICE(nr).EQ.10.OR.
     '                          ANAL_CHOICE(nr).EQ.11) THEN
                              IF(DABS(XP(1,1,1,np)).LE.1.0d-1.OR.
     '                          DABS(XP(1,1,2,np)).LE.1.0d-1.OR.
     '                          DABS(XP(1,1,2,np)-ANISO_H).LE.1.0d-1.OR.
     '                          DABS(XP(1,1,1,np)-ANISO_L).LE.1.0d-1)
     '                          THEN
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ELSE IF(NJT.EQ.3) THEN
                            IF(ANAL_CHOICE(nr).GE.12.AND.
     '                        ANAL_CHOICE(nr).LE.15) THEN
C                             Do Nothing
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE
                      ny1=NYNP(nk,nv,nh,np,0,2,nr)
C                      IF(NJT.EQ.3.AND.ANAL_CHOICE(nr).GE.5.AND.
C     '                  ANAL_CHOICE(nr).LE.8) THEN
CC***                    Dipole solutions - set flux to zero on outer
CC***                    boundary
C                        R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
C     '                    XP(1,nv,3,np)**2)
C                        IF(DABS(R-MESH3_RAD(NSPHERES))/
C     '                    MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
CC***                      On outer surface
C                          ny2=NYNP(nk,nv,nh,np,0,2,nr)
C                          FIX(ny2,1)=.TRUE.
C                          ne=NENP(np,1)
C                          CALL DIPOLE_EVALUATE(nk+3,np,NP_INTERFACE,
C     '                      nr,nx,CE(1,ne),
C     '                      XP(1,nv,1,np),
C     '                      YP(ny2,1),XP(1,nv,1,np),XP(1,nv,2,np),
C     '                      XP(1,nv,3,np),ERROR,*9999)
C                        ENDIF
                      IF(NJT.EQ.2) THEN
                        IF(ANAL_CHOICE(nr).EQ.5) THEN !single circle
                          IF(FIX_ZERO) THEN
                            YP(ny1,1)=YP(ny1,7)
                            FIX(ny1,1)=.TRUE.
                          ELSE
                            IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1) THEN
                              IF(POT_BC_TYPE.EQ.2) THEN !fix deriv,val
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ELSE IF(ANAL_CHOICE(nr).EQ.6) THEN !mul circles
                          R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2)
                          IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                      MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
                            IF(FIX_ZERO) THEN
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ELSE
                              IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1)
     '                          THEN
                                IF(POT_BC_TYPE.EQ.2) THEN !fix deri,val
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ELSE
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ENDIF
                          ENDIF
                        ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                          R=XP(1,nv,1,np)
                          IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                      MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
C***                        On outer surface
                            ny2=NYNP(nk,nv,nh,np,0,2,nr)
                            YP(ny2,1)=0.0d0
                            FIX(ny2,1)=.TRUE.
                          ENDIF
                        ELSE IF(ANAL_CHOICE(nr).EQ.10.OR.
     '                      ANAL_CHOICE(nr).EQ.11) THEN
C                          IF(DABS(XP(1,1,1,np)).LE.1.0d-1) THEN
C                            ny2=NYNP(nk,nv,nh,np,0,2,nr)
C                            YP(ny2,1)=0.0d0
C                            FIX(ny2,1)=.TRUE.
C                          ENDIF
                        ENDIF
                      ELSE IF(NJT.EQ.3) THEN
                        IF(ANAL_CHOICE(nr).EQ.12.OR.
     '                    ANAL_CHOICE(nr).EQ.14) THEN
                          IF(FIX_ZERO) THEN
                            YP(ny1,1)=YP(ny1,7)
                            FIX(ny1,1)=.TRUE.
                          ELSE
                            IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1) THEN
                              IF(POT_BC_TYPE.EQ.2) THEN !fix deriv,val
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ELSE
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ENDIF
                          ENDIF
                        ELSE IF(ANAL_CHOICE(nr).EQ.13.OR.
     '                      ANAL_CHOICE(nr).EQ.15) THEN
                          R=DSQRT(XP(1,nv,1,np)**2+XP(1,nv,2,np)**2+
     '                      XP(1,nv,3,np)**2)
                          IF(DABS(R-MESH3_RAD(NSPHERES))/
     '                      MESH3_RAD(NSPHERES).LE.SPHERE_RAD_TOL) THEN
                            IF(FIX_ZERO) THEN
                              YP(ny1,1)=YP(ny1,7)
                              FIX(ny1,1)=.TRUE.
                            ELSE
                              IF(np.EQ.ANAL_FIXEDNODE.AND.nk.EQ.1)
     '                          THEN
                                IF(POT_BC_TYPE.EQ.2) THEN !fix der,val
                                  YP(ny1,1)=YP(ny1,7)
                                  FIX(ny1,1)=.TRUE.
                                ENDIF
                              ELSE
                                YP(ny1,1)=YP(ny1,7)
                                FIX(ny1,1)=.TRUE.
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF !nc
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nh
            ENDDO !nonode (np)
          ENDDO !nc

        ELSE !Not generating initial conditions

c cpb 1/5/95 Replacing Fourier analysis with Quasi-static analysis
c        IF(ITYP5(nr,nx).EQ.2.OR.ITYP5(nr,nx).EQ.4.OR.
          IF((ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).NE.11).
     '      OR.ITYP5(nr,nx).EQ.5.
     '      OR.(ITYP5(nr,nx).EQ.1.AND.ITYP6(nr,nx).EQ.2)) THEN
            FORMAT='('' Specify whether initial solution is [1]:'''//
     '        '/''   (1) Zero'''//
     '        '/''   (2) Read in'''//
     '        '/''  *(3) Restart from previous solution'''//
     '        '/$,''    '',I1)'
            IF(IOTYPE.EQ.3) THEN !always write out initial conditions
              KTYP5=2
              IDATA(1)=KTYP5
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP5=IDATA(1)
            IF(KTYP5.EQ.3) THEN
              RESTAR=.TRUE.
              CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
            ELSE
              RESTAR=.FALSE.
            ENDIF
          ELSE IF(ITYP2(nr,nx).EQ.11.AND.ITYP3(nr,nx).EQ.2)THEN
            !water vapour and heat transfer/transport
            FORMAT='('' Specify whether initial temperature is [1]:'''//
     '        '/''   (1) Zero'''//
     '        '/''   (2) Equal to body temperature'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=1
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP5=IDATA(1)
            RESTAR=.FALSE.
          ELSE
            KTYP5=1
          ENDIF

C KAT 26Aug97:  Everything is already initiallized to zero
C          IF(KTYP5.EQ.1) THEN !Initial solution is zero
C            DO nc=1,NCM
C              DO no_nynr=1,NYNR(0,0,nc) !loop over global variables
C                ny=NYNR(no_nynr,0,nc) !global variable number
C                YP(ny,3)=0.0d0
C              ENDDO
C            ENDDO
C          ELSE IF(KTYP5.EQ.3) THEN
CC***        Restart from previous solution
CC??? cpb 6/12/94 Dump YP(ny,1) into YP(ny,4) ???
C          ENDIF

          IF(KTYP5.EQ.2) THEN !read in initial conditions
            LOOPT=4
          ELSE
            IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
              LOOPT=3
            ELSE
              LOOPT=2
C AJP 3/2/94 Only need two question for linear elliptic probs
            ENDIF
          ENDIF

          IF((ITYP6(nr,nx).EQ.1.AND.(ITYP2(nr,nx).NE.5.
     '      OR.ITYP3(nr,nx).NE.2).AND.ITYP2(nr,nx).NE.11).
     '      OR.ITYP5(nr,nx).EQ.5) THEN
C***        Linear and wavefront path problems
            RDEFLT(1)=1.1d6
            NPMC=0

            DO loop=1,LOOPT !    ------- begin loop -------
              nonode=0
              FLUXBC=.FALSE.
              IF(loop.EQ.1) THEN
                FORMAT='('' Essential boundary conditions'//
     '            ' defined at nodes:'')'
                nc=1
              ELSE IF(loop.EQ.2) THEN
                FORMAT='(/'' Flux boundary conditions'//
     '            ' defined at nodes:'')'
                FLUXBC=.TRUE.
                nc=2
              ELSE IF(loop.EQ.3) THEN
C KAT 28Jul97:  Incremented boundary conditions not used
                GOTO 710
C                IF(ITYP6(nr,nx).EQ.1) THEN !linear
C                  GOTO 710
C                ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
C                  IF(ITYP7(nr,nx).EQ.1) THEN !nonlinear equation
C                    FORMAT='(/'' Incremented boundary conditions:'')'
C                    nc=1
C                  ENDIF
C                ENDIF !ityp6
              ELSE IF(loop.EQ.4)THEN
                FORMAT='(/'' Initial conditions defined at nodes:'')'
                nc=1
              ENDIF
              CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)

              IF(IOTYPE.EQ.3) KTYPMBC=1
 6000         IF(FLUXBC) THEN
C PJH 13/9/94 Adding point values for flux boundary conditions
                FORMAT='('' Specify whether the flux values are '
     '            //'[Exit]:'''//
     '            '/''   (1) Integrated (W/m^2)'''//
     '            '/''   (2) Point values (Watts)'''//
     '            '/''   (3) Proportional to u'''//
     '            '/$,''    '',I1)'
CC cpb 12/2/98 Only write out integrated fluxes
               IF(IOTYPE.EQ.3) IDATA(1)=KTYPMBC
C                IF(IOTYPE.EQ.3) IDATA(1)=1
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '            0,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '            *9999)
                IF(IOTYPE.NE.3) KTYPMBC=IDATA(1)
              ENDIF !fluxbc

              IF(.NOT.FLUXBC.OR.(FLUXBC.AND.KTYPMBC.NE.0)) THEN !ask for nodes - otherwise exit
 6100           FORMAT='($,'' Enter node #s/name [EXIT]: '',I5)'
                NO_DATA=1
                IF(IOTYPE.EQ.3) THEN
                  NEXTNODE=.TRUE.
                  DO WHILE(NEXTNODE)
                    nonode=nonode+1
                    IF(nonode.LE.NPNODE(0,nr)) THEN
                      np=NPNODE(nonode,nr)
C KAT 27Aug97:        Check that there is info to be output for node.
                      IF(loop.EQ.4.AND.ITYP6(nr,nx).EQ.1) THEN
                        NEXTNODE=.FALSE. !output all initial conditions
                      ELSE
                        DO nhx=1,NHP(np) !variables
                          nh=nh_loc(nhx,nx)
                          DO nv=1,NVHP(nh,np,nc) !versions
                            DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                              ny=NYNP(nk,nv,nh,np,0,nc,nr) !derivatives
                              IF(loop.LT.3) THEN !boundary conds.
                                IF(FIX(ny,1)) NEXTNODE=.FALSE.
                              ELSE !IF(loop.EQ.4) THEN !initial conds.
                                IF(.NOT.FIX(ny,1)) NEXTNODE=.FALSE.
                              ENDIF
                            ENDDO !nk
                          ENDDO !nv
                        ENDDO !nh
                      ENDIF
                      IDATA(1)=np
                      NPLIST(0)=1
                      NPLIST(1)=np
                    ELSE
                      NEXTNODE=.FALSE.
                      NO_DATA=0
                      IDATA(1)=0 !for default exit condition
                    ENDIF
                  ENDDO !WHILE(NEXTNODE)
                ENDIF !iotype=3

 6500           CDATA(1)='NODES' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '            FILEIP,FORMAT,NO_DATA,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '            IDATA,IZERO,0,NPT(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '            RMAX,INFO,ERROR,*9999)
                IF(IDATA(1).NE.0) THEN !not default exit
                  IF(IOTYPE.NE.3) THEN
                    NPLIST(0)=IDATA(0)
                    DO n=1,IDATA(0)
                      NPLIST(n)=IDATA(n)
                      np=IDATA(n)
                      IF(.NOT.INLIST(np,NPNODE(1,nr),
     '                  NPNODE(0,nr),N1)) THEN
                        WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '                    //'belong to the current region'')') np
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                        GOTO 6500
                      ENDIF
                    ENDDO !n
                  ENDIF !iotype.NE.3

C***              Define bdry condition for first node in group
                  np=NPLIST(1) !rest of group is filled at end of nh loop
                  DO nhx=1,NHP(np)
                    nh=nh_loc(nhx,nx)
                    WRITE(CHAR1,'(I1)') nhx
                    FORMAT='('' Dependent variable number '//CHAR1(1:1)
     '                //' :'')'
                    CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &                RMIN,RMAX,INFO,ERROR,*9999)
                    DO nv=1,NVHP(nh,np,nc)
                      IF(NVHP(nh,np,nc).GT.1) THEN
                        WRITE(CHAR1,'(I2)') nv
                        FORMAT='('' For version number '//CHAR1(1:2)
     '                    //':'')'
                        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,
     '                    IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '                    RMIN,RMAX,INFO,ERROR,*9999)
                      ENDIF
                      DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
                        IF(nk.EQ.1) THEN
                          FORMAT='($,'' The value of the dependent '
     '                     //'variable is [no b.c.]: '',G12.5,1X,G12.5)'
                        ELSE IF(nk.GT.1) THEN
                          WRITE(CHAR1,'(I1)') nk-1
                          FORMAT='($,'' The value of derivative number '
     '                      //CHAR1(1:1)//' is [no b.c.]: '',G12.5,1X,'
     '                      //'G12.5)'
                        ENDIF
                        ny=NYNP(nk,nv,nh,np,0,nc,nr)
                        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPINI3_2)
                          WRITE(OP_STRING,'('' GETNYP: ny='',I6)') ny
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(IPINI3_2)
                        ENDIF
C AJP 18-1-95 The above used to be ny=NYNP(nk,nv,nh,np,2,nc,nr)
C which I think is wrong.
                        IF(IOTYPE.EQ.3) THEN
                          IF(loop.LT.3) THEN
                            IF(FIX(ny,1)) THEN
                              RDATA(1)=YP(ny,1)
                              NO_DATA=1
                            ELSE
                              NO_DATA=0
                            ENDIF
C KAT 28Jul97:             loop is never 3 in this structure
C                          ELSE IF(loop.EQ.3) THEN
C                            RDATA(1)=RDEFLT(1)
C                            RDATA(2)=RDEFLT(1)
C                            RDATA(3)=RDEFLT(1)
                          ELSE IF(loop.EQ.4) THEN
C KAT 26Aug97:  Changed output from residual to current solution
                            RDATA(1)=YP(ny,1)
                          ENDIF
                        ENDIF !iotype=3

                        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                    FILEIP,FORMAT,NO_DATA,ADATA,ADEFLT,CDATA,
     '                    CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,
     '                    LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     '                    *9999)

                        IF(IOTYPE.NE.3) THEN
                          IF(RDATA(1).NE.RDEFLT(1)) THEN
                            IF(loop.EQ.1) THEN !essential bcs
                              YP(ny,1)=RDATA(1)
                              FIX(ny,1)=.TRUE.
                              IF(DOP) THEN
C$OMP CRITICAL(IPINI3_3)
                                WRITE(*,'('' FIX('',I4,'',1)='',L1)')
     '                            ny,FIX(ny,1)
C$OMP END CRITICAL(IPINI3_3)
                              ENDIF
C KAT 27Jan99: FIX(ny,5) not used
C                              IF(FILEIP) THEN
C                                FIX(ny,5)=.TRUE.
C                              ELSE IF(.NOT.FILEIP) THEN
C                                FIX(ny,5)=.FALSE.
C                              ENDIF
C KAT 28Jul97:                 loop is never 3 in this structure
C                            ELSE IF(loop.EQ.3) THEN !nonlinear increms
C                              YP(ny,2)=RDATA(1)
C                              FIX(ny,2)=.TRUE.
C                              IF(FILEIP) THEN
C                                FIX(ny,5)=.TRUE.
C                              ELSE IF(.NOT.FILEIP) THEN
C                                FIX(ny,5)=.FALSE.
C                              ENDIF
                            ELSE IF(loop.EQ.2) THEN !flux conditions (Note:nc=2)
C PJH 19May99 Add bc for flux=au
                              IF(KTYPMBC.LE.2) THEN !flux value is specified
                                FIX(ny,1)=.TRUE.   !nc=2 ny FIX'd for niy=1
                                FIX(ny,2)=.FALSE.  !nc=2 ny not FIX'd for niy=2
                                YP(ny,1)=RDATA(1)
                              ELSE IF(KTYPMBC.EQ.3) THEN !flux proportional to u
                                FIX(ny,1)=.TRUE.   !nc=2 ny FIX'd for niy=1
                                FIX(ny,2)=.TRUE.   !nc=2 ny FIX'd for niy=2
                                ny1=NYNP(nk,nv,nh,np,0,1,nr) !is ny# for nc=1
                                YP(ny1,2)=RDATA(1) !is coeff for u dependence
                              ENDIF !KTYPMBC

                            ELSE IF(loop.EQ.4) THEN !initial conditions
                              YP(ny,3)=RDATA(1)
C                             Copy initial condition to current solution
                              IF(.NOT.FIX(ny,1)) YP(ny,1)=RDATA(1)
                            ENDIF !loop
                          ELSE IF(loop.LE.2) THEN !essential or flux bcs
                            FIX(ny,1)=.FALSE.
                          ENDIF !rdata.ne.rdeflt/loop<=2
                        ENDIF !iotype.ne.3
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !nh

C***              Apply bdry conditions to rest of group
                  DO n=2,NPLIST(0)
                    np=NPLIST(n)
                    IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPINI3_4)
                      WRITE(OP_STRING,'('' Group data: n='',I3,'' np='''
     '                  //',I5)') n,np
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(IPINI3_4)
                    ENDIF !DOP
                    DO nhx=1,MIN(NHP(np),NHP(NPLIST(1)))
                      nh=nh_loc(nhx,nx)
                      DO nv=1,MIN(NVHP(nh,np,nc),NVHP(nh,NPLIST(1),nc))
                        DO nk=1,MIN(MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1),
     '                    NKH(nh,NPLIST(1),nc))
                          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPINI3_5)
                          WRITE(OP_STRING,'('' nh='',I1,'' nv='',I1,'
     &                      //''' nk='',I1,'' np='',I3)') nh,nv,nk,np
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(IPINI3_5)
                          ENDIF !DOP
                          ny=NYNP(nk,nv,nh,np,0,nc,nr)  !Note:nc=2
                          ny_first=NYNP(nk,nv,nh,NPLIST(1),0,nc,nr)
!                       WRITE(OP_STRING,'(''1 ny='',I4,'' 1st='',I4)')
!      '                  ny,ny_first
!                       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          FIX(ny,1)=FIX(ny_first,1)
                          YP(ny,1)=YP(ny_first,1)
C PJH 19May99 Add bc for flux=au
                          IF(KTYPMBC.EQ.3) THEN  !flux proportional to u
                            FIX(ny,2)=FIX(ny_first,2)    !nc=2 ny FIX'd
                            ny1=NYNP(nk,nv,nh,np,0,1,nr) !is ny# for nc=1
                            ny1_first=NYNP(nk,nv,nh,NPLIST(1),0,1,nr) !is ny_first# for nc=1
                            YP(ny1,2)=YP(ny1_first,2)
                          ENDIF !KTYPMBC
                          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(IPINI3_6)
                             WRITE(OP_STRING,'('' FIX('',I4,'
     '                        //''',1)='',L1)') ny,FIX(ny,1)
                             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(IPINI3_6)
                          ENDIF
                        !ENDDO !nv2
                        ENDDO !nk
!                           WRITE(OP_STRING,'(''1a nh='',I1,'' nv='',I1,'
!      &                      //''' nk='',I1,'' np='',I3)') nh,nv,nk,np
!                           CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDDO !nv
                      ! Copy to all the other versions in the target
                      ! these will get overwritten on the next loop
                      ! if they are defined in the first node, but
                      ! if they are not defined in the first node then
                      ! they will get the value of the previous version
                      ! that was defined.
                      nv_last=nv-1 ! last version that was defined for the group
                      DO nv=nv,NVHP(nh,np,nc)
                      WRITE(OP_STRING,'(''>>Warning: for node '',I4,'
     &                  //''' copying version '',I1,'' to version '','
     &                  //'I1)') np,nv_last,nv
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        DO nk=1,MIN(MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1),
     &                    NKH(nh,NPLIST(1),nc))
                          ny_first=NYNP(nk,nv_last,
     &                      nh,NPLIST(1),0,nc,nr)
                          ny=NYNP(nk,nv,nh,np,0,nc,nr)  !Note:nc=2
!                       WRITE(OP_STRING,'(''2 ny='',I4,'' 1st='',I4)')
!      '                  ny,ny_first
!                       CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          FIX(ny,1)=FIX(ny_first,1)
                          YP(ny,1)=YP(ny_first,1)
                          IF(KTYPMBC.EQ.3) THEN  !flux proportional to u
                            FIX(ny,2)=FIX(ny_first,2)    !nc=2 ny FIX'd
                            ny1=NYNP(nk,nv,nh,np,0,1,nr) !is ny# for nc=1
                            ny1_first=NYNP(nk,nv_last,
     &                        nh,NPLIST(1),0,1,nr) !is ny_first# for nc=1
                            YP(ny1,2)=YP(ny1_first,2)
                          ENDIF !KTYPMBC
                        ENDDO !nk
                      ENDDO !nv
                    
                    ENDDO !nh
                  ENDDO !n

                  GO TO 6100 !for more nodes
                ENDIF !idata(1).NE.0

                IF(loop.EQ.1) THEN !essential bdry conditions
C CPB ??? 6/12/94
                  DO no_nynr=1,NYNR(0,0,nc) !loop over global variables
                    ny=NYNR(no_nynr,0,nc) !global variable number
                    YP(ny,4)=YP(ny,1) !to inc any essential bcs in YP(4)
                  ENDDO

                ELSE IF(loop.EQ.2.AND.IOTYPE.NE.3.AND.KTYPMBC.EQ.2) THEN
C***              Compute integrated flux values
                  DO no_nynr=1,NYNR(0,0,nc) !loop over global variables
                    ny=NYNR(no_nynr,0,nc) !global variable number
                    YP(ny,5)=YP(ny,1) !temporarily store fluxes in y(5)
                    YP(ny,1)=0.d0
                  ENDDO
                  nh=NH_LOC(1,nx)
                  CALL YPZP(5,NBH,NEELEM,NHE,NHP,NKH,NPNODE,nr,
     '              NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    nb=NBH(nh,2,ne)
C DPN 08 February 1999 - adding integrated point fluxes for 2D
c                    CALL ASSERT(NIT(nb).EQ.3,
c     '                '>>Point fluxes not implemented for 2D',ERROR,*9999)
                    IF(NIT(nb).EQ.3) THEN !3D
                      DO nef=1,NFE(nb)
                        nf=NFF(nef,ne)
                        IF(nf.NE.0) THEN
                          FACE_FIXED=.TRUE.
                          CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     '                      NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,
     '                      NPNE(1,1,ne),NPNF,nr,NVJE(1,1,1,ne),NVJF,
     '                      SE(1,1,ne),SF,ERROR,*9999)
                          nbf=NBHF(nh,2,nf)
                          nv=1 !? PJH 30Mar95 nv was undefined
                          DO nn=1,NNT(nbf)
                            np=NPNF(nn,nbf)
                            DO nk=1,NKT(nn,nbf)
                              ny=NYNP(nk,nv,nh,np,0,2,nr)
                              IF(.NOT.FIX(ny,1)) FACE_FIXED=.FALSE.
                            ENDDO !nk
                          ENDDO !nn
                          IF(FACE_FIXED) THEN !update fluxes for these nodes
                            CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,
     '                        NVJF,SF,XA(1,1,1),XE,XP,ERROR,*9999)
                            CALL CALC_FACE_INFORMATION_DEP(NBH(1,2,ne),
     '                        NBHF(1,2,nf),nef,NHE(ne),
     '                        NKHE(1,1,1,ne),NKEF,NKHF,NNF,NPNE(1,1,ne),
     '                        NPNF,NVHE(1,1,1,ne),NVHF,nx,SE(1,1,ne),SF,
     '                        ERROR,*9999)
                            CALL ZPZE(NBHF(1,1,nf),2,NHE(ne),
     '                        NKHF,NPF(1,nf),NPNF,nr,NVHF,NW(ne,1),nx,
     '                        CURVCORRECT(1,1,1,ne),SF,ZA(1,1,1,ne),
     '                        ZE,ZP,ERROR,*9999)
                            ns=0
                            DO nn=1,NNT(nbf)
                              DO nk=1,NKT(nn,nbf)
                                ns=ns+1
                                SUM=0.0d0
                                DO ng=1,NGT(nbf)
                                  CALL XEXG(NBJF(1,nf),ng,nr,PG,XE,XG,
     '                              ERROR,*9999)
                                  CALL XGMG(0,NIT(nbf),NBJF(1,nf),nr,
     '                              DXIX,GL,GU,RG(ng),XG,ERROR,*9999)
                                  SUM=SUM+PG(ns,1,ng,nbf)*SF(ns,nbf)*
     '                              PXI(IBT(1,1,nbf),IDO(1,1,0,nbf),
     '                              INP(1,1,nbf),nbf,1,XIG(1,ng,nbf),
     '                              ZE(1,nh))*RG(ng)*WG(ng,nbf)
                                ENDDO !ng
                                np=NPNF(nn,nbf)
                                ny=NYNP(nk,nv,nh,np,0,2,nr)
                                YP(ny,1)=YP(ny,1)+SUM
                              ENDDO !nk
                            ENDDO !nn
                          ENDIF !FACE_FIXED
                        ENDIF !nf.NE.0
                      ENDDO !nef
                    ELSE IF(NIT(nb).EQ.2) THEN !2D
                      nv = 1 ! DPN ???????
                      !loop over local lines of current element
                      DO kl=1,NLE(nb)
                        nl=NLL(kl,ne) !set the global line number
                        IF(nl.NE.0) THEN
                          LINE_FIXED = .TRUE.
                          !loop over all nodes on the line, checking if
                          !the ny's are fixed
                          DO nn=1,NNL(0,kl,nb)
                            np=NPL(nn+1,1,nl)
                            DO nk=1,NKT(nn,nb)
                              ny = NYNP(nk,nv,nh,np,0,nc,nr)
                              IF (.NOT.FIX(ny,1)) LINE_FIXED = .FALSE.
                            ENDDO !nk
                          ENDDO !nn
                          IF(LINE_FIXED) THEN !update the fluxes for
                                              !these nodes
                            DO nn=1,NNL(0,kl,nb)
                              DO nk=1,NKT(nn,nb)
                                nt=2
                                CALL CALC_CONTRIB_COEFF(COEFF,nb,nb,nt,
     '                            NPL(1,0,nl),PG,WG,ERROR,*9999)
                                np=NPL(nn+1,1,nl)
                                ny=NYNP(nk,nv,nh,np,0,nc,nr)
                                SUM=COEFF(nn)*DL(3,nl)*YP(ny,5)
                                YP(ny,1)=YP(ny,1)+SUM
                              ENDDO !nk
                            ENDDO !nn
                          ENDIF !LINE_FIXED
                        ENDIF !nl.NE.0
                      ENDDO !kl
                    ELSE
                      CALL ASSERT(.FALSE.,
     '                  '>>Point fluxes not implemented for this case',
     '                  ERROR,*9999)
                    ENDIF !NIT(nb)
                  ENDDO !noelem
                  DO no_nynr=1,NYNR(0,0,nc) !loop over global variables
                    ny=NYNR(no_nynr,0,nc) !global variable number
                    YP(ny,5)=0.0d0
                  ENDDO

                ENDIF !loop=1 or 2
                IF(IOTYPE.EQ.3) KTYPMBC=KTYPMBC-1
                IF(FLUXBC) GOTO 6000 !for more flux bdry conditions
              ENDIF !FLUXBC.AND.KTYPMBC.NE.0
 710          CONTINUE
            ENDDO !loop       ------- end loop -------

          ELSE IF(ITYP6(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.11)THEN
            IF(ITYP3(nr,nx).LE.2)THEN
              WRITE(OP_STRING,'('' Use define initial;c '')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6.OR.
     &        ITYP3(nr,nx).EQ.7)THEN
C... Pulmonary capillaries (KSB:11/2001) & arterioles/venules (KSB:02/2004) &
C... gas exchange (AJS 8/2007)
              FORMAT='(/$,'' Enter the initial hematocrit value '
     '          //' [0.4]: '',D12.4)'
              RDEFLT(1)=0.4d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.4d0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,1.d0,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) INITIAL_HD=RDATA(1)
            ENDIF !ITYP3
            IF(ITYP3(nr,nx).EQ.3.OR.ITYP3(nr,nx).EQ.6)THEN
C... Pulmonary capillaries & arterioles/venules
              FORMAT='('' Set lung pressure conditions at? [1]:'''//
     '          '/''   (1) Rest '''//
     '          '/''   (2) Inspiration '''//
     '          '/''   (3) Passive expiration '''//
     '          '/$,''    '',I1)'
              IDEFLT(1)=1
              IF(IOTYPE.EQ.3) IDATA(1)=1
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.NE.3)THEN !Setting default values
                IF(IDATA(1).EQ.1) THEN !At rest
                  RDEFLT(1)=-5.d0 !P_pleural
                  PALV_DEFLT=0.d0 !default alveolar pressure
                ELSE IF(IDATA(1).EQ.2) THEN !Inspiration
                  RDEFLT(1)=-7.5d0
                  PALV_DEFLT=-1.d0
                ELSE IF(IDATA(1).EQ.3) THEN !Expiration
                  RDEFLT(1)=-5.d0
                  PALV_DEFLT=1.d0
                ENDIF
              ENDIF
              WRITE(CHAR2,'(D12.4)') RDEFLT(1)
              FORMAT='(/$,'' Enter the pleural pressure (cm H2O) '
     '          //' ['//CHAR2//']: '',D12.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-20.d0,5.d0,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) P_pleural=RDATA(1)
              RDEFLT(1)=PALV_DEFLT
              WRITE(CHAR2,'(D12.4)') RDEFLT(1)
              FORMAT='(/$,'' Enter the alveolar pressure (cm H2O) '
     '          //' ['//CHAR2//']: '',D12.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-10.d0,5.d0,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) P_alveolar=RDATA(1)
              FORMAT='(/$,'' Enter the % of segments blocked by WBCs: '
     '          //' [0.0]: '',D12.4)'
              RDEFLT(1)=0.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=0.0d0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,100.d0,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) BLOCKED=RDATA(1)/100.d0 !% -> fraction
              MANUAL=.FALSE.
              FORMAT='('' Specify boundary element by: [1]:'''//
     '          '/''   (1) Automatically finding ne #'''//
     '          '/''   (2) Manually defining ne #'''//
     '          '/$,''    '',I1)'
              IDEFLT(1)=1
              IF(IOTYPE.EQ.3) IDATA(1)=1
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '          1,2,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IOTYPE.EQ.3) IDATA(1)=1
              IF(IDATA(1).EQ.2) MANUAL=.TRUE. !so can manually define ne #
              NUMBC=0
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBH(1,1,ne) !nh & nc=1
                CE(nm_Hd,ne)=INITIAL_HD !initial hematocrit value for ne
                IF(NENP(NPNE(1,nb,ne),0,nr).EQ.1.OR.
     '            NENP(NPNE(2,nb,ne),0,nr).EQ.1) THEN !node not connected, needs BC
C... also define initial diameter for feed vessels, larger than
C... capillaries
                  NUMBC=NUMBC+1 !number BC's required for solution
                  IF(ITYP3(nr,nx).EQ.3) THEN
                    CE(nm_C0,ne)=0.04d0 !mm measurement approx from Huang et al
                    CE(nm_a0,ne)=0.023d0 !mm Huang - diameter of arteriole=0.02mm
                    perimeter=CE(nm_C0,ne)+CE(nm_a0,ne)
                    CE(nm_b,ne)=CE(nm_a0,ne)/2.d0
                    CE(nm_a,ne)=DSQRT(DABS(perimeter**2.d0/
     &                (2.d0*PI**2.d0)-CE(nm_b,ne)**2.d0))
                    CE(nm_Dh,ne)=4.d0*PI*CE(nm_a,ne)*CE(nm_b,ne)
     &                /perimeter!(mm) 4Ac/perimeter
                    CE(nm_length,ne)=CE(nm_length,ne)/scale_factor
C                    !NB/ this undo the previous *scale_factor in CAP_NE!
                  ENDIF
                  IF(.NOT.MANUAL) THEN
                    WRITE(CHAR,'(I5)') ne !specify boundary condition
                    FORMAT=
     '                '('' Enter the type of boundary condition for ne'
     '                //CHAR(1:5)//' [1]:'''//
     '                '/''   (1) Inlet pressure (cm H2O)'''//
     '                '/''   (2) Outlet pressure (cm H2O)'''//
     '                '/''   (3) Flow (mm**3/s)'''//'/$,''    '',I1)'
                    IDEFLT(1)=1
                    IF(IOTYPE.EQ.3) IDATA(1)=1
                    CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)
                    IF(IOTYPE.NE.3)THEN
                      IF(IDATA(1).EQ.1.OR.IDATA(1).EQ.2)THEN !pressure BC
                        CDATA(1)='pressure (cm H2O)'
                        IF(IDATA(1).EQ.1) THEN !inlet pressure boundary
                          RDEFLT(1)=10.0d0 !cm H2O (Huang:2001)
C... also check that feed nodes are around correct way
C... for inlet feed node should be NPNE(1)                        
C                        np1=NPNE(1,nb,ne)
C                        np2=NPNE(2,nb,ne)
                          CAP_INLET=ne
                        ELSE IF(IDATA(1).EQ.2) THEN !outlet pressure boundary
                          RDEFLT(1)=2.0d0
C... for outlet feed node should be NPNE(2)
                          CAP_OUTLET=ne
                        ENDIF
                      ELSE !flow BC
                        CDATA(1)='flow (mm**3/s)'
                        RDEFLT(1)=5.69d-4 !mm**3/s
                        CAP_INLET=ne
                      ENDIF
                    ENDIF
                    WRITE(CHAR2,'(D12.4)') RDEFLT(1)
                    CALL STRING_TRIM(CDATA(1),IBEG,IEND)
                    WRITE(CHAR3,'(A)') CDATA(1)(IBEG:IEND)
                    FORMAT='(/$,'' Enter '//CHAR3(IBEG:IEND)//
     '                ' value for '//' ne '//CHAR(1:5)//' ['//CHAR2//']:
     '                '', D12.4)'
                    IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &                -RMAX,RMAX,INFO,ERROR,*9999)
                    IF(IOTYPE.NE.3) bc=RDATA(1)
                    IF(CDATA(1)(IBEG:IEND).EQ.'pressure (cm H2O)') THEN
                      IF(NENP(NPNE(1,nb,ne),0,nr).EQ.1) THEN
                        np=NPNE(1,nb,ne)
                      ELSE IF(NENP(NPNE(2,nb,ne),0,nr).EQ.1) THEN
                        np=NPNE(2,nb,ne)
                      ENDIF
                      ny=NYNP(1,1,1,np,0,1,nr)
                      YP(ny,1)=bc*98.06d0 !cm H2O -> Pa
                      YP(ny,3)=bc*98.06d0
                      FIX(ny,1)=.TRUE.
                    ELSE IF(CDATA(1)(IBEG:IEND).EQ.'flow (mm**3/s)')THEN
                      ny=NYNE(1,1,0,1,ne) !inlet flow
                      YP(ny,1)=bc !only need to specify inlet flow
                      YP(ny,3)=bc
                      FIX(ny,1)=.TRUE.
                    ENDIF
                  ENDIF !.NOT.MANUAL
                ENDIF !node not connected
              ENDDO !noelem
              IF(MANUAL) THEN
                WRITE(OP_STRING,'('' Enter '',I3,'' boundary '
     '            //'conditions'')') NUMBC
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ne=1
                WRITE(CHAR,'(I5)') ne
                DO i=1,NUMBC
                  FORMAT='(/$,'' Enter the element # ['//CHAR//']:'',
     '              I5)'
                  IDEFLT(1)=ne
                  IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '              RMIN,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) ne=IDATA(1)
                  WRITE(CHAR,'(I5)') ne !specify boundary condition
                  FORMAT='('' Enter the type of boundary condition '
     '              //'for ne'//CHAR(1:5)//' [1]:'''//
     '              '/''   (1) Inlet pressure (cm H2O)'''//
     '               '/''   (2) Outlet pressure (cm H2O)'''//
     '              '/''   (3) Flow (mm**3/s)'''//'/$,''    '',I1)'
                  IDEFLT(1)=1
                  IF(IOTYPE.EQ.3) IDATA(1)=1
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     '              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '              IDATA,IDEFLT,1,3,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '              RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3)THEN
                    IF(IDATA(1).EQ.1.OR.IDATA(1).EQ.2)THEN !pressure BC
                      CDATA(1)='pressure (cm H2O)'
                      IF(IDATA(1).EQ.1) THEN
                        RDEFLT(1)=10.0d0 !cm H2O (Huang:2001)
                        CAP_INLET=ne
                      ELSEIF(IDATA(1).EQ.2) THEN
                        RDEFLT(1)=2.0d0
                        CAP_OUTLET=ne
                      ENDIF
                    ELSE !flow BC
                      CDATA(1)='flow (mm**3/s)'
                      RDEFLT(1)=5.69d-4 !mm**3/s
                    ENDIF
                  ENDIF
                  WRITE(CHAR2,'(D12.4)') RDEFLT(1)
                  CALL STRING_TRIM(CDATA(1),IBEG,IEND)
                  WRITE(CHAR3,'(A)') CDATA(1)(IBEG:IEND)
                  FORMAT='(/$,'' Enter '//CHAR3(IBEG:IEND)//
     '              ' value for '//' ne '//CHAR(1:5)//' ['//CHAR2//']:
     '              '', D12.4)'
                  IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     '              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     '              RMIN,RMAX,INFO,ERROR,*9999)
                  IF(IOTYPE.NE.3) bc=RDATA(1)
                  IF(CDATA(1)(IBEG:IEND).EQ.'pressure (cm H2O)') THEN
                    IF(NENP(NPNE(1,nb,ne),0,nr).EQ.1) THEN
                      np=NPNE(1,nb,ne)
                    ELSE IF(NENP(NPNE(2,nb,ne),0,nr).EQ.1) THEN
                      np=NPNE(2,nb,ne)
                    ELSE
                      np=NPNE(1,nb,ne) !if node is connected both ends
                    ENDIF
                    ny=NYNP(1,1,1,np,0,1,nr)
                    YP(ny,1)=bc*98.06d0 !cm H2O -> Pa
                    YP(ny,3)=bc*98.06d0
                    FIX(ny,1)=.TRUE.
                  ELSE IF(CDATA(1)(IBEG:IEND).EQ.'flow (mm**3/s)') THEN
                    ny=NYNE(1,1,0,1,ne) ! flow
                    YP(ny,1)=bc !only need to specify inlet flow
                    YP(ny,3)=bc
                    FIX(ny,1)=.TRUE.
                  ENDIF
                ENDDO !NUMBC
              ENDIF !MANUAL
              CALL ASSERT(NUMBC.GT.0,
     &          'No boundary conditions set, check IPINIT file',ERROR,
     &          *9999) !checks boundary conditions are set
            ENDIF !ITYP3.EQ.3 OR ITYP3.EQ.6
C...        Initial conditions for gas exchange AJS 9/2007
            IF(ITYP3(nr,nx).EQ.3.AND.ITYP7(nr,nx).GE.2.OR.
     &        ITYP3(nr,nx).EQ.6.OR.ITYP3(nr,nx).EQ.7)THEN
	      FORMAT='(/$,'' Enter the initial blood oxygen partial'
     '            //' pressure (mmHg) [40.0]: '',D12.4)'
              RDEFLT(1)=40.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=40.0d0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,200.d0,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) INITIAL_POB=RDATA(1)
              FORMAT='(/$,'' Enter the initial blood carbon dioxide'
     '            //' partial pressure (mmHg) [46.0]: '',D12.4)'
              RDEFLT(1)=46.0d0
              IF(IOTYPE.EQ.3) RDATA(1)=46.0d0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,200.d0,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) INITIAL_PCB=RDATA(1)
	      FORMAT='(/$,'' Enter the initial alveolar oxygen'
     '            //' partial pressure (mmHg) [97.5]: '',D12.4)'
              RDEFLT(1)=97.5d0
              IF(IOTYPE.EQ.3) RDATA(1)=97.5d0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,200.d0,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) INITIAL_POA=RDATA(1)
              FORMAT='(/$,'' Enter the initial alveolar carbon dioxide'
     '            //' partial pressure (mmHg) [37.5]: '',D12.4)'
              RDEFLT(1)=37.5d0
              IF(IOTYPE.EQ.3) RDATA(1)=37.5d0
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.d0,200.d0,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) INITIAL_PCA=RDATA(1)
	    ENDIF !ITYP7 or ITYP3

          ELSE !IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear

C*** KAT 28Jul97:  The following code does not work as nc is not set
C***               correctly.  Also loop is never 5 so several IFs
C***               should not be here.

            CALL ASSERT(.FALSE.,'>>Non-linear not implemented',ERROR,
     '        *9999)

CC***        Nonlinear problems
C            nv=1 !temporary
CC***        Initialize YP and FIX arrays
C            DO nonode=1,NPNODE(0,nr)
C              np=NPNODE(nonode,nr)
C              DO nhx=1,NHP(np)
C                nh=nh_loc(nhx,nx)
C                DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
C                  DO nrc=0,2
C                    ny=NYNP(nk,nv,nh,np,nrc,nc,nr)
C                    DO loop=1,LOOPT
C                      YP(ny,loop) = 0.0d0
C                      FIX(ny,loop)=.FALSE.
C                    ENDDO
C                  ENDDO !nrc
C                ENDDO
C              ENDDO
C            ENDDO
C            DO loop=2,LOOPT
C              IF(IOTYPE.EQ.3) THEN
C                CALL YPZP(loop,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
C     '            nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C              ELSE IF(iotype.ne.3) THEN
C                DO nonode=1,NPNODE(0,nr)
C                  np=NPNODE(nonode,nr)
C                  DO nhx=1,NHP(np)
C                    nh=nh_loc(nhx,nx)
C                    DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
C                      ZP(nk,nv,nh,np,nc)=0.0d0
C                    ENDDO
C                  ENDDO
C                ENDDO
C                DO noelem=1,NEELEM(0,nr)
C                  ne=NEELEM(noelem,nr)
C                  DO nhx=1,NHE(ne)
C                    nh=nh_loc(nhx,nx)
C                    DO na=1,NAT(NBH(nh,nc,ne))
C                      ZA(na,nh,nc,ne)=0.0d0
C                    ENDDO
C                  ENDDO
C                ENDDO
C              ENDIF
C              IF(loop.EQ.2) THEN
C                FORMAT='('' Flux boundary conditions:'')'
C                CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '            ERROR,*9999)
Cc               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
C              ELSE IF(loop.EQ.3) THEN
C                FORMAT='(/'' Dependent variable boundary conditions:'')'
C                CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '            ERROR,*9999)
Cc               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
C              ELSE IF(loop.EQ.4) THEN
C                FORMAT='(/'' Initial conditions:'')'
C                CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '            ERROR,*9999)
Cc               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
C              ELSE IF(loop.EQ.5) THEN
C                FORMAT='(/'' Radiation/convection boundary '
C     '            //'conditions:'')'
C                CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '            ERROR,*9999)
Cc               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
C                FORMAT='($,'' Enter constants a,b,u0 '
C     '            //'in du/dn=a(u^4-u0^4)+b(u-u0): '',3E12.4)'
C                IF(IOTYPE.EQ.3) THEN
C                  RDATA(1)=DBLE(RADIATION_A)
C                  RDATA(2)=DBLE(RADIATION_B)
C                  RDATA(3)=DBLE(RADIATION_U0)
C                ENDIF
C                CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,
C     '            FILEIP,FORMAT,3,ADATA,ADEFLT,CDATA,CDEFLT,
C     '            ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
C     '            RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
C                IF(iotype.ne.3) THEN
C                  RADIATION_A =SNGL(RDATA(1))
C                  RADIATION_B =SNGL(RDATA(2))
C                  RADIATION_U0=SNGL(RDATA(3))
C                ENDIF
C              ENDIF
C
C              DO nhx=1,NH_LOC(0,nx)
C                nh=nh_loc(nhx,nx)
C                WRITE(CHAR1,'(I1)') nhx
C                FORMAT='(/'' Dependent variable/equation number '
C     '            //CHAR1(1:1)//' : '')'
C                CALL GINOUT(IOTYPE,0,IVDU,IFILE,1,0,NOQUES,FILEIP,
C     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
C     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
C     '            ERROR,*9999)
Cc               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
C                nb=NBH(nh,1,NEELEM(1,nr))
C                IF(NNT(nb).GT.0) THEN
C                  IF(IOTYPE.EQ.3) THEN
C                    N1NODE=1
C !                 ny=0
C                    IF(loop.EQ.4) nonode=0 !to initialize node loop
C                                           !for writing i.c.s
C                  ENDIF
C 7100             FORMAT='($,'' Enter node number [EXIT]: '',I3)'
C                  IF(IOTYPE.EQ.3) THEN
C                    IF(loop.LT.4) THEN
C                      DO nonode=N1NODE,NPNODE(0,nr)
C                        np=NPNODE(nonode,nr)
C                        N2NODE=nonode
C !                     DO n1h=1,NHP(np)
C !                       DO nk=1,MAX(NKH(nc,n1h,np)-KTYP93(nc,nr),1)
C !                         ny=ny+1
C !                         IF((n1h.eq.nh).AND.FIX(ny,loop)) GOTO 6200
C !                        ENDDO
C !                     ENDDO
C                        DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
C                          IF(FIX(NYNP(nk,nv,nh,np,0,nc,nr),loop))
C     '                      GOTO 7200
C                        ENDDO
C                      ENDDO
C                      N2NODE=0
C 7200                 IF(N2NODE.EQ.0) THEN
C                        IDATA(1)=0
C                      ELSE
C                        IDATA(1)=np
C                        N1NODE=N2NODE+1
C                      ENDIF
C                    ELSE IF(loop.EQ.4) THEN !initial conditions
C                      nonode=nonode+1
C                      IF(nonode.LE.NPNODE(0,nr)) THEN
C                        np=NPNODE(nonode,nr)
C                        IDATA(1)=np
C                      ELSE
C                        IDATA(1)=0 !to terminate node loop
C                      ENDIF
C                    ELSE IF(loop.EQ.5) THEN !radiation b.c.s
C                      nonode=nonode+1
C                      IF(nonode.LE.NPNODE(0,nr)) THEN
C                        np=NPNODE(nonode,nr)
C                        IDATA(1)=np
C                      ELSE
C                        IDATA(1)=0 !to terminate node loop
C                      ENDIF
C                    ENDIF !loop
C                  ENDIF !iotype=3
C                  CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IZERO,1,NPT(nr), LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
C     '              RMAX,INFO,ERROR,*9999)
C                  IF(IDATA(1).ne.0) THEN
C                    np=IDATA(1)
C                    DO nk=1,MAX(NKH(nh,np,nc)-KTYP93(nc,nr),1)
C                      IF(loop.EQ.2) THEN
C                        IF(nk.EQ.1) THEN
C                          FORMAT='($,'' Do you want to prescribe the'//
C     '                      ' flux component [Y]? '',A)'
C                          ADEFLT(1)='Y'
C                        ELSE IF(nk.GT.1) THEN
C                          WRITE(CHAR1,'(I1)') nk
C                          FORMAT='($,'' Do you want to prescribe moment'
C     '                      //' number '//CHAR1(1:1)//' [N]? '',A)'
C                          ADEFLT(1)='N'
C                        ENDIF
C                      ELSE IF(loop.EQ.3) THEN
C                        IF(nk.EQ.1) THEN
C                          FORMAT='($,'' Do you want to prescribe the'//
C     '                      ' dependent variable [Y]? '',A)'
C                          ADEFLT(1)='Y'
C                        ELSE IF(nk.GT.1) THEN
C                          WRITE(CHAR1,'(I1)') nk
C                          FORMAT='($,'' Do you want to prescribe '
C     '                      //'derivative '
C     '                      //'number '//CHAR1(1:1)//' [N]? '',A)'
C                          ADEFLT(1)='N'
C                        ENDIF
C                      ENDIF
C                      IF(loop.EQ.2.OR.loop.EQ.3) THEN !flux or dep.var.
C                        IF(IOTYPE.EQ.3) THEN
C                          IF(FIX(NYNP(nk,nv,nh,np,0,nc,nr),loop)) THEN
C                            ADATA(1)='Y'
C                          ELSE
C                            ADATA(1)='N'
C                          ENDIF
C                        ENDIF
C                        CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,
C     '                    NOQUES,FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,
C     '                    CDEFLT,
C     '                    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
C     '                    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C                        IF(iotype.ne.3) THEN
C                          IF(ADATA(1).EQ.'Y') THEN
C                            FIX(NYNP(nk,nv,nh,np,0,nc,nr),loop)=.TRUE.
C                          ELSE IF(ADATA(1).EQ.'N') THEN
C                            FIX(NYNP(nk,nv,nh,np,0,nc,nr),loop)=.FALSE.
C                          ENDIF
C                        ENDIF
C                        IF(ADATA(1).EQ.'Y') THEN
C                          FORMAT='($,'' The increment is [0.0]: '','
C     '                      //'G12.5)'
C                          IF(IOTYPE.EQ.3) RDATA(1)=ZP(nk,nv,nh,np,nc)
C                          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,
C     '                      FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
C     '                      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
C     '                      RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
C                          IF(iotype.ne.3) ZP(nk,nv,nh,np,nc)=RDATA(1)
C                        ENDIF
C                      ELSE IF(loop.EQ.4) THEN !initial conditions
C                        RDEFLT(1)=XP(nk,nc,nh,np)
C                        WRITE(CHAR1,'(E12.5)') RDEFLT(1)
C                        IF(nk.EQ.1) THEN
C                          FORMAT='($,'' The value of the dependent '
C     '                      //'variable is ['//CHAR1(1:12)//']: '','
C     '                      //'G12.5)'
C                        ELSE IF(nk.GT.1) THEN
C                          WRITE(CHAR2,'(I1)') nk
C                          FORMAT='($,'' The value of derivative number '
C     '                      //CHAR2(1:1)//' is ['//CHAR1(1:12)//']: '','
C     '                      //'G12.5)'
C                        ENDIF
C                        IF(IOTYPE.EQ.3) RDATA(1)=ZP(nk,nv,nh,np,nc)
C                        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,
C     '                    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
C     '                    ICHAR,
C     '                    IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,
C     '                    RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C                        IF(iotype.ne.3) ZP(nk,nv,nh,np,nc)=RDATA(1)
C                      ELSE IF(loop.EQ.5) THEN !radiation b.c.s
C                        IF(iotype.ne.3) THEN
C                          FIX(NYNP(nk,nv,nh,np,0,nc,nr),1)=.TRUE.
C                          FIX(NYNP(nk,nv,nh,np,0,nc,nr),2)=.TRUE.
C                        ENDIF
C                      ENDIF !loop
C                    ENDDO !nk
C                    GO TO 7100
C                  ENDIF !idata>0
C                ENDIF !nnt>0
C
C                IF(NAT(nb).GT.0) THEN
C                  IF(IOTYPE.EQ.3) THEN
C                    N1ELEM=1
C                    IF(loop.EQ.4) noelem=0 !to initialize elem loop
C !                                         for writing i.c.s
C                  ENDIF
C 7300             FORMAT='($,'' Enter element number [EXIT]: '',I3)'
C                  IF(IOTYPE.EQ.3) THEN
C                    IF(loop.LT.4) THEN
C                      DO noelem=N1ELEM,NEELEM(0,nr)
C                        ne=NEELEM(noelem,nr)
C                        N2ELEM=noelem
C !                     DO n1h=1,NHE(ne)
C !                       DO na=1,NAT(NBH(n1h,nc,ne))
C !                         ny=ny+1
C !                         IF((n1h.eq.nh).AND.FIX(ny,loop)) GOTO 7400
C !                       ENDDO
C !                     ENDDO
C                        DO na=1,NAT(NBH(nh,nc,ne))
C                          IF(FIX(NYNE(na,nh,0,nc,ne),loop)) GOTO 7400
C                        ENDDO !na
C                      ENDDO
C                      N2ELEM=0
C 7400                 IF(N2ELEM.EQ.0) THEN
C                        IDATA(1)=0
C                      ELSE
C                        IDATA(1)=ne
C                        N1ELEM=N2ELEM+1
C                      ENDIF
C                    ELSE IF(loop.EQ.4) THEN !initial conditions
C                      noelem=noelem+1
C                      IF(noelem.LE.NEELEM(0,nr)) THEN
C                        ne=NEELEM(noelem,nr)
C                        IDATA(1)=ne
C                      ELSE
C                        IDATA(1)=0 !to terminate element loop
C                      ENDIF
C                    ENDIF
C                  ENDIF
C                  CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,
C     '              FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
C     '              IZERO,1,NET(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
C     '              RMAX,INFO,ERROR,*9999)
C                  IF(IDATA(1).ne.0) THEN
C                    ne=IDATA(1)
C                    DO na=1,NAT(NBH(nh,nc,ne))
C                      WRITE(CHAR1,'(I1)') na
C                      IF(loop.EQ.2.OR.loop.EQ.3) THEN
C                        FORMAT='($,'' Do you want to prescribe '
C     '                    //'auxiliary variable/rhs number '//CHAR1(1:1)
C     '                    //' [N]? '',A)'
C                        IF(IOTYPE.EQ.3) THEN
C                          IF(FIX(NYNE(na,nh,0,nc,ne),loop)) THEN
C                            ADATA(1)='Y'
C                          ELSE
C                            ADATA(1)='N'
C                          ENDIF
C                        ENDIF
C                        CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,
C     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,
C     '                    IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,
C     '                    RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C                        IF(iotype.ne.3) THEN
C                          IF(ADATA(1).EQ.'Y') THEN
C                            FIX(NYNE(na,nh,0,nc,ne),loop)=.TRUE.
C                            FORMAT='($,'' The increment is [0.0]: '','
C     '                        //'G12.5)'
C                            IF(IOTYPE.EQ.3) RDATA(1)=ZA(na,nh,nc,ne)
C                            CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,
C     '                        NOQUES,FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,
C     '                        CDEFLT,
C     '                        ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
C     '                        RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
C                            IF(iotype.ne.3) ZA(na,nh,nc,ne)=RDATA(1)
C                          ELSE IF(ADATA(1).EQ.'N') THEN
C                            FIX(NYNE(na,nh,0,nc,ne),loop)=.FALSE.
C                          ENDIF
C                        ENDIF
C                      ELSE IF(loop.EQ.4) THEN
C                        FORMAT='($,'' The prescribed value of auxiliary'
C     '                    //' variable number '//CHAR1(1:1)
C     '                    //' is [0.0]: '',G12.5)'
C                        IF(IOTYPE.EQ.3) RDATA(1)=ZA(na,nh,nc,ne)
C                        CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,
C     '                    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
C     '                    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
C     '                    RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
C                        IF(iotype.ne.3) ZA(na,nh,nc,ne)=RDATA(1)
C                      ENDIF !loop
C                    ENDDO !na
C                    GO TO 7300
C                  ENDIF !idata>0
C                ENDIF !nat>0
C              ENDDO !nh
C              IF(iotype.ne.3) THEN
C                CALL ZPYP(loop,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
C     '            nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
C              ENDIF
C            ENDDO !loop
          ENDIF !linear/nonlinear
        ENDIF !generating initial conditions
      ENDIF !equation type

      IF(KTYPMBC.EQ.2.AND.NPMC.GT.0) THEN
!       CPB 28/10/92 Calculate the effective integrated mixed bc
!       values from the point values specified on the boundary.
        DO ny=1,NYT(2,1,nx) !should use NYNR here
          IF(FIX(ny,3)) THEN
C***        Current dof has a mixed bc so find the corresponding np,nh,
C***        nk
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO nhx=1,NHP(np)
                nh=nh_loc(nhx,nx)
                DO nk=1,MAX(NKH(nh,np,1)-KTYP93(1,nr),1)
                  IF(ny.EQ.NYNP(nk,nv,nh,np,0,nc,nr)) GOTO 7000
                ENDDO
              ENDDO
            ENDDO
 7000       CONTINUE
C***        Have found the np,nh,nk appropriate to the ny so now find
C***        all elements that contain that np.
            DO noelem=1,NEELEM(0,nr) !Loop over elements
              ne=NEELEM(noelem,nr)
              nb=NBH(NH_LOC(1,nx),1,ne)
              DO nn=1,NNT(nb) !Loop over dep. variable element nodes
                IF(np.EQ.NPNE(nn,nb,ne))THEN !Node is local node nn in
! element ne
C This array needs NPNE(np,0,nr)=NPNE(np,0,nr)+1
C to be passed     NENP(ne,0,nr)=NENP(ne,0,nr)+1
C                  NPNE(np,NPNE(np,0,nr),nr)=ne
C                  NENP(ne,NENP(ne,0,nr),nr)=np
                ENDIF
              ENDDO
            ENDDO
C            IF(NENP(np,0).GT.1) THEN !BC node is in more than one elem
C            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF(FILEIP.AND..NOT.ALL_REGIONS) CALL CLOSEF(IFILE,ERROR,*9999)

      CALL EXITS('IPINI3')
      RETURN
 9999 CALL ERRORS('IPINI3',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPINI3')
      RETURN 1
      END


