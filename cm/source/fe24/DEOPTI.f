      SUBROUTINE DEOPTI(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '  IDO,INP,ISIZE_MFI,ISIZE_PHI,LD,LDR,NBH,NBJ,NDIPOLES,
     '  NEELEM,NENP,NHP,
     '  NKB,NKHE,NKH,NKJ,NLNO,NMNO,NNB,NONL,NONM,NONY,NP1OPT,
     '  NP2OPT,NP3OPT,NPL,NPLIST4,
     '  NPNE,NPNODE,NPNY,NRLIST,NVHE,NVHP,
     '  NVJP,NXI,NXLIST,NYNE,NYNO,NYNP,
     '  NYNR,NYNY,PAOPTY,
     '  CONY,CYNO,CYNY,DIPOLE_CEN,DIPOLE_DIR,DL,PAOPTI,
     '  PBOPTI,PMAX,PMIN,XP,YP,STRING,FIX,ERROR,*)

C#### Subroutine: DEOPTI
C###  Description:
C###    DEOPTI defines optimisation parameters.
C###    Note: This command deals with all currently defined regions
C###    since often the optimisation involves parameters from multiple
C###    regions.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER
     '  DIPOLE_CEN_NTIME(NDIPOLEM,NRM),DIPOLE_DIR_NTIME(NDIPOLEM,NRM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),ISIZE_MFI(3,NSSM),
     '  ISIZE_PHI(2),LD(NDM),LDR(0:NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDIPOLES(NRM,NXM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NHP(NPM,0:NRM,NXM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NLNO(NOPM,NXM),NMNO(1:2,0:NOPM,NXM),NNB(4,4,4,NBFM),
     '  NONL(NLM,NXM),NONM(NMM,NPM,NXM),NONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  NP1OPT(NOPM),NP2OPT(NOPM),NP3OPT(NOPM),NPL(5,0:3,NLM),
     '  NPLIST4(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM),NYNY(0:NYYM,NYM,NRM,NXM),
     '  PAOPTY(NOPM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM,NXM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM),CYNY(0:NYYM,NYM,NRM,NXM),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DL(3,NLM),PAOPTI(*),PBOPTI(*),PMAX(*),PMIN(*),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)

!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,
     '  IPFILE,N3CO,no,no_nrlist,no_nynr,nr,nx_opt,nx_sol,nxc,
     '  NX_ACTION,ny
      REAL*8 TIME
      CHARACTER FILE*(MXCH),STATUS*3
      LOGICAL ALL_REGIONS,CALCU,CBBREV,FILIO,GENER,MOUSE

!     Function
      REAL*8 RFROMC

      CALL ENTERS('DEOPTI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define optimise;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Defines optimisation parameters. Parameters are read from or
C###    written to a file FILENAME.ipopti.
C###  Parameter:      <time (#)[0.0]>
C###    Specify the time to use for setting the initial values
C###    for static optimisation problems
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to define the optimisation for.
C###    The all value specifies all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number to define
C###    the optimisation for.
C###  Parameter:      <(lock/nolock)[lock]>
C###    The lock option locks the optimisation problem, this prevents
C###    other problems being given the same nx (problem type number)
C###    and overwriting the optimisation problem's arrays.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<time (#)[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #1>[1]'
        OP_STRING(5)=BLANK(1:15)//'<(lock/nolock)[lock]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEOPTI',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 4-Feb-1994

        CALL PARSE_QUALIFIERS('DLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C GBS 6/10/94 Locking for NX_LOC
        IF(CBBREV(CO,'NOLOCK',2,noco+1,NTCO,N3CO)) THEN
          NX_ACTION=NX_ALLOCATE
        ELSE
          NX_ACTION=NX_ALLOCATE_AND_LOCK
        ENDIF

        IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
          TIME=RFROMC(CO(N3CO+1))
        ELSE
          TIME=0.D0
        ENDIF

C CPB 8/6/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
        IF(nx_opt.EQ.0) THEN
          CALL NX_LOC(NX_ACTION,nxc,nx_opt,NX_OPTI,ERROR,*9999)
        ENDIF
        CALL NX_LOC(NX_INQUIRE,nxc,nx_sol,NX_SOLVE,ERROR,*9999)

        IF(FILIO) THEN
          CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'opti',STATUS,
     '      ERR,ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
c cpb 20/3/95 Need to prompt for regions here
CC AJPs Don't drop the arrays down by nx_opt
C            CALL IPOPTI(LD,LDR,NEELEM,NHP(1,nr,nx_opt),
C     '        NKH(1,1,1,nr),
C     '        NLNO(1,nx_opt),NMNO(1,0,nx_opt),
C     '        NONL(1,nx_opt),NONM(1,1,nx_opt),NONY(0,1,1,nr,nx_opt),
C     '        NP1OPT,NP2OPT,NP3OPT,
C     '        NPL,NPNODE,NPNY(0,1,0,nx_opt),nr,NVHP(1,1,1,nr),NVJP,
C     '        nx_opt,nx_sol,NYNE,NYNO(0,1,1,nr,nx_opt),NYNP,
C     '        NYNR(0,0,1,nr,nx_opt),PAOPTY,CONY(0,1,1,nr,nx_opt),
C     '        CYNO(0,1,1,nr,nx_opt),DL,PAOPTI,PBOPTI,PMAX,PMIN,
C     '        XP,YP(1,1,nx_opt),FIX(1,1,nx_opt),ERROR,*9999)

            CALL IPOPTI(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '        IDO,INP,ISIZE_MFI,ISIZE_PHI,LD,LDR,NBH,NBJ,
     '        NDIPOLES(1,nxc),NEELEM,NENP,NHP,NKB,NKHE,
     '        NKH(1,1,1,nr),NKJ,NLNO,NMNO,NNB,NONL,NONM,NONY,NP1OPT,
     '        NP2OPT,NP3OPT,NPL,NPLIST4,NPNE,NPNODE,NPNY,nr,NVHE,NVHP,
     '        NVJP,nx_opt,nx_sol,NXI,NYNE,NYNO,NYNP,
     '        NYNR,NYNY,PAOPTY,CONY,CYNO,CYNY,DIPOLE_CEN,DIPOLE_DIR,DL,
     '        PAOPTI,PBOPTI,PMAX,PMIN,TIME,XP,YP,FIX,ERROR,*9999)
CC AJPe
cPJH 9Mar96 Copy mapping info from soln arrays to optim arrays
            IF(KTYP26.EQ.1) THEN
              IF(ITYP4(nr,nx_sol).NE.4) THEN
                !material parameter optimisation, but not a grid problem
                NOT(1,1,nr,nx_opt)=NOT(1,1,nr,nx_sol)
                DO no=1,NOT(1,1,nr,nx_opt)
                  NYNO(0,no,1,nr,nx_opt)=NYNO(0,no,1,nr,nx_sol)
                  NYNO(1,no,1,nr,nx_opt)=NYNO(1,no,1,nr,nx_sol)
                  CYNO(0,no,1,nr,nx_opt)=CYNO(0,no,1,nr,nx_sol)
                  CYNO(1,no,1,nr,nx_opt)=CYNO(1,no,1,nr,nx_sol)
                ENDDO !no
                DO no_nynr=1,NYNR(0,1,1,nr,nx_sol)
                  ny=NYNR(no_nynr,1,1,nr,nx_sol)
                  NONY(0,ny,1,nr,nx_opt)=NONY(0,ny,1,nr,nx_sol)
                  NONY(1,ny,1,nr,nx_opt)=NONY(1,ny,1,nr,nx_sol)
                  CONY(0,ny,1,nr,nx_opt)=CONY(0,ny,1,nr,nx_sol)
                  CONY(1,ny,1,nr,nx_opt)=CONY(1,ny,1,nr,nx_sol)
                  DO i=0,6
                    NPNY(i,ny,1,nx_opt)=NPNY(i,ny,1,nx_sol)
                  ENDDO !i
                ENDDO !no_nynr
              ENDIF !not grid
            ENDIF !ktyp26

          ENDDO !no_nrlist
          CALL CLOSEF(IFILE,ERROR,*9999)

C LKC 30-NOV-1999 Adding a warning that define inverse
C should be called for activation problems
          IF(.NOT.CALL_INVE.AND.KTYP27.EQ.12) THEN
            OP_STRING(1)='WARNING: Should define inverse first'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

        ENDIF
      ENDIF

      CALL EXITS('DEOPTI')
      RETURN
 9999 CALL ERRORS('DEOPTI',ERROR)
      IF(nx_opt.GT.0) THEN
        CALL NX_LOC(NX_FREE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
      ENDIF
      CALL EXITS('DEOPTI')
      RETURN 1
      END


