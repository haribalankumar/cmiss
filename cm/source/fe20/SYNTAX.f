      SUBROUTINE SYNTAX(ISEG,CSEG,END,ERROR,*)

C#### Subroutine: SYNTAX
C###  Description:
C###    Determines whether CO is a valid command. If CO is known the
C###    command qualifiers are checked. When the qualifiers are
C###    understood the command is executed.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi02.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'diag00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks000.cmn'
      INCLUDE 'gtstr00.cmn'
      INCLUDE 'learn00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'ofst00.cmn'
C$    INCLUDE 'openmp.inc'
      INCLUDE 'opti00.cmn'
C$    INCLUDE 'mp00.cmn'
      INCLUDE 'time02.cmn'

      INCLUDE 'parameters.inc'

!     Parameter List
      INTEGER ISEG(*)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL END
!     Work array pointers
      INTEGER*4 INTWORK_PTR,REALWORK_PTR
      SAVE INTWORK_PTR,REALWORK_PTR

!     Local Variables
      INTEGER CODE,i,IBEG,IBEG1,IEND,IEND1,
     '  IUNIT,iw,N3CO,nb,nc,NMAX,NUM_THREADS,noiw,nosg,nr,
     '  nrc,NTCH1,nx,SQUID_CONFIG
      REAL CPUTIME,REALTIME
      REAL*8 BASELINE,TIME
      CHARACTER C_IW_NJT*1,C1*(MXCH),COMFILENAME*(MXCH),FILE*(MXCH),
     '  OPTION_1(40)*80,STRING*(MXCH)
      LOGICAL ABBREV,CBBREV,COMMAND,ECHO,LIST_1

!     Function
      INTEGER IFROMC
      REAL*8 RFROMC

      CALL ENTERS('SYNTAX',*9999)

      IF(CMGUI_LINK) THEN
C       Update the front end
        CALL WH_OUTPUT_F_UPDATE(CMGUI_DATA_O,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='Could not update data output'
          GOTO 9999
        ENDIF
      ENDIF

      IF(FIRST_SYNTAX.OR.REALLOCATE_SYNTAX) THEN
        IF(FIRST_SYNTAX) THEN
          MEM_INIT=.FALSE.
          FIRST_FEM=.TRUE.
          IO1=1
          IO2=2
          IO3=IOIP
          IO4=IOOP
          IVDU=IOIP
          IFILE=10
          IF(.NOT.USEPARAMFILE) THEN
            NAM=NAMX
            NBM=NBMX
            NCM=NCMX
            NDM=NDMX
            NEM=NEMX
            NEGM=NEGMX
            NEIM=NEIMX
            NELM=NELMX
            NEPM=NEPMX
            NNEPM=NNEPMX
            NNEPM=NNEPMX
            NE_R_M=NE_R_MX
            NFM=NFMX
            NF_R_M=NF_R_MX
            NFVM=NFVMX
            NFVCM=NFVCMX
            NGM=NGMX
            NGFM=IDNINT(NGM**(2.0d0/3))
            NGRSEGM=NGRSEGMX
            NHM=NHMX
            NIM=NIMX
            NISC_GDM=NISC_GDMX
            NISR_GDM=NISR_GDMX
            NISC_GKM=NISC_GKMX
            NISR_GKM=NISR_GKMX
            NISC_GKKM=NISC_GKKMX
            NISR_GKKM=NISR_GKKMX
            NISC_GMM=NISC_GMMX
            NISR_GMM=NISR_GMMX
            NISC_GMMM=NISC_GMMMX
            NISR_GMMM=NISR_GMMMX
            NISC_GQM=NISC_GQMX
            NISR_GQM=NISR_GQMX
            NJM=NJMX
            NKM=NKMX
            NLM=NLMX
            NL_R_M=NL_R_MX
            NMAQM=NMAQMX
            NMM=NMMX
            NNM=NNMX
            NOM=NOMX
            NPM=NPMX
            NP_R_M=NP_R_MX
            NQIM=NQIMX
            NQISVM=NQISVMX
            NQRM=NQRMX
            NQRSVM=NQRSVMX
            NQM=NQMX
            NRM=NRMX
            NSM=NSMX
            NSSM=NSSMX
            NSFM=NSFMX
            NTM=NTMX
            NTSM=NTSMX
            NUM=NUMX
            NVCBM=NVCBMX
            NVCM=NVCMX
            NVM=NVMX
            NWM=NWMX
            NXM=NXMX
            NYM=NYMX
            NYQM=NYQMX
            NY_R_M=NY_R_MX
C PJH 14/6/98 NZM=NZMX
            NZ_FNY_M  =NZ_FNY_MX
            NZ_GD_M   =NZ_GD_MX
            NZ_GK_M   =NZ_GK_MX
            NZ_GKK_M  =NZ_GKK_MX
            NZ_GM_M   =NZ_GM_MX
            NZ_GMM_M  =NZ_GMM_MX
            NZ_GQ_M   =NZ_GQ_MX
            NBFM=NBFMX
            NCOM=NCOMX
            NDEM=NDEMX
            NDIPOLEM=NDIPOLEMX
            NDIPTIMM=NDIPTIMMX
            NIQM=NIQMX
            NIFEXTM=NIFEXTMX
            NIYM=NIYMX
            NIYFIXM=NIYFIXMX
            NIYGM=NIYGMX
            NIYGFM=NIYGFMX
            NLCM=NLCMX
            NOPM=NOPMX
            NORM=NORMX
            NOYM=NOYMX
            NPDM=NPDMX
            NQEM=NQEMX
            NQGM=NQGMX
            NQSCM=NQSCMX
            NRCM=NRCMX
            NREM=NREMX
            NTIMEPOINTSM=NTIMEPOINTSMX
            NTIMEVARSM=NTIMEVARSMX
            NYOM=NYOMX
            NLISTM=NLISTMX
            NOOPM=NOOPMX
            NZ_MINOSM=NZ_MINOS
            NIMAGEM  =NIMAGEMX
            NYROWM=NYROWMX
            NY_TRANSFER_M=NY_TRANSFER_MX
            NYYM=NYYMX
            USE_BEM      =USE_BEMX
            USE_DATA     =USE_DATAX
            USE_DIPOLE   =USE_DIPOLEX
            USE_GAUSS_PT_MATERIALS =USE_GAUSS_PT_MATERIALSX
            USE_GRAPHICS =USE_GRAPHICSX
            USE_GRID     =USE_GRIDX
            USE_LUNG     =USE_LUNGX
            USE_MAGNETIC =USE_MAGNETICX
            USE_MINOS    =USE_MINOSX
            USE_NLSTRIPE =USE_NLSTRIPEX
            USE_NONLIN   =USE_NONLINX
            USE_NPSOL    =USE_NPSOLX
C DBs 25/3/01: F77 does not allow functions in parameter statements
            USE_OPTI     =MAX(USE_NPSOLX,USE_MINOSX)
C DBe        USE_OPTI     =USE_OPTIX
            USE_SPARSE   =USE_SPARSEX
            USE_TIME     =USE_TIMEX
            USE_TRANSFER =USE_TRANSFERX
            USE_VORONOI  =USE_VORONOIX
          ENDIF
          INTWORK_PTR=0
          REALWORK_PTR=0
        ENDIF
        IF(REALLOCATE_SYNTAX) THEN
C***      NAM
          NMAX=0
          DO nb=1,NBT
            IF(NAT(nb).GT.NMAX) NMAX=NAT(nb)
          ENDDO !nb
          NAM=MAX(NAM,NMAX)
C***      NBM
          NBM=MAX(NBM,NBT)
C***      NBFM
          NBFM=MAX(NBFM,NBFT)
C***      NCM
          NMAX=0
          DO nx=1,NXM
            DO nr=1,NRT
              IF(NCT(nr,nx).GT.NMAX) NMAX=NCT(nr,nx)
            ENDDO !nr
          ENDDO !nr
          NCM=MAX(NCM,NMAX)
C***      NCOM
          NCOM=MAX(NCOM,NTCNTR)
C***      NDM
          NDM=MAX(NDM,NDT)
C***      NDEM
C          NDEM=MAX(NDEM,NDEM) !TO DO
C***      NEM
          NEM=MAX(NEM,NET(0))
C***      NE_R_M
C KAT 30Oct97:  This doesn't do what it should.
C          NMAX=0
C          DO nr=1,NRT
C            IF(NET(nr).GT.NMAX) NMAX=NET(nr)
C          ENDDO !nr
C          NE_R_M=MAX(NE_R_M,NMAX)
           DO nr=1,NRT
             NET(nr)=0
           ENDDO !nr
C***      NFM
          NFM=MAX(NFM,NFT)
C***      NF_R_M
C          NF_R_M=MAX(NF_R_M,NF_R_M) !TO DO
C***      NGM
          NMAX=0
          DO nb=1,NBT
            IF(NGT(nb).GT.NMAX) NMAX=NGT(nb)
          ENDDO !nb
          NGM=MAX(NGM,NMAX)
          NGFM=IDNINT(NGM**(2.0d0/3))
          NEGM=1+NEM*(NGM+2)
C***      NHM
          NHM=MAX(NHM,NH_LOC(0,0))
C***      NIM
C          NIM=MAX(NIM,NIM) !TO DO
C***      NIMAGEM
C          NIMAGEM=MAX(NIMAGEM,NIMAGEM) !TO DO
C***      NIQM
C          NIQM=MAX(NIQM,NIQM) !TO DO
C***      NIYM
C          NIYM=MAX(NIYM,NIYM) !TO DO
C***      NIYGM
C***      NIYGFM
C***      NJM
          NJM=MAX(NJM,NJ_LOC(0,0,0))
C***      NKM
          NMAX=0
          DO nb=1,NBT
            IF(NKT(0,nb).GT.NMAX) NMAX=NKT(0,nb)
          ENDDO !nb
          NKM=MAX(NKM,NMAX)
C***      NJM
          NJM=MAX(NJM,NJ_LOC(0,0,0))
C***      NKM
          NMAX=0
          DO nb=1,NBT
            IF(NKT(0,nb).GT.NMAX) NMAX=NKT(0,nb)
          ENDDO !nb
          NKM=MAX(NKM,NMAX)
C***      NLM
          NLM=MAX(NLM,NLT)
C***      NLCM
C          NLCM=MAX(NLCM,NLCM) !TO DO
C***      NL_R_M
C          NL_R_M=MAX(NL_R_M,NL_R_M) !TO DO
C***      NMAQM
C          NMAQM=MAX(NMAQM,NMAQM) !TO DO
C***      NMM
C          NMM=MAX(NMM,NMM) !TO DO
C***      NNM
          NMAX=0
          DO nb=1,NBT
            IF(NNT(nb).GT.NMAX) NMAX=NNT(nb)
          ENDDO !nb
          NNM=MAX(NNM,NMAX)
C***      NOM
          NMAX=0
          DO nx=1,NXM
            DO nr=1,NRT
              DO nrc=1,2
                IF(NOT(nrc,1,nr,nx).GT.NMAX) NMAX=NOT(nrc,1,nr,nx)
              ENDDO !nrc
            ENDDO !nr
          ENDDO !nx
          NOM=MAX(NOM,NMAX)
C***      NOPM
          NOPM=MAX(NOPM,NTOPTI)
C***      NORM
C          NORM=MAX(NORM,NORM) !TO DO
C***      NOYM
C          NOYM=MAX(NOYM,NOYM) !TO DO
C***      NPM
          NPM=MAX(NPM,NPT(0))
C***      NPDM
C          NPDM=MAX(NPDM,NPDM) !TO DO
C***      NP_R_M
          NMAX=0
          DO nr=1,NRT
            IF(NPT(nr).GT.NMAX) NMAX=NPT(nr)
          ENDDO !nr
          NP_R_M=MAX(NP_R_M,NMAX)
C***      NQIM
          NQIM=MAX(NQIM,NQIT)
C***      NQISVM
          NQISVM=MAX(NQISVM,NQISVT)
C***      NQRM
          NQRM=MAX(NQRM,NQRT)
C***      NQRSVM
          NQRSVM=MAX(NQRSVM,NQRSVT)
C***      NQM
          NQM=MAX(NQM,NQT)
C***      NRM
          NRM=MAX(NRM,NRT)
C***      NRCM
C          NRCM=MAX(NRCM,NRCM) !TO DO
C***      NREM
          NREM=MAX(NREM,NT_RES)
C***      NSM
          NMAX=0
          DO nb=1,NBT
            IF(NST(nb).GT.NMAX) NMAX=NST(nb)
          ENDDO !nb
          NSM=MAX(NSM,NMAX)
C***      NSFM
          NMAX=MIN(16,NSM-NSM/2)
          NSFM=MAX(NSFM,NMAX)
C***      NTM
C          NTM=MAX(NTM,NTM) !TO DO
C***      NUM
          NMAX=0
          DO nb=1,NBT
            IF(NUT(nb).GT.NMAX) NMAX=NUT(nb)
          ENDDO !nb
          NUM=MAX(NUM,NMAX)
C***      NVM
C          NVM=MAX(NVM,NVM) !TO DO
C***      NWM
C          NWM=MAX(NWM,NWM) !TO DO
C***      NXM
          NXM=MAX(NXM,NX_LIST(0))
C***      NYM
C          NYM=MAX(NYM,NYM) !TO DO !!!!!!!!
C***      NYOM
C          NYOM=MAX(NYOM,NYOM) !TO DO
C          NYQM=MAX(NYQM,NYQM) !TO DO !!!!!!!!
C***      NY_R_M
C          NY_R_M=MAX(NY_R_M,NY_R_M) !TO DO  !!!!!!!!
C***      NYROWM
          NMAX=0
          DO nx=1,NXM
            DO nr=1,NRT
              DO nc=1,NCT(nr,nx)
                DO nrc=1,2
                  IF(NYT(nrc,nc,nx).GT.NMAX) NMAX=NYT(nrc,nc,nx)
                ENDDO !nrc
              ENDDO !nc
            ENDDO !nr
          ENDDO !nx
          NYROWM=MAX(NYROWM,NMAX)
C***      NY_TRANSFER_M
C          NY_TRANSFER_M=MAX(NY_TRANSFER_M,NY_TRANSFER_M) !TO DO
C***      NZ_FNY_M
C          NZ_FNY_M=MAX(NZ_FNY_M,NZ_FNY_M) !TO DO
C***      NZ_GD_M
          NMAX=0
          DO nx=1,NXM
            IF(NZT(3,nx).GT.NMAX) NMAX=NZT(3,nx)
          ENDDO !nx
          NZ_GD_M=MAX(NZ_GD_M,NMAX)
C***      NZ_GK_M
          NMAX=0
          DO nx=1,NXM
            IF(NZT(1,nx).GT.NMAX) NMAX=NZT(1,nx)
          ENDDO !nx
          NZ_GK_M=MAX(NZ_GK_M,NMAX)
C***      NZ_GKK_M
          NMAX=0
          DO nx=1,NXM
            DO nr=1,NRT
              IF(NZZT(1,nr,nx).GT.NMAX) NMAX=NZZT(1,nr,nx)
            ENDDO !nr
          ENDDO !nx
          NZ_GKK_M=MAX(NZ_GKK_M,NMAX)
C***      NZ_GM_M
          NMAX=0
          DO nx=1,NXM
            IF(NZT(4,nx).GT.NMAX) NMAX=NZT(4,nx)
          ENDDO !nx
          NZ_GM_M=MAX(NZ_GM_M,NMAX)
C***      NZ_MM_M
          NMAX=0
          DO nx=1,NXM
            DO nr=1,NRT
              IF(NZZT(4,nr,nx).GT.NMAX) NMAX=NZZT(4,nr,nx)
            ENDDO !nr
          ENDDO !nx
          NZ_GMM_M=MAX(NZ_GMM_M,NMAX)
C***      NZ_GQ_M
          NMAX=0
          DO nx=1,NXM
            IF(NZT(2,nx).GT.NMAX) NMAX=NZT(2,nx)
          ENDDO !nx
          NZ_GQ_M=MAX(NZ_GQ_M,NMAX)
C***      NISC_GDM
C          NISC_GDM=MAX(NISC_GDM,NISC_GDM) !TO DO
C***      NISR_GDM
C          NISR_GDM=MAX(NISR_GDM,NISR_GDM) !TO DO
C***      NISC_GKM
C          NISC_GKM=MAX(NISC_GKM,NISC_GKM) !TO DO
C***      NISR_GKM
C          NISR_GKM=MAX(NISR_GKM,NISR_GKM) !TO DO
C***      NISC_GKKM
C          NISC_GKKM=MAX(NISC_GKKM,NISC_GKKM) !TO DO
C***      NISR_GKKM
C          NISR_GKKM=MAX(NISR_GKKM,NISR_GKKM) !TO DO
C***      NISC_GM
C          NISC_GMM=MAX(NISC_GMM,NISC_GMM) !TO DO
C***      NISR_GM
C          NISR_GMM=MAX(NISR_GMM,NISR_GMM) !TO DO
C***      NISC_GMM
C          NISC_GMMM=MAX(NISC_GMMM,NISC_GMMM) !TO DO
C***      NISR_GMM
C          NISR_GMMM=MAX(NISR_GMMM,NISR_GMMM) !TO DO
C***      NISC_GQM
C          NISC_GQM=MAX(NISC_GQM,NISC_GQM) !TO DO
C***      NISR_GQM
C          NISR_GQM=MAX(NISR_GQM,NISR_GQM) !TO DO
C***      NZ_MINOSM
C          NZ_MINOSM=MAX(NZ_MINOSM,NZ_MINOSM) !TO DO
C***      NELM
C          NELM=MAX(4,NVM)
C          NELM=MAX(NELM,NELM) !TO DO
C***      NNEPM
          NNEPM=MAX(NEM,NPM)
C***      NLISTM
          NLISTM=MAX(NEM,NPM,NLM,NDM,NQM)
C***      NOOPM
          NOOPM=MAX(NOM,NOPM)

        ENDIF

C       Allocate memory used by the optimisers, and calculate the
C       memory offset addresses.
        CALL OPTIMISER_ALLOC(INTWORK_PTR,REALWORK_PTR,ERROR,*9999)

        IF(FIRST_SYNTAX) FIRST_SYNTAX=.FALSE.

        REALLOCATE_SYNTAX=.FALSE.

      ENDIF


      STRING=' '
      DO noco=1,NTCO-1
        CALL STRING_TRIM(STRING,IBEG,IEND)
        STRING=STRING(1:IEND)//' '//CO(noco)
      ENDDO
      noco=1

      IF(ABBREV(CO(1),'?',1)) THEN
        COMMAND=.FALSE.
        LIST_1=.TRUE.
      ELSE
        COMMAND=.TRUE.
        LIST_1=.FALSE.
      ENDIF

      IF(LIST_1) THEN
        OPTION_1( 1)='Angiography'
        OPTION_1( 2)='Assign'
        OPTION_1( 3)='Beads'
        OPTION_1( 4)='Close'
        OPTION_1( 5)='Constants'
        OPTION_1( 6)='Deassign'
        OPTION_1( 7)='Define'
        OPTION_1( 8)='Display'
        OPTION_1( 9)='Evaluate'
        OPTION_1(10)='FEM'
        OPTION_1(11)='GEN'
        OPTION_1(12)='Help'
        OPTION_1(13)='IMP'
        OPTION_1(14)='Inquire'
        OPTION_1(15)='Iterate'
        OPTION_1(16)='List'
        OPTION_1(17)='Menu'
        OPTION_1(18)='Open'
        OPTION_1(19)='Optimise'
C MPN 26Mar2002: changed >print to >fem gxprint
C        OPTION_1(20)='Print'
        OPTION_1(20)='Quit'
        OPTION_1(21)='Read'
        OPTION_1(22)='Refresh'
        OPTION_1(23)='Set'
        OPTION_1(24)='System'
        OPTION_1(25)='Units'
        OPTION_1(26)='Return'
        NTCH1=26
        CALL LIST_COMMANDS(1,NTCH1,OPTION_1,ERROR,*9999)
      ENDIF

      IF(ABBREV(CO(noco),'ANGIOGRAPHY',2)) THEN
        CALL ANGIO(ISEG,CSEG,STRING,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'ASSIGN',2)) THEN
         CALL ASSIGN(STRING,noco,CO,COQU,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'BEADS',1)) THEN
        CALL BEADS(ISEG,CSEG,STRING,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'CLOSE',2)) THEN
        noco=noco+1
        IF(ABBREV(CO(noco),'?',1)) THEN
          OPTION_1( 1)='Binary'
          OPTION_1( 2)='Files'
          OPTION_1( 3)='Metafile'
          OPTION_1( 4)='Postcript'
          OPTION_1( 5)='Window'
          OPTION_1( 6)='Exit'
          NTCH1=6
          CALL LIST_COMMANDS(1,NTCH1,OPTION_1,ERROR,*9999)
        ENDIF
        IF(ABBREV(CO(noco),'BINARY',1)) THEN
          IF(CO(noco+1).EQ.'?') THEN
            CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: close binary <#[1]>
C###  Description:
C###    Closes the binary file.
            OP_STRING(1)=STRING(1:IEND)//' <#[1]>'
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          ELSE !close file NUMBER
            IF(NTCO.EQ.3) THEN
              IUNIT=IFROMC(CO(3))
            ELSE
              IUNIT=1
            ENDIF
            CALL BINCLOSEFILE(IUNIT,ERROR,*9999)
          ENDIF
        ELSE IF(ABBREV(CO(noco),'FILES',1)) THEN
          IF(CO(noco+1).EQ.'?') THEN
            CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: close files <(#/all)[all]>
C###  Description:
C###    Closes a file(s).
            OP_STRING(1)=STRING(1:IEND)//' <(#/all)[all]>'
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          ELSE IF(NTCO.EQ.3) THEN !close file NUMBER
            IUNIT=IFROMC(CO(3))
            CALL CLOSEF(IUNIT,ERROR,*9999)
          ELSE IF(NTCO.EQ.2) THEN !close file all
            DO i=1,100
              IF(i.ne.ioip.and.i.ne.ioop.and.i.ne.IO4.AND
     '          .i.ne.ioer.and.i.ne.IOTR) THEN
                CALL CLOSEF(i,ERROR,*9999)
              ENDIF
            ENDDO
          ENDIF
        ELSE IF(ABBREV(CO(noco),'METAFILE',1)) THEN
          IF(CO(noco+1).EQ.'?') THEN
            CALL STRING_TRIM(STRING,IBEG,IEND)
            OP_STRING(1)=STRING(1:IEND)
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          ELSE !close metafile
c           CALL CLOSE_PRINT_FILE('METAFILE',ERROR,*9999)
          ENDIF
        ELSE IF(ABBREV(CO(noco),'POSTSCRIPT',1)) THEN
          IF(CO(noco+1).EQ.'?') THEN
            CALL STRING_TRIM(STRING,IBEG,IEND)
            OP_STRING(1)=STRING(1:IEND)
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          ELSE !close postscript
c           CALL CLOSE_PRINT_FILE('POSTSCRIPT',ERROR,*9999)
          ENDIF
        ELSE IF(ABBREV(CO(noco),'WINDOW',1)) THEN
          IF(CO(noco+1).EQ.'?') THEN
            CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: close window <NUMBER>[all]
C###  Description:
C###    Closes the graphics window.
            OP_STRING(1)=STRING(1:IEND)//' <NUMBER>[all]'
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          ELSE
            IF(NTCO.EQ.3) THEN !close window NUMBER
              iw=IFROMC(CO(3))
              CALL DAWK(iw,1,ERROR,*9999)
!new MPN 7-Apr-95: deleting GKS stuff
              CALL CLOSE_WS(iw,ERROR,*9999)
!old
!              CALL GKS_CLWK(iw,ERROR,*9999)
!              IWKS(iw)=0
            ELSE IF(NTCO.EQ.2) THEN !close window all
!new MPN 7-Apr-95: deleting GKS stuff
              CALL CLWS(ERROR,*9999)
!old
!              DO iw=1,99
!                CALL DAWK(iw,1,ERROR,*9999)
!                CALL CLOSE_WS(iw,ERROR,*9999)
!                CALL GKS_CLWK(iw,ERROR,*9999)
!                IWKS(iw)=0
!              ENDDO
            ENDIF
          ENDIF
        ELSE IF(ABBREV(CO(noco),'EXIT',1)) THEN
        ENDIF
        noiw=0
        DO iw=1,99      !to record defined workstations in IWKDEF
          IF(IWKS(iw).GT.0) THEN !iw is defined
            noiw=noiw+1
            IWKDEF(noiw)=iw
          ENDIF
        ENDDO
        IWKDEF(0)=noiw

      ELSE IF(ABBREV(CO(noco),'CONSTANTS',2)) THEN
        noco=noco+1
c       CALL CONSTANTS(CO,noco,STRING,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'DEFINE',2)) THEN
        noco=noco+1
        IF(ABBREV(CO(noco),'?',1)) THEN
          OPTION_1( 1)='Iterate'
          OPTION_1( 2)='Exit'
          NTCH1=2
          CALL LIST_COMMANDS(2,NTCH1,OPTION_1,ERROR,*9999)
        ELSE
          IF(ABBREV(CO(noco),'ITERATE',1)) THEN
            CALL DEITER(STRING,ERROR,*9999)
          ELSE
            CALL STAND(END,STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ABBREV(CO(noco),'DEASSIGN',3)) THEN
        CALL DEASSI(STRING,noco,CO,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'DISPLAY',1)) THEN
c       CALL DISPLAY(CO,ERROR,*9999)

C cpb 7/2/96 This should be at FEM level not ROOT level
C      ELSE IF(ABBREV(CO(noco),'EVALUATE',1)) THEN
C        noco=noco+1
C        IF(ABBREV(CO(noco),'?',1)) THEN
C          OPTION_1( 1)='Residuals'
C          OPTION_1( 2)='Exit'
C          NTCH1=2
C          CALL LIST_COMMANDS(2,NTCH1,OPTION_1,ERROR,*9999)
C        ELSE
C          IF(ABBREV(CO(noco),'RESIDUALS',1)) THEN
C            IF(CO(noco+1).EQ.'?') THEN
C              CALL STRING_TRIM(STRING,IBEG,IEND)
C              WRITE(OP_STRING,*) STRING(1:IEND)
C              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C            ELSE
C!             Note offset to PAOPTI in place of XC
C              CALL FUNCT2(0,NT_RES,NTOPTI,NREM,REALWORK(OS_PAOPTI),
C     '          REALWORK(OS_F),REALWORK(OS_FJAC),0,INTWORK,REALWORK)
C            ENDIF
C          ENDIF
C        ENDIF
C
      ELSE IF(ABBREV(CO(noco),'FEM',1)) THEN
        CALL FEM(ISEG,CSEG,END,STRING,%VAL(INTWORK_PTR),
     '    %VAL(REALWORK_PTR),ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'GEN',1)) THEN
C        CALL GEN(ISEG,CSEG,END,STRING,ERROR,*9999)
        CALL ASSERT(.FALSE.,'>>Obsolete',ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'HELP',1)) THEN
c       CALL HELP(STRING,noco,NTCO,CO,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'IMP',2)) THEN
C        CALL IMP(ISEG,CSEG,END,STRING,ERROR,*9999)
        CALL ASSERT(.FALSE.,'>>Obsolete',ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'INQUIRE',2)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          OP_STRING(1)=STRING(1:IEND)//' not used'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        ELSE
        ENDIF

      ELSE IF(ABBREV(CO(noco),'ITERATE',2)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
C#### Command: iterate <(echo/noecho)[noecho]>
C###  Description:
C###
          OP_STRING(1)=' iterate <(echo/noecho)[noecho]>'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        ELSE
          IF(CBBREV(CO,'ECHO',1,noco,NTCO,N3CO)) THEN
            ECHO=.TRUE.
          ELSE IF(CBBREV(CO,'NOECHO',1,noco,NTCO,N3CO)) THEN
            ECHO=.FALSE.
          ELSE
            ECHO=.FALSE.
          ENDIF
C cpb 7/2/96 Buffering call to iterate
          CALL ITERATE(ISEG,CSEG,ECHO,END,STRING,%VAL(INTWORK_PTR),
     '      %VAL(REALWORK_PTR),ERROR,*9999)
        ENDIF

      ELSE IF(ABBREV(CO(noco),'LIST',1)) THEN
        noco=noco+1
        IF(ABBREV(CO(noco),'?',1)) THEN
C#### Command: list
C###  Parameter:  Assign
C###  Parameter:  Iterate
C###  Parameter:  Num_threads
C###  Parameter:  Time
C###  Parameter:  Dtime
C###  Parameter:  Exit
C###  Description:
C###
          OPTION_1( 1)='Assign'
          OPTION_1( 2)='Iterate'
          OPTION_1( 3)='Num_threads'
          OPTION_1( 4)='Time'
          OPTION_1( 5)='Dtime'
          OPTION_1( 6)='Exit'
          NTCH1=6
          CALL LIST_COMMANDS(2,NTCH1,OPTION_1,ERROR,*9999)
        ELSE
          IF(ABBREV(CO(noco),'ASSIGN',1)) THEN
            CALL LIASSI(STRING,noco,CO,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco),'ITERATE',1)) THEN
            CALL LIITER(STRING,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco),'NUM_THREADS',1)) THEN
            NUM_THREADS=1
C$          NUM_THREADS=OMP_GET_MAX_THREADS()
            WRITE(OP_STRING(1),'('' NUM_THREADS = '',I3)') NUM_THREADS
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco),'TIME',1)) THEN
            CALL CPU_TIMER(CPU_TOTAL,CPUTIME)
            CALL REAL_TIMER(REAL_TOTAL,REALTIME)
            WRITE(OP_STRING,'(A,6X,F16.5,A)') ' CPU Time:      ',
     '        CPUTIME,' s'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(A,6X,F16.5,A)') ' Wallclock Time:',
     '        REALTIME,' s'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco),'DTIME',1)) THEN
            CALL DCPU_TIMER(CPUTIME)
            CALL DREAL_TIMER(REALTIME)
            WRITE(OP_STRING,'(A,F16.5,A)') ' Delta CPU Time:      ',
     '        CPUTIME,' s'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(A,F16.5,A)') ' Delta Wallclock Time:',
     '        REALTIME,' s'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          ELSE
            CALL STAND(END,STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ABBREV(CO(noco),'MENU',1)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          OP_STRING(1)=STRING(1:IEND)
     '      //' <(2D/3D/mapping/angiography)[2D]>'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        ELSE
          IF(CBBREV(CO,'MAPPING',1,noco,NTCO,N3CO)) THEN
            CALL ACWK(95,1,ERROR,*9999)
            CALL DAWK(95,1,ERROR,*9999)
          ELSE IF(CBBREV(CO,'ANGIOGRAPHY',1,noco,NTCO,N3CO)) THEN
            CALL ANGIO(ISEG,CSEG,STRING,ERROR,*9999)
          ELSE
            IF(CBBREV(CO,'2D',1,noco,NTCO,N3CO)) THEN
              NJT=2
              C_IW_NJT='1'
            ELSE IF(CBBREV(CO,'3D',1,noco,NTCO,N3CO)) THEN
              NJT=3
              C_IW_NJT='3'
            ELSE
              IF(NJT.EQ.3) THEN
                C_IW_NJT='3'
              ELSE
                C_IW_NJT='1'
              ENDIF
            ENDIF
            CO(1)='FEM'
            CO(2)='define'
            CO(3)='window'
            CO(4)='on'
            CO(5)=C_IW_NJT
            NTCO=5
            CALL FEM(ISEG,CSEG,END,STRING,%VAL(INTWORK_PTR),
     '        %VAL(REALWORK_PTR),ERROR,*9999)
            CO(1)='FEM'
            CO(2)='define'
            CO(3)='axes'
            CO(4)='on'
            CO(5)=C_IW_NJT
            NTCO=5
            COQU(3,1)='s'
            NTCOQU(3)=1
            CALL FEM(ISEG,CSEG,END,STRING,%VAL(INTWORK_PTR),
     '        %VAL(REALWORK_PTR),ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ABBREV(CO(noco),'OPEN',3)) THEN
        noco=noco+1
 101    IF(CO(noco).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          CALL STRING_TRIM(FILE00,IBEG1,IEND1)
C           CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C#### Command: open commands<;<PATH/>FILENAME>
C###  Description:
C###    Opens the command file FILENAME.com.

          OP_STRING(1)=STRING(1:IEND)//' commands'//
     '    '<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
C     '    //'<;(example)['//PATH00(IBEG2:IEND2)//']>'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        ELSE
          IF(NTCOQU(noco).EQ.0) THEN
            CO(noco)='?'
            STRING='>>Reenter: '//CO(noco-1)
            GO TO 101
          ENDIF
          CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*101)
C KAT 1Feb00: Closing any comfile already opened.
          IF(IREC_COMFILE(0).NE.-1) THEN
            CALL CLOSEF(OPEN_COM_UNIT,ERROR,*9999)
            IREC_COMFILE(0)=-1 !in case new com doesn't exist
          ENDIF
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPEN_COM_UNIT=FIRST_COM_UNIT
C KAT 25/2/00: If OPEN_COM_UNIT has been set to a `read' comfile
C         previously and there are no `read' comfiles open, com_unit can
C         be reset.
          IF(NEST.EQ.0) COM_UNIT=FIRST_COM_UNIT
          CALL OPENF(OPEN_COM_UNIT,'DISK',FILE(IBEG:IEND)//'.com','OLD',
     '      'DIRECT','FORMATTED',132,ERROR,*9999)
          IREC_COMFILE(0)=0
        ENDIF

      ELSE IF(ABBREV(CO(noco),'OPTIMISE',3)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: optimise <update> <warm/cold>[cold]
C###  Parameter:      <time #[0.0]>
C###    Specify the current time for a time dependent problem
C###  Parameter:      <region #[1]>
C###    Specify the region which the parameters are being optimised
C###  Parameter:      <squid_config #[0]>
C###    Specify the configuration of the SQUID being used. Refer to
C###    fem evalulate solution for more info.
C###  Parameter:      <baseline #[0.0]>
C###    Specify the baseline of the SQUID sensors
C###    Refer to fem evalulate solution for more info.
C###  Description:
C###    Invokes an optimiser (NPSOL or MINOS) to optimise a previously-setup problem.
C###    Using the warm qualifier means that the previous values of the Jacobian etc will be
C###    used (for use when an optimisation process has terminated due to eg. too many
C###    iterations and you want to continue).
C###    The update qualifier is used when gks graphics windows are to be updated during the
C###    optimisation.

          OP_STRING(1)=STRING(1:IEND)//
     '      ' <update> <warm/cold>[cold]'
          OP_STRING(2)=STRING(1:IEND)//' <time #[0.0]>'
          OP_STRING(3)=STRING(1:IEND)//' <region #[1]>'
          OP_STRING(4)=STRING(1:IEND)//' <squid_config #[0]>'
          OP_STRING(5)=STRING(1:IEND)//' <baseline #[0.0]>'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
          
        ELSE
          IF(CBBREV(CO,'UPDATE',1,noco,NTCO,N3CO)) THEN
            UPVUOP=.TRUE.
          ELSE
            UPVUOP=.FALSE.
          ENDIF
          IF(CBBREV(CO,'WARM',1,noco,NTCO,N3CO)) THEN
            WARMST=.TRUE.
          ELSE
            WARMST=.FALSE.
          ENDIF

          IF(CBBREV(CO,'TIME',3,noco,NTCO,N3CO)) THEN
            TIME=RFROMC(CO(N3CO+1))
          ELSE
            TIME=0.D0
          ENDIF

          IF(CBBREV(CO,'REGION',3,noco,NTCO,N3CO)) THEN
            nr=IFROMC(CO(N3CO+1))
          ELSE
            nr=1
          ENDIF

          IF(CBBREV(CO,'SQUID_CONFIG',3,noco,NTCO,N3CO)) THEN
            SQUID_CONFIG=IFROMC(CO(N3CO+1))
          ELSE
            SQUID_CONFIG=0
          ENDIF

          IF(CBBREV(CO,'BASELINE',3,noco,NTCO,N3CO)) THEN
            BASELINE=RFROMC(CO(N3CO+1))
          ELSE
            BASELINE=0.d0
          ENDIF

c         CALL FRPRMN(ISEG,ITER,NTOPTI,CSEG,END,PAOPTI,TOL,FMIN,
c    '      STRING,ERROR,*9999)
c         CALL POWELL(ISEG,ITER,NTOPTI,CSEG,END,PAOPTI,TOL,FMIN,
c    '      STRING,ERROR,*9999)
          !Put ISEG & CSEG into temporary storage & retrieve in FUNCT
          DO nosg=1,NTSG
            ISEG_TEMP(nosg)=ISEG(nosg)
            CSEG_TEMP(nosg)=CSEG(nosg)
          ENDDO

C cpb 7/6/96 Buffering optimistation calls

          IF(KTYP26.LE.4) THEN
C CPB 1/6/94 Adding MINOS optimisation package
            IF(KTYP29.EQ.1) THEN !NPSOL routines
              CALL NAGMIN(%VAL(INTWORK_PTR),nr,SQUID_CONFIG,
     '          BASELINE,%VAL(REALWORK_PTR),TIME,WARMST,ERROR,*9999)
            ELSE IF(KTYP29.EQ.2) THEN !MINOS routines
              IF(KTYP26.LE.2) THEN

               CALL MINOSOPTI(%VAL(INTWORK_PTR),%VAL(REALWORK_PTR),
     '           ERROR,*9999)

               ENDIF
             ENDIF
c          ELSE IF(KTYP26.EQ.2) THEN
c            CALL GEOMIN(INTWORK,REALWORK,ERROR,*9999)
          ENDIF
        ENDIF

C MPN 26Mar2002: changed >print to >fem gxprint
C      ELSE IF(ABBREV(CO(noco),'PRINT',1)) THEN
C        CALL PRSCRN(STRING,ERROR,*9999)

      ELSE IF(ABBREV(CO(noco),'READ',3)) THEN
        noco=noco+1
        IF(ABBREV(CO(noco),'?',1)) THEN
          OPTION_1( 1)='Commands'
          NTCH1=1
          CALL LIST_COMMANDS(noco,NTCH1,OPTION_1,ERROR,*9999)

        ELSEIF(ABBREV(CO(noco),'COMMANDS',1)) THEN

 201      IF(CO(noco+1).EQ.'?') THEN
            CALL STRING_TRIM(STRING,IBEG,IEND)
            CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C#### Command: read commands<;FILENAME<;example>[$current]>
C###  Description:
C###    Read commands from the file FILENAME.com.

            OP_STRING(1)=STRING(1:IEND)//
     '        '<;FILENAME<;example>['//
     '        FILE00(IBEG1:IEND1)//']>'
            CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

          ELSE
            IOTYPE=2
            CALL CHECKF(1,noco,NTCOQU,CO,COQU,FILE,STRING,*201)

            CALL STRING_TRIM(FILE,IBEG,IEND)
            IF(FILE(MAX(1,IEND-3):IEND).NE.'.com') THEN
              FILE(IEND+1:)='.com'
              IEND=IEND+4
            ENDIF
            CALL READCOM(FILE(IBEG:IEND),END,ERROR,*9999)
          ENDIF

        ELSE
          CALL STAND(END,STRING,ERROR,*9999)
        ENDIF !?/command

      ELSE IF(ABBREV(CO(noco),'REFRESH',3)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
          OP_STRING(1)='   refresh graphics'
          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
        ELSE
          IF(ABBREV(CO(noco),'GRAPHICS',1)) THEN
            CALL REGRAP(STRING,ERROR,*9999)
          ELSE
            CALL STAND(END,STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(ABBREV(CO(noco),'SET',2)) THEN
        noco=noco+1
        IF(ABBREV(CO(noco),'?',1)) THEN
          OPTION_1( 1)='Array_processor'
          OPTION_1( 2)='Diagnostics'
          OPTION_1( 3)='Directory'
          OPTION_1( 4)='Echo'
          OPTION_1( 5)='Fatal_error'
          OPTION_1( 6)='Learn'
          OPTION_1( 7)='Numerical_library'
          OPTION_1( 8)='Num_threads'
          OPTION_1( 9)='Output'
          OPTION_1(10)='Trace'
          OPTION_1(11)='Return'
          NTCH1=11
          CALL LIST_COMMANDS(2,NTCH1,OPTION_1,ERROR,*9999)

        ELSE
          IF(ABBREV(CO(noco),'ARRAY_PROCESSOR',1)) THEN

C LC 24/2/97 archived section : old MPN   7-Apr-95: unused?
C CS 2/6/98 This command is nolonger available with the current GX
C          ELSE IF(ABBREV(CO(noco),'BACKGROUND',2)) THEN
C            IF(CO(noco+1).EQ.'?') THEN
C              CALL STRING_TRIM(STRING,IBEG,IEND)
CC#### Command: set Background <colour R,G,B[1.0,1.0,1.0]> <on (WS_LIST/all)[all]>
CC###  Description:
CC###
C              WRITE(OP_STRING,*) STRING(1:IEND)
C     '          //' <colour R,G,B[1.0,1.0,1.0]>'
C     '          //' <on (WS_LIST/all)[all]>'
C              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C            ELSE
C              IF(ABBREV(CO(noco+1),'COLOUR',2)) THEN
C                CALL PARSRL(CO(noco+2),3,NTCL,COLOUR8,ERROR,*9999)
Cc                COLOUR(1)=REAL(COLOUR8(1))
Cc                COLOUR(2)=REAL(COLOUR8(2))
Cc                COLOUR(3)=REAL(COLOUR8(3))
C              ELSE
Cc                COLOUR(1)=1.0
Cc                COLOUR(2)=1.0
Cc                COLOUR(3)=1.0
C              ENDIF
C              IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
C                CALL PARSIL(CO(N3CO+1),6,NTIW,IWK,ERROR,*9999)
C              ELSE
C                NTIW=2*NJT-3+IMAP
C                DO noiw=1,NTIW
C                  IWK(noiw)=noiw
C                ENDDO
C              ENDIF
C              DO noiw=1,NTIW
C                iw=IWK(noiw)
C                CALL ACWK(iw,1,ERROR,*9999)
Cc               CALL SET_COLOUR_ONE(iw,0,COLOUR,ERROR,*9999)
C                CALL DAWK(iw,1,ERROR,*9999)
C              ENDDO
C            ENDIF

          ELSE IF(ABBREV(CO(noco),'DIRECTORY',3)) THEN
            IF(CO(noco+1).EQ.'?') THEN
              CALL STRING_TRIM(STRING,IBEG,IEND)

C#### Command: set directory [current] UNIX_PATHNAME
C###  Parameter:             example EXAMPLE_NUMBER
C###  Description:
C###    Sets the current/example directory path

              OP_STRING(1)=STRING(1:IEND)
     '          //' <(current|example)[current]> UNIX_PATHNAME'
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ELSE
              IF(CBBREV(CO,'EXAMPLE',2,noco,NTCO,N3CO)) THEN
                CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
                CALL SETEXAMPLEDIR(CO(N3CO+1)(IBEG:IEND),COMFILENAME,
     '            ERROR,*9999)
C KAT 31/3/00: Echoing example name instead of directory
                WRITE(OP_STRING,'('' >>Example is '',A)')
     '            CO(N3CO+1)(IBEG:IEND)
C                CALL STRING_TRIM(Example_subdir,IBEG,IEND)
C                WRITE(OP_STRING,'('' >>Example subdirectory is '',A)')
C     '            Example_subdir(IBEG:IEND)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ELSE IF(CBBREV(CO,'DOC',3,noco,NTCO,N3CO)) THEN
                WRITE(OP_STRING,'('' >>WARNING: doc is obsolescent. '
     '            //'Change doc to example'')')
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
                CALL SETEXAMPLEDIR(CO(N3CO+1)(IBEG:IEND),COMFILENAME,
     '            ERROR,*9999)
                WRITE(OP_STRING,'('' >>Example is '',A)')
     '            CO(N3CO+1)(IBEG:IEND)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ELSE
                IF(CBBREV(CO,'CURRENT',7,noco,NTCO,N3CO)) THEN
                  CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
                  PATH00=CO(N3CO+1)(IBEG:IEND)
                ELSE
                  CALL STRING_TRIM(CO(noco+1),IBEG,IEND)
                  PATH00=CO(noco+1)(IBEG:IEND)
                ENDIF
C KAT 1/5/00: prepending current directory
                CALL STRING_TRIM(PATH00,IBEG,IEND)
                IF(PATH00.EQ.' ') THEN
                  PATH00='./'
                ELSEIF(PATH00(IEND:IEND).NE.'/') THEN
                  IEND=IEND+1
                  PATH00(IEND:IEND)='/'
                ENDIF
                WRITE(OP_STRING,'('' >>Default directory is '',A)')
     '            PATH00(IBEG:IEND)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(noco),'DIAGNOSTICS',1)) THEN
            IF(CO(noco+1).EQ.'?') THEN
              CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: set Diagnostics<;FILENAME> ( on / off / in SUBROUTINES / from SUBROUTINES)
C###  Parameter: <thread (THREAD_NUM#/all)[all]>
C###  Description:
C###    Set diagnostic output.
              OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'//
     '          ' (on / off / in SUBROUTINES / from SUBROUTINES)'
C$            OP_STRING(2)=BLANK(1:IEND+11)//
C$   &          ' <thread (THREAD_NUM#/all)[all]>'
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ELSE
              IODI=6
              IF(NTCOQU(noco).GT.0) THEN
                FILE=COQU(noco,1)
                CALL STRING_TRIM(FILE,IBEG,IEND)
C CS 17/10/2000 CHANGING to USE IODI and setting to 18
                IODI=18
C                IO4=4
C                CALL OPENF(IO4,'DISK',FILE(IBEG:IEND)//'.opdop','NEW',
C     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
                CALL OPENF(IODI,'DISK',FILE(IBEG:IEND)//'.opdop','NEW',
     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
              ENDIF
              IF(ABBREV(CO(noco+1),'ON',2)) THEN
                DIAGNO=.TRUE.
                DOP=.TRUE.
                ALLSUB=.TRUE.
                FROMSUB=.FALSE.
C$              IF(CBBREV(CO,'THREAD',1,noco+1,NTCO,N3CO)) THEN
C$                IF(ABBREV(CO(N3CO+1),'ALL',1)) THEN
C$                  THREAD_NUM=-1
C$                ELSE IF(IFROMC(CO(N3CO+1)).GT.-1) THEN
C$                  THREAD_NUM=IFROMC(CO(N3CO+1))
C$                ENDIF
C$              ELSE
C$                THREAD_NUM=-1
C$              ENDIF
              ELSE IF(ABBREV(CO(noco+1),'IN',1)) THEN
                DIAGNO=.TRUE.
                DOP=.FALSE.
                ALLSUB=.FALSE.
                FROMSUB=.FALSE.
                CALL CUPPER(CO(noco+2),C1)
                CALL PARSSL(C1,MXSUB,NT_SUB,SUBNAM,ERROR,*9999)
C                CALL PARSSL(CUPPER(CO(noco+2)),MXSUB,NT_SUB,SUBNAM,
C     '            ERROR,*9999)
              ELSE IF(ABBREV(CO(noco+1),'FROM',1)) THEN
                DIAGNO=.TRUE.
                DOP=.FALSE.
                ALLSUB=.FALSE.
                FROMSUB=.TRUE.
                CALL CUPPER(CO(noco+2),C1)
                CALL PARSSL(C1,MXSUB,NT_SUB,SUBNAM,ERROR,*9999)
C                CALL PARSSL(CUPPER(CO(noco+2)),MXSUB,NT_SUB,SUBNAM,
C     '            ERROR,*9999)
              ELSE IF(ABBREV(CO(noco+1),'OFF',2)) THEN
                DIAGNO=.FALSE.
                DOP=.FALSE.
                ALLSUB=.FALSE.
                FROMSUB=.FALSE.
                IF(IO4.EQ.4) THEN
                  CALL CLOSEF(4,ERROR,*9999)
                  IO4=IOOP
                ENDIF
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(noco),'ECHO',1)) THEN
            ECHO_RAW_COM=ECHO_INTERP_COM
            ECHO_INTERP_COM=.NOT.ECHO_INTERP_COM
C            CALL SET_ECHO(STRING,noco,NTCO,CO,ERROR,*9999)

          ELSE IF(ABBREV(CO(noco),'FATAL_ERROR',1)) THEN
            IF(CO(noco+1).EQ.'?') THEN
              CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: set fatal (on/off)
C###  Description:
C###    Sets fatal error handler on or off
              WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)//' (on/off)'
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ELSE
c cpb 26/10/95 changing condition handlers
              IF(ABBREV(CO(noco+1),'ON',2)) THEN
                CALL SET_HANDLER(ERROR,*9999)
              ELSE IF(ABBREV(CO(noco+1),'OFF',2)) THEN
                CALL RESET_HANDLER(ERROR,*9999)
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(noco),'LEARN',1)) THEN
            IF(CO(noco+1).EQ.'?') THEN
              CALL STRING_TRIM(STRING,IBEG,IEND)
              IF(LEARN) THEN
                OP_STRING(1)=STRING(1:IEND)//' off'
              ELSE IF(.NOT.LEARN) THEN
                OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> on'
              ENDIF
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ELSE
              IF(NTCOQU(noco).GT.0) THEN
                FILE=COQU(noco,1)
                CALL STRING_TRIM(FILE,IBEG,IEND)
                FILE00=FILE(IBEG:IEND)
                COM_UNIT=COM_UNIT+1
                NEST=NEST+1
                CALL OPENF(COM_UNIT,
     '            'DISK',FILE(IBEG:IEND)//'.com','NEW',
     '            'SEQUEN','FORMATTED',132,ERROR,*9999)
              ENDIF
              IF(ABBREV(CO(noco+1),'ON',2)) THEN
                LEARN=.TRUE.
              ELSE IF(ABBREV(CO(noco+1),'OFF',2)) THEN
                LEARN=.FALSE.
                CALL CLOSEF(COM_UNIT,ERROR,*9999)
                COM_UNIT=COM_UNIT-1
                NEST=NEST-1
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(noco),'NUMERICAL_LIBRARY',4)) THEN
            IF(CO(noco+1).EQ.'?') THEN
              CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: set numerical_library <(nag/lapack)[nag]>
C###  Description:
C###    Allows choice of numerical routine, obsolescent
              OP_STRING(1)=STRING(1:IEND)//' <(nag/lapack)[nag]>'
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ELSE
              IF(ABBREV(CO(noco+1),'NAG',1)) THEN
                NUM_LIBRARY=0 !NAG
              ELSE IF(ABBREV(CO(noco+1),'LAPACK',1)) THEN
                NUM_LIBRARY=1 !LAPACK
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(noco),'NUM_THREADS',4)) THEN
            CALL SET_NUM_THREADS(STRING,ERROR,*9999)

          ELSE IF(ABBREV(CO(noco),'OUTPUT',1)) THEN
 301        IF(CO(noco+1).EQ.'?') THEN
              CALL STRING_TRIM(STRING,IBEG,IEND)
              CALL STRING_TRIM(FILE00,IBEG1,IEND1)
C#### Command:  set output<;FILENAME>[file] on
C###  Description:
C###    Sets output to file FILENAME.out
              OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>['
     '          //FILE00(IBEG1:IEND1)//'] on'
C#### Command:  set output off
C###  Description:
C###    Stops output to file
              OP_STRING(2)=STRING(1:IEND)//' off'
              CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
            ELSE
              IF(ABBREV(CO(noco+1),'ON',2)) THEN
                ECHO_OUTPUT=.TRUE.
                CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*301)
                CALL STRING_TRIM(FILE,IBEG,IEND)
                CALL OPENF(IOOUT,'DISK',FILE(IBEG:IEND)//'.out','NEW',
     '            'SEQUEN','FORMATTED',0,ERROR,*9999)
              ELSE IF(ABBREV(CO(noco+1),'OFF',2)) THEN
                ECHO_OUTPUT=.FALSE.
                CALL CLOSEF(IOOUT,ERROR,*9999)
              ENDIF
            ENDIF

          ELSE IF(ABBREV(CO(noco),'TRACE',1)) THEN
            CALL TRACE(STRING,noco,NTCO,CO,ERROR,*9999)

          ELSE
            CALL STAND(END,STRING,ERROR,*9999)
          ENDIF
        ENDIF !?/command

      ELSE IF(ABBREV(CO(noco),'SYSTEM',2)) THEN
        noco=noco+1
        IF(CO(noco).EQ.'?') THEN
          CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command:  system
C###  Description:
C###    Execute a shell comand
        ELSE
            CALL SYSTEM(CO(noco))
        ENDIF

      ELSE IF(ABBREV(CO(noco),'UNITS',1)) THEN
        noco=noco+1
c       CALL UNITS(CO,noco,STRING,ERROR,*9999)

C rgb 21/09/1999
C      ELSE IF(ABBREV(CO(noco),'VALIDATE',1)) THEN
C        noco=noco+1
C        IF(CO(noco).EQ.'?') THEN
C          CALL STRING_TRIM(STRING,IBEG,IEND)
CC#### Command: validate <(commands/current)>
CC###  Description:
CC###
C          WRITE(OP_STRING,*) STRING(1:IEND)//' <(commands/current)>'
C          CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C        ELSE
C          IF(CBBREV(CO,'COMMANDS',2,noco,NTCO,N3CO)) THEN
CC#### Command: validate commands commands<;FILENAME<;example>[file]>
CC###  Description:
CC###
C            CO(1)='read'
C            CO(2)='commands'
C            NTCO=2
C            COQU(2,1)='com'
C            COQU(2,2)='example'
C            NTCOQU(2)=2
C            CALL READC(ISEG,CSEG,STRING,END,ERROR,*9999)
C          ENDIF
C        ENDIF

      ELSE
        IF(COMMAND) THEN
          CALL STAND(END,STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(CMGUI_LINK) THEN
C       Update the front end
        CALL WH_INPUT_F_UPDATE(CMGUI_COMMAND_I,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='Could not update command input'
          GOTO 9999
        ENDIF
        CALL WH_INPUT_F_UPDATE(CMGUI_DATA_I,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='Could not update data input'
          GOTO 9999
        ENDIF
      ENDIF

      CALL EXITS('SYNTAX')
      RETURN
 9999 CALL ERRORS('SYNTAX',ERROR)
      CALL EXITS('SYNTAX')
      RETURN 1
      END


