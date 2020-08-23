      SUBROUTINE DEMOTI(IBT,NBH,NBJ,NEELEM,NELIST,NELIST2,NENP,NHE,NHP,
     &  NKH,NKJ,NKJE,NORD,NPNE,NPNODE,NRLIST,NVJE,NVJP,NXI,NXLIST,
     &  BBM,CE,XAB,XP,YP,FIX,STRING,ERROR,*)

C#### Subroutine: DEMOTI
C###  Description:
C###    DEMOTI defines motion parameters with prompted input or from
C###    filename.ipmoti.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'moti00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBJ(NJM,NEM),NBH(NHM,NCM,NEM),
     &  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NELIST2(0:NEM),
     &  NENP(NPM,0:NEPM,0:NRM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     &  NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     &  NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM,NXM),XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(*)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER ERR,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,IPFILE, 
     &  NJJ_COEFF,NJJ_FLOW,NJJ_PRESSURE,NJJ_RADIUS,NJJ_PBO2,NJJ_PAO2,
     &  NJJ_SOURCE,nr,nx,nxc,N3CO,ne,noelem,r_index
      REAL*8 comp,FRC,MINPCNT,RFROMC,TLC 
      CHARACTER FILE*(MXCH),SCALE*(10),STATUS*3 
      LOGICAL ALL_REGIONS,BELOW,CALCU,CBBREV,FILIO,FIRST_TIME,GENER, 
     &  LUMPED_PARAMETER,MOUSE,SET_FRC,SET_TLC,VMIN_STATE 

      CALL ENTERS('DEMOTI',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CALL STRING_TRIM(PATH00,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM define motion;d/l/p/r/w<;FILENAME[$current]><;(PATH/example)[$current]>
C###  Description:
C###    Define motion parameters. Motion parameters can be read from or written
C###    to the file FILENAME.ipmoti in the directory specified by PATH, with
C###    $current specifing the current default file.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';d/l/p/r/w'
     '    //'<;FILENAME['//FILE00(IBEG1:IEND1)//']>'
     '    //'<;(PATH/example)['//PATH00(IBEG2:IEND2)//']>'
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM define motion;c lung
C###  Description:
C###    Define motion parameters for pulmonary tidal breathing.
C###  Parameter:      <tidal [VT]>
C###    Specify the tidal volume in litres.
C###  Parameter:      <Tinsp [Ti]>
C###    Specify the time for inspiration in seconds.
C###  Parameter:      <Bhold [Tbh]>
C###    Specify the time for breath hold in seconds.
C###  Parameter:      <Texpn [Te]>
C###    Specify the time for expiration in seconds.
C###  Parameter:      <constant/sine [constant]>
C###    Specify whether flow is constant or fit to sinusoid.
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//';c lung'
        OP_STRING(2)=BLANK(1:15)//'<tidal [VT]>'
        OP_STRING(3)=BLANK(1:15)//'<Tinsp [Ti]>'
        OP_STRING(4)=BLANK(1:15)//'<Bhold [Tbh]>'
        OP_STRING(5)=BLANK(1:15)//'<Texpn [Te]>'
        OP_STRING(6)=BLANK(1:15)//'<FRC [FRC]>'
        OP_STRING(7)=BLANK(1:15)//'<TLC [TLC]>'
        OP_STRING(8)=BLANK(1:15)//'<constant/sine/pv [constant]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(10)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DEMOTI',ERROR,*9999)
      ELSE
        IPFILE=1 !is input file version number on 24-Jan-1990

        CALL PARSE_QUALIFIERS('CDLPRW',noco,1,CO,COQU,
     '    CALCU,FILIO,GENER,MOUSE,STATUS,ERROR,*1)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CBBREV(CO,'ELEMENTS',3,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,ERROR,
     &      *9999)
        ELSE
          NELIST(0)=0
        ENDIF
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        IF(FILIO) CALL CHECKF(2,noco,NTCOQU,CO,COQU,FILE,STRING,*1)

C CPB 8/6/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
c        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
c     '    ERROR,*9999)

c cpb 20/3/95 Need to loop over regions here

        nr=NRLIST(1)
        IF(FILIO) THEN
          FIRST_TIME=.TRUE.
          BACKUP=.FALSE.
          DO WHILE(FIRST_TIME.OR.BACKUP)
            FIRST_TIME=.FALSE.
            CALL OPEN_FILE(ALL_REGIONS,IPFILE,FILE,'moti',
     '        STATUS,ERR,ERROR,*9999)
            CALL IPMOTI(IBT,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NKJ,nr,YP(1,1,nx),FIX(1,1,nx),ERROR,*9999)
            CALL CLOSEF(IFILE,ERROR,*9999)
            IF(BACKUP) IOTYPE=2
          ENDDO
        ELSE IF(CALCU)THEN
          IF(CBBREV(CO,'LUNG',3,noco+1,NTCO,N3CO)) THEN
            KTYP58(nr)=4
            IF(CBBREV(CO,'UNIFORM',4,noco+1,NTCO,N3CO)) THEN
              UNIFORM_DISTRIBUTION=.TRUE.
            ELSE
              UNIFORM_DISTRIBUTION=.FALSE.
            ENDIF
            IF(CBBREV(CO,'LUMPED_PARAMETER',4,noco+1,NTCO,N3CO)) THEN
              LUMPED_PARAMETER=.TRUE.
c              NELIST(0)=0
            ELSE
              LUMPED_PARAMETER=.FALSE.
            ENDIF
            VMIN_STATE=.FALSE.
            IF(CBBREV(CO,'INLET_FLOW',2,noco+1,NTCO,N3CO)) 
     &        INLET_FLOW(nx)=RFROMC(CO(N3CO+1)) !*1.d6 !(L to mm^3)
            IF(CBBREV(CO,'FLOW_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_FLOW=IFROMC(CO(N3CO+1))
            ELSE
              NJJ_FLOW=0
            ENDIF
            IF(CBBREV(CO,'PRESSURE_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_PRESSURE=IFROMC(CO(N3CO+1))
            ELSE
              NJJ_PRESSURE=0
            ENDIF
            IF(CBBREV(CO,'RADIUS_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_RADIUS=IFROMC(CO(N3CO+1))
            ELSE
              NJJ_RADIUS=0
            ENDIF
            IF(CBBREV(CO,'COEFFICIENT_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_COEFF=IFROMC(CO(N3CO+1))
              nj_coeff=NJ_LOC(NJL_FIEL,NJJ_COEFF,nr)
            ELSE
              NJJ_COEFF=0
            ENDIF
            IF(CBBREV(CO,'FLOW_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_FLOW=IFROMC(CO(N3CO+1))
              nj_flow=NJ_LOC(NJL_FIEL,NJJ_FLOW,nr)
            ENDIF
            IF(CBBREV(CO,'PBO2_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_PBO2=IFROMC(CO(N3CO+1))
            ELSE
              NJJ_PBO2=0
            ENDIF
            IF(CBBREV(CO,'PAO2_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_PAO2=IFROMC(CO(N3CO+1))
            ELSE
              NJJ_PAO2=0
            ENDIF
            IF(CBBREV(CO,'SOURCE_FIELD',4,noco+1,NTCO,N3CO)) THEN
              NJJ_SOURCE=IFROMC(CO(N3CO+1))
            ELSE
              NJJ_SOURCE=0
            ENDIF

            IF(CBBREV(CO,'UNIFORM',2,noco+1,NTCO,N3CO)) 
     &        UNIFORM_DISTRIBUTION=.TRUE.
              
            MIXING=1 !perfect gas mixing efficiency
            IF(CBBREV(CO,'MIXING',2,noco+1,NTCO,N3CO))THEN
              noco=N3CO
              IF(CBBREV(CO,'PERFECT',3,noco+1,NTCO,N3CO))THEN
                MIXING=1
              ELSE IF(CBBREV(CO,'REGRESSION',3,noco+1,NTCO,N3CO))THEN
                MIXING=2
              ENDIF
            ENDIF

            IF(CBBREV(CO,'FRC',3,noco+1,NTCO,N3CO)) THEN
              SET_FRC=.TRUE.
              FRC=RFROMC(CO(N3CO+1)) !*1.d6 !entered in litres
              IF(CBBREV(CO,'DIAMETER',3,noco+1,NTCO,N3CO)) THEN
                SCALE='DIAMETER'
              ELSE IF(CBBREV(CO,'LENGTH',3,noco+1,NTCO,N3CO)) THEN
                SCALE='LENGTH'
              ELSE IF(CBBREV(CO,'VOLUME',3,noco+1,NTCO,N3CO)) THEN
                SCALE='VOLUME'
              ELSE
                SCALE='VOLUME'
              ENDIF
c              IF(LUMPED_PARAMETER)THEN
                IF(CBBREV(CO,'MINIMUM_PERCENT',3,noco+1,NTCO,N3CO)) THEN
                  VMIN_STATE=.TRUE.
                  MINPCNT=RFROMC(CO(N3CO+1))
                ENDIF
c              ENDIF
            ELSE
              SET_FRC=.FALSE.
            ENDIF
            IF(CBBREV(CO,'TLC',3,noco+1,NTCO,N3CO)) THEN
              SET_TLC=.TRUE.
              TLC=RFROMC(CO(N3CO+1))*1.d6 !TLC entered in litres
            ELSE
              SET_TLC=.FALSE.
            ENDIF

 ! true if the element group defines the boundary between elements
 ! that are fixed and elements that are changing volume
            IF(CBBREV(CO,'BELOW',3,noco+1,NTCO,N3CO)) THEN
              BELOW=.TRUE.
            ELSE
              BELOW=.FALSE.
            ENDIF
            
            IF(PRESSURE_DISTRIBUTION)THEN
              IF(CBBREV(CO,'VFILE',3,noco+1,NTCO,N3CO)) THEN
                FB_EXPN_TIME=IFROMC(CO(N3CO+1))
c	      write(*,*)'FB_EXPN_TME',FB_EXPN_TIME
              ENDIF
              IF(CBBREV(CO,'RM_a',3,noco+1,NTCO,N3CO)) THEN
                REMOVE_a=.TRUE.
              ENDIF
              IF(CBBREV(CO,'NOG',3,noco+1,NTCO,N3CO)) THEN
                NO_GRAV=.TRUE.
              ENDIF
              IF(CBBREV(CO,'ASET',3,noco+1,NTCO,N3CO)) THEN
                A_VEN=RFROMC(CO(N3CO+1))
                ASET=.TRUE.
              ENDIF
              IF(CBBREV(CO,'BSET',3,noco+1,NTCO,N3CO)) THEN
                B_VEN=RFROMC(CO(N3CO+1))
                BSET=.TRUE.
              ENDIF
              IF(CBBREV(CO,'CSET',3,noco+1,NTCO,N3CO)) THEN
                C_VEN=RFROMC(CO(N3CO+1))
              ENDIF
              IF(CBBREV(CO,'DSET',3,noco+1,NTCO,N3CO)) THEN
                D_VEN=RFROMC(CO(N3CO+1))
              ENDIF
              IF(CBBREV(CO,'RSET',3,noco+1,NTCO,N3CO)) THEN
                C_RND=RFROMC(CO(N3CO+1))
              ENDIF
              IF(CBBREV(CO,'CSPRD',3,noco+1,NTCO,N3CO)) THEN
                C_SPREAD=RFROMC(CO(N3CO+1))
              ENDIF
              IF(CBBREV(CO,'CPL',3,noco+1,NTCO,N3CO)) THEN
                comp=RFROMC(CO(N3CO+1))!amount to change compliance by
                CHANGE_COMPLIANCE=comp
              ENDIF

c.............Use to set up an end inspiratory breath hold
c.............Note T_insp=Breath_hold+T_insp              
              IF(CBBREV(CO,'BHO',3,noco+1,NTCO,N3CO))THEN
                BREATH_HOLD=RFROMC(CO(N3CO+1)) !TLC entered in cmH20
c               write(*,*)BREATH_HOLD
              ENDIF

              IF(CBBREV(CO,'RCMP',3,noco+1,NTCO,N3CO))THEN
                RNDCOMP=.TRUE.
              ENDIF

              IF(CBBREV(CO,'INDEX',3,noco+1,NTCO,N3CO))THEN
                r_index=IFROMC(CO(n3CO+1))
              ENDIF

              IF(CBBREV(CO,'VALUE',3,noco+1,NTCO,N3CO))THEN
                R_VALUE(r_index)=RFROMC(CO(n3CO+1))
              ENDIF
              
c.............Required to constrict the airways dynamically in MESH_FLOW
              IF(CBBREV(CO,'CONSTRICTION',3,noco+1,NTCO,N3CO))THEN
                CONSTRICTION=.TRUE.
                CONSTRICT=RFROMC(CO(N3CO+1)) !Narrowing as percent
              ENDIF
              
              IF(CBBREV(CO,'ZED',3,noco+1,NTCO,N3CO))THEN
                ZED=RFROMC(CO(N3CO+1)) !normalised height
                ZED_DISTRIBUTION=.TRUE.
c                write(*,*)'ZED=',ZED
c                write(*,*)'ZED_DISTRIBUTION=',ZED_DISTRIBUTION
              ENDIF
              
              IF(CBBREV(CO,'SQUISHED',3,noco+1,NTCO,N3CO))THEN
                CNC=IFROMC(CO(N3CO+1)) !index of array   
              ENDIF

              IF(CBBREV(CO,'PATCHY',3,noco+1,NTCO,N3CO))THEN
                PATCHY_PERCENT(CNC)=IFROMC(CO(N3CO+1)) !mask
c              DO nn=1,CNC
c                 write(*,*)'index is',nn,'of MASK',PATCHY_PERCENT(nn)
c              ENDDO
              ENDIF
              IF(CBBREV(CO,'ATELEM',3,noco+1,NTCO,N3CO))THEN
                USE_ATELEM=.TRUE.
                ATELEM(CNC)=IFROMC(CO(N3CO+1)) !TERMINAL_ELEMENT NO
c              DO nn=1,CNC
c                write(*,*)'CNC is',nn,'at Element',ATELEM(nn)
c              ENDDO
              ENDIF
              
              IF(CBBREV(CO,'TIDAL_VOL',3,noco+1,NTCO,N3CO))THEN
                TIDAL_V=RFROMC(CO(N3CO+1))!Tidal volume in litres
                TIDAL_V=TIDAL_V*1.0d6
c		write(*,*)TIDAL_V
              ENDIF
            ENDIF  

            FLOW_CHANGE(nxc)=.TRUE.
            IF(PRESSURE_DISTRIBUTION)THEN
              CALL IPMOTI_LUNG_PV(NBJ,NEELEM(0,nr),NELIST,NELIST2,
     &          NENP(1,0,nr),NJJ_COEFF,NJJ_FLOW,NJJ_PAO2,NJJ_PBO2,
     &          NJJ_PRESSURE,NJJ_RADIUS,NJJ_SOURCE,NKJ,NKJE,NORD,
     &          NPNE,NPNODE(0,nr),nr,NVJE,NVJP,NXI,BBM,
     &          CE(1,1,nxc),FRC,MINPCNT,TLC,XAB,XP,BELOW,
     &          LUMPED_PARAMETER,SET_FRC,SET_TLC,SCALE,VMIN_STATE,
     &          ERROR,*9999)
            ELSE
              CALL IPMOTI_LUNG(NBJ,NEELEM(0,nr),NELIST,NELIST2,
     &          NENP(1,0,nr),NJJ_COEFF,NJJ_FLOW,NJJ_PAO2,NJJ_PBO2,
     &          NJJ_PRESSURE,NJJ_RADIUS,NJJ_SOURCE,NKJ,NKJE,NORD,NPNE,
     &          NPNODE(0,nr),nr,NVJE,NVJP,NXI,BBM,CE(1,1,nxc),
     &          FRC,MINPCNT,TLC,XAB,XP,BELOW,LUMPED_PARAMETER,SET_FRC,
     &          SET_TLC,SCALE,VMIN_STATE,ERROR,*9999)
            ENDIF
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              CE(1,ne,2)=CE(1,ne,1)
            ENDDO
            INLET_FLOW(2)=-INLET_FLOW(1)
            
          ELSE
            WRITE(OP_STRING,'('' DEFINE MOTION;C only '
     '        //'implemented for pulmonary problems'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(MOUSE) THEN
        ENDIF
        CALL_MOTI=.TRUE.
      ENDIF

      CALL EXITS('DEMOTI')
      RETURN
 9999 CALL ERRORS('DEMOTI',ERROR)
      CALL EXITS('DEMOTI')
      RETURN 1
      END


