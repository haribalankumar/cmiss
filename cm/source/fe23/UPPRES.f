      SUBROUTINE UPPRES(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NELIST,NGAP,NHE,
     '  NHP,NKH,NKHE,NKJE,NMNO,NPF,NPNE,NPNODE,NPNY,NRE,NRLIST,
     '  NVHE,NVHP,NVJE,
     '  NW,NXI,NYNE,NYNP,NYNR,CE,CP,CURVCORRECT,FEXT,PAOPTI,PF,PG,SE,WG,
     '  XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,STRING,FIX,ERROR,*)

C#### Subroutine: UPPRES
C###  Description:
C###    UPPRES updates finite element auxiliary (pressure) variable
C###    vector ZA (and YP) after undeformed and deformed nodal
C###    coordinates have been entered. Fluid velocity is related to
C###    pressure thru Darcy's law.

      IMPLICIT NONE
      INCLUDE 'acti00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'load00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ptr00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NGAP(NIM,NBM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NMNO(1:2,0:NOPM,NXM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CE(NMM,NEM,NXM),CP(NMM,NPM,NXM),CURVCORRECT(2,2,NNM,NEM),
     '  FEXT(NIFEXTM,NGM,NEM),PAOPTI(*),PF(2,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER aux,I,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,iface,IFROMC,
     '  J,K,mi,mj,mz,N1LIST,N3CO,na,na1,na2,na3,
     '  NAUX,nc,NCW,ne,ng,ngi1,ngi2,ngi3,
     '  ng_near,nh,nhx,ni,ni2,NITB,nj,nm,noelem,no_nynr,noopti,np,nr,
     '  NTLIST,nxc,nx_opt,nx,ny,nz
      PARAMETER (NCW=35) !CW must be dimen.d the same size as CE array
      REAL*8 AXU(3,3),AZ,AZL(3,3),AZU(3,3),CHTOFF(3,3,3),CONDUCTK,
     '  CW(NCW),DBM(3,3,3),delsqPGP2,delsqPGP3,DETERM,DET_F_NU,
     '  DNUREFDZ(3,3),DXIX(3,3),DXIXN(3,3),DXIZN(3,3),
     '  DXNZN(3,3),DZDNU(3,3),DZDNUREF(3,3),DZNXI(3,3),DZNXN(3,3),
     '  EG(3,3),GXL(3,3),GXU(3,3),GZ,GZL(3,3),GZU(3,3),P(5),PALIST(20),
     '  RFROMC,RGX,RI1,RI2,RI3,RZWG,SUM,SUM1,SUM2,SUM3,TC(3,3),
     '  TCNUREF(3,3),TCNUREF33,TCRC(3,3),TEMP(20),TERM1(3,3),TERM2(3,3),
     '  TERM3(3,3),TERM33,TNA,TG(3,3),XI(3),X3G(4,3),ACTIVE_STRESS
      CHARACTER CHAR1*144,CHAR2*11,STRESSTYPE*17
      LOGICAL ALL_REGIONS,BOUNDARY,CBBREV,INLIST,USE_OPT_PARAMS,NONZERO

      CALL ENTERS('UPPRES',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        nx_opt=1 !temporary
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR1,'(E11.4)') PAOPTI(1)
        DO noopti=2,NMNO(1,0,nx_opt)
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          WRITE(CHAR2,'(E11.4)') PAOPTI(noopti)
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          CHAR1=CHAR1(IBEG1:IEND1)//','//CHAR2(IBEG2:IEND2)
        ENDDO
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM update pressure
C###  Parameter:      <with PARAM_VALUES[0.0]>
C###  Specify the parameter value
C###  Parameter:      <increment FACTOR[0.0]>
C###  Specify the increment value
C###  Parameter:      <region #[1]>
C###    Specify the region numbers to update.
C###  Description: updates finite element auxiliary (pressure) variable
C###    vector ZA (and YP) after undeformed and deformed nodal
C###    coordinates have been entered. Fluid velocity is related to
C###    pressure thru Darcy's law.
C###

        OP_STRING(1)=STRING(1:IEND)
     '    //' <with PARAM_VALUES['//CHAR1(IBEG1:IEND1)//']>'
        OP_STRING(2)=BLANK(1:15)//'<increment FACTOR[0.0]>'
        OP_STRING(3)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update pressure optimisation_parameters
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//' optimisation_parameters'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------
C#### Command: FEM update pressure boundary
C###  Parameter:      <element (#s/all)[all]>
C###  Specify the parameter value
C###  Parameter:      <auxillary (1/2)[1]>
C###  Specify the parameter value
C###  Parameter:      <increment FACTOR[0.0]>
C###  Specify the increment value
C###  Description: updates finite element auxiliary (pressure)
C###  variable in YP
C###

        OP_STRING(1)=STRING(1:IEND)//' boundary'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<aux (1/2)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<increment FACTOR[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------


      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPPRES',ERROR,*9999)
      ELSE
        nxc=1 !temporary
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL ASSERT(NCW.EQ.NMM,'>>Dimension of CW array (NCW)'
     '    //' must equal dimension of CE array (NMM)',ERROR,*9999)

!Initialise DXIX,X3G
        DO i=1,3
          DO j=1,3
            DXIX(i,j)=0.0d0
          ENDDO
          DO j=1,4
            X3G(j,i)=0.0d0
          ENDDO
        ENDDO

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'BOUNDARY',2,noco+1,NTCO,N3CO)) THEN
          BOUNDARY=.TRUE.
        ELSE
          BOUNDARY=.FALSE.
        CALL ASSERT(KTYP51(nr).EQ.3,'Must be a 3D problem',ERROR,*9999)
        CALL ASSERT(KTYP53(nr).EQ.2,'Stresses must be referred to nu '//
     '    'coordinates in the constitutive law (passive stress only).',
     '    ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'INCREMENT',1,noco+1,NTCO,N3CO)) THEN
          FACTOR=RFROMC(CO(N3CO+1))
        ELSE
          FACTOR=0.0d0
        ENDIF


        IF(CBBREV(CO,'WITH',2,noco+1,NTCO,N3CO)) THEN
          USE_OPT_PARAMS=.TRUE.
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',ERROR,*9999)
          CALL PARSRL(CO(N3CO+1),20,NTLIST,PALIST,ERROR,*9999)
        ELSE IF(CBBREV(CO,'OPTIMISATION_PARAMETERS',1,noco+1,
     '      NTCO,N3CO))THEN
          USE_OPT_PARAMS=.TRUE.
          CALL NX_LOC(NX_INQUIRE,nxc,nx_opt,NX_OPTI,ERROR,*9999)
          CALL ASSERT(nx_opt.GT.0,
     '      '>>No nx defined for this optimisation class',ERROR,*9999)
          CALL ASSERT(NMNO(1,0,nx_opt).LE.20,
     '      '>>Increase dimension of PALIST',ERROR,*9999)
          DO noopti=1,NMNO(1,0,nx_opt)
            PALIST(noopti)=PAOPTI(noopti)
          ENDDO
        ELSE
          USE_OPT_PARAMS=.FALSE.
        ENDIF

        IF(BOUNDARY) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          IF(CBBREV(CO,'AUXILLARY',2,noco+1,NTCO,N3CO)) THEN
            aux=IFROMC(CO(N3CO+1))
          ELSE
            aux=1
          ENDIF
          DO nc=1,NCT(nr,nx) !loop over RHS(displ) and LHS(force) vars
            DO no_nynr=1,NYNR(0,0,nc,nr,nx) !loop over global variables
              ny=NYNR(no_nynr,0,nc,nr,nx) !is global variable number
              IF(FIX(ny,1,nx)) THEN !ny has incremented press bdry cond
                IF(NPNY(0,ny,nc,nx).NE.1) THEN ! auxillary variable
                  IF(INLIST(NPNY(4,ny,nc,nx),NELIST(1),
     '              NELIST(0),N1LIST)) THEN
C new MPN/VY/IK 16Dec2008
                    IF((NYNE(aux,4,1,nc,NPNY(4,ny,nc,nx)).EQ.ny))
     '                THEN
                      YP(ny,1,nx)=YP(ny,1,nx)+YP(ny,2,nx)*FACTOR
                    ENDIF
C old
C                    IF((NYNE(1,4,1,nc,NPNY(4,ny,nc,nx)).EQ.ny)
C     '                .AND.(aux.EQ.1)) THEN
C                       YP(ny,1,nx)=YP(ny,1,nx)+YP(ny,2,nx)*FACTOR
C                    ELSE IF ((NYNE(2,4,1,nc,NPNY(4,ny,nc,nx)).EQ.ny)
C     '                .AND.(aux.EQ.2)) THEN
C                       YP(ny,1,nx)=YP(ny,1,nx)+YP(ny,2,nx)*FACTOR
C                    ENDIF
C end old
                  ENDIF
                ENDIF
              ENDIF
            ENDDO !no_nynr (ny)
          ENDDO !nc

C ***   Calculate the auxiliary parameters
        ELSE IF(KTYP52(nr).EQ.2) THEN      !Incompressibility constraint
          WRITE(OP_STRING,'('' >>Cannot calculate auxiliary'
     '      //' pressure parameters for incompressible (only)'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GOTO 9999

        ELSE IF(KTYP52(nr).EQ.3) THEN !Incompressible + Fluid constraint

C ***     Apply bc increments
          DO nc=1,2 !LHS then RHS vars
            DO no_nynr=1,NYNR(0,0,nc,nr,nx) !loop over global vars
              ny=NYNR(no_nynr,0,nc,nr,nx) !is global variable number
              IF(FIX(ny,1,nx)) THEN !ny has increm'ted ess bdry cond
                IF(NPNY(0,ny,0,nx).EQ.1) THEN
                  np=NPNY(4,ny,0,nx)
C GMH 8/1/97 Update cmgui link
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                ENDIF
                YP(ny,1,nx)=YP(ny,1,nx)+YP(ny,2,nx)*FACTOR
              ENDIF
            ENDDO !no_nynr (ny)
          ENDDO !nc
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '      ERROR,*9999)

C ***     Main element loop
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nhx=NHE(ne,nx)
            nh=NH_LOC(nhx,nx)
            NITB=NIT(NBJ(1,ne))
            NAUX=NAT(NBH(nh,1,ne))
            CALL ASSERT(NAUX.GE.3,'Must have 3 aux pressure params',
     '        ERROR,*9999)
            DO na=1,NAUX
              IF(NAN(3,na,NBH(nh,1,ne)).EQ.-1) THEN
                PF(1,ne)=ZA(na,nh,1,ne)
              ELSE IF(NAN(3,na,NBH(nh,1,ne)).EQ.-2) THEN
                PF(2,ne)=ZA(na,nh,1,ne)
              ELSE IF(NAN(3,na,NBH(nh,1,ne)).GE.0) THEN
                P(na)=0.0d0
                ZA(na,nh,1,ne)=0.0d0
                IF(NAN(3,na,NBH(nh,1,ne)).EQ.0) THEN
                  na1=na
                ELSE IF(NAN(3,na,NBH(nh,1,ne)).EQ.1.OR.
     '              NAN(3,na,NBH(nh,1,ne)).EQ.3) THEN
                  na2=na
                ELSE IF(NAN(3,na,NBH(nh,1,ne)).EQ.2.OR.
     '              NAN(3,na,NBH(nh,1,ne)).EQ.4.OR.
     '              NAN(3,na,NBH(nh,1,ne)).EQ.5) THEN
                  na3=na
                ENDIF
              ENDIF
            ENDDO
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Element: ne='',I2)')ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''   nh='',I2,'' nb='',I2)')
     '          nh,NBH(nh,1,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

            IF(USE_OPT_PARAMS) THEN
C             Store material params temporarily for the current elem
              CALL ASSERT(NMNO(1,0,nx_opt).LE.20,
     '          '>>Increase dim of TEMP',ERROR,*9999)
              DO noopti=1,NMNO(1,0,nx_opt)
                nm=NMNO(1,noopti,nx_opt)
C!!!            NOTE: this needs updating for nodally based
C!!!            material parameter interpolation
                TEMP(noopti)=CE(nm,ne,nx)
                CE(nm,ne,nx)=PALIST(noopti)
              ENDDO
            ENDIF

            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),
     '        NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)

C ***       Calculate the first two pressure parameters using the bdy
C ***       conditions applied at the internal and external faces.
            DO mi=1,NITB
              XI(mi)=0.5d0
            ENDDO
            DO iface=1,2

C BAD CODING
C              IF(iface.EQ.1.AND.(NXI(-3,1,ne).EQ.0.OR.
C     '          nr.NE.NRE(NXI(-3,1,ne)))
C     '          .AND.(NW(ne,1,nx).EQ.2.OR.NW(ne,1,nx).EQ.4).OR.
C     '          iface.EQ.2.AND.(NXI( 3,1,ne).EQ.0.OR.
C     '          nr.NE.NRE(NXI( 3,1,ne)))
C     '          .AND.(NW(ne,1,nx).EQ.3.OR.NW(ne,1,nx).EQ.4)) THEN

              IF(iface.EQ.1) THEN
                ni2=-3
              ELSE
                ni2=3
              ENDIF
              NONZERO=.FALSE.
              IF(NXI(ni2,1,ne).EQ.0) THEN
                NONZERO=.TRUE.
              ELSE
                IF(nr.NE.NRE(NXI(ni2,1,ne))) NONZERO=.TRUE.
              ENDIF
              IF(iface.EQ.1.AND.NONZERO.AND.
     '          (NW(ne,1,nx).EQ.2.OR.NW(ne,1,nx).EQ.4).OR.
     '          iface.EQ.2.AND.NONZERO.AND.
     '          (NW(ne,1,nx).EQ.3.OR.NW(ne,1,nx).EQ.4)) THEN

C                 ext press bc on Xi3=0 face or
C                 ext press bc on Xi3=1 face

C old
C              IF(NW(ne,1,nx).NE.(4-iface) !press bc on current face
C     '          .OR.iface.EQ.1.AND.NXI(-3,1,ne).NE.0
CC                Xi3=0 face in 2 adj elems
C     '          .OR.iface.EQ.2.AND.NXI( 3,1,ne).NE.0) THEN
CC                Xi3=1 face in 2 adj elems
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' Face '',I2,''  PF='',D12.4)')
     '              iface,PF(iface,ne)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
C               IFE=iface+4

C!!! The variable nf has been put into the 'set but not
C!!! used' .omit file to prevent it showing in the parsed fortran
C!!! check output, because the line that uses it has been comment out
C!!! When the line is uncomment remove nf from the .omit file
C!!! CS 12-FEB-1997
C               nf=NFF(IFE,ne)
C                NBJF=NPF(1,nf) !basis fn for first geometric var
                XI(3)=DBLE(iface-1)

                CALL CPXI(1,IBT,IDO,INP,NPNE(1,1,ne),nr,nx,
     '            CE(1,1,nx),CP(1,1,nx),CW,XI,ERROR,*9999)
C old
C                DO il=1,ILT(1,nr,nx)
C                  IF(ILP(il,1,nr,nx).EQ.3) THEN
C                    CW(il)=0.d0
C                    DO nnbf=1,NNT(NBFF)
C                      CW(il)=CW(il)
C     '                  +CP(il,NPNE(nnbf,NBH(NH_LOC(1,nx),1,ne),ne),nx)
C                    ENDDO
C                    CW(il)=CW(il)/DBLE(NNT(NBFF))
C                  ELSE IF(ILP(il,1,nr,nx).ne.3) THEN
C                    CW(il)=CE(il,ne,nx)
C                  ENDIF
C                ENDDO
C               Interpolate midwall geometric vars XG and derivs wrt Xi
                CALL XEXW(0,IBT,IDO,INP,NAN,NBJ(1,ne),nr,XE,XG,XI,
     '            ERROR,*9999)
C               Calc undeformed metric tensors wrt Xi (GXL,GXU) and
C               derivs of Xi wrt undeformed nu coords (DXIXN)
                CALL XGMG(1,NITB,NBJ(1,ne),nr,DXIXN,GXL,GXU,RGX,XG,
     '            ERROR,*9999)
C               Interpolate dep var.s ZG and derivs wrt undeformed nu
                CALL ZEZW(0,1,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NHE(ne,nx),nr,nx,DXIXN,ZE,ZG,XI,ERROR,*9999)
C               Calc deformed metric tensors wrt undef nu (AZL,AZU)
                CALL ZGMG(NBH(NH_LOC(1,nx),1,ne),nr,AZ,AZL,AZU,ZG,
     '            ERROR,*9999)
C news VJ 10Dec2003
C calculate ng_near
C!!!            Pick Gauss pt nearest to centre of current face in
C!!!            current element. To be strictly correct need to either
C!!!            store FEXT at face basis Gauss pts or somehow interp
C!!!            elem Gauss pt values to the central point on the face
                ngi1=NGAP(1,NBH(NH_LOC(1,nx),1,ne))
                ngi2=NGAP(2,NBH(NH_LOC(1,nx),1,ne))
                ngi3=NGAP(3,NBH(NH_LOC(1,nx),1,ne))
                ng_near=(ngi1+1)/2 + ((ngi2-1)/2)*ngi1
     '            + (iface-1)*(ngi3-1)*ngi2*ngi1 !need int division
C newe VJ 10Dec2003
C               Get contravariant cpts of 2nd Piola-Kirchhoff stress
C               tensor (TG) wrt undeformed nu coordinates
                IF(KTYP51(nr).EQ.1) THEN
                  CALL ZGTG51(NBH(NH_LOC(1,nx),1,ne),nr,nx,AXU,AZ,AZL,
     '              AZU,CE(1,ne,nx),RI1,RI2,
     '              RI3,TG,YG(1,ng_near,ne),ZG,ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.2) THEN
                  CALL ZGTG52(NBH(NH_LOC(1,nx),1,ne),nr,nx,AXU,AZ,AZL,
     '              AZU,CE(1,ne,nx),RI1,RI2,RI3,TG,
     '              YG(1,ng_near,ne),ZG,ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.3) THEN
                  STRESSTYPE=' '
                  CALL ZGTG53(STRESSTYPE,NBH(NH_LOC(1,nx),1,ne),
     '              nr,nx,AXU,AZ,AZL,AZU,CE(1,ne,nx),EG,RI1,RI2,RI3,
     '              TG,XG,YG(1,ng_near,ne),ZG,ERROR,*9999)
                ELSE IF(KTYP51(nr).EQ.4) THEN
                  CALL ZGTG54(NBH(NH_LOC(1,nx),1,ne),nr,AXU,AZ,AZL,
     '              AZU,CE(1,ne,nx),EG,RI1,RI2,RI3,TG,
     '              YG(1,ng_near,ne),ERROR,*9999)
                ENDIF
C new MPN 4-May-96: new way of handling sheets
C               Get derivs of Xi wrt deformed Nu coords, DXIZN
                CALL DXIDZM(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),
     '            0,NHE(ne,nx),nr,nx,
     '            DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Fibre',ERROR,*9999)
C               Calculate derivs of deformed Nu wrt undeformed Nu
C               (DZNXN) and inverse (DXNZN)
                DO ni=1,NITB
                  DO mi=1,NITB
                    SUM=0.0d0
                    DO k=1,NITB
                      SUM=SUM+DZNXI(ni,k)*DXIXN(k,mi)
                    ENDDO
                    DZNXN(ni,mi)=SUM
                  ENDDO
                ENDDO
                CALL INVERT(NITB,DZNXN,DXNZN,DET_F_NU)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(CHAR1,'(I1)') NITB
                  WRITE(OP_STRING,'(''  DZNXN:'','//CHAR1(1:1)//'D12.4,'
     '              //'''  DXNZN:'','//CHAR1(1:1)//'D12.4,'
     '              //'/(8X,'//CHAR1(1:1)//'D12.4,8X,'//CHAR1(1:1)
     '              //'D12.4))')
     '              ((DZNXN(mi,ni),ni=1,NITB),(DXNZN(mi,ni),ni=1,NITB),
     '              mi=1,NITB)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                IF(KTYP53(nr).EQ.3) THEN !Active stress cmpt included
C!!!              Pick Gauss pt nearest to centre of current face in
C!!!              current element. To be strictly correct need to either
C!!!              store FEXT at face basis Gauss pts or somehow interp
C!!!              elem Gauss pt values to the central point on the face
C                  ngi1=NGAP(1,NBH(NH_LOC(1,nx),1,ne))
C                  ngi2=NGAP(2,NBH(NH_LOC(1,nx),1,ne))
C                  ngi3=NGAP(3,NBH(NH_LOC(1,nx),1,ne))
C                  ng_near=(ngi1+1)/2 + ((ngi2-1)/2)*ngi1
C     '              + (iface-1)*(ngi3-1)*ngi2*ngi1 !need int division
C     
C     OR 15-08-06
C     
C     Changes to determine the proper active stress value depending 
C     on the definitions in IPACTI. KTYP59==3 got introduced. It allows
C     the user to specify a particular CellML variable to be added
C     to a specified component of the stress tensor
C     
                  IF (KTYP59(nr).EQ.3) THEN !Active stress component defined
                                ! within CellML
                    CALL EVALASC(ne,ng_near,%VAL(NQNE_PTR),%VAL(YQS_PTR)
     &                   ,%VAL(RCQS_SPATIAL_PTR),%VAL(ICQS_SPATIAL_PTR)
     &                   ,ASC_ARRAYNAME(nr) ,ASC_CELLVARINDEX(nr)
     &                   ,ACTIVE_STRESS,ERROR,*9999)
                  ELSE
                    ACTIVE_STRESS = YG(1,ng_near,ne)
                  ENDIF                  
                  CALL ZGTG5A(NBH(NH_LOC(1,nx),1,ne),nr, FEXT(1,ng_near
     &                 ,ne),DXNZN,DZNXN,DET_F_NU,TG,TNA,ACTIVE_STRESS
     &                 ,ERROR,*9999)
C                  CALL ZGTG5A(NBH(NH_LOC(1,nx),1,ne),nr,
C     '                 FEXT(1,ng_near,ne),DXNZN,DZNXN,DET_F_NU,
C     '                 TG,TNA,YG(1,ng_near,ne),ERROR,*9999)
                ENDIF
C old
CC               Interpolate dependent var.s ZG and derivs wrt Xi
C                CALL ZEZW(0,0,IBT,IDO,INP,NAN,NBH(1,1,ne),
C     '            NHE(ne,nx),nr,nx,DXIX,ZE,ZG,XI,ERROR,*9999)
CC               Calculate deformed metric tensors wrt Xi (GZL,GZU)
C                CALL ZGMG(NBH(NH_LOC(1,nx),1,ne),nr,GZ,GZL,GZU,ZG,
C     '            ERROR,*9999)
CC               Get derivs of Xi wrt deformed nu coords, DXIZN
C                CALL dXidNu(NBH(1,1,ne),nr,DXIZN,DZNXI,GZL,GZU,XG,
C     '            ERROR,*9999)
C               Calculate derivs of deformed nu wrt undeformed nu
C                ni=3
C                DO mi=1,NITB
C                  DZNXN(ni,mi)=0.0d0
C                  DO k=1,NITB
C                    DZNXN(ni,mi)=DZNXN(ni,mi)+DZNXI(ni,K)*DXIXN(K,mi)
C                  ENDDO !k
C                ENDDO !mi
CC               Get Physical Cauchy stress TC33 wrt def Nu_3 from TG
C                TC33=0.0d0
C                SUM1=0.0d0
C                DO nj=1,NITB
C                  DO mj=1,NITB
C                    TC33=TC33+DZNXN(3,nj)*DZNXN(3,mj)*TG(nj,mj)
C                    SUM1=SUM1+DZNXN(3,nj)*DZNXN(3,mj)*AZU(nj,mj)
C                  ENDDO !mj
C                ENDDO !nj
C                TC33=TC33/DSQRT(RI3)
C                SUM1=SUM1/DSQRT(RI3)
C end old
C new MPN 5-May-96: rotate TW from def material coordinates to
C                   deformed fibre reference (wall) coordinates
C               Get Physical Cauchy stress tensor TC wrt deformed
C               nu-material coordinates from TW
                DO mz=1,3
                  DO nz=1,3
                    SUM=0.0d0
                    SUM1=0.0d0
                    DO mj=1,NITB
                      DO nj=1,NITB
                        SUM=SUM+DZNXN(mz,mj)*TG(mj,nj)*DZNXN(nz,nj)
                        SUM1=SUM1+DZNXN(mz,mj)*AZU(mj,nj)*DZNXN(nz,nj)
                      ENDDO !nj
                    ENDDO !mj
                    TC(mz,nz)=SUM/DET_F_NU
                    TERM1(mz,nz)=SUM1/DET_F_NU
                  ENDDO !nz
                ENDDO !mz
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' TC:'',12X,3D12.4,/(16X,3D12.4))')
     '              ((TC(mz,nz),nz=1,3),mz=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                IF(NJ_LOC(NJL_FIBR,0,nr).GE.2) THEN
C                 imbric (+ sheet) angles defined
C                 Compute def anatomical fibre vects wrt rc coords at XI
                  CALL MAT_VEC_DEF(IBT,IDO,INP,NAN,
     '              NBH,NBJ,0,NHE(1,1),nr,nx,
     '              DZDNU(1,1),DZDNU(1,2),DZDNU(1,3),
     '              PG,XE,XG,XI,ZE,ZG,.TRUE.,ERROR,*9999)
C                 Compute cmpts of Cauchy stress tensor wrt rc coords
C                 by rotating def material coord system into rc coords
                  DO mi=1,3
                    DO ni=1,3
                      SUM=0.0d0
                      SUM1=0.0d0
                      DO mj=1,3
                        DO nj=1,3
                          SUM=SUM+DZDNU(mi,mj)*TC(mj,nj)*DZDNU(ni,nj)
                          SUM1=SUM1+DZDNU(mi,mj)*TERM1(mj,nj)*
     '                      DZDNU(ni,nj)
                        ENDDO !nj
                      ENDDO !mj
                      TCRC(mi,ni)=SUM
                      TERM2(mi,ni)=SUM1
                    ENDDO !ni
                  ENDDO !mi
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' TCRC:'',12X,3D12.4,'
     '                //'/(18X,3D12.4))') ((TCRC(mi,ni),ni=1,3),mi=1,3)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
C                 Compute deformed fibre ref vectors wrt rc coords at XI
                  CALL FIBRE_REF_VECS_DEF(IBT,IDO,INP,NAN,NBH(1,1,ne),
     '              0,NHE(ne,nx),NITB,nr,nx,DZDNUREF(1,1),DZDNUREF(1,2),
     '              DZDNUREF(1,3),PG,XI,ZE,ZG,ERROR,*9999)
                  CALL INVERT(NITB,DZDNUREF,DNUREFDZ,DETERM)
C                 Compute cmpts of Cauchy stress tensor wrt deformed
C                 fibre reference coords by rotating rc coord system
C                 into deformed fibre ref coords
C                 NOTE: only use TCNUREF(3,3); could delete outer loops
                  DO mi=1,3
                    DO ni=1,3
                      SUM=0.0d0
                      SUM1=0.0d0
                      DO mj=1,3
                        DO nj=1,3
                          SUM=SUM+DNUREFDZ(mi,mj)*TCRC(mj,nj)*
     '                      DZDNUREF(nj,ni)
                          SUM1=SUM1+DNUREFDZ(mi,mj)*TERM2(mj,nj)*
     '                      DZDNUREF(nj,ni)
                        ENDDO !nj
                      ENDDO !mj
                      TCNUREF(mi,ni)=SUM
                      TERM3(mi,ni)=SUM1
                    ENDDO !ni
                  ENDDO !mi
                  IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                  call mp_setlock()
                    WRITE(OP_STRING,'('' TCNUREF:'',12X,3D12.4,'
     '                //'/(21X,3D12.4))')
     '                ((TCNUREF(mi,ni),ni=1,3),mi=1,3)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                  call mp_unsetlock()
                  ENDIF
                ELSE !at most fibres defined
C                 TCNUREF33 is aligned with TC(3,3) for
C                 zero imbric/sheet angles
                  TCNUREF(3,3)=TC(3,3)
                  TERM3(3,3)=TERM1(3,3)
                ENDIF
                TCNUREF33=TCNUREF(3,3)
                TERM33=TERM3(3,3)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' TCNUREF33='',D12.4)') TCNUREF33
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' TERM33='',D12.4)') TERM33
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF

                IF(NW(ne,1,nx).NE.(4-iface)) THEN !press bc on current face
                  IF(iface.EQ.1) THEN                 !Internal face
                    P(na1)=(-PF(1,ne)-TCNUREF33)/(2.0d0*TERM33)
                  ELSE IF(iface.EQ.2) THEN            !External face
                    P(na2)=(-PF(2,ne)-TCNUREF33)/(2.0d0*TERM33)-P(na1)
                  ENDIF
                ELSE IF(iface.EQ.1.AND.NXI(-3,1,ne).NE.0) THEN
C                 no pressure bc on Xi3=0 face and
C                 adjacent element exists in -Xi3 dirn
                  ERROR='Cannot calculate aux elem params'
                  GO TO 9999
                ELSE IF(iface.EQ.2.AND.NXI( 3,1,ne).NE.0) THEN
C                 no pressure bc on Xi3=1 face and
C                 adjacent element exists in +Xi3 dirn
                  ERROR='Cannot calculate aux elem params'
                  GO TO 9999
                ENDIF
              ENDIF
            ENDDO !iface loop

C ***       Calculate the third pressure parameter using Darcy's law
C ***       constraint.
            SUM1=0.0d0
            SUM2=0.0d0
            SUM3=0.0d0
C           Gauss point loop to calculate integrals
            DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
C             Interpolate Gauss pt geometric vars XG and derivs wrt Xi
              CALL XEXG(NBJ(1,ne),ng,nr,PG,
     '          XE,XG,ERROR,*9999)
C             Calculate undeformed metric tensors wrt Xi (GXL,GXU) and
C             derivs (DXIXN) of Xi wrt undeformed nu coords.
              CALL XGMG(1,NITB,NBJ(1,ne),nr,DXIXN,
     '          GXL,GXU,RGX,XG,ERROR,*9999)
C             Interpolate dependent var.s ZG and derivs wrt undef nu
              CALL ZEZG(1,NBH(1,1,ne),ng,NHE(ne,nx),nx,DXIXN,PG,
     '          ZE,ZG,ERROR,*9999)
C             Calculate deformed metric tensors wrt undef nu (AZL,AZU)
              CALL ZGMG(NBH(NH_LOC(1,nx),1,ne),nr,AZ,AZL,AZU,ZG,
     '          ERROR,*9999)
              RI3=AZ

C             Calc deformed Christoffel symbols wrt undef nu coords
C             NOTE: ZG needs derivs wrt nu not Xi !
              CALL TOFFEL(ITYP10(nr),NBJ(1,ne),nr,CHTOFF,DBM,AZU,ZG,X3G,
     '          .FALSE.,ERROR,*9999)
              SUM=0.0d0
              DO J=1,NITB
                DO K=1,NITB
                  SUM=SUM+CHTOFF(3,J,K)*AZU(J,K)
                ENDDO
              ENDDO
              delsqPGP2=PG(2,8,ng,NBH(nh,1,ne))*AZU(3,3)-
     '        PG(2,7,ng,NBH(nh,1,ne))*SUM
              delsqPGP3=PG(3,8,ng,NBH(nh,1,ne))*AZU(3,3)-
     '        PG(3,7,ng,NBH(nh,1,ne))*SUM

C             Calculate the Jacobian (wrt def coords) for integration:
C             Interpolate dependent var.s ZG and derivs wrt Xi
              CALL ZEZG(0,NBH(1,1,ne),ng,NHE(ne,nx),nx,DXIX,PG,
     '        ZE,ZG,ERROR,*9999)
C             Calculate deformed metric tensors wrt Xi (GZL,GZU)
              CALL ZGMG(NBH(NH_LOC(1,nx),1,ne),nr,GZ,GZL,GZU,ZG,
     '          ERROR,*9999)
              RZWG=DSQRT(GZ)*WG(ng,NBH(NH_LOC(1,nx),1,ne))
              IF(JTYP4.EQ.2) RZWG=RZWG*2.d0*PI*ZG(1,1) !cyl sym about x
              IF(JTYP4.EQ.3) RZWG=RZWG*2.d0*PI*ZG(2,1) !cyl sym about y
              IF(JTYP4.EQ.4) RZWG=RZWG*4.d0*PI*ZG(1,1)**2.d0 !sph sym

C new MPN 2Apr98: Fixing fluid shift constraint
              SUM1=SUM1+(DSQRT(RI3)-1.d0)*
     '          PG(3,1,ng,NBH(nh,1,ne))*RZWG
C old              SUM1=SUM1+(DSQRT(RI3)-1.d0)/DSQRT(RI3)*
C old     '          PG(3,1,ng,NBH(nh,1,ne))*RZWG
              SUM2=SUM2+delsqPGP2*PG(3,1,ng,NBH(nh,1,ne))*RZWG
              SUM3=SUM3+delsqPGP3*PG(3,1,ng,NBH(nh,1,ne))*RZWG

              IF(DOP)THEN
              WRITE(OP_STRING,'('' PG(2,7,ng,nb)='',D12.4,'
     '          //'''  PG(2,8,ng,nb)='',D12.4)')
     '          PG(2,7,ng,NBH(nh,1,ne)),PG(2,8,ng,NBH(nh,1,ne))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' PG(3,7,ng,nb)='',D12.4,'
     '          //'''  PG(3,8,ng,nb)='',D12.4)')
     '          PG(3,7,ng,NBH(nh,1,ne)),PG(3,8,ng,NBH(nh,1,ne))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' PG(3,1,ng,nb)='',D12.4)')
     '          PG(3,1,ng,NBH(nh,1,ne))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' delsqPGP2='',D12.4,'
     '          //'''  delsqPGP3='',D12.4)')
     '          delsqPGP2,delsqPGP3
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' SQRT(RI3)= '',D12.4,'
     '          //'''  RZWG= '',D12.4)') DSQRT(RI3),RZWG
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(/'' AZU:'',11X,3D12.4,/(16X,3D12.4))')
     '          ((AZU(I,J),J=1,3),I=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(/'' CHTOFF:'',8X,3D12.4,'
     '          //'2(/16X,3D12.4),2(/3(/16X,3D12.4)))')
     '          (((CHTOFF(I,J,K),K=1,3),J=1,3),I=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !Gauss pt loop

            CONDUCTK=CE(IL_fluid_conductivity,ne,nx)
            CALL ASSERT(DABS(CONDUCTK).GE.1.0D-9,'>>3rd aux param is '
     '        //'undetermined for thru wall conductivity=0',ERROR,*9999)
            CALL ASSERT(DT.GE.1.0D-9,'>>DT=0',ERROR,*9999)
            P(na3)=(-SUM1/(DT*CONDUCTK) - P(na2)*SUM2)/SUM3

            IF(DOP)THEN
              WRITE(OP_STRING,'('' sum1           : '',D12.4)')SUM1
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' sum2           : '',D12.4)')SUM2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' sum3           : '',D12.4)')SUM3
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Thru wall cond.: '',D12.4)')CONDUCTK
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Time increment : '',D12.4)')DT
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

C ***       Put pressure parameters,P(i), into ZA for element, ne.
            DO na=1,NAUX
              IF(NAN(3,na,NBH(nh,1,ne)).GE.0) THEN
                ZA(na,nh,1,ne)=P(na)
              ENDIF
            ENDDO !na loop
            IF(DOP)THEN
              WRITE(OP_STRING,'('' ZA(1..,nh='',I2,'',1,ne)= '','
     '          //'5D13.5)')
     '          nh,(ZA(na,nh,1,ne),na=1,NAUX)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(USE_OPT_PARAMS) THEN
              !Return mat params to orig values for the current elem
              DO noopti=1,NMNO(1,0,nx_opt)
                nm=NMNO(1,noopti,nx_opt)
                CE(nm,ne,nx)=TEMP(noopti)
              ENDDO
            ENDIF

          ENDDO !element loop
          CALL ZPYP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '      ERROR,*9999)

        ENDIF
      ENDIF

      CALL EXITS('UPPRES')
      RETURN
 9999 CALL ERRORS('UPPRES',ERROR)
      CALL EXITS('UPPRES')
      RETURN 1
      END


