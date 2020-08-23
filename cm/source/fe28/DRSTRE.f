      SUBROUTINE DRSTRE(IBT,IDO,INP,ISEG,ISSTRE,NAN,NBH,NBJ,
     '  NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '  NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RG,SE,XA,XE,XG,XIG,XP,YG,YP,
     '  ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*)

C#### Subroutine: DRSTRE
C###  Description:
C###    DRSTRE draws principal stresses for vector plots.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISSTRE(NEM,NGRSEGM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),CG(NMM,NGM),CGE(NMM,NGM,NEM,NXM),
     '  CP(NMM,NPM,NXM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,INDEX,INDEX_POLYLINE,iw,IWK(6),
     '  N3CO,ne,noelem,noiw,no_nelist,no_nrlist,nr,nr1,NTIW,nx,nxc
      REAL*8 RFROMC
      LOGICAL ALL_REGIONS,CBBREV

      CALL ENTERS('DRSTRE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw stress
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specifies the elements in which to draw the stress vectors.
C###    The "all" keyword specifies all currently defined elements.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the worksation(s). The all keyword specifies all
C###    currently defined workstations.
C###  Parameter:    <scale FACTOR#[1]>
C###    Specify the scale factor to use when drawing the stress
C###    vectors.
C###  Parameter:    <rgb=RGB[yellow]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description:
C###    Draws principal stress vectors.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<scale FACTOR#[1]>'
        OP_STRING(5)=BLANK(1:15)//'<rgb=RGB[yellow]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRSTRE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=0
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            DO noelem=1,NEELEM(0,nr)
              NELIST(noelem+NELIST(0))=NEELEM(noelem,nr)
            ENDDO !noelem
            NELIST(0)=NELIST(0)+NEELEM(0,nr)
          ENDDO !no_nrlist (nr)
        ENDIF
        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLACK')
        ENDIF

        IF(ADD) THEN
          NTSTRE=NTSTRE+1
        ELSE IF(NTSTRE.EQ.0) THEN
          NTSTRE=1
        ENDIF
        CALL ASSERT(NTSTRE.LE.NRM,'>>NRM too small',ERROR,*9999)

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          SCALE=RFROMC(CO(N3CO+1))
        ELSE !eval largest (PRSTMAX) + smallest (PRSTMIN) princ stresses
          SCALE=1.0D0
        ENDIF

        nr1=0
        DO no_nelist=1,NELIST(0)
          ne=NELIST(no_nelist)
          nr=NRE(ne)
          IF(nr.NE.nr1) THEN
            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),
     '        NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,
     '        NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)
            nr1=nr
          ENDIF
          DO noiw=1,NTIW
            iw=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            CALL SGSTRE(INDEX,IBT,IDO,INP,ISEG,ISSTRE,iw,
     '        NAN,NBH(1,1,ne),NBJ(1,ne),ne,NHE(ne,nx),
     '        NKHE(1,1,1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '        NVHE(1,1,1,ne),NVJE(1,1,1,ne),NW(ne,1,nx),nx,
     '        CE(1,ne,nx),CG,CGE(1,1,ne,nx),
     '        CP(1,1,nx),CSEG,CURVCORRECT(1,1,1,ne),
     '        FEXT(1,1,ne),PG,RG,SE(1,1,ne),XA,XE,XG,
     '        XIG,XP,YG(1,1,ne),ZA(1,1,1,ne),ZE,ZG,ZP,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO !noiw (iw)
        ENDDO !no_nelist (ne)
      ENDIF !nr.NE.nr1

      CALL EXITS('DRSTRE')
      RETURN
 9999 CALL ERRORS('DRSTRE',ERROR)
      CALL EXITS('DRSTRE')
      RETURN 1
      END


