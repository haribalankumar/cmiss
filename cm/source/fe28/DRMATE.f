      SUBROUTINE DRMATE(IBT,IDO,INP,ISEG,ISMATE,
     '  NAN,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRE,NVJE,NXLIST,
     '  CE,CQ,SE,XA,XE,XG,XP,XQ,CSEG,STRING,ERROR,*)

C#### Subroutine: DRMATE
C###  Description:
C###    DRMATE defines materials with prompted input or from
C###    filename.mat.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b13.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gks001.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISMATE(NWM,NEM),NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXLIST(0:NXM)
      REAL*8 CE(NMM,NEM,NXM),CQ(NMM,NQM,NXM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  XQ(NJM,NQM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER i,IBEG,IBEG1,IBEG3,ID_AXES(0:3),IEND,IEND1,IEND3,
     '  IFROMC,iw,IWK(6),
     '  N3CO,ne,nm,noelem,nolist,nq,nr,NTIW,NT_LENGTHS,
     '  NT_XI1,NT_XI2,NT_XI3,nx,nxc
      REAL*8 AXES_LENGTHS(3),AXES_XI(3,3)
      CHARACTER CHAR3*9,TYPE*5
      LOGICAL CBBREV

      CALL ENTERS('DRMATE',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        CALL STRING_TRIM(FILE00,IBEG1,IEND1)

C---------------------------------------------------------------------

C#### Command: FEM draw materials
C###  Parameter:    <axes AXES#s[1,2,3]>
C###   Specifies the material axes to be drawn
C###  Parameter:    <lengths XI#s[0.1,0.1,0.1]>
C###   Specifies the XI lengths
C###  Parameter:    <xi_1 RANGE#s[0.0,1.0,0.2]>
C###   Specifies the xi1 range
C###  Parameter:    <xi_2 RANGE#s[0.0,1.0,0.2]>
C###   Specifies the xi1 range
C###  Parameter:    <xi_3 RANGE#s[0.0,1.0,0.2]>
C###   Specifies the xi3 range
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###   Specifies the element numbers
C###  Parameter:    <on WS_ID#[1]>
C###  Description: draws the material axes
C###

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<axes AXES#s[1,2,3]>'
        OP_STRING(3)=BLANK(1:15)//'<lengths XI#s[0.1,0.1,0.1]>'
        OP_STRING(4)=BLANK(1:15)//'<xi_1 RANGE#s[0.0,1.0,0.2]>'
        OP_STRING(5)=BLANK(1:15)//'<xi_2 RANGE#s[0.0,1.0,0.2]>'
        OP_STRING(6)=BLANK(1:15)//'<xi_3 RANGE#s[0.0,1.0,0.2]>'
        OP_STRING(7)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(8)=BLANK(1:15)//'<on WS_ID#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        IF(NINDICES.GE.150) THEN
          CHAR3='colour'
        ELSE
          CHAR3='greyscale'
        ENDIF
        CALL STRING_TRIM(CHAR3,IBEG3,IEND3)

C---------------------------------------------------------------------

C#### Command: FEM draw materials field
C###  Parameter:    <number MATERIAL_ID#[1]>
C###   Specifies the material number
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###   Specifies the element numbers
C###  Parameter:    <region #[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <on WS_ID#[1]>
C###   Specifies the work station number
C###  Parameter:    <(greyscale/colour)[greyscale]>
C###    Draws field in greyscale or colour.
C###  Description: draws the material field
C###

        OP_STRING(1)=STRING(1:IEND)//' field'
        OP_STRING(2)=BLANK(1:15)//'<number MATERIAL_ID#[1]>'
        OP_STRING(3)=BLANK(1:15)//'<in (all/ELEMENT#s)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<region #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(6)=BLANK(1:15)//'<on WS_ID#[1]>'
        OP_STRING(7)=BLANK(1:15)
     '    //'<(greyscale/colour)['//CHAR3(IBEG3:IEND3)//']>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRMATE',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)

        CALL ASSERT(CALL_EQUA,'>>Problem type not defined',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)

        IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
        ELSE IF(CBBREV(CO,'AXES',2,noco+1,NTCO,N3CO)) THEN
          TYPE='AXES'
          CALL PARSIL(CO(N3CO+1),3,ID_AXES(0),ID_AXES(1),ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
          IF(nr.GT.NRT) NRT=nr
        ELSE
          nr=1
        ENDIF

        IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE
          NELIST(0)=NEELEM(0,nr)
          DO noelem=1,NEELEM(0,nr)
            NELIST(noelem)=NEELEM(noelem,nr)
          ENDDO
        ENDIF

        IF(TYPE(1:5).EQ.'FIELD') THEN
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
        ELSE
          nx=1
        ENDIF !type

        IF(TYPE(1:4).EQ.'AXES') THEN
          IF(CBBREV(CO,'LENGTHS',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NT_LENGTHS,AXES_LENGTHS,ERROR,
     '        *9999)
          ELSE
            DO I=1,3
              AXES_LENGTHS(I)=0.1D0
            ENDDO
          ENDIF
          IF(CBBREV(CO,'XI_1',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NT_XI1,AXES_XI(1,1),ERROR,*9999)
          ELSE
            AXES_XI(1,1)=0.0D0
            AXES_XI(2,1)=1.0D0
            AXES_XI(3,1)=0.2D0
          ENDIF
          IF(CBBREV(CO,'XI_2',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NT_XI2,AXES_XI(1,2),ERROR,*9999)
          ELSE
            AXES_XI(1,2)=0.0D0
            AXES_XI(2,2)=1.0D0
            AXES_XI(3,2)=0.2D0
          ENDIF
          IF(CBBREV(CO,'XI_3',4,noco+1,NTCO,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NT_XI3,AXES_XI(1,3),ERROR,*9999)
          ELSE
            AXES_XI(1,3)=0.0D0
            AXES_XI(2,3)=1.0D0
            AXES_XI(3,3)=0.2D0
          ENDIF

        ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
          IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
            nm=IFROMC(CO(N3CO+1))
          ELSE
            nm=1
          ENDIF
C MHT 24-03-00 COLOUR not used, removed from parameter list
C          IF(CBBREV(CO,'GREYSCALE',1,noco+1,NTCO,N3CO)) THEN
C            COLOUR=.FALSE.
C          ELSE
C            IF(NINDICES.GE.150) THEN
C              COLOUR=.TRUE.
C            ELSE
C              COLOUR=.FALSE.
C            ENDIF
C          ENDIF
          IF(ILP(nm,1,nr,nx).LE.2) THEN      !nm uses CE array
            ZMINI=CE(nm,NEELEM(1,nr),nx)
            ZMAXI=CE(nm,NEELEM(1,nr),nx)
            DO noelem=2,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(CE(nm,ne,nx).LT.ZMINI) ZMINI=CE(nm,ne,nx)
              IF(CE(nm,ne,nx).GT.ZMAXI) ZMAXI=CE(nm,ne,nx)
            ENDDO !noelem
          ELSE IF(ILP(nm,1,nr,nx).EQ.4) THEN !nm uses CQ array
            ZMINI=CQ(nm,1,nx)
            ZMAXI=CQ(nm,1,nx)
            DO nq=2,NQT
              IF(CQ(nm,nq,nx).LT.ZMINI) ZMINI=CQ(nm,nq,nx)
              IF(CQ(nm,nq,nx).GT.ZMAXI) ZMAXI=CQ(nm,nq,nx)
            ENDDO !nq
            WRITE(*,'('' ZMINI='',E12.3)') ZMINI
            WRITE(*,'('' ZMAXI='',E12.3)') ZMAXI
          ENDIF !ILP
          ZDIFF=ZMAXI-ZMINI
        ENDIF !type

        iw=IWK(1)
        CALL ACWK(iw,1,ERROR,*9999)
        IF(TYPE(1:4).EQ.'AXES') THEN
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL SGMATE(IBT,IDO,INP,ISEG,ISMATE(iw,ne),iw,
     '        NAN,NBJ(1,ne),ne,nm,NRE(ne),AXES_LENGTHS,AXES_XI,
     '        CE(nm,ne,nx),CQ(1,1,nx),CSEG,TYPE,XE,XG,XQ,ERROR,*9999)
          ENDDO !nolist
        ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
          IF(ILP(nm,1,nr,nx).LE.2) THEN      !nm uses CE array
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              CALL SGMATE(IBT,IDO,INP,ISEG,ISMATE(iw,ne),iw,
     '          NAN,NBJ(1,ne),ne,nm,NRE(ne),AXES_LENGTHS,AXES_XI,
     '          CE(nm,ne,nx),CQ(1,1,nx),CSEG,TYPE,XE,XG,XQ,ERROR,*9999)
            ENDDO !nolist
          ELSE IF(ILP(nm,1,nr,nx).EQ.4) THEN !nm uses CQ array
            CALL SGMATE(IBT,IDO,INP,ISEG,ISMATE(iw,1),iw,
     '        NAN,NBJ(1,1),1,nm,NRE(1),AXES_LENGTHS,AXES_XI,
     '        CE(nm,1,nx),CQ(1,1,nx),CSEG,TYPE,XE,XG,XQ,ERROR,*9999)
          ENDIF !ILP
        ENDIF !type

C        IF(TYPE(1:5).EQ.'FIELD') THEN
C          IF(NTSCAL.EQ.0) NTSCAL=1
C          CALL SGSCAL(1,ISEG,ISSCAL(iw,NTSCAL),iw,NTSCAL,COLOUR,
C     '      CSEG,ERROR,*9999)
C        ENDIF
        CALL DAWK(iw,1,ERROR,*9999)
      ENDIF

      CALL EXITS('DRMATE')
      RETURN
 9999 CALL ERRORS('DRMATE',ERROR)
      CALL EXITS('DRMATE')
      RETURN 1
      END


