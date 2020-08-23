      SUBROUTINE DRGRID(ISEG,ISGRID,NAQ,NELIST,NLQ,NQET,NQLIST,NQNE,NQS,
     '  NWQ,NXLIST,NXQ,DXDXIQ,DXDXIQ2,XQ,YQ,CSEG,STRING,ERROR,*)

C#### Subroutine: DRGRID
C###  Description:
C###    DRGRID draws finite difference collocation point grid.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISGRID(NWM),NAQ(NQM,NAM),NELIST(0:NEM),
     '  NLQ(NQM),NQET(NQSCM),NQLIST(0:NQM),NQNE(NEQM,NQEM),NQS(NEQM),
     '  NWQ(8,0:NQM,NAM),NXLIST(0:NXM),NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),
     '  XQ(NJM,NQM),YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INDEX,INDEX_POLYMARKER,iw,IWK(6),
     '  na,ne,N3CO,noiw,nolist,nq,nqq,nqa,NqEnd,NqStart,
     '  nr,NTIW,nx,nxc,ILISTMBR,ngr,nqgroup
      REAL*8 RFROMC,SCALE
      LOGICAL ADAPTIVE,CBBREV,LGRMATCH
      CHARACTER TYPE*7,CHARGR1*(GR_MAXNAME),CHARGR2*(GR_MAXNAME)

      CALL ENTERS('DRGRID',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw grid
C###  Parameter:    <(points/numbers/values)[points]>
C###   Specify whether points, numbers or values are displayed.
C###  Parameter:    <(all/adaptive)[all]>
C###   Specify whether all or only adaptive points are displayed.
C###  Parameter:    <level #[1]>
C###   Specify the level of grid points.
C###  Parameter:    <rgb=RGB[black]>
C###    Specify the colour.  The options are:
C###    BLACK,RED,GREEN,BLUE,CYAN,YELLOW,WHITE,LTBLUE,GREY.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify the workstation (GX window) to draw the grid on.
C###  Parameter:    <from NQ_START#[1]>
C###   Specify the grid point to start from.
C###  Parameter:    <to NQ_END#[NQT]>
C###   Specify the grid point to end with.
C###  Parameter:    <in ELEMENT#s>
C###    Specify the elements withn which to draw the grid points in.
C###  Parameter:    <region #[1]>
C###    Draw grid points only in one region.
C###    This prevents the use of the 'from' and 'to' parameters.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:    <group GRIDGROUP>
C###    Draws only those gridpoints from a given grid group.
C###  Description:
C###    Draws finite difference collocation grids.

        OP_STRING(1)=STRING(1:IEND)
     '    //' <(points/numbers/values)[points]>'
        OP_STRING(2)=BLANK(1:15)//'<(all/adaptive)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<level #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<rgb=RGB[black]>'
        OP_STRING(5)=BLANK(1:15)//'<on (all/WS#s)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<from NQ_START#[1]>'
        OP_STRING(7)=BLANK(1:15)//'<to NQ_END#[NQT]>'
        OP_STRING(8)=BLANK(1:15)//'<in ELEMENT#s>'
        OP_STRING(9)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(10)=BLANK(1:15)//'<group GRIDGROUP>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM draw grid normals
C###  Parameter:    <scale #[1.0]>
C###   Scale the normal vectors
C###  Parameter:    <from NQ_START#[1]>
C###   Specify the grid point to start from.
C###  Parameter:    <to NQ_END#[NQT]>
C###   Specify the grid point to end with.
C###  Parameter:    <in ELEMENT#s>
C###    Specify the elements withn which to draw the grid points in.
C###  Parameter:    <region #[1]>
C###    Draw grid points only in one region.
C###    This prevents the use of the 'from' and 'to' parameters.
C###  Description:
C###    Draws normal vectors on boundary collocation points

        OP_STRING(1)=STRING(1:IEND)
     '    //' normals'
        OP_STRING(2)=BLANK(1:15)//'<scale #[1.0]>'
        OP_STRING(3)=BLANK(1:15)//'<from NQ_START#[1]>'
        OP_STRING(4)=BLANK(1:15)//'<to NQ_END#[NQT]>'
        OP_STRING(5)=BLANK(1:15)//'<in ELEMENT#s>'
        OP_STRING(6)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRGRID',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        IW=2*NJT-3

        CALL ASSERT(NQT.GT.0,'No grid points to draw',ERROR,*9999)

        IF(CBBREV(CO,'NORMALS',1,noco+1,NTCO,N3CO)) THEN
          TYPE='NORMALS'
        ENDIF

        IF(TYPE(1:7).NE.'NORMALS') THEN
          IF(CBBREV(CO,'POINTS',1,noco+1,NTCO,N3CO)) THEN
            TYPE='POINTS'
          ELSE IF(CBBREV(CO,'NUMBERS',1,noco+1,NTCO,N3CO)) THEN
            TYPE='NUMBERS'
          ELSE IF(CBBREV(CO,'VALUES',1,noco+1,NTCO,N3CO)) THEN
            TYPE='VALUES'
          ELSE
            TYPE='POINTS'
          ENDIF

          IF(TYPE(1:6).EQ.'VALUES') THEN
            CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
            nxc=NXLIST(1)
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '        ERROR,*9999)
          ELSE
            nx=1 !since an nx is required in call sggrid
          ENDIF !values
        ENDIF !normals

        IF(CBBREV(CO,'ADAPTIVE',1,noco+1,NTCO,N3CO)) THEN
          ADAPTIVE=.TRUE.
        ELSE
          ADAPTIVE=.FALSE.
        ENDIF !adaptive

        IF(CBBREV(CO,'REGION',2,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
          NqStart=NQR(1,nr)
          NqEnd=NQR(2,nr)
        ELSE
          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            NqStart=IFROMC(CO(N3CO+1))
          ELSE
            NqStart=1
          ENDIF
          IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
            NqEnd=IFROMC(CO(N3CO+1))
          ELSE
            NqEnd=NQT
          ENDIF
        ENDIF !region

        IF(TYPE(1:7).NE.'NORMALS') THEN
          IF(CBBREV(CO,'LEVEL',1,noco+1,NTCO,N3CO)) THEN
            na=IFROMC(CO(N3CO+1))
          ELSE
            na=1
          ENDIF !levels
        ELSE
          IF(CBBREV(CO,'SCALE',3,noco+1,NTCO,N3CO)) THEN
            SCALE=RFROMC(CO(N3CO+1))
          ELSE
            SCALE=1.0d0
          ENDIF !levels
          na=1
        ENDIF !normals

        IF(CBBREV(CO,'RGB',3,noco+1,NTCO,N3CO)) THEN
          INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1',CO(N3CO+1))
        ELSE
          INDEX=INDEX_POLYMARKER(0,'PLUS','SIZE1','BLACK')
        ENDIF !rgb

        IF(CBBREV(CO,'IN',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NE_R_M,NELIST(0),NELIST(1),ERROR,*9999)
          nqa=0
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            DO nqq=1,NQET(NQS(ne))
              nq=NQNE(ne,nqq)
              IF(ADAPTIVE) THEN !adaptive grid pts only
                IF(NLQ(nq).NE.0) THEN
                  nqa=nqa+1
                  NQLIST(nqa)=nq
                ENDIF !NLQ
              ELSE !all points
                nqa=nqa+1
                NQLIST(nqa)=nq
              ENDIF !adaptive/all
            ENDDO !nqq
            NQLIST(0)=nqa
          ENDDO !ne
        ELSE IF (CBBREV(CO,'GROUPS',2,noco+1,NTCO,N3CO))  THEN !Draw groups
          !Get the group name the user wanted
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHARGR1)

          !Try and match the group name
          LGRMATCH=.FALSE.
          DO nqgroup=1,NTGRGR !Groups names
            CALL STRING_TRIM(LAGRGR(nqgroup),IBEG,IEND)
            CALL CUPPER(LAGRGR(nqgroup)(IBEG:IEND),CHARGR2)
            IF (CHARGR1.EQ.CHARGR2) THEN
              LGRMATCH=.TRUE.
              ngr=nqgroup !Store grid group number
            ENDIF
          ENDDO
          CALL ASSERT(LGRMATCH,'>>No grid groups of that name',
     '      ERROR,*9999)
          !Copy array across
          NQLIST(0)=NLIGRGR(ngr)
          DO nq=1,NLIGRGR(ngr)
            NQLIST(nq)=ILISTMBR(%VAL(LIGRGR_PTR(ngr)),nq)
          ENDDO
        ELSE IF(na.EQ.1) THEN !fine grid level na=1
          nqa=0
          DO nq=NqStart,NqEnd
            IF(ADAPTIVE) THEN !adaptive grid pts only
              IF(NLQ(nq).NE.0) THEN
                nqa=nqa+1
                NQLIST(nqa)=nq
              ENDIF !NLQ
            ELSE !all points
              nqa=nqa+1
              NQLIST(nqa)=nq
            ENDIF !adaptive/all
          ENDDO !nq
          NQLIST(0)=nqa
        ELSE !coarse grid level na
          nqa=0
          DO nq=1,NQT
            IF(NAQ(nq,na).EQ.0) THEN !nq is in grid na
              IF(ADAPTIVE) THEN !adaptive grid pts only
                IF(NLQ(nq).NE.0) THEN
                  nqa=nqa+1
                  NQLIST(nqa)=nq
                ENDIF !NLQ
              ELSE !all points
                nqa=nqa+1
                NQLIST(nqa)=nq
              ENDIF !adaptive/all
            ENDIF !NAQ
          ENDDO !nq
          NQLIST(0)=nqa
        ENDIF !na/in element

        IF(DOP) THEN
          WRITE(OP_STRING,'('' NQLIST:'',/(10I8))')
     '      (NQLIST(nqa),nqa=1,NQLIST(0))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !DOP

        IF(TYPE(1:6).EQ.'VALUES') THEN !find range
          IF(CBBREV(CO,'ZMINI',3,noco+1,NTCO,N3CO)) THEN
            ZMINI=IFROMC(CO(N3CO+1))
          ELSE
            ZMINI=YQ(1,1,na,nx)
            DO nq=2,NQT
              IF(YQ(nq,1,na,nx).LT.ZMINI) ZMINI=YQ(nq,1,na,nx)
            ENDDO !nq
          ENDIF !ZMINI
          IF(CBBREV(CO,'ZMAXI',3,noco+1,NTCO,N3CO)) THEN
            ZMAXI=IFROMC(CO(N3CO+1))
          ELSE
            ZMAXI=YQ(1,1,na,nx)
            DO nq=2,NQT
              IF(YQ(nq,1,na,nx).GT.ZMAXI) ZMAXI=YQ(nq,1,na,nx)
            ENDDO !nq
          ENDIF !ZMAXI
          ZDIFF=ZMAXI-ZMINI
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ZMINI='',E12.3,'' ZMAXI='',E12.3)')
     '        ZMINI,ZMAXI
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF !dop
        ENDIF !values

        IF(TYPE(1:7).EQ.'NORMALS') THEN
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)
            CALL SGGRID_NORM(ISEG,iw,NQLIST,NWQ(1,0,1),
     '        NXQ(-NIM,0,0,1),DXDXIQ,DXDXIQ2,SCALE,XQ,CSEG,ERROR,*9999)
            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO !noiw
        ELSE
          DO noiw=1,NTIW
            IW=IWK(noiw)
            CALL ACWK(iw,1,ERROR,*9999)

            CALL SGGRID(INDEX,ISEG,ISGRID(iw),iw,NLQ,NQLIST,
     '        CSEG,ADAPTIVE,TYPE,XQ,YQ(1,1,na,nx),ERROR,*9999)

            CALL DAWK(iw,1,ERROR,*9999)
          ENDDO !noiw
        ENDIF
      ENDIF !?

      CALL EXITS('DRGRID')
      RETURN
 9999 CALL ERRORS('DRGRID',ERROR)
      CALL EXITS('DRGRID')
      RETURN 1
      END


