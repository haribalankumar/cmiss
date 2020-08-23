      SUBROUTINE EVFIEL(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NELIST,
     '  NGLIST,NHE,NHP,NKH,NKHE,NKJE,NLL,NPF,NPL,NPNE,NPNODE,NRLIST,
     '  NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '  CURVCORRECT,DL,SE,XA,XE,XG,XIG,XP,YP,ZA,ZE,ZP,STRING,ERROR,*)

C#### Subroutine: EVFIEL
C###  Description:
C###    EVFIEL evaluates geometric, field or dependent variables
C###    at specified Xi increments or Gauss points.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NGLIST(0:NGM),NHE(NEM,NXM),NHP(NPM,0:NRM,NXM),
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),
     '  SE(NSM,NBFM,NEM),XA(NAM,NJM),XE(NSM,NJM),
     '  XG(NJM,NUM),
     '  XIG(NIM,NGM,NBM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER eval_nr,i,IBEG,IEND,LD,loop_num,nae,
     '  N3CO,nb,nbb,nc,ndata,ne,ne_eval,ne_last,ng1,ng2,nh,nhx,
     '  ni,NITB,nj,njj1,njj2,nl,noelem,
     '  no_nelist,no_nrlist,NT_dXi,NT_dX,nr,nu(1),
     '  NUM_LOOPS,nx,nxc,NUMX(3),nXi(3),TMP(1)
      REAL*8 dX(3),dXi(3),PXI,Xi(3),Xi_last(3),
     '  X(20),XE_eval(NSM,NJM),XRC(20),Z(6),
     '  A(3),B(3),C(3) ,LILEN(3)
      CHARACTER FILE*(MXCH),NUMDEP*1,NUMFIBR*1,NUMFIEL*2,
     '  NUMGEOFIEL*2,NUMGEOM*1,OUTPUTFORMAT*7,TYPE*9,
     '  DATASTRING*500
      LOGICAL ALL_REGIONS,CBBREV,FOUND,
     '  FROM_REGION,GAUSS,OPFILE,USEDX

      CALL ENTERS('EVFIEL',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate field<;FILENAME>
C###  Parameter:        <(geometry/field/dependent)[geometry]>
C###    Specify which field to evaluate
C###  Parameter:        <element (#s/all)[all]>
C###    Specify in which elements the field is to be evaluated.
C###    The element can be a specified by a single number
C###    or a comma separated list, by default the field is
C###    evaluated in all elements.
C###  Parameter:        <(xi STEP#s/gauss (GAUSS_PT#s/all))[xi 0.2,0.2,0.2]>
C###    Specify the increments of xi at which the field is to
C###    be evaluated at, or the Gauss point numbers
C###  Parameter:        <nu #[1]>
C###    <HTML>
C###    Specify the derivative number of the field to evaluate.
C###    <PRE>
C###    nu=1 function
C###    nu=2 derivative wrt Xi1
C###    nu=3 derivative wrt Xi1^2
C###    nu=4 derivative wrt Xi2
C###    nu=5 derivative wrt Xi2^2
C###    nu=6 derivative wrt Xi1 Xi2
C###    nu=7 derivative wrt Xi3
C###    nu=8 derivative wrt Xi3^2
C###    nu=9 derivative wrt Xi1 Xi3
C###    nu=10 derivative wrt Xi2 Xi3
C###    nu=11 derivative wrt Xi1 XI2 Xi3
C###    </PRE>
C###    </HTML>
C###  Parameter:        <(full/reduced/data)[full]>
C###    Specify full, a reduced general format output or a 
C###    ipdata file output.
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers
C###  Parameter:        <class #[1]>
C###    Specify the problem class number
C###  Parameter:        <dx [0,0,0]>
C###  Parameter:        <from_region (#)[1]>
C###    Specify the region to be sampled. The xi locations can
C###    be specifid by a mesh in one region while the fields
C###    to sampled can be from a mesh in another
C###  Description:
C###

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
     '    //' <geometry/field/dependent>[geometry]'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(xi STEP#s/gauss (GAUSS_PT#s/all))'
     '   //'[xi 0.2,0.2,0.2]>'
        OP_STRING(4)=BLANK(1:15)//'<nu #[1]>'
        OP_STRING(5)=BLANK(1:15)//'<(full/reduced/data)[full]>'
        OP_STRING(6)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(7)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(8)=BLANK(1:15)//'<dx [0,0,0]>'
        OP_STRING(9)=BLANK(1:15)//'<from_region (#)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVFIEL',ERROR,*9999)
      ELSE

        nc=1 !temporary
        ndata = 0
C DMAL 06 JUNE 2003 : Adding ability to output to ipdata file
        IF(CBBREV(CO,'FULL',2,noco+1,NTCO,N3CO)) THEN
          OUTPUTFORMAT='FULL'
        ELSE IF(CBBREV(CO,'REDUCED',2,noco+1,NTCO,N3CO)) THEN
          OUTPUTFORMAT='REDUCED'
        ELSE IF(CBBREV(CO,'DATA',3,noco+1,NTCO,N3CO)) THEN
          OUTPUTFORMAT='DATA'
        ELSE
          OUTPUTFORMAT='FULL'
        ENDIF

        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          IF (OUTPUTFORMAT(1:4).EQ.'DATA') THEN
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.ipdata','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
            WRITE(OP_STRING,'(''Data File'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opfiel','NEW',
     '        'SEQUEN','FORMATTED',132,ERROR,*9999)
          ENDIF
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'GAUSS',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NGM,NGLIST(0),NGLIST(1),ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'NU',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),1,i,nu,ERROR,*9999)
        ELSE
          nu(1)=1
        ENDIF

        IF(CBBREV(CO,'FROM_REGION',4,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),1,i,TMP,ERROR,*9999)
          eval_nr=TMP(1)
          FROM_REGION=.TRUE.
        ELSE
          FROM_REGION=.FALSE.
          FOUND=.FALSE.
        ENDIF

C       Set up character strings for use in format statements
        WRITE(NUMGEOM,'(I1)') NJ_LOC(NJL_GEOM,0,0)
        WRITE(NUMFIBR,'(I1)') NJ_LOC(NJL_FIBR,0,0)
        WRITE(NUMFIEL,'(I2)') NJ_LOC(NJL_FIEL,0,0)
        IF(CBBREV(CO,'GEOMETRIC',2,noco+1,NTCO,N3CO)) THEN
          TYPE='GEOMETRIC'
          CALL ASSERT(NJ_LOC(NJL_GEOM,0,0)+NJ_LOC(NJL_FIBR,0,0).LE.20,
     '      '>>Increase size of X and XRC arrays',ERROR,*9999)
        ELSE IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
          TYPE='FIELD'
          CALL ASSERT(NJ_LOC(NJL_GEOM,0,0)+NJ_LOC(NJL_FIEL,0,0).LE.20,
     '      '>>Increase size of X and XRC arrays',ERROR,*9999)
        ELSE IF(CBBREV(CO,'DEPENDENT',2,noco+1,NTCO,N3CO)) THEN
          TYPE='DEPENDENT'
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)
          CALL ASSERT(NHE(1,nx).LE.6,
     '      '>>Increase size of Z array',ERROR,*9999)
          WRITE(NUMDEP,'(I1)') NHE(1,nx)
        ELSE
          TYPE='GEOMETRIC'
          CALL ASSERT(NJ_LOC(NJL_GEOM,0,0)+NJ_LOC(NJL_FIBR,0,0).LE.12,
     '      '>>Increase size of X and XRC arrays',ERROR,*9999)
        ENDIF

        GAUSS=.FALSE.
        NITB=NIT(NBJ(1,NEELEM(1,NRLIST(1))))
        IF(CBBREV(CO,'XI',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),NITB,NT_dXi,dXi,ERROR,*9999)
          DO i=1,NT_dXi
            CALL ASSERT(dXi(i).GT.0.d0.AND.dXi(i).LE.1.d0,
     '        '>> must specify 0 < dXi <= 1',ERROR,*9999)
            NUMX(i)=INT(1.d0/dXi(i))
          ENDDO !i
          DO i=NT_dXi+1,3
            NUMX(i)=1
          ENDDO !i
        ELSE IF(CBBREV(CO,'GAUSS',1,noco+1,NTCO,N3CO)) THEN
          GAUSS=.TRUE.
        ELSE
          DO i=1,3
            dXi(i)=0.20d0
            NUMX(i)=1+INT(1.d0/dXi(i))
          ENDDO !i
        ENDIF

        IF(CBBREV(CO,'DX',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),NITB,NT_dX,dX,ERROR,*9999)
          USEDX=.TRUE.
        ELSE
          USEDX=.FALSE.
        ENDIF

        IF(NITB.LE.2) THEN
          dXi(3)=1.1d0 !initialize so only one step in Xi(3) if 1D or 2D
          dX(3)=0.d0 !initialize
        ENDIF
        IF(NITB.EQ.1) THEN
          dXi(2)=1.1d0 !initialize so only one step in Xi(2) if 1D
          dX(2)=0.d0 !initialize
        ENDIF

        IF(TYPE(1:9).EQ.'GEOMETRIC') THEN
          IF(OUTPUTFORMAT(1:7).EQ.'REDUCED') THEN
            WRITE(OP_STRING,'(5X,''x'',9X,''y'',9X,'
     '        //'''z       Fibre vector:'',12X,''Sheet vector:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !full

        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)

          WRITE(NUMGEOFIEL,'(I2)') NJ_LOC(NJL_GEOM,0,0)+
     '      NJ_LOC(NJL_FIEL,0,nr)

          IF(TYPE(1:9).EQ.'DEPENDENT') THEN
            CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '        NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,
     '        ERROR,*9999)
          ENDIF

          DO no_nelist=1,NELIST(0)
            ne=NELIST(no_nelist)

            IF(USEDX) THEN
              nbb=NBJ(1,NEELEM(1,NRLIST(1)))
              DO i=1,3
                LILEN(i)=0.0d0
              ENDDO
C LKC 11-OCT-97 Used before set, changed NLE(nb) -> NLE(nbb)
              DO nae=1,NLE(nbb)
                nl=NLL(nae,ne)
                IF(nl.NE.0) THEN
                  IF(DL(3,nl).GT.LILEN(NPL(1,0,nl)))
     '              LILEN(NPL(1,0,nl))=DL(3,nl)
                ENDIF
              ENDDO
              DO i=1,3
                IF(dX(i).GT.ZERO_TOL) THEN
C new MPN 22Sep2006: round rather than truncate
                  NUMX(i)=INT(LILEN(i)/dX(i)+0.5d0)
C old                  NUMX(i)=INT(LILEN(i)/dX(i))
                ELSE
                  NUMX(i)=1
                ENDIF
                IF(NUMX(i).GT.0) THEN
                  dXi(i)=1.0d0/DBLE(NUMX(i))
                ELSE
                  dXi(i)=1.0d0
                ENDIF
              ENDDO
            ENDIF

            IF((OUTPUTFORMAT(1:4).EQ.'FULL'.OR.
     '        TYPE(1:9).NE.'GEOMETRIC').AND.
     '        (OUTPUTFORMAT(1:4).NE.'DATA').AND.
     '        .NOT.FROM_REGION) THEN
              WRITE(OP_STRING,'(/'' Element '',I3)') ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA,XE,XP,ERROR,*9999)
            IF(TYPE(1:9).EQ.'DEPENDENT') THEN
              CALL ZPZE(NBH(1,1,ne),nc,NHE(ne,nx),NKHE(1,1,1,ne),
     '          NPF(1,1),NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '          CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '          ERROR,*9999)
            ENDIF

C CS 11/6/97 Changed to add the ability to evaluate fields
C at Gauss points
C SEN 28/01/03 Removed floating do loops.
            IF(.NOT.GAUSS) THEN
C             determine the number of points we evaluate the field at
C MPN 4Apr2003: fixed NUM_LOOPS for varying number of coords
C               and whether Xi or dX increments were specified
              NUM_LOOPS=1
              DO ni=1,NITB
                NUM_LOOPS=NUM_LOOPS*(1+NUMX(ni))
              ENDDO !ni
C old              NUM_LOOPS=(1+NUMX(1))*(1+NUMX(2))*(1+NUMX(3))
            ELSE
              IF(NGLIST(0).GE.1) THEN
                NUM_LOOPS=NGLIST(0)
              ELSE
                ng1=0
                IF(TYPE(1:9).EQ.'GEOMETRIC') THEN
                  DO njj1=1,2    !geom + fibres + sheets
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      nb=NBJ(nj,ne)
                      ng2=NGT(nb)
                      IF(ng1.EQ.0) THEN
                        ng1=ng2
                      ENDIF
                      CALL ASSERT(ng1.EQ.ng2,
     '                  '>>Number of Gauss points for fields differs',
     '                  ERROR,*9999)
                    ENDDO
                  ENDDO
                ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
                  DO njj1=1,3,2    !geom + field
                    DO njj2=1,NJ_LOC(njj1,0,nr)
                      nj=NJ_LOC(njj1,njj2,nr)
                      nb=NBJ(nj,ne)
                      ng2=NGT(nb)
                      IF(ng1.EQ.0) THEN
                        ng1=ng2
                      ENDIF
                      CALL ASSERT(ng1.EQ.ng2,
     '                  '>>Number of Gauss points for fields differs',
     '                  ERROR,*9999)
                    ENDDO
                  ENDDO
                ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                  DO nhx=1,NHE(ne,nx)
                    nh=NH_LOC(nhx,nx)
                    nb=NBH(nh,1,ne)
                    ng2=NGT(nb)
                    IF(ng1.EQ.0) THEN
                      ng1=ng2
                    ENDIF
                    CALL ASSERT(ng1.EQ.ng2,
     '                '>>Number of Gauss points for fields differs',
     '                ERROR,*9999)
                  ENDDO !nhx (nh)
                ENDIF
                NUM_LOOPS=ng1
              ENDIF
            ENDIF

            IF(FROM_REGION) THEN
              ne_last=NEELEM(1,eval_nr)
              Xi_last(1)=0.5d0
              Xi_last(2)=0.5d0
              Xi_last(3)=0.5d0
            ENDIF

            nXi(1)=0
            nXi(2)=0
            nXi(3)=0
            DO loop_num=1,NUM_LOOPS
              IF(.NOT.GAUSS) THEN
C!!! CS 24/9/2003 If input xi increments, expect them to be used
C rather than the even spacing recalculated here
C                Xi(1)=DBLE(nXi(1))/DBLE(NUMX(1))
C                Xi(2)=DBLE(nXi(2))/DBLE(NUMX(2))
C                Xi(3)=DBLE(nXi(3))/DBLE(NUMX(3))
                Xi(1)=nXi(1)*dXi(1)
                Xi(2)=nXi(2)*dXi(2)
                Xi(3)=nXi(3)*dXi(3)
                nXi(1)=nXi(1)+1
                IF(nXi(1).GT.NUMX(1)) THEN
                  nXi(1)=0
                  nXi(2)=nXi(2)+1
                  IF(nXi(2).GT.NUMX(2)) THEN
                    nXi(2)=0
                    nXi(3)=nXi(3)+1
                    IF(nXi(3).GT.NUMX(3)) THEN
                      nXi(3)=0
C??? CS 24/9/2003 Don't see why this is necessary, infact it is 
C??? always encountered on the last loop, so removing
C???                      ERROR='>> nXi(3) > NUMX(3)'
C???                      GOTO 9999
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF

              IF(TYPE(1:9).EQ.'GEOMETRIC') THEN
                DO njj1=1,2    !geom + fibres + sheets
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    nb=NBJ(nj,ne)
                    IF(GAUSS) THEN
                      IF(NGLIST(0).GE.1) THEN
                        Xi(1)=XIG(1,NGLIST(loop_num),nb)
                        Xi(2)=XIG(2,NGLIST(loop_num),nb)
                        Xi(3)=XIG(3,NGLIST(loop_num),nb)
                      ELSE
                        Xi(1)=XIG(1,loop_num,nb)
                        Xi(2)=XIG(2,loop_num,nb)
                        Xi(3)=XIG(3,loop_num,nb)
                      ENDIF
                    ENDIF
                    IF(nb.GT.0) THEN
                      X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,nu(1),Xi,XE(1,nj))
                    ENDIF
                  ENDDO
                ENDDO

                IF(FROM_REGION) THEN
                  FOUND=.FALSE.
                  noelem=-1
                  LD=0
                  DO WHILE (.NOT.FOUND
     '              .AND.noelem.LT.NEELEM(0,eval_nr))
                    noelem=noelem+1
                    IF(noelem.EQ.0) THEN
                      ne_eval=ne_last
                    ELSE
                      ne_eval=NEELEM(noelem,eval_nr)
                    ENDIF

                    CALL XPXE(NBJ(1,ne_eval),NKJE(1,1,1,ne_eval),
     '                NPF(1,1),NPNE(1,1,ne_eval),eval_nr,
     '                NVJE(1,1,1,ne_eval),SE(1,1,ne_eval),XA,
     '                XE_eval,XP,ERROR,*9999)
                    NITB=NIT(NBJ(1,ne_eval))

                    IF(ne_eval.EQ.ne_last) THEN
                      Xi(1)=Xi_last(1)
                      Xi(2)=Xi_last(2)
                      Xi(3)=Xi_last(3)
                    ELSE
                      Xi(1)=0.5d0
                      Xi(2)=0.5d0
                      Xi(3)=0.5d0
                    ENDIF
                    CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne_eval,NITB,
     '                eval_nr,0.d0,XE_eval,
     '                XI,XI,X,.FALSE.,ERROR,*9999)

                    IF(LD.NE.0) THEN
                      FOUND=.TRUE.
                      Xi_last(1)=Xi(1)
                      Xi_last(2)=Xi(2)
                      Xi_last(3)=Xi(3)
                      ne_last=LD
                    ENDIF
                  ENDDO

                 DO njj1=1,2    !geom + fibres + sheets
                  DO njj2=1,NJ_LOC(njj1,0,eval_nr)
                    nj=NJ_LOC(njj1,njj2,eval_nr)
                    nb=NBJ(nj,ne_eval)
                      X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,nu(1),Xi,XE_eval(1,nj))
                  ENDDO
                 ENDDO

                ENDIF

                IF(.NOT.FROM_REGION.OR.FOUND) THEN

                IF(OUTPUTFORMAT(1:4).EQ.'FULL') THEN
                  IF(GAUSS) THEN
                    IF(NGLIST(0).GE.1) THEN
                      WRITE(OP_STRING,'(''ng= '',I2,'' Xi:'',
     '                  3F5.2,''  X:'','
     '                  //NUMGEOM//'D12.4,'//NUMFIBR//'F7.3)')
     '                  (NGLIST(loop_num),Xi(ni),ni=1,3),
     '                  ((X(NJ_LOC(njj1,njj2,nr)),
     '                  njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,2)
                    ELSE
                      WRITE(OP_STRING,'(''ng= '',I2,'' Xi:'',
     '                  3F5.2,''  X:'','
     '                  //NUMGEOM//'D12.4,'//NUMFIBR//'F7.3)')
     '                  (loop_num,Xi(ni),ni=1,3),
     '                  ((X(NJ_LOC(njj1,njj2,nr)),
     '                  njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,2)
                    ENDIF
                  ELSE
                    IF(FROM_REGION) THEN
                      WRITE(OP_STRING,'('' Xi:'',3F5.2,''  X:'','
     '                  //NUMGEOM//'D12.4,'//NUMFIBR//'F7.3)')
     '                  (Xi(ni),ni=1,3),
     '                  ((X(NJ_LOC(njj1,njj2,eval_nr)),
     '                  njj2=1,NJ_LOC(njj1,0,eval_nr)),njj1=1,2)
                    ELSE
                      WRITE(OP_STRING,'('' Xi:'',3F5.2,''  X:'','
     '                  //NUMGEOM//'D12.4,'//NUMFIBR//'F7.3)')
     '                  (Xi(ni),ni=1,3),
     '                  ((X(NJ_LOC(njj1,njj2,nr)),
     '                  njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,2)
                    ENDIF
                  ENDIF
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
                !print r.c. coords
                IF(FROM_REGION) THEN
                  IF (ITYP10(eval_nr).GT.1) THEN
                    CALL XZ(ITYP10(eval_nr),X,XRC)
                    IF(OUTPUTFORMAT(1:4).EQ.'FULL') THEN
                      WRITE(OP_STRING,'(11X,''R.C. coords:'','
     '                  //NUMGEOM//'D12.4)')
     '                  (XRC(NJ_LOC(NJL_GEOM,njj2,eval_nr)),
     '                njj2=1,NJ_LOC(NJL_GEOM,0,eval_nr))
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDIF
                  ENDIF
                ELSE IF(ITYP10(nr).GT.1) THEN
                  CALL XZ(ITYP10(nr),X,XRC)
                  IF(OUTPUTFORMAT(1:4).EQ.'FULL') THEN
                    WRITE(OP_STRING,'(11X,''R.C. coords:'','
     '                //NUMGEOM//'D12.4)')
     '                (XRC(NJ_LOC(NJL_GEOM,njj2,nr)),
     '              njj2=1,NJ_LOC(NJL_GEOM,0,nr))
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                IF((NJ_LOC(NJL_FIBR,0,nr).EQ.3).OR.FROM_REGION) THEN
C                 print fibres and sheets
                  IF(FROM_REGION) THEN
                    CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne_eval)
     '                ,eval_nr,A,B,C,XE_eval,XG,Xi,.TRUE.,
     '                ERROR,*9999)
                  ELSE
                    CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne)
     '                ,nr,A,B,C,XE,XG,Xi,.TRUE.,ERROR,*9999)
                  ENDIF
                  IF(OUTPUTFORMAT(1:4).EQ.'FULL') THEN
                    WRITE(OP_STRING,'('' A:'',3F7.3,'
     '               //'''  B:'',3F7.3,''  C:'',3F7.3)')
     '               (A(i),i=1,3),(B(i),i=1,3),(C(i),i=1,3)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ELSEIF(OUTPUTFORMAT(1:4).EQ.'DATA') THEN
                    IF(FROM_REGION) THEN
                      CALL ASSERT(.FALSE.,'>>Not implemented',
     '                  ERROR,*9999)
                    ELSE
                      CALL ASSERT(.FALSE.,'>>Not implemented',
     '                  ERROR,*9999)
                    ENDIF
                  ELSE !REDUCED output
C new MPN 22Sep2006: XRC was all zero for 'reduced' output - need to CALL XZ
                    CALL XZ(ITYP10(nr),X,XRC)
C end new
                    IF(FROM_REGION) THEN
                      WRITE(OP_STRING,'('//NUMGEOM//'F10.4,1X,3F8.4,'
     '                  //'1X,3F8.4)')
     '                  (XRC(NJ_LOC(NJL_GEOM,njj2,eval_nr)),
     '                  njj2=1,NJ_LOC(NJL_GEOM,0,eval_nr)),
     '                  (A(i),i=1,3),(B(i),i=1,3)
                    ELSE
                      WRITE(OP_STRING,'('//NUMGEOM//'F10.4,1X,3F8.4,'
     '                  //'1X,3F8.4)')
     '                  (XRC(NJ_LOC(NJL_GEOM,njj2,nr)),
     '                  njj2=1,NJ_LOC(NJL_GEOM,0,nr)),
     '                  (A(i),i=1,3),(B(i),i=1,3)
                    ENDIF
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF !full
                ENDIF

                ENDIF

              ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
                DO njj1=1,3,2    !geom + field
                  DO njj2=1,NJ_LOC(njj1,0,nr)
                    nj=NJ_LOC(njj1,njj2,nr)
                    nb=NBJ(nj,ne)
                    IF(GAUSS) THEN
                      IF(NGLIST(0).GE.1) THEN
                        Xi(1)=XIG(1,NGLIST(loop_num),nb)
                        Xi(2)=XIG(2,NGLIST(loop_num),nb)
                        Xi(3)=XIG(3,NGLIST(loop_num),nb)
                      ELSE
                        Xi(1)=XIG(1,loop_num,nb)
                        Xi(2)=XIG(2,loop_num,nb)
                        Xi(3)=XIG(3,loop_num,nb)
                      ENDIF
                    ENDIF
                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,nu(1),Xi,XE(1,nj))
                  ENDDO
                ENDDO
                IF(OUTPUTFORMAT(1:4).EQ.'FULL') THEN
                  WRITE(OP_STRING,'('' Xi:'',3F5.2,'
     '              //'''  Ref coords:'','//NUMGEOM//'D12.4,'
     '              //'''  Field var(s):'','//NUMFIEL//'D12.4)')
     '              (Xi(ni),ni=1,3),
     '              (X(NJ_LOC(NJL_GEOM,njj2,nr)),
     '              njj2=1,NJ_LOC(NJL_GEOM,0,nr)),
     '              (X(NJ_LOC(NJL_FIEL,njj2,nr)),
     '              njj2=1,NJ_LOC(NJL_FIEL,0,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  IF(ITYP10(nr).GT.1) THEN !print r.c. coords
                    CALL XZ(ITYP10(nr),X,XRC)
                    WRITE(OP_STRING,'(11X,''R.C. coords:'','
     '                //NUMGEOM//'D12.4)')
     '                (XRC(NJ_LOC(NJL_GEOM,njj2,nr)),
     '                njj2=1,NJ_LOC(NJL_GEOM,0,nr))
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSEIF(OUTPUTFORMAT(1:4).EQ.'DATA') THEN
                  ndata=ndata+1
                  WRITE(DATASTRING,'('//'I8,'//NUMGEOFIEL//'D12.4,'
     '              //NUMGEOFIEL//'F5.1)') ndata,
     '              ((X(NJ_LOC(njj1,njj2,nr)),
     '              njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3,2),
     '              ((1.0D0,njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3,2)
                  CALL STRING_TRIM(DATASTRING,IBEG,IEND)
                  CALL WRITE_LINE(IOFI,DATASTRING(IBEG:IEND),
     '              ERROR,*9999)
                ELSE !REDUCED output
                  WRITE(OP_STRING,'('' Xi:'',3F6.3,'
     '              //'''  X:'','//NUMGEOFIEL//'D12.4)')
     '              (Xi(ni),ni=1,3),((X(NJ_LOC(njj1,njj2,nr)),
     '              njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,3,2)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF

              ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
C MPN 15Apr2003: rewritten to correct logic!
                DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj2,nr)
                  nb=NBJ(nj,ne)
                  IF(GAUSS) THEN
                    IF(NGLIST(0).GE.1) THEN
                      Xi(1)=XIG(1,NGLIST(loop_num),nb)
                      Xi(2)=XIG(2,NGLIST(loop_num),nb)
                      Xi(3)=XIG(3,NGLIST(loop_num),nb)
                    ELSE
                      Xi(1)=XIG(1,loop_num,nb)
                      Xi(2)=XIG(2,loop_num,nb)
                      Xi(3)=XIG(3,loop_num,nb)
                    ENDIF
                  ENDIF
                  IF(nb.GT.0) THEN
                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),nb,nu(1),Xi,XE(1,nj))
                  ENDIF
                ENDDO !njj2
                DO nhx=1,NHE(ne,nx)
                  nh=NH_LOC(nhx,nx)
                  nb=NBH(nh,1,ne)
                  Z(nh)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,nu(1),Xi,ZE(1,nhx))
                ENDDO !nhx (nh)
                IF(OUTPUTFORMAT(1:4).EQ.'FULL') THEN
                  WRITE(OP_STRING,'(/'' Xi:'',3F5.2,'
     '              //'''  Dep var(s):'','//NUMDEP//'D14.6)')
     '              (Xi(ni),ni=1,3),
     '              (Z(NH_LOC(nhx,nx)),nhx=1,NHE(ne,nx))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  IF(ITYP10(nr).EQ.1) THEN !print r.c. coords only
                    WRITE(OP_STRING,'('' Ref coords:'','
     '                //NUMGEOM//'D12.4)')
     '                (X(NJ_LOC(NJL_GEOM,njj2,nr)),
     '                njj2=1,NJ_LOC(NJL_GEOM,0,nr))
                        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ELSE !ref coords not RC so print r.c. coords also
                    CALL XZ(ITYP10(nr),X,XRC)
                    WRITE(OP_STRING,'('' Ref coords:'','
     '                //NUMGEOM//'D12.4,'//'''  R.C. coords:'','
     '                //NUMGEOM//'D12.4)')
     '                (X(NJ_LOC(NJL_GEOM,njj2,nr)),
     '                njj2=1,NJ_LOC(NJL_GEOM,0,nr)),
     '                (XRC(NJ_LOC(NJL_GEOM,njj2,nr)),
     '                njj2=1,NJ_LOC(NJL_GEOM,0,nr))
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ELSEIF(OUTPUTFORMAT(1:4).EQ.'DATA') THEN
                   CALL ASSERT(.FALSE.,'>>Not implemented',ERROR,*9999)
                ELSE ! not FULL output
                  WRITE(OP_STRING,'('' Xi: '',3F5.2,'
     '              //'''  Z:'','//NUMDEP//'D12.4)')
     '              (Xi(ni),ni=1,3),
     '              (Z(NH_LOC(nhx,nx)),nhx=1,NHE(ne,nx))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
C                ENDDO !Xi(1)
C              ENDDO !Xi(2)
C            ENDDO !Xi(3)
            ENDDO !loop_num
          ENDDO !no_nelist
        ENDDO !no_nrlist

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('EVFIEL')
      RETURN
 9999 CALL ERRORS('EVFIEL',ERROR)
      CALL EXITS('EVFIEL')
      RETURN 1
      END


