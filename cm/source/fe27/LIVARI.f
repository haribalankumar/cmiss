      SUBROUTINE LIVARI(NBH,NEELEM,NENQ,NHE,NHP,
     '  NKH,NPNODE,NP_INTERFACE,NQNP,NQS,NQXI,NRLIST,
     '  NVHP,NXLIST,NXQ,NYNE,NYNP,NYNR,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '  XQ,YP,YQ,ZA,ZP,STRING,FIX,FIXQ,ERROR,*)

C#### Subroutine: LIVARI
C###  Description:
C###    LIVARI lists variables.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),NQNP(NPM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NXLIST(0:NXM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM,NXM),DNUDXQ(3,3,NQM),
     &  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),XQ(NJM,NQM),YP(NYM,NIYM,NXM),
     &  YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),
     &  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM),FIXQ(NYQM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IY,N3CO,nc,NCLIST(0:4),no_nc,nolist,
     '  nr,nrb,nrg,nx,nxb,nxc,nxg,nx_upd
      CHARACTER FILE*100,TYPE*4
      LOGICAL ABBREV,ALL_NCs,ALL_REGIONS,CBBREV,OPFILE

      CALL ENTERS('LIVARI',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list variables<;FILENAME>
C###  Parameter:    <number IY#[1]{>0}>
C###  Parameter:    <array (yp/zp)[yp]>
C###  Parameter:    <nc (#s/all)[all]>
C###  Parameter:    <region (#s/all)[1]>
C###  Parameter:    <using (solve/fit/optimise)[solve]>
C###  Parameter:    <class #[1]>
C###  Description:
C###    Lists specified dependent variables to screen or file
C###    FILENAME.opvari if qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> <number IY#[1]{>0}>'
        OP_STRING(2)=BLANK(1:15)//'<array (yp/zp)[yp]>'
        OP_STRING(3)=BLANK(1:15)//'<nc (#s/all)[all]>'
        OP_STRING(4)=BLANK(1:15)//'<using (fit/optimise/solve)[solve]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM list variables<;FILENAME> grid
C###  Description:
C###    List the differences in the solutions of an extracellular
C###    potential problem and a boundary element problem when
C###    iterating to match potentials and fluxes on epicardial
C###    and endocardial boundaries.
C###    Output values are written to the file FILENAME
C###    (with extension .opvari) in the directory specified by PATH
C###    if qualifier present.
C###  Parameter: <grid_region #>[2]
C###    Specify the region which contains the extracellular
C###    solution.
C###  Parameter: <bem_region #>[1]
C###    Specify the region which contains the boundary element
C###    solution.
C###  Parameter: <grid_class #s>[2]
C###    Specify the class of the extracellular potential problem
C###  Parameter: <grid_update_class #s>[3]
C###    Specify the class of the extracellular update problem
C###  Parameter: <bem_class #>[1]
C###    Specify the class of the boundary element problem.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME> grid'
        OP_STRING(2)=BLANK(1:15)//'<grid_region #>[2]'
        OP_STRING(3)=BLANK(1:15)//'<bem_region #>[1]'
        OP_STRING(4)=BLANK(1:15)//'<grid_class #>[2]'
        OP_STRING(5)=BLANK(1:15)//'<grid_update_class #>[3]'
        OP_STRING(6)=BLANK(1:15)//'<bem_class #>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C----------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','LIVARI',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN !file output
          OPFILE=.TRUE.
          CALL CHECKF(3,noco,NTCOQU,CO,COQU,FILE,STRING,*1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.opvari','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'GRID',4,noco+1,NTCO,N3CO)) THEN
          TYPE='GRID'
          !parse class information
          IF(CBBREV(CO,'GRID_CLASS',6,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=2
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nxg,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nxg.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'BEM_CLASS',5,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=1
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nxb,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nxb.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

          IF(CBBREV(CO,'GRID_UPDATE_CLASS',5,noco+1,NTCO,N3CO)) THEN
            nxc=IFROMC(CO(N3CO+1))
          ELSE
            nxc=3
          ENDIF
          CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx_upd.GT.0,
     '      '>>No nx defined for this solve class',ERROR,*9999)

          !parse region information
          IF(CBBREV(CO,'GRID_REGION',6,noco+1,NTCO,N3CO)) THEN
            nrg=IFROMC(CO(N3CO+1))
          ELSE
            nrg=2
          ENDIF
          CALL ASSERT(nrg.GT.0,'>>Invalid grid region number',
     '      ERROR,*9999)

          IF(CBBREV(CO,'BEM_REGION',5,noco+1,NTCO,N3CO)) THEN
            nrb=IFROMC(CO(N3CO+1))
          ELSE
            nrb=1
          ENDIF
          CALL ASSERT(nrb.GT.0,'>>Invalid bem region number',
     '      ERROR,*9999)

          !some initialisations - values should not be needed.
          IY=1
          NCLIST(0)=0

          CALL OPVARI(IY,NBH,NCLIST,NEELEM,NENQ,NHE(1,nxb),
     '      NHP(1,nrb,nxb),NKH(1,1,1,nrb),NPNODE,NP_INTERFACE,NQNP,NQS,
     '      NQXI,nrb,nrg,NVHP(1,1,1,nrb),nxb,nxg,NXQ,NYNE,NYNP,
     '      NYNR(0,0,1,nrb,nxb),AQ,CQ(1,1,nxg),DNUDXQ,DXDXIQ,DXDXIQ2,XQ,
     '      YP(1,1,nxb),YQ(1,1,1,nxg),ZA,ZP,FIX(1,1,nxb),
     '      FIXQ(1,1,nx_upd),TYPE,ERROR,*9999)
        ELSE
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          IF(CBBREV(CO,'NUMBER',1,noco+1,NTCO,N3CO)) THEN
            IY=IFROMC(CO(N3CO+1))
          ELSE
            IY=1
          ENDIF

          IF(CBBREV(CO,'ARRAY',1,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'ZP',1)) THEN
              TYPE='ZP'
            ELSE
              TYPE='YP'
            ENDIF
          ELSE
            TYPE='YP'
          ENDIF

          IF(CBBREV(CO,'NC',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'ALL',2)) THEN
              ALL_NCs=.TRUE.
            ELSE
              ALL_NCs=.FALSE.
              CALL PARSIL(CO(N3CO+1),4,NCLIST(0),NCLIST(1),ERROR,*9999)
              DO no_nc=1,NCLIST(0)
                nc=NCLIST(no_nc)
                CALL ASSERT(nc.GT.0.AND.nc.LE.NCM,
     '            '>>Invalid nc number',ERROR,*9999)
              ENDDO !no_nc
            ENDIF
          ELSE
            ALL_NCs=.TRUE.
          ENDIF

          IF(CBBREV(CO,'USING',2,noco+1,NTCO,N3CO)) THEN
            IF(ABBREV(CO(N3CO+1),'FIT',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_FIT,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this fit class',ERROR,*9999)
            ELSE IF(ABBREV(CO(N3CO+1),'OPTIMISE',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_OPTI,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this optimise class',ERROR,*9999)
            ELSE IF(ABBREV(CO(N3CO+1),'SOLVE',2)) THEN
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ELSE
              CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
              CALL ASSERT(nx.GT.0,
     '          '>>No nx defined for this solve class',ERROR,*9999)
            ENDIF
          ELSE
            CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
            CALL ASSERT(nx.GT.0,
     '        '>>No nx defined for this solve class',ERROR,*9999)
          ENDIF

          DO nolist=1,NRLIST(0)
            nr=NRLIST(nolist)
            IF(ALL_NCs) THEN
              NCLIST(0)=NCT(nr,nx)
              DO no_nc=1,NCLIST(0)
                NCLIST(no_nc)=no_nc
              ENDDO !no_nc
            ENDIF
            CALL OPVARI(IY,NBH,NCLIST,NEELEM,NENQ,NHE(1,nx),
     '        NHP(1,nr,nx),NKH(1,1,1,nr),NPNODE,NP_INTERFACE,NQNP,NQS,
     '        NQXI,nr,nr,NVHP(1,1,1,nr),nx,nx,NXQ,NYNE,NYNP,
     '        NYNR(0,0,1,nr,nx),AQ,CQ(1,1,nx),DNUDXQ,DXDXIQ,DXDXIQ2,
     '        XQ,YP(1,1,nx),
     '        YQ(1,1,1,nx),ZA,ZP,FIX(1,1,nx),FIXQ(1,1,nx),TYPE,
     '        ERROR,*9999)
          ENDDO
        ENDIF

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF

      CALL EXITS('LIVARI')
      RETURN
 9999 CALL ERRORS('LIVARI',ERROR)
      CALL EXITS('LIVARI')
      IF(OPFILE) THEN
        CALL CLOSEF(IOFI,ERROR,*9999)
        IOFI=IOOP
      ENDIF
      RETURN 1
      END
