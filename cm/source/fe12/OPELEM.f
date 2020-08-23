      SUBROUTINE OPELEM(IBT,NBJ,NBH,NEELEM,NELIST,NFF,NHE,
     '  NHP,nh_update,NKH,NKHE,NKJE,NLL,NP_average,NPF,
     '  NPLIST,NPNE,NPNODE,NRE,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,VOLTC,
     '  ELEM_VOL,nj_volume,CURVCORRECT,PG,RG,SE,VOL,VOLT,WG,XA,XAB,XE,
     '  XG,XP,YP,ZA,ZE,ZG,ZP,AVGENODE,OUTPUT_TOT_ONLY,TYPE,ERROR,*)

C#### Subroutine: OPELEM
C###  Description:
C###    OPELEM outputs element topology.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NFF(6,NEM),
     '  NHE(NEM),NHP(NPM,0:NRM),nh_update,nj_volume,
     '  NKH(NHM,NPM,NCM,0:NRM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NP_average,NPF(9,NFM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),VOLTC(NBFM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),VOL(NBFM),VOLT(NBFM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XAB(NORM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL AVGENODE,ELEM_VOL,OUTPUT_TOT_ONLY
!     Local Variables
      INTEGER ERR,nb,nc,ne,nh,nhx,noelem,nonode,np
      REAL*8 DUMMY(1),SUM,VOLAVG,Z(3),ZRC(3),ZRC_AVGE(3)
      CHARACTER NAME1(3)*6

      DATA NAME1(1)/'Length'/,NAME1(2)/'  Area'/,NAME1(3)/'Volume'/

      CALL ENTERS('OPELEM',*9999)

C LKC 2-MAY-1999 init dummy
      DUMMY(1)=0.0d0

      WRITE(OP_STRING,'(/'' Listing elements:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(NELIST(0).GT.0) THEN
        CALL WRITE_LONG(INTTYPE,1,1,IOFI,NELIST(0),10,10,NELIST(1),
     '    DUMMY,'(1X,10I6)','(1X,10I6)',ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'(/'' No elements defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        GO TO 9998
      ENDIF

      DO nb=1,NBFT
          VOLT(nb)=0.0d0
          VOLTC(nb)=0
      ENDDO
      VOLAVG=0.d0

      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
        nr=NRE(ne) !changed AJP 14/4/97
        IF(TYPE(1:10).EQ.'UNDEFORMED')THEN
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          CALL OPELEM1(ELEM_VOL,IBT,NBJ(1,ne),ne,
     '      NFF(1,ne),nj_volume,NKJE(1,1,1,ne),
     '      NLL(1,ne),NPNE(1,1,ne),nr,NRE,NVJE(1,1,1,ne),VOLTC,
     '      PG,RG,SE(1,1,ne),VOL,VOLT,WG,XAB,XE,XG,
     '      OUTPUT_TOT_ONLY,ERROR,*9999)
        ELSE IF(TYPE(1:8).EQ.'DEFORMED')THEN
          IF(TYPE(1:15).EQ.'DEFORMED_CAVITY')THEN
            nr=1
            nx=1
          ENDIF
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,
     '      ERROR,*9999)

          IF(AVGENODE) THEN
C new MPN 17Apr97: average RC coords instead of curv. coords.
C           Calc RC coords of node to be changed
            DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
              nh=NH_LOC(nhx,nx)
              Z(nhx)=ZP(1,1,nh,NP_average,1)
            ENDDO !nhx (nh)
            CALL XZ(ITYP11(nr),Z,ZRC_AVGE)
C           Calc the average nh_update RC coord for nodes in NPLIST
            SUM=0.0d0
            DO nonode=1,NPLIST(0)
              np=NPLIST(nonode)
              DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
                nh=NH_LOC(nhx,nx)
                Z(nhx)=ZP(1,1,nh,np,1)
              ENDDO !nhx (nh)
C             Convert curv. coords to rc coords
              CALL XZ(ITYP11(nr),Z,ZRC)
              SUM=SUM+ZRC(nh_update)
            ENDDO !nonode (np)
            ZRC_AVGE(nh_update)=SUM/DBLE(NPLIST(0))
C           Convert averaged RC coords to curv. coords
            CALL ZX(ITYP11(nr),ZRC_AVGE,Z)
C           Put averaged curv. coords into node NP_average
C           (unless there is more that one version)
            DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
              nh=NH_LOC(nhx,nx)
              IF(NVHP(nh,NP_average,1,nr).EQ.1)
     '          ZP(1,1,nh,NP_average,1)=Z(nhx)
            ENDDO !nhx
C old MPN 17Apr97
CC           Calc the average nj_update value for nodes in NPLIST
C            SUM=0.0d0
C            DO nonode=1,NPLIST(0)
C              np=NPLIST(nonode)
C              SUM=SUM+ZP(1,1,nj_update,np,1)
C            ENDDO !np
C            ZP(1,1,nj_update,NP_average,1)=SUM/DBLE(NPLIST(0))
          ENDIF

          IF(TYPE(1:15).EQ.'DEFORMED_CAVITY')THEN
            DO nh=1,NHM
              DO nc=1,NCM
                NBH(nh,nc,ne)=3
              ENDDO
            ENDDO
            NHE(ne)=3
            CALL ZPZE_CAVITY(NBH(1,1,ne),1,
     '        NKJE(1,1,1,ne),NPNE(1,1,ne),NPNODE,nr,
     '        NVJE(1,1,1,ne),SE(1,1,ne),XP,ZE,ZP,ERROR,*9999)
            CALL OPELEMD(NBH(1,1,ne),NBJ(1,ne),
     '        ne,NFF(1,ne),NHE(ne),NKJE(1,1,1,ne),
     '        NLL(1,ne),NPNE(1,1,ne),nr,NRE,NVHE(1,1,1,ne),nx,VOLTC,
     '        PG,RG,SE(1,1,ne),VOL,VOLT,WG,ZE,ZG,
     '        OUTPUT_TOT_ONLY,ERROR,*9999)
          ELSE
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),
     '        NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
     '        NVHE(1,1,1,ne),NW(ne,1),nx,CURVCORRECT(1,1,1,ne),
     '        SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,ERROR,*9999)
            CALL OPELEMD(NBH(1,1,ne),NBJ(1,ne),
     '        ne,NFF(1,ne),NHE(ne),NKJE(1,1,1,ne),
     '        NLL(1,ne),NPNE(1,1,ne),nr,NRE,NVHE(1,1,1,ne),nx,VOLTC,
     '        PG,RG,SE(1,1,ne),VOL,VOLT,WG,ZE,ZG,
     '        OUTPUT_TOT_ONLY,ERROR,*9999)
          ENDIF


        ENDIF
      ENDDO

C LKC 3-FEB-1999 Modify volume calculations
C      DO nb=1,NBFT
C        NITOT=NIT(nb)
C        IF(DABS(VOLT(NITOT,nb)).GT.1.0d-6) THEN
C          IF(ITYP1(nr,1).EQ.6) NITOT=NJ_LOC(NJL_GEOM,0,nr)
C          WRITE(OP_STRING,'(/(4X,''Total   '',A,'' = '',D13.5,'
C     '      //'''  (Basis type '',I2,'')''))')
C     '      (NAME1(ni+JTYP4-1),VOLT(ni,nb),nb,ni=NIT(nb),NITOT)
C          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C        ENDIF
C      ENDDO

      IF(TYPE(1:15).EQ.'DEFORMED_CAVITY') THEN
        IF (NRE(NELIST(1)).EQ.3) THEN
          CALL SET_USER_DOUBLE("LV_CAVITY_VOLUME",VOLT(3),ERR)
        ELSE IF(NRE(NELIST(1)).EQ.4) THEN
          CALL SET_USER_DOUBLE("RV_CAVITY_VOLUME",VOLT(3),ERR)
        ENDIF
      ENDIF

      ne=0 ! use as a total element counter
      DO nb=1,NBFT
        WRITE(OP_STRING,'(1P/4X,''Total     '',A,'' = '',E13.5,'
     '    //'''  (Basis type '',I2,'')'')')
     '    NAME1(NIT(nb)+JTYP4-1),VOLT(nb),nb
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(NIT(nb).EQ.3)THEN
          CALL SET_USER_DOUBLE('MESH_VOLUME',VOLT(nb),ERROR)
        ENDIF

        WRITE(OP_STRING,'(4X,''Num. Elements    = '',I13)')
     '    VOLTC(nb)
        ne=ne+VOLTC(nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(VOLTC(nb).NE.0) THEN
          WRITE(OP_STRING,'(1P,4X,''Average   '',A,'' = '',E13.5)')
     '      NAME1(NIT(nb)+JTYP4-1),VOLT(nb)/DBLE(VOLTC(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          VOLAVG=VOLAVG+VOLT(nb)
        ENDIF
      ENDDO
      WRITE(OP_STRING,'(1P/4X,''Tot. Avg. '',A,'' = '',E13.5)')
     '    NAME1(NIT(1)+JTYP4-1),VOLAVG/NELIST(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(4X,''Tot. Elements    = '',I13)') ne
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' '')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

 9998 CALL EXITS('OPELEM')
      RETURN
 9999 CALL ERRORS('OPELEM',ERROR)
      CALL EXITS('OPELEM')
      RETURN 1
      END


