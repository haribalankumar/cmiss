      SUBROUTINE OPSTRA(IBT,IDO,INDEX,INP,
     '  IPOINTTYP,IXI,NAN,NBH,NBJ,NDDL,NDLT,
     '  NEELEM,NELIST,NGLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '  NPNODE,nr,NTPOIN,NVHE,NVHP,NVJE,
     '  NW,nx,nxc,NYNE,NYNP,CE,CG,CGE,CP,CURVCORRECT,PG,RG,SE,XA,XE,XG,
     '  XID,XIG,XIPOS,XP,YP,ZA,ZE,ZG,ZP,COORDS,STATIC,HIGH_PRECISION,
     '  EXTREMA_ONLY,FULL,TYPE,ERROR,*)

C#### Subroutine: OPSTRA
C###  Description:
C###    OPSTRA outputs strain tensors, invariants and principal strains
C###    at Gauss pts or at equally spaced points along an arbitrary
C###    xi coordinate line.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'moti00.cmn'
      INCLUDE 'stra00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IPOINTTYP,IXI,NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NDDL(NEM,NDEM),NDLT(NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NGLIST(0:NGM),
     '  NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,NTPOIN,
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,nxc,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  curvcorrect(2,2,NNM,NEM),DZDX(3,3),EG(3,3),
     '  PG(NSM,NUM,NGM,NBM),PHI(3),PST(3),
     '  R(3,3),RG(NGM),RM(3,3),SE(NSM,NBFM,NEM),U(3,3),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XID(NIM,NDM),XIG(NIM,NGM,NBM),XIPOS(3),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER COORDS*(*),TYPE*3,ERROR*(*)
      LOGICAL HIGH_PRECISION,EXTREMA_ONLY,FULL,STATIC
!     Local Variables
      INTEGER i,IBEG,IEND,ig,j,na,nb,NBFF,NBG,nd,nde,ne,ng,nh,nhx,
     '  nj,njj,nk,nn,nolist,nopoin,NSF,NSG
      REAL*8 PF1,RI1,RI2,RI3,TEMP,TEMP2

      CALL ENTERS('OPSTRA',*9999)

      DO i=1,3
        DO j=1,3
          DZDX(i,j)=0.0d0
        ENDDO
      ENDDO

      WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '  //''') Region '',I1,'':'')') nxc,nx,nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(.NOT.STATIC) THEN
        CALL YPZP(MOTION_IY,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
      ENDIF

      CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '  nr,NVHP,nx,NYNE,NYNP,YP(1,1,nx),ZA,ZP,ERROR,*9999)

      IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
        IF(TABLE.AND..NOT.EXTREMA_ONLY) THEN
          CALL STRING_TRIM(COORDS,IBEG,IEND)
          FORMAT='('' '//COORDS(IBEG:IEND)//' strains: '//
     '      'np,((EG(mi,ni),mi=1,ni),ni=1,NIT),'//
     '      '(PST(ni),ni=1,NIT),(Phi(ni),ni=1,NIT):'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='('' ne='',30(I5,'',''),:/(4X,30(I5,'','')))'
          WRITE(OP_STRING,FORMAT) (NELIST(nolist),nolist=1,NELIST(0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        nb=NBJ(1,ne)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '    NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '    CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '    ZE,ZP,ERROR,*9999)

        IF(.NOT.STATIC) THEN
C         output strain at gauss points only at pres.
          CALL ASSERT(IPOINTTYP.EQ.1,
     '      '>>Only Gauss pt strains implemented for non-static',
     '      ERROR,*9999)
          DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
            nj=NH_LOC(NJL_GEOM,njj)
            nhx=njj
            nh=NH_LOC(nhx,nx)
            NBFF=NBH(nh,1,ne)  !time basis
            NBG=NBJ(nj,ne)  !corresponding spatial basis
            NSG=0
            DO nn=1,NNT(NBG)
              DO nk=1,NKT(nn,NBG)
                NSG=NSG+1
                TEMP=0.0D0
                TEMP2=0.0D0
                DO na=1,IBT(2,NIT(NBFF),NBFF) !#terms in series
                  NSF=(nn-1)*NKT(0,NBFF)+(na-1)*NKT(0,NBG)+nk
                  TEMP=TEMP+PF1(na,1,TIME)*ZE(NSF,nhx)
                  TEMP2=TEMP2+PF1(na,1,TIME_REF)*ZE(NSF,nhx)
                ENDDO !na
                ZE(NSG,nh)=TEMP+XE(NSG,nj)
                XE(NSG,nj)=TEMP2+XE(NSG,nj)
              ENDDO !nk
            ENDDO !nn
          ENDDO !nj
          CALL STRING_TRIM(COORDS,IBEG,IEND)
          FORMAT='(/'' Element'',I3,'' Gauss point strains'
     '      //' wrt '//COORDS(IBEG:IEND)//' coords:'')'
          WRITE(OP_STRING,FORMAT) ne
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO ig=1,NGLIST(0)
            ng=NGLIST(ig)
            IF(KTYP51(nr).EQ.5) THEN !string theory
              CALL OPEG55(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),
     '          ng,NHE(ne),nr,nx,EG,PG,RG(ng),
     '          XE,XG,XIG(1,ng,nb),ZE,ZG,ERROR,*9999)
            ELSE !plane stress,strain/3D/membrane/shell
C KAT 11Feb98: Don't need CG
C              CALL CPCG(NW(ne,1),NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
C     '          CE(1,ne),CG,CP,PG,ERROR,*9999)
              CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBJ(1,ne),
     '          NBJ(1,ne),ng,NHE(ne),NPNE(1,1,ne),nr,nx,
     '          DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '          R,RG(ng),RI1,RI2,RI3,RM,U,
     '          XE,XG,XIG(1,ng,nb),ZE,ZG,ERROR,*9999)
              IF(.NOT.EXTREMA_ONLY) THEN
                CALL OPEG50(COORDS,IPOINTTYP,NBH(1,1,ne),ng,nr,nx,
     '            DZDX,EG,PHI,PST,R,RI1,RI2,RI3,RM,U,
     '            XG,XIG(1,ng,nb),ZG,FULL,.FALSE.,ERROR,*9999)
              ENDIF
            ENDIF
          ENDDO !ig

        ELSE IF(ITYP1(nr,nx).EQ.5) THEN !finite elasticity
          IF(IPOINTTYP.EQ.1) THEN !strain at Gauss points
C          IF(IXI.EQ.0) THEN !strain at Gauss points
            IF(.NOT.TABLE.AND..NOT.EXTREMA_ONLY) THEN
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              FORMAT='(/'' Element'',I3,'' Gauss point strains'
     '          //' wrt '//COORDS(IBEG:IEND)//' coords:'')'
              WRITE(OP_STRING,FORMAT) ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            IF(KTYP55(nr).EQ.3.AND.KTYP56(nr).EQ.3) THEN !pole zero law
              CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
     '          CE(1,ne),CG,CGE(1,1,ne),CP,PG,ERROR,*9999)
            ENDIF
            DO ig=1,NGLIST(0)
              ng=NGLIST(ig)
              IF(KTYP51(nr).EQ.5) THEN !string theory
                CALL OPEG55(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),
     '            ng,NHE(ne),nr,nx,EG,PG,RG(ng),
     '            XE,XG,XIG(1,ng,nb),ZE,ZG,ERROR,*9999)
              ELSE !plane stress,strain/3D/membrane/shell
C KAT 11Feb98: Only need to calculate CG once per element
C                CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
C     '            CE(1,ne),CG,CP,PG,ERROR,*9999)
                CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NBJ(1,ne),ng,NHE(ne),NPNE(1,1,ne),nr,nx,
     '            DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '            R,RG(ng),RI1,RI2,RI3,RM,U,
     '            XE,XG,XIG(1,ng,nb),ZE,ZG,ERROR,*9999)
                IF(.NOT.EXTREMA_ONLY) THEN
C TVK 16Nov99: Option for higher precision output
                  CALL OPEG50(COORDS,IPOINTTYP,NBH(1,1,ne),ng,nr,nx,
     '              DZDX,EG,PHI,PST,R,RI1,RI2,RI3,RM,U,
     '              XG,XIG(1,ng,nb),ZG,FULL,
     '              HIGH_PRECISION,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO !ig

          ELSE IF(IPOINTTYP.EQ.2) THEN !strain at data points
            IF(.NOT.TABLE.AND..NOT.EXTREMA_ONLY) THEN
              CALL STRING_TRIM(COORDS,IBEG,IEND)
              FORMAT='(/'' Element'',I3,'' Data point strains'
     '          //' wrt '//COORDS(IBEG:IEND)//' coords:'')'
              WRITE(OP_STRING,FORMAT) ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nde=1,NDLT(ne)
              nd=NDDL(ne,nde)
              IF(KTYP51(nr).EQ.5) THEN !string theory
                CALL OPEG55(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),
     '            0,NHE(ne),nr,nx,EG,PG,RG(ng),
     '            XE,XG,XID(1,nd),ZE,ZG,ERROR,*9999)
              ELSE !plane stress,strain/3D/membrane/shell
                CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,nx,
     '            DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST, !(Don't need CG)
     '            R,RG(1),RI1,RI2,RI3,RM,U,
     '            XE,XG,XID(1,nd),ZE,ZG,ERROR,*9999)
                IF(.NOT.EXTREMA_ONLY) THEN
                  CALL OPEG50(COORDS,IPOINTTYP,NBH(1,1,ne),nd,nr,nx,
     '              DZDX,EG,PHI,PST,R,RI1,RI2,RI3,RM,U,
     '              XG,XID(1,nd),ZG,FULL,.FALSE.,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO !nde

          ELSE IF(IXI.GE.1.AND.IXI.LE.4) THEN !strain at Xi point
            IF(.NOT.TABLE.AND..NOT.EXTREMA_ONLY) THEN
              FORMAT='(/'' Element'',I3,'' strain solutions:'')'
              WRITE(OP_STRING,FORMAT) ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nopoin=1,NTPOIN
              IF((NTPOIN.LE.1).AND.(IXI.NE.4)) THEN
                XIPOS(IXI)=0.5D0
              ELSE IF(IXI.NE.4) THEN
                XIPOS(IXI)=DBLE(nopoin-1)/DBLE(NTPOIN-1)
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' IXI='',I2,'' XI_1='',F5.3,'
     '            //''' XI_2='',F5.3,'' XI_3='',F5.3)')
     '            IXI,XIPOS(1),XIPOS(2),XIPOS(3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(KTYP51(nr).EQ.5) THEN !string theory
                CALL OPEG55(IBT,IDO,INP,NAN,NBH(1,1,ne),NBJ(1,ne),
     '            0,NHE(ne),nr,nx,EG,PG,RG(ng),
     '            XE,XG,XIPOS,ZE,ZG,ERROR,*9999)
              ELSE !plane stress,strain/3D/membrane/shell
C KAT 11Feb98: Don't need CG
C                CALL CPCG(1,NBH(1,1,ne),NPNE(1,1,ne),nr,nx,
C     '            CE(1,ne),CG,CP,PG,ERROR,*9999)
                CALL ZEEX50(COORDS,IBT,IDO,INP,NAN,NBH(1,1,ne),
     '            NBJ(1,ne),0,NHE(ne),NPNE(1,1,ne),nr,nx,
     '            DZDX,CE(1,ne),CG,CP,EG,PG,PHI,PST,
     '            R,RG(1),RI1,RI2,RI3,RM,U,
     '            XE,XG,XIPOS,ZE,ZG,ERROR,*9999)
                IF(.NOT.EXTREMA_ONLY) THEN
                  CALL OPEG50(COORDS,IPOINTTYP,NBH(1,1,ne),nopoin,nr,nx,
     '              DZDX,EG,PHI,PST,R,RI1,RI2,RI3,RM,U,
     '              XG,XIPOS,ZG,FULL,.FALSE.,ERROR,*9999)

                  IF(TYPE.EQ.'E11') THEN
                    YP(INDEX,7,1)=EG(1,1)
                  ELSE IF(TYPE.EQ.'E22') THEN
                    YP(INDEX,7,1)=EG(2,2)
                  ELSE IF(TYPE.EQ.'E33') THEN
                    YP(INDEX,7,1)=EG(3,3)
                  ELSE IF(TYPE.EQ.'E12') THEN
                    YP(INDEX,7,1)=EG(1,2)
                  ELSE IF(TYPE.EQ.'E13') THEN
                    YP(INDEX,7,1)=EG(1,3)
                  ELSE IF(TYPE.EQ.'E23') THEN
                    YP(INDEX,7,1)=EG(2,3)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO !nopoin
          ENDIF
        ENDIF
      ENDDO !nolist (ne)

      IF(FULL) THEN
        WRITE(OP_STRING,'(/'' Min princ strain for listed points = '','
     '    //'D12.4)') PRSTMIN
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Max princ strain for listed points = '','
     '    //'D12.4)') PRSTMAX
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPSTRA')
      RETURN
 9999 CALL ERRORS('OPSTRA',ERROR)
      CALL EXITS('OPSTRA')
      RETURN 1
      END


