      SUBROUTINE OPFUNC(IBT,IDO,ID_TYPE,INP,LD,NBH,NBJ,
     '  NEELEM,NHE,NHP,NKH,NKHE,NKJE,NO_OBJECT,
     '  NPF,NPNE,NPNODE,nr,NVHE,NVHP,NVJE,NW,NYNE,NYNP,
     '  CURVCORRECT,SE,XA,XE,XID,XP,YP,ZA,ZE,ZP,
     '  TYPE,ERROR,*)

C CPB 2/6/94 I have changed the call to objfun but I can't find
C where this subroutine is called from so I will just comment it out

C#### Subroutine: OPFUNC
C###  Description:
C###    OPFUNC outputs objective function.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'obje00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),ID_TYPE,
     '  INP(NNM,NIM,NBFM),LD(NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NO_OBJECT,
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*9
!     Local Variables
      INTEGER nb,nd,ND1,ND2,ne,ni,NITB,nj,NU1(0:3),nx
      REAL*8 DD,DXDXI(3,3),DXIDX(3,3),PXI,SUM,X(3),XI(3),ZVAL1(0:3),
     '  ZVAL2(3)
      DATA NU1/1,2,4,7/

      CALL ENTERS('OPFUNC',*9999)
      nx=1 !temporary

      IF(TYPE(1:5).EQ.'FIELD'.OR.TYPE(1:9).EQ.'DEPENDENT') THEN
!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
!old
c        IF(ITYP6(nr,nx).EQ.1) THEN      !linear
c          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
c          CALL YPZP(4,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ENDIF

        ND1=NSOBJE(3,NO_OBJECT) !is 1st nd in object
        ND2=NSOBJE(4,NO_OBJECT) !is 2nd nd in object
        DO nd=ND1,ND2
          ne=LD(nd)
          NITB=NIT(NBJ(1,ne))
          DO ni=1,NITB
            XI(ni)=XID(ni,nd)
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' nd='',I4,'' ne='',I4,'' Xi(ni):'',3E12.3)')
     '        nd,ne,(XI(ni),ni=1,NITB)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(TYPE(1:5).EQ.'FIELD') THEN
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            WRITE(OP_STRING,'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nj=1,NJT
              nb=NBJ(nj,ne)
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '          XE(1,nj))
              DO ni=1,NIT(nb)
                DXDXI(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,NU1(ni),XI,XE(1,nj))
              ENDDO
              WRITE(OP_STRING,'('' X('',I1,'')='',E13.5,'
     '          //''' dX('',I1,'')/dXi:'','
     '          //'3E13.5)') nj,X(nj),nj,(DXDXI(nj,ni),ni=1,NIT(nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
            nb=NBJ(NJ_LOC(NJL_GEOM,0,nr)+ID_TYPE,ne)
            NITB=NIT(nb)
            DO ni=0,NITB
              ZVAL1(ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          NU1(ni),XI,XE(1,NJ_LOC(NJL_GEOM,0,nr)+ID_TYPE))
            ENDDO

          ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            WRITE(OP_STRING,'('' '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nj=1,NJT
              nb=NBJ(nj,ne)
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '          XE(1,nj))
              DO ni=1,NIT(nb)
                DXDXI(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,NU1(ni),XI,XE(1,nj))
              ENDDO
              WRITE(OP_STRING,'('' X('',I1,'')='',E13.5,'
     '          //''' dX('',I1,'')/dXi:'','
     '          //'3E13.5)') nj,X(nj),nj,(DXDXI(nj,ni),ni=1,NIT(nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
            nb=NBH(ID_TYPE,1,ne)
            NITB=NIT(nb)
            DO ni=0,NITB
              ZVAL1(ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '          NU1(ni),XI,ZE(1,ID_TYPE))
            ENDDO
          ENDIF
          WRITE(OP_STRING,'('' Function f='',E13.5,'' d(f)/d(Xi):'','
     '      //'3E13.5)') (ZVAL1(ni),ni=0,NITB)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(NJT.EQ.2.AND.NITB.EQ.2) THEN
            CALL INVERT(2,DXDXI,DXIDX,DD)
            DO nj=1,NJT
              SUM=0.d0
              DO ni=1,NITB
                SUM=SUM+ZVAL1(ni)*DXIDX(ni,nj)
              ENDDO
              ZVAL2(nj)=SUM
            ENDDO
            WRITE(OP_STRING,'('' Gradient d(f)/d(Xj):'',3E13.5)')
     '        (ZVAL2(nj),nj=1,NJT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO

      ELSE IF(TYPE(1:9).EQ.'OBJECTIVE') THEN
C        CALL OBJFUN(IBT,IDO,INP,LD,LGE,NAN,NBH,NBJ,
C     '    NEELEM,NFF,NGAP,NHE,NHP,NKE,NKEF,NKH,
C     '    NLL,NLNO,NNF,NNL,NONL,NONY,NP_INTERFACE,NP1OPT,NPF,NPNE,
C     '    NPL,NPNODE,NPNY,nr,NRE,NVHE,NVHP,NVJE,NVJP,NW,
C     '   NYNE,NYNO,NYNP,PAOPTY,
C     '    CE,CG,CP,DL,FEXT,FGRAD,PAOPTI,PF,PG,RE1,RESID,RESIDM,
C     '    RESJAC,RG,SE,WG,WU,XA,XE,XF,XG,XID,XIG,XP,YG,YP,
C     '    ZA,ZD,ZE,ZF,ZG,ZP,FIX,ERROR,*9999)
C        WRITE(OP_STRING,'('' Objective function = '',D12.4)') FUNC
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPFUNC')
      RETURN
 9999 CALL ERRORS('OPFUNC',ERROR)
      CALL EXITS('OPFUNC')
      RETURN 1
      END


