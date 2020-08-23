      SUBROUTINE D_ZPRP(PARAMTYPE,IBT,IDO,INP,LGE,NAN,NBH,NBJ,NBJF,
     '  NEELEM,NFF,NGAP,NHE,NKEF,NKHE,NKJE,NMNO,NNF,NPF,NPNE,NPNY,
     '  nr,NRE,NVHE,NVJE,NW,nx,NXI,NYNE,NYNP,NYNR,
     '  CE,CG,CGE,CP,CURVCORRECT,D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,
     '  FIX,PG,RE1,RE2,RG,SE,WG,XA,XE,XG,XP,YG,
     '  ZA,ZE,ZE1,ZG,ZG1,ZP,ERROR,*)

C#### Subroutine: D_ZPRP
C###  Description:
C###    D_ZPRP calculates derivatives of global residual vector
C###    YP(ny,4) at current solution, wrt the material params to be
C###    optimised.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  LGE(NHM*NSM,NRCM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),NFF(6,NEM),
     '  NGAP(NIM,NBM),NHE(NEM),NKEF(0:4,16,6,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NMNO(1:2,0:NOPM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),
     '  nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),NYNR(0:NY_R_M,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     '  CURVCORRECT(2,2,NNM,NEM),D_RE(NSM,NHM,NOPM),
     '  D_RI3(NHM*NSM),D_RP(NYM,NYM),D_TG(3,3,NHM*NSM),
     '  D_ZG(NHM,NUM,NHM*NSM),ES(NHM*NSM,NHM*NSM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),
     '  RE1(NSM,NHM),RE2(NSM,NHM),RG(NGM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZE1(NSM,NHM),
     '  ZG(NHM,NUM),ZG1(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER PARAMTYPE*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER nb,nc,ne,nh,nh1,nhs,nhs1,NHST(2),nhx,nhx1,noelem,noopti,
     '  no_nynr,no_nynr1,no_nynr2,ns,ns1,ny,ny1,ny2

      CALL ENTERS('D_ZPRP',*9999)

      nc=1 !LHS

      IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
        CALL ASSERT(NTOPTI.LE.NYM,'NYM too small',ERROR,*9999)
        DO no_nynr=1,NYNR(0,1,1) !loop over rows
          ny=NYNR(no_nynr,1,1) !is local[global] row number
          DO noopti=1,NTOPTI
            D_RP(ny,noopti)=0.0d0
          ENDDO !noopti
        ENDDO !no_nynr (ny)
      ELSE IF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
        DO no_nynr1=1,NYNR(0,1,1) !loop over rows
          ny1=NYNR(no_nynr1,1,1) !is local[global] row number
          DO no_nynr2=1,NYNR(0,0,1) !loop over global columns
            ny2=NYNR(no_nynr2,0,1) !is global column number
            D_RP(ny1,ny2)=0.0d0
          ENDDO !no_nynr2 (ny2)
        ENDDO !no_nynr1 (ny1)
      ENDIF

      DO noelem=1,NEELEM(0,nr) !is main element loop
        ne=NEELEM(noelem,nr)
        IF(NW(ne,1).GE.0) THEN

          CALL MELGE(LGE,NBH(1,1,ne),1,ne,NHE(ne),NHST,
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nx,NYNE,NYNP,ERROR,*9999)

          IF(PARAMTYPE(1:19).EQ.'MATERIAL_PARAMETERS') THEN
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '        ZE,ZP,ERROR,*9999)
            nb=NBH(NH_LOC(1,nx),1,ne)
            IF(ITYP1(nr,nx).EQ.3) THEN
              ERROR=' Analytic derivs not implemented for PDE''s'
              GO TO 9999
            ELSE IF(ITYP1(nr,nx).EQ.4) THEN
              ERROR=' Analytic derivs not implemented for lin. elast.'
              GO TO 9999
            ELSE IF(ITYP1(nr,nx).EQ.5) THEN
              CALL CPCG(1,nb,NPNE(1,1,ne),nr,nx,CE(1,ne),CG,
     '          CGE(1,1,ne),CP,PG,
     '          ERROR,*9999)
              CALL D_ZERE50(PARAMTYPE,IBT,IDO,INP,NAN,NBH(1,nc,ne),
     '          NBJ(1,ne),NBJF,NGAP,ne,NFF(1,ne),NHE(ne),
     '          NKEF,NMNO,NNF,NPNE,nr,NRE,NW(ne,1),nx,NXI,
     '          CE(1,ne),CG,CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT(1,1,ne),
     '          PG,RG,SE,WG,XE,XG,YG(1,1,ne),ZE,ZE1,ZG,ZG1,
     '          ERROR,*9999)
            ENDIF
            nhs=0
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              DO ns=1,NST(NBH(nh,1,ne))+NAT(NBH(nh,1,ne))
                nhs=nhs+1
                ny=ABS(LGE(nhs,1))
                DO noopti=1,NTOPTI
                  D_RP(ny,noopti)=D_RP(ny,noopti)+D_RE(ns,nh,noopti)
                ENDDO
              ENDDO
            ENDDO

          ELSEIF(PARAMTYPE(1:20).EQ.'GEOMETRIC_PARAMETERS') THEN
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,
     '        XP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,
     '        ZP,ERROR,*9999)
            CALL ZEES(IBT,IDO,INP,LGE,NAN,NBH(1,nc,ne),NBJ(1,ne),
     '        NBJF,ne,NFF(1,ne),NGAP,NHE(ne),NKEF,NMNO,
     '        NNF,NPNE,NPNY,nr,NRE,NW(ne,1),nx,NXI,NYNE,NYNP,
     '        CE(1,ne),CG,CGE(1,1,ne),
     '        CP,D_RE,D_RI3,D_TG,D_ZG,ES,FEXT(1,1,ne),
     '        PG,RE1,RE2,RG,SE,WG,
     '        XE,XG,YG(1,1,ne),ZE,ZE1,%VAL(0),ZG,ZG1,FIX,ERROR,*9999)
            nhs=0
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              DO ns=1,NST(NBH(nh,1,ne))+NAT(NBH(nh,1,ne))
                nhs=nhs+1
                ny=ABS(LGE(nhs,1))
                nhs1=0
                DO nhx1=1,NH_LOC(0,nx)
                  nh1=NH_LOC(nhx1,nx)
                  DO ns1=1,NST(NBH(nh1,1,ne))+NAT(NBH(nh1,1,ne))
                    nhs1=nhs1+1
                    ny1=ABS(LGE(nhs1,2))
                    D_RP(ny,ny1)=D_RP(ny,ny1)
     '                +ES(nhs,nhs1)*SE(ns1,NBH(nh1,1,ne),ne)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('D_ZPRP')
      RETURN
 9999 CALL ERRORS('D_ZPRP',ERROR)
      CALL EXITS('D_ZPRP')
      RETURN 1
      END


