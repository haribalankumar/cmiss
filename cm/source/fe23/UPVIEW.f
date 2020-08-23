      SUBROUTINE UPVIEW(IBT,IDO,INP,ISEG,ISELNO,ISFIBR,
     '  ISFIEL,ISLINE,ISLINO,ISNONO,IWK,MXI,NAN,NBH,NBJ,NEELEM,
     '  NGAP,NHE,NHP,NKH,NKHE,NKJE,NLATNE,NLL,NLLIST,NPF,NPL,NPNE,
     '  NPNODE,NQNLAT,NQNE,NQS,NQXI,NRE,NTIW,NVHE,NVHP,NVJE,NW,nx,
     '  NYNE,NYNP,CURVCORRECT,DL,SE,XA,XE,XG,XP,XQ,YG,YP,YQ,ZA,ZE,ZP,
     '  CSEG,FIX,ERROR,*)

C#### Subroutine: UPVIEW
C###  Description:
C###    UPVIEW updates mesh on workstation viewports.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISELNO(NWM,NEM),ISFIBR(NWM,NEM,NGRSEGM),ISFIEL(NWM,NEM),
     '  ISLINE(NWM,2*NGRSEGM),ISLINO(NWM),
     '  ISNONO(NWM,NPM),IWK(6),MXI(2,NEM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NGAP(NIM,NBM),NHE(NEM),NHP(NPM,0:NRM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NLATNE(NEQM+1),
     '  NLL(12,NEM),NLLIST(0:NLM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NQNE(NEQM,NQEM),
     '  NQNLAT(NEQM*NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),NRE(NEM),NTIW,
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),DL(3,NLM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),YG(NIYGM,NGM,NEM),
     '  YP(NYM,NIYM),YQ(NYQM,NIQM,NAM,NXM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER ID_TYPE,INDEX,INDEX_POLYLINE,INDEX_TEXT,iw,nb,nc,ne,nj,nn,
     '  noelem,nofibr,noiw,NOLINE,nonode,NOSEGM,np,nr
      REAL*8 RFROMC,X(3),XI(3),Z(3),ZCENTR(3)
      LOGICAL FINITE,STATIC_VIEW

      CALL ENTERS('UPVIEW',*9999)

      nc=1 !temporary cpb 22/11/94

C cpb 4/3/96 Check to see if nx is not zero
      IF(nx.EQ.0) THEN
        FINITE=.FALSE.
      ELSE
        DO nr=1,NRT
          IF((ITYP2(nr,nx).EQ.1.OR.ITYP2(nr,nx).EQ.2).
     '      AND.ITYP6(nr,nx).EQ.2) THEN
            FINITE=.TRUE.
            CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '        nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
        ELSE
          FINITE=.FALSE.
        ENDIF
      ENDDO
      ENDIF

      DO noiw=1,NTIW
        iw=IWK(noiw)
        IF(IWKS(iw).GT.0) THEN
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*)' update on iw=',iw
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          CALL ACWK(iw,1,ERROR,*9999)
          IF(iw.EQ.5.OR.iw.EQ.6) THEN
C old MPN unused?
            CALL ASSERT(.FALSE.,'>>Not Implemented',ERROR,*9999)
c           CALL FBCLEAR(iw)
          ENDIF
          IF(ISNONO(iw,NPNODE(1,1)).GT.0) THEN !node segment is defined
            INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLACK')
            DO nr=1,NRT
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                NOSEGM=ISNONO(iw,np)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' np='',I4,'' NOSEGM='',I4)')
     '              np,NOSEGM
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                IF(ISEG(NOSEGM).EQ.2) THEN !..and visible
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    IF(FINITE) THEN
                      X(nj)=ZP(1,1,nj,np,nc)
                    ELSE
                      X(nj)=XP(1,1,nj,np)
                    ENDIF
                  ENDDO
                  CALL XZ(ITYP10(nr),X,Z)
                  CALL SGNODE(INDEX,ISEG,NOSEGM,iw,NHP(1,nr),
     '              NKH(1,1,1,nr),np,nr,nx,NYNP,FIX,Z(1),Z(2),Z(3),CSEG,
     '              ERROR,*9999)
                ENDIF
              ENDDO !nonode
            ENDDO !nr
          ENDIF

          nr=1 !temporarily
          DO NOLINE=1,NTLINE
            NOSEGM=ISLINE(iw,NOLINE)
            IF(NOSEGM.GT.0) THEN !line segment is defined
              IF(ISEG(NOSEGM).EQ.2) THEN !..and visible
                IF(CSEG(ISLINE(iw,NOLINE))(43:48).ne.'      ') THEN
                  STATIC_VIEW=.FALSE.
                  TIME=RFROMC(CSEG(ISLINE(iw,NOLINE))(43:52))
                ELSE
                  STATIC_VIEW=.TRUE.
                ENDIF
                IF(FINITE) THEN
                  INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
                  CALL SGLINE(INDEX,ISEG,ISLINE(iw,nr),ISLINO(iw),iw,
     '              NLLIST,NOLINE,NPL,nr,nx,CSEG,'DEFORMED',DL,'SOLID',
     '              STATIC_VIEW,XP,ZP,ERROR,*9999)
                ELSE
                  INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
                  CALL SGLINE(INDEX,ISEG,ISLINE(iw,nr),ISLINO(iw),iw,
     '              NLLIST,NOLINE,NPL,nr,nx,CSEG,'UNDEFORMED',DL,
     '              'SOLID',STATIC_VIEW,XP,ZP,ERROR,*9999)
                ENDIF
              ENDIF
            ENDIF
          ENDDO

          INDEX=INDEX_TEXT(0,'WIDTH1','FONT1','BLUE')
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NOSEGM=ISELNO(iw,ne)
              IF(NOSEGM.GT.0) THEN !element segment is defined
                IF(ISEG(NOSEGM).EQ.2) THEN !..and visible
                  DO nj=1,NJT
                    ZCENTR(nj)=0.d0
                  ENDDO
                  nb=NBJ(1,ne)
                  DO nn=1,NNT(nb)
                    np=NPNE(nn,nb,ne)
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      IF(FINITE) THEN
                        X(nj)=ZP(1,1,nj,np,nc)
                      ELSE
                        X(nj)=XP(1,1,nj,np)
                      ENDIF
                    ENDDO
                    CALL XZ(ITYP10(nr),X,Z)
                    DO nj=1,NJT
                      ZCENTR(nj)=ZCENTR(nj)+Z(nj)
                    ENDDO
                  ENDDO
                  DO nj=1,NJT
                    ZCENTR(nj)=ZCENTR(nj)/NNT(nb)
                  ENDDO
                  CALL SGELEM(INDEX,ISEG,NOSEGM,iw,MXI(1,ne),
     '              NBJ(1,ne),ne,NLL(1,ne),NPL,nr,.FALSE.,.FALSE.,
     '              CSEG,DL,XP,ZCENTR,ERROR,*9999)
                ENDIF
              ENDIF
            ENDDO
          ENDDO

          INDEX=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
          DO nofibr=1,NTFIBR
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                NOSEGM=ISFIBR(iw,ne,nofibr)
                IF(NOSEGM.GT.0) THEN !fibre segment is defined
                  IF(ISEG(NOSEGM).EQ.2) THEN !..and visible
                    IF(ITYP2(nr,nx).EQ.14.OR.ITYP2(nr,nx).EQ.15) THEN
C!!! this needs checking
!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
                      CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),
     '                  NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,
     '                  NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
!old
c                      CALL YPZP(4,NBH,NEELEM,NHE,NHP(1,nr),
C     '                  NKH(1,1,1,nr),NPNODE,nr,NVHP(1,1,1,nr),nx,
C     '                  NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
                    ENDIF
                    CALL SGFIBR(INDEX,IBT,IDO,INP,ISEG,
     '                ISFIBR(iw,ne,nr),iw,
     '                MXI(1,ne),NAN,NBJ(1,ne),ne,
     '                NKHE(1,1,1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '                NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NVJE(1,1,1,ne),
     '                NW(ne,1),nx,CSEG,CURVCORRECT(1,1,1,ne),
     '                SE(1,1,ne),'AXIS',XA,XE,XG,XI,XP,ZA(1,1,1,ne),
     '                ZE,ZP,.FALSE.,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              NOSEGM=ISFIEL(iw,ne)
              IF(NOSEGM.GT.0) THEN !field segment is defined
                IF(ISEG(NOSEGM).EQ.2) THEN !..and visible
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA,XE,XP,ERROR,*9999)
                  ID_TYPE=1
                  nb=NBJ(1,ne)
                  IF(CALL_FIEL.AND..NOT.CALL_SOLV) THEN
                    CALL SGFIEL(IBT,IDO,INP,ISEG,NOSEGM,iw,MXI(1,ne),
     '                NBJ(NJ_LOC(NJL_FIEL,ID_TYPE,nr),ne),NBJ(1,ne),ne,
     '                NGAP(1,nb),NLATNE,NQNE,NQNLAT,NQS,NQXI,CSEG,
     '                'FIELD',XE,
     '                XQ,YG(1,1,ne),YQ(1,1,1,nx),
     '                XE(1,NJ_LOC(NJL_FIEL,ID_TYPE,nr)),ERROR,*9999)
                  ELSE IF(CALL_SOLV) THEN
                    CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),
     '                NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '                YP,ZA,ZP,ERROR,*9999)
                    CALL ZPZE(NBH(1,1,ne),nc,NHE(ne),NKHE(1,1,1,ne),
     '                NPF(1,1),NPNE(1,1,ne),NRE(ne),
     '                NVHE(1,1,1,ne),NW(ne,1),nx,
     '                CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '                ZE,ZP,ERROR,*9999)
                    CALL SGFIEL(IBT,IDO,INP,ISEG,NOSEGM,iw,MXI(1,ne),
     '                NBH(ID_TYPE,1,ne),NBJ(1,ne),ne,NGAP(1,nb),
     '                NLATNE,NQNE,NQNLAT,
     '                NQS,NQXI,CSEG,'DEPENDENT',XE,XQ,YG(1,1,ne),
     '                YQ(1,1,1,nx),ZE(1,ID_TYPE),ERROR,*9999)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO

          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDDO

      CALL EXITS('UPVIEW')
      RETURN
 9999 CALL ERRORS('UPVIEW',ERROR)
      CALL EXITS('UPVIEW')
      RETURN 1
      END

      
