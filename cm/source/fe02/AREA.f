      SUBROUTINE AREA(NBJ,NBJF,NKJE,NKEF,NNF,NPF,NPNE,NPNF,NRE,
     '  NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*)

C#### Subroutine: AREA
C###  Description:
C###    AREA calculates area DF(nf) of face segments nf using
C###    Gaussian quadrature.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBJF(NJM,NFM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM)
      REAL*8 DF(NFM),PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  SF(NSM,NBFM),WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,mn,nb,ne,nef,nf,ng,nj,NKJF(NKM,NNM,NJM),nn,nr,ns,nu,
     '  NODES_REP(0:4)
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3)
      LOGICAL ALL_ZERO,FOUNDnn,ZERO_AREA

      CALL ENTERS('AREA',*9999)
      DO nf=1,NFT
        DF(nf)=0.0d0
        ne=NPF(6,nf)
        nef=NPF(8,nf)
        nr=NRE(ne)
        CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),nef,
     '    NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
     '    NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
C       Determine the number of repeated nodes for the current face
C       and their numbers
        nb=NBJF(1,nf)
        NODES_REP(0)=0
        DO nn=1,NNT(nb)
          mn=nn+1
          FOUNDnn=.FALSE.
          DO WHILE(mn.LE.NNT(nb).AND..NOT.FOUNDnn)
            IF(NPNF(nn,nb).EQ.NPNF(mn,nb)) THEN
              NODES_REP(0)=NODES_REP(0)+1
              NODES_REP(NODES_REP(0))=NPNF(nn,nb)
              FOUNDnn=.TRUE.
            ENDIF
            mn=mn+1
          ENDDO !mn
        ENDDO !nn
        ZERO_AREA=.FALSE.
        IF(NNT(nb).EQ.4.AND.NODES_REP(0).GE.2.OR.
     '     NNT(nb).EQ.3.AND.NODES_REP(0).GE.1) THEN
C         Checks need to be done on radial/angular coords for each
C         coord system to avoid calculating areas for for zero
C         area faces (do this by setting ZERO_AREA to .TRUE.)
          IF(ITYP10(nr).EQ.2.OR. !cyl polar
     '      ITYP10(nr).EQ.3) THEN !sph polar - UNCHECKED maybe check phi
C           Check if radius is zero for all repeated nodes
            nj=NJ_LOC(NJL_GEOM,1,nr) !radius
            ALL_ZERO=.TRUE.
            DO i=1,NODES_REP(0)
              IF(XP(1,1,nj,NODES_REP(i)).GT.ZERO_TOL) ALL_ZERO=.FALSE.
            ENDDO
            IF(ALL_ZERO) ZERO_AREA=.TRUE.
          ELSE IF(ITYP10(nr).EQ.4) THEN !prolate sph
C           Check if Lambda is zero for all repeated nodes
            nj=NJ_LOC(NJL_GEOM,1,nr) !lambda
            ALL_ZERO=.TRUE.
            DO i=1,NODES_REP(0)
              IF(XP(1,1,nj,NODES_REP(i)).GT.ZERO_TOL) ALL_ZERO=.FALSE.
            ENDDO
            IF(ALL_ZERO) ZERO_AREA=.TRUE.
C           Check if mu is zero for all repeated nodes
            nj=NJ_LOC(NJL_GEOM,2,nr) !mu
            ALL_ZERO=.TRUE.
            DO i=1,NODES_REP(0)
              IF(XP(1,1,nj,NODES_REP(i)).GT.ZERO_TOL) ALL_ZERO=.FALSE.
            ENDDO
C CS 3/10/00 New checks to ensure that face is not infact collapsed
C twice and is now line like
            IF((NNT(nb).EQ.4).AND.(NPF(8,nb).EQ.6).AND.
     '        (NPNF(1,nb).EQ.NPNF(3,nb)).AND.
     '        (NPNF(2,nb).EQ.NPNF(4,nb))) THEN
              ZERO_AREA=.TRUE.
            ENDIF
            IF((NNT(nb).EQ.4).AND.(NPF(8,nb).EQ.2).AND.
     '        (NPNF(1,nb).EQ.NPNF(2,nb)).AND.
     '        (NPNF(3,nb).EQ.NPNF(4,nb))) THEN
              ZERO_AREA=.TRUE.
            ENDIF
            IF(ALL_ZERO) ZERO_AREA=.TRUE.
          ELSE IF(ITYP10(nr).EQ.5) THEN !oblate sph
          ENDIF
        ENDIF
        IF(.NOT.ZERO_AREA) THEN !face area is non-zero
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Face #nf= '',I6,'' nb= '',I2,'
     '        //''' SF:'',(/10D10.3))') nf,nb,(SF(ns,nb),ns=1,NST(nb))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          CALL XPXE(NBJF(1,nf),NKJF,NPF(1,nf),NPNF,nr,NVJF,SF,XA(1,1,1),
     '      XE,XP,ERROR,*9999)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' XE:'',10D11.4)')
     '        ((XE(ns,nj),ns=1,NST(nb)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DO ng=1,NGT(nb)
            CALL XEXG(NBJF(1,nf),ng,nr,PG,XE,XG,ERROR,*9999)
            IF(DOP) THEN
C KAT 14May01: mp_setlock not OPENMP.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' XG:'',6D11.3,/(4X,6D11.3))')
     '          ((XG(nj,nu),nu=1,NUT(NBJF(nj,nf))),
     '          nj=1,NJ_LOC(NJL_GEOM,0,nr))
CC$            call mp_unsetlock()
            ENDIF
            CALL XGMG(0,2,nb,nr,DXIX,GL,GU,RG(ng),XG,ERROR,*9999)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' ng='',I2,'', RG='',D12.4)') ng,RG(ng)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            DF(nf)=DF(nf)+RG(ng)*WG(ng,nb)
          ENDDO !ng
        ENDIF !.NOT.ZERO_AREA
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' DF('',I6,'')='',D11.4)') nf,DF(nf)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDDO !nf

      CALL EXITS('AREA')
      RETURN
 9999 CALL ERRORS('AREA',ERROR)
      CALL EXITS('AREA')
      RETURN 1
      END


