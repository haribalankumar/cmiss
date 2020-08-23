      SUBROUTINE LINCAL(IBT,IDO,INP,JDER,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '  NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,DL,SE,XP,ERROR,*)

C#### Subroutine: LINCAL
C###  Description:
C###    LINCAL calculates global line and face parameters.
C**** If JDER=1 nodal derivs are updated to be wrt arclength

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),JDER,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNL(0:4,12,NBFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

      CALL ENTERS('LINCAL',*9999)

      CALL LINSEG(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,NKJE,NLL,NLLINE,
     '  NNL,NPL,NPNE,NPNODE,NVJE,NVJL,XP,ERROR,*9999)

      CALL LINSCA(IBT,IDO,JDER,0,NBJ,NEELEM,NEL,
     '  NLL,NLLINE,NNL,NPL,NPNE,0,NRE,NVJL,DL,SE,XP,ERROR,*9999)

C KAT 31Aug98:  Old.  Doesn't handle averaging properly.
C      DO nl=1,NLT !loop over global lines
C        TYPE=' '
C        DO N=1,NEL(0,nl) !loop over elements adjoining line nl
C          ne=NEL(N,nl)
C          DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(nr)) !loop over geom variables
C            nb=NBJ(nj,ne)
C            IF(NBI(nb).EQ.1.AND.TYPE(1:1).EQ.' ') THEN
C              !Note: for a linear element adjacent to a cubic element
C              !the derivative type is that associated with the cubic
C              TYPE='Unit derivs'
C            ELSE IF(NBI(nb).EQ.2) THEN
C              TYPE='Element derivs'
C            ELSE IF(NBI(nb).EQ.3) THEN
C              TYPE='Global derivs'
Cc cpb 5/12/96 re-adding arc length scaling
Cc cpb 8/10/94 Swaping over nbi(nb) = 4 / 5
CC            ELSE IF(NBI(nb).EQ.4) THEN
CC              TYPE='Arc-length derivs'
CC            ELSE IF(NBI(nb).EQ.5) THEN
CC              TYPE='Angle-change derivs'
C            ELSE IF(NBI(nb).EQ.4) THEN
C              TYPE='Angle-change derivs'
C            ELSE IF(NBI(nb).EQ.5) THEN
C              TYPE='Arc-length derivs'
C            ELSE IF(NBI(nb).EQ.6.OR.NBI(nb).EQ.7) THEN
C              TYPE='Ave. Arc-length derivs'
C            ENDIF
C          ENDDO
C        ENDDO
C
C        IF(DOP) THEN
CC$        call mp_setlock()
C          WRITE(OP_STRING,'('' For global line '',I6,'' nbj(1,1)='',I2,'
C     '      //''' type='',A)') nl,NBJ(1,1),TYPE
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C        ENDIF
C
C        IF(TYPE(1:11).EQ.'Unit derivs') THEN
CC ***     Set DL array derivatives to unity
C          DL(1,nl)=1.0D0
C          DL(2,nl)=1.0D0
C
C        ELSE IF(TYPE(1:14).EQ.'Element derivs') THEN
CC ***     The SE array has already been read in
C
C        ELSE IF(TYPE(1:13).EQ.'Global derivs') THEN
CC ***     Derivs have been read into DL(1,nl) and DL(2,nl) already
C
C        ELSE IF(TYPE(1:19).EQ.'Angle-change derivs') THEN
CC ***     Calculate DL array derivatives from angle change
C          CALL ANGSCA(NPL(1,0,nl),DL(1,nl),XP,ERROR,*9999)
C
C        ELSE IF(TYPE(1:17).EQ.'Arc-length derivs') THEN
CC ***     Calculate DL array derivatives from arclength
C          CALL ARCSCA(IDO,JDER,0,0,NBJ,NEL(0,nl),nl,
C     '      NPL(1,0,nl),NPNE,DL,1.0d-6,XP,ERROR,*9999)
C
C        ELSE IF(TYPE(1:22).EQ.'Ave. Arc-length derivs') THEN
CC ***     Calculate DL array derivatives from arclength
C          CALL ARCSCA(IDO,JDER,0,0,NBJ,NEL(0,nl),nl,
C     '      NPL(1,0,nl),NPNE,DL,1.0d-6,XP,ERROR,*9999)
C
C        ENDIF
C
CC ***   Calculate line segment lengths DL(3,nl)
C        CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,NPL(1,0,nl),NPNE,
C     '    DL,XP,ERROR,*9999)
C
CC ***   Transfer DL(3,nl) to DL(1,nl) & DL(2,nl) for lines which have
CC ***   linear basis - for use by nongeom. params with cubic variation
C        LINEAR=.TRUE.
C        DO nj=1,NJT
C          IF(NPL(1,nj,nl).NE.1) LINEAR=.FALSE.
C        ENDDO
C        IF(LINEAR) THEN
C          DL(1,nl)=DL(3,nl)
C          DL(2,nl)=DL(3,nl)
C        ENDIF
C        IF(DOP) THEN
CC$        call mp_setlock()
C          WRITE(OP_STRING,'('' DL(n,'',I4,''):'',3E12.4)')
C     '      nl,(DL(N,nl),N=1,3)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C        ENDIF
C      ENDDO !nl
C
C      IF(TYPE(1:14).NE.'Element derivs') THEN
CC ***   Calculate SE from DL
C        DO nb=1,NBFT
C          CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,SE,ERROR,*9999)
C        ENDDO
C      ENDIF

      CALL_LINE=.TRUE.

      CALL EXITS('LINCAL')
      RETURN
 9999 CALL ERRORS('LINCAL',ERROR)
      CALL EXITS('LINCAL')
      RETURN 1
      END


