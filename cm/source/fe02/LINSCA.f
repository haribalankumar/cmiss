      SUBROUTINE LINSCA(IBT,IDO,JDER,JEST,NBJ,NEELEM,NEL,
     '  NLL,NLLINE,NNL,NPL,NPNE,nr,NRE,NVJL,DL,SE,XP,ERROR,*)

C#### Subroutine: LINSCA
C###  Description:
C###    LINSCA calculates line lengths and scale factors from global
C###    line parameters.
C###    If JEST=1 current scale factors are used as estimates.
C***  JDER=1 probably won't do what you want it to

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'linc00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  JDER,JEST,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),
     '  NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),
     '  NNL(0:4,12,NBFM),NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),
     '  nr,NRE(NEM),NVJL(4,NJM,NLM)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ITER,N,NEXTNL,nn,nb,ne,nj,nl,NL_T,no_nl,TYPE
      REAL*8 DIFF,LAST_MAX_DIFF,MAX_DIFF,OLDDL,TEMP
      LOGICAL AVERAGE,ITERATE,LINEAR

      CALL ENTERS('LINSCA',*9999)

      IF(nr.EQ.0) THEN
        NL_T=NLT
      ELSE
        NL_T=NLLINE(0,nr)
      ENDIF
      ITER=0
      ITERATE=.TRUE.

      DO WHILE(ITERATE.AND.ITER.LE.LINC_ITMAX)

        MAX_DIFF=0.0d0
C LKC 9-FEB-2001 Unneeded variable?
C        HERMITE_ELEM=.FALSE.

        IF(ITER.EQ.0) ITERATE=.FALSE. !assume for now only 1 iteration needed

        DO no_nl=1,NL_T !loop over global lines
          IF(nr.EQ.0) THEN
            nl=no_nl
          ELSE
            nl=NLLINE(no_nl,nr)
          ENDIF
C         Note that on lines, it is assumed that only one scale
C         factor type is used for all geometric variables.
C         If NBI is 1 then this may correspond to a Lagrange
C         interpolation so if there is an NBI that is not 1 the first
C         one found is chosen.
          TYPE=1
          DO N=1,NEL(0,nl) !loop over elements adjoining line nl
            ne=NEL(N,nl)
            DO nj=1,NJ_LOC(NJL_GEOM,0,NRE(ne)) !loop over geom variables
              nb=NBJ(nj,ne)
              IF(NBI(nb).NE.1) THEN
                TYPE=NBI(nb)
                GOTO 100
              ENDIF
            ENDDO
          ENDDO
 100      CONTINUE

          IF(ITER.EQ.0) THEN
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,
     '          '('' For global line '',I6,'' nbj(1,1)='',I2,'
     '          //''' type='',I1)') nl,NBJ(1,1),TYPE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

            IF(TYPE.LE.4) THEN
C             DL(1) and DL(2) are either known or easily calculated.
C             For TYPE=2,3 scale factors have already been read in.
C!!!          Let's hope line numbers haven't changed. i.e. no new elements.
              IF(TYPE.EQ.1) THEN !unit scale factors
C               Set DL array derivatives to unity
                DL(1,nl)=1.0D0
                DL(2,nl)=1.0D0
              ELSE IF(TYPE.EQ.4) THEN !Angle-change derivs
C               Calculate DL(1) and DL(2) from angle change
                CALL ANGSCA(NPL(1,0,nl),DL(1,nl),XP,ERROR,*9999)
              ENDIF
C             Calculate line segment lengths, DL(3)
              CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,NPL(1,0,nl),NPNE,
     '          NVJL(1,1,nl),DL,XP,ERROR,*9999)
C             Transfer DL(3,nl) to DL(1,nl) & DL(2,nl) for lines which have
C             linear basis - for use by nongeom. params with cubic variation
              LINEAR=.TRUE.
              DO nj=1,NJT
                IF(NPL(1,nj,nl).NE.1) LINEAR=.FALSE.
              ENDDO
              IF(LINEAR) THEN
                DL(1,nl)=DL(3,nl)
                DL(2,nl)=DL(3,nl)
              ENDIF
            ELSE IF(TYPE.EQ.5) THEN !arc length based scale factors
C             DL(1) and DL(2) depend on DL(3) but only on this line
              CALL ARCSCA(IDO,JDER,0,0,NBJ,NEL(0,nl),nl,NPL(1,0,nl),
     '          NPNE,NVJL(1,1,nl),DL,LOOSE_TOL,XP,ERROR,*9999)
            ELSE !IF(TYPE.GE.6) THEN !average arc length scale factors
C             Scale factors depend also on surrounding lines so
C             iteration is necessary.
C             ARCSCA provides an initial estimate.
              ITERATE=.TRUE.
C              HERMITE_ELEM=.TRUE.
              IF(JEST.NE.1) THEN
                CALL ARCSCA(IDO,JDER,0,0,NBJ,NEL(0,nl),nl,NPL(1,0,nl),
     '            NPNE,NVJL(1,1,nl),DL,0.1d0,XP,ERROR,*9999)
              ELSE
                OLDDL=DL(3,nl)
                CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,NPL(1,0,nl),NPNE,
     '            NVJL(1,1,nl),DL,XP,ERROR,*9999)
                DIFF=DABS(DL(3,nl)-OLDDL)/(1.0d0+OLDDL)
                IF(DIFF.GT.MAX_DIFF) THEN
                  MAX_DIFF=DIFF
                ENDIF
              ENDIF
            ENDIF

            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' DL(n,'',I4,''):'',3E12.4)')
     '          nl,(DL(N,nl),N=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF

          ELSE IF(TYPE.GE.6) THEN !and iterating
C           Update DL(1,nl),DL(2,nl) from current estimate of line
C           segment lengths DL(3,*).
            IF(JTYP2B.NE.1.OR.NPL(4,0,nl).GE.0) THEN
              DO nn=1,2
C               First find out if there is a neighbouring line
                NEXTNL=NPL(1+nn,0,nl)
                IF(NEXTNL.NE.0) THEN
                  IF(JTYP2B.EQ.1.AND.NPL(4,0,nl).NE.0) THEN !Nonstd line mapping
C                   only average if the 2 lines have consistent xi dirn
                    AVERAGE=NPL(4,0,NEXTNL).GT.0
                  ELSE !Standard line mapping
                    AVERAGE=.TRUE.
                  ENDIF
                ELSE !no neighbouring line
                  AVERAGE=.FALSE.
                ENDIF
C               Calculate the average
                IF(.NOT.AVERAGE) THEN
                  DL(nn,nl)=DL(3,nl)
                ELSE IF(NPL(4-nn,0,NEXTNL).NE.nl) THEN
C                 Next line does not know it is connected to this line.
C                 This may happen if a line connects to more than one
C                 line.  In order to keep continuity of scale factors,
C                 this scale factor here is simply set to the value of
C                 the adjacent scale factor on next line.  This value
C                 will converge to the average arc length of NEXTNL and
C                 its adjoining line.
                  DL(nn,nl)=DL(3-nn,NEXTNL)
                ELSE IF(NBI(nb).EQ.6) THEN !arithmetic mean
                  DL(nn,nl)=(DL(3,nl)+DL(3,NEXTNL))/2.0d0
                ELSE !IF(NBI(nb).EQ.7) THEN !harmonic mean
C                 ensure denominator is not zero
                  TEMP=DL(3,nl)*DL(3,NEXTNL)
                  IF(TEMP.EQ.0.0d0) THEN
                    DL(nn,nl)=0.0d0
                  ELSE
                    DL(nn,nl)=2.0d0*TEMP/(DL(3,nl)+DL(3,NEXTNL))
                  ENDIF
                ENDIF
              ENDDO !n
C             Calculate line segment length, DL(3,nl) using current
C             estimate of DL(1,nl) & DL(2,nl).
              OLDDL=DL(3,nl)
              CALL ARCLEN(IDO,NBJ,NEL(0,nl),nl,NPL(1,0,nl),NPNE,
     '          NVJL(1,1,nl),DL,XP,ERROR,*9999)
              DIFF=DABS(DL(3,nl)-OLDDL)/(1.0d0+OLDDL)
              IF(DIFF.GT.MAX_DIFF) THEN
                MAX_DIFF=DIFF
              ENDIF
              IF(JTYP2B.EQ.1.AND.NPL(4,0,nl).GT.0) THEN !non-standard mapping
C               Update equivalent mapped line
                DL(1,NPL(4,0,nl))=-1.0d0*DL(2,nl)
                DL(2,NPL(4,0,nl))=-1.0d0*DL(1,nl)
                DL(3,NPL(4,0,nl))=DL(3,nl)
              ENDIF
            ENDIF !JTYP2B.NE.1.OR.NPL(4,0,nl).GE.0
          ENDIF !ITER=0
        ENDDO !nl

        IF(ITER.EQ.0.AND.JEST.EQ.1) ITER=1
        IF(ITER.GT.0) THEN
          ITERATE=MAX_DIFF.GT.LINC_TOL
          IF(ITERATE) THEN
            IF(ITER.EQ.1) THEN
              LAST_MAX_DIFF=MAX_DIFF
            ELSE IF(MAX_DIFF.LT.LOOSE_TOL
     '        .AND.MAX_DIFF.GE.LAST_MAX_DIFF) THEN
C             seem to be at numerical limit
              ITERATE=.FALSE.
            ELSE
              LAST_MAX_DIFF=MAX_DIFF
C LKC 9-FEB-2001 Unneeded variable?
C              HERMITE_ELEM=.TRUE.
            ENDIF
          ENDIF
        ENDIF
        ITER=ITER+1
C LKC 9-FEB-2001 Unneeded variable?
C        IF(HERMITE_ELEM) ITERATE=.TRUE.
      ENDDO !iterate

      IF(ITERATE) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>WARNING!!! Iteration in LINSCA has not'
     '    //' converged.'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(14X,''Estimate of maximum error:'',D9.2,''.'')') MAX_DIFF
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IF(.NOT.((JTYP2C.EQ.1).AND.IN_REFINE))THEN
        DO nb=1,NBFT
          IF(NBI(nb).NE.2.AND.NBI(nb).NE.3) THEN
            !SE hasn't been set up from values read in from a file.
C ***       Calculate SE from DL
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,SE,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('LINSCA')
      RETURN
 9999 CALL ERRORS('LINSCA',ERROR)
      CALL EXITS('LINSCA')
      RETURN 1
      END


