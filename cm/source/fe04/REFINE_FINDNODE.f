      SUBROUTINE REFINE_FINDNODE(IBT,IDO,INP,nb,NBJ,ne,NKJ,NODE,
     '  NPLIST,NPNODE,NPTNEW,NPTOLD,nr,NUNK,NVJP,XE,XI,XP,EXISTS,
     '  ERROR,*)

C#### Subroutine: REFINE_FINDNODE
C###  Description:
C###    REFINE_FINDNODE returns a node number. It takes a XI location
C###    within an element and works out if there is already a node at
C###    that XI location. If so it will return that node number
C###    otherwise it will return a new node number.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),nb,NBJ(NJM,NEM),ne,NKJ(NJM,NPM),NODE,
     '  NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),NPTNEW,NPTOLD,nr,
     '  NUNK(NKM,NJM,NPM),NVJP(NJM,NPM)
      REAL*8 XE(NSM,NJM),XI(3),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL EXISTS
!     Local Variables
      INTEGER nb1,ni,nj,njj1,njj2,nonode,np,NP_TEST,nr2
      REAL*8 DIFF,LLOOSE_TOL,PXI,XS(3),XTEMP
      LOGICAL FOUND,MATCH

      CALL ENTERS('REFINE_FINDNODE',*9999)

      LLOOSE_TOL=DSQRT(LOOSE_TOL)

C*** Find coordinates of proposed new node position
      DO nj=1,NJT
        nb1=NBJ(nj,ne)
        XS(nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),
     '    nb1,1,XI,XE(1,nj))
        IF((nj.EQ.2.AND.ITYP10(nr).GE.2).OR.
     '    (nj.EQ.3.AND.ITYP10(nr).GE.3)) THEN
C KAT 18Oct00:  XPXE often modifies nodal angles usually by adding 2pi.
C               Here we make a half hearted attempt to keep them
C               within [0,2pi).
          XTEMP=XS(nj)-2d0*PI
          IF(XTEMP.GE.0d0) XS(nj)=XTEMP
        ENDIF
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Proposed node:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' XI:'',3(1X,D12.5))') (XI(ni),ni=1,NIT(nb))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' XS:'',3(1X,D12.5))') (XS(nj),nj=1,NJT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF


C*** Search for existing node to see if this can be used instead
C***  of creating a new node

      EXISTS=.FALSE.

C KAT 22Jul99: Trying to speed this up
      np=NPTNEW
      DO WHILE(.NOT.EXISTS.AND.np.GT.0)

C LKC 6-APR-2001 This assumes that we have continuous
C  nodes starting from 1. Adding an IF statement
C  to by pass unitialised nodes
C!!! rgb initialising nj
        nj=1
        IF(NVJP(nj,np).GT.0) THEN

          MATCH=.TRUE.

          nj=1
          DO WHILE(MATCH.AND.nj.LE.NJT)
            DIFF=DABS(XS(nj)-XP(1,1,nj,np)) !all vers should be at same place
            IF(DIFF.GT.(LLOOSE_TOL*DABS(XS(nj))).AND.DIFF.GT.LOOSE_TOL)
     '        THEN
              MATCH=.FALSE. !unless the follow special cases
              IF(ITYP10(nr).NE.1) THEN
C             angles of 0 and 2pi are the same and the angle may
C             have no affect on the node position.
                IF((ITYP10(nr).EQ.2.AND.nj.EQ.2).OR.
     '            (ITYP10(nr).EQ.3.AND.nj.GE.2)) THEN
                  MATCH=DABS(2d0*PI-DIFF).LT.LOOSE_TOL.OR.
     '              DABS(XS(1)).LT.LOOSE_TOL
                ELSE IF(ITYP10(nr).EQ.4.AND.nj.EQ.3) THEN
                  MATCH=DABS(2d0*PI-DIFF).LT.LOOSE_TOL.OR.
     '              DABS(XS(2)).LT.LOOSE_TOL.OR.
     '              DABS(XS(2)-PI).LT.LOOSE_TOL
                ENDIF
              ENDIF
            ENDIF
            nj=nj+1
          ENDDO !nj
C      np=1
C      DO WHILE(.NOT.EXISTS.AND.np.LE.NPTNEW)
C        MATCH=.TRUE.
C        DO nj=1,NJT
C          MATCHNJ=.FALSE.
C          DO nv=1,NVJP(nj,np)
C            DIFF=DABS(XS(nj)-XP(1,nv,nj,np))
C            IF((DIFF.LT.(LLOOSE_TOL*DABS(XS(nj))).OR.
C     '        DIFF.LT.LOOSE_TOL).OR.
C     '        (ITYP10(nr).EQ.2.AND.nj.EQ.2.AND.
C     '        (DABS(XP(1,nv,1,np)).LT.LOOSE_TOL.OR.
C     '        DABS(DIFF-2.0d0*PI).LT.LOOSE_TOL)).OR.
C     '        (ITYP10(nr).EQ.3.AND.nj.EQ.2.AND.
C     '        (DABS(XP(1,nv,1,np)).LT.LOOSE_TOL.OR.
C     '        DABS(DIFF-2.0d0*PI).LT.LOOSE_TOL)).OR.
C     '        (ITYP10(nr).EQ.3.AND.nj.EQ.3.AND.
C     '        (DABS(XP(1,nv,1,np)).LT.LOOSE_TOL.OR.
C     '        DABS(DIFF-2.0d0*PI).LT.LOOSE_TOL)).OR.
C     '        (ITYP10(nr).EQ.4.AND.nj.EQ.3.AND.
C     '        (DABS(XP(1,nv,2,np)).LT.LOOSE_TOL.OR.
C     '        DABS(XP(1,nv,2,np)-PI).LT.LOOSE_TOL.OR.
C     '        DABS(DIFF-2.0d0*PI).LT.LOOSE_TOL))) THEN
C              MATCHNJ=.TRUE.
C            ENDIF
C          ENDDO !nv
C          MATCH=MATCH.AND.MATCHNJ
C        ENDDO !nj

          IF(MATCH) THEN
            EXISTS=.TRUE.
            NODE=np
          ELSE
            np=np-1
          ENDIF !match
        ELSE
          np=np-1
        ENDIF !NVJP if node exists
      ENDDO !not found


      IF(.NOT.EXISTS) THEN !no existing node found
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Proposed node does not '
     '      //'exist'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        NPTNEW=NPTNEW+1
        NPNODE(0,0)=NPNODE(0,0)+1
        NPNODE(0,nr)=NPNODE(0,nr)+1
        IF(NPNODE(0,0).GT.NPM) THEN
          ERROR='>>Increase NPM'
          GOTO 9999
        ELSE IF(NPNODE(0,nr).GT.NP_R_M) THEN
          ERROR='>>Increase NP_R_M'
          GOTO 9999
        ELSE
          NPNODE(NPNODE(0,nr),nr)=NPTNEW
          NODE=NPTNEW
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(NODE,.TRUE.,ERROR,*9999)
          NPT(nr)=NPTNEW
          NPT( 0)=NPTNEW
          NPTOLD=NPNODE(NPNODE(0,nr)-1,nr)
          NPLIST(0)=NPLIST(0)+1
          NPLIST(NPLIST(0))=NODE
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' New node ='',I5)') NODE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
C KAT 22Jul99: Giving the node a position
C       Initialize for all nj
        DO nj=1,NJT
          NVJP(nj,NODE)=1
          XP(1,1,nj,NODE)=XS(nj)
          NKJ(nj,NODE)=1
          NUNK(1,nj,NODE)=1
        ENDDO
        DO njj1=2,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            NKJ(nj,NODE)=0
            NVJP(nj,NODE)=0
          ENDDO !njj2
        ENDDO !njj1

      ELSE !node found, but possible only in another region
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' Proposed node does exist'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        FOUND=.FALSE.
        nonode=1
        DO WHILE((nonode.LE.NPNODE(0,nr)).AND.(.NOT.FOUND))
          NP_TEST=NPNODE(nonode,nr)
          IF(NP_TEST.eq.np)THEN !Node np is in current region
            FOUND=.TRUE.
          ELSE
            nonode=nonode+1
          ENDIF
        ENDDO
        IF(.NOT.FOUND)THEN !Update arrays to include np
          NPNODE(0,nr)=NPNODE(0,nr)+1
          NPNODE(NPNODE(0,nr),nr)=np
          IF(np.GT.NPT(nr))THEN
            NPT(nr)=np
          ENDIF
C         Node exists but in another region. Update current
C         region arrays. Firstly find other region number of
C         node np
          FOUND=.FALSE.
          nr2=1
          DO WHILE ((nr2.LE.NRT).AND.(.NOT.FOUND))
            IF(nr2.ne.nr)THEN
              nonode=1
              DO WHILE ((nonode.LE.NPNODE(0,nr2)).AND.
     '          (.NOT.FOUND))
                NP_TEST=NPNODE(nonode,nr2)
                IF(np.EQ.NP_TEST) THEN
                  FOUND=.TRUE.
                ELSE
                  nonode=nonode+1
                ENDIF
              ENDDO
              IF(.NOT.FOUND) nr2=nr2+1
            ELSE
              nr2=nr2+1
            ENDIF
          ENDDO
          IF(nr2.GT.NRT)THEN
            ERROR='>>Can not find node in other regions?'
            GOTO 9999
          ENDIF
        ENDIF !.not.found
      ENDIF

      CALL EXITS('REFINE_FINDNODE')
      RETURN
 9999 CALL ERRORS('REFINE_FINDNODE',ERROR)
      CALL EXITS('REFINE_FINDNODE')
      RETURN 1
      END


