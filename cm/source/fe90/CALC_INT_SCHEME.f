      SUBROUTINE CALC_INT_SCHEME(NBJ,NLL,NNSPLIT,NPB,NPNE,NPNODE,nr,
     '  NW,DL,XP,ERROR,*)

C#### Subroutine: CALC_INT_SCHEME
C###  Description:
C###    CALC_INT_SCHEME determines the integration scheme to use in the
C###    current element for each node in the boundary element region
C###    and sets up NPB accordingly.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM),NLL(12),NNSPLIT(*),NPB(0:NP_R_M,5),
     '  NPNE(NNM,NBFM),NPNODE(0:NP_R_M),nr,NW
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER intscheme,nb,nj,NNMIN,nonode,np
      REAL*8 MINDIST,XPFP(3)

      CALL ENTERS('CALC_INT_SCHEME',*9999)

      nb=NBJ(1) !Geometric basis function number
      DO intscheme=1,5
        NPB(0,intscheme)=0
      ENDDO
      DO nonode=1,NPNODE(0)
        np=NPNODE(nonode)
        DO nj=1,NJT
          XPFP(nj)=XP(1,1,nj,np)
        ENDDO !nj
        CALL DIST(intscheme,nb,NLL,nnmin,NPNE,nr,NW,DL,MINDIST,
     '    XP,XPFP,ERROR,*9999)
        NPB(0,intscheme)=NPB(0,intscheme)+1
        NPB(NPB(0,intscheme),intscheme)=np
        IF(intscheme.EQ.1) THEN !element splitting
          NNSPLIT(NPB(0,1))=nnmin
        ENDIF
      ENDDO !nonode
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CALC_INT_SCHEME_1)
        WRITE(OP_STRING,'(/'' Integration scheme node lists:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO intscheme=1,5
          WRITE(OP_STRING,'('' Integration scheme '',I1,'', Number of '
     '     //'nodes in scheme '',I5)') intscheme,NPB(0,intscheme)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Nodes: '',10I6,/:(8X,10I6))')
     '      (NPB(nonode,intscheme),nonode=1,NPB(0,intscheme))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          IF(intscheme.EQ.1) THEN
            WRITE(OP_STRING,'('' Local node numbers: '',8I6,'
     '        //'/:(21X,8I6))') (NNSPLIT(nonode),nonode=1,NPB(0,1))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !intscheme
CC$OMP END CRITICAL(CALC_INT_SCHEME_1)
      ENDIF

      CALL EXITS('CALC_INT_SCHEME')
      RETURN
9999  CALL ERRORS('CALC_INT_SCHEME',ERROR)
      CALL EXITS('CALC_INT_SCHEME')
      RETURN 1
      END


