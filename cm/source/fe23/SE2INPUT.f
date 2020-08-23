      SUBROUTINE SE2INPUT(nb,ne,ner,nl,ns,NLL,NPL,NPNE,NSB,SE,SE2,
     '  VTOE,VTOE2d,ERROR,*)

C#### Subroutine: SE2INPUT
C###  Description:
C###    SE2INPUT enters scalefactors into array SE2

C#### Variable: SE2
C###  Type: REAL*8
C###  Set_up: UPLINE and SE2INPUT
C###  Description:
C###    SE2 contains the scale-factors for two nodes along
C###    a line shared by two elements (or more). The line
C###    may have inconsitent Xi-directions and have reversed
C###    nodes for the two elements. Called from UPLINE.
C###    The two rows of SE2(c,r) contains info about respective
C###    elements.
C###    c=1,7   -> Scale-factors for node 1 of elem r
C###    c=8,14  -> Scale-factors for node 2 of elem r
C###    c=15    -> Xi direction of line in elem r. If 0, then nl is not
C###               a local line in element.
C###    c=16,17 -> Global node 1 and 2 of line for elem r
C###    c=18    -> Xi-direction in elem r for other line also in common
C###               with the other element. 0 if no other line in common.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,ne,ner,nl,ns,NLL(12,NEM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NSB(NKM,NNM,NBFM),VTOE(3,12),VTOE2d(3,4)
      REAL*8 SE(NSM,NBFM,NEM),SE2(18,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nae,nk,nn1,nn2,NP1,NP2
      LOGICAL PROGRESS

      CALL ENTERS('SE2INPUT',*9999)
      SE2(18,ner)=0
      PROGRESS=.TRUE.
      nae=0
      SE2(15,ner)=0
      NP1=NPL(2,1,nl) !is 1st node#
      NP2=NPL(3,1,nl) !is 2nd node#
      DO WHILE(nae.LE.(NLE(nb)-1).AND.PROGRESS) ! numb of lines for basis nb
        nae=nae+1

        IF(NIT(nb).EQ.3) THEN
          nn1=VTOE(1,nae)
          nn2=VTOE(2,nae)
        ELSE
          nn1=VTOE2d(1,nae)
          nn2=VTOE2d(2,nae)
        ENDIF

        IF(nl.EQ.NLL(nae,ne).OR.
     '     NP1.EQ.NPNE(nn1,nb,ne).AND.
     '     NP2.EQ.NPNE(nn2,nb,ne).OR.
     '     NP1.EQ.NPNE(nn2,nb,ne).AND.
     '     NP2.EQ.NPNE(nn1,nb,ne)) THEN
          SE2(15,ner)=VTOE(3,nae)
          SE2(16,ner)=NPNE(nn1,nb,ne)
          SE2(17,ner)=NPNE(nn2,nb,ne)
          DO nk=2,NKT(nn1,nb)
            ns=NSB(nk,nn1,nb)
            SE2(nk-1,ner)=SE(ns,nb,ne)
          ENDDO !nk
          DO nk=2,NKT(nn2,nb)
            ns=NSB(nk,nn2,nb)
            SE2(6+nk,ner)=SE(ns,nb,ne)
          ENDDO !nk
          PROGRESS=.FALSE.
        ENDIF !nl
      ENDDO !nae

      CALL EXITS('SE2INPUT')
      RETURN
 9999 CALL ERRORS('SE2INPUT',ERROR)
      CALL EXITS('SE2INPUT')
      RETURN 1
      END


