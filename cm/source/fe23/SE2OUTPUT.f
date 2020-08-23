      SUBROUTINE SE2OUTPUT(nb,ne,ner,nl,ns,NLL,NSB,SE,SE2,
     '  VTOE,VTOE2d,ERROR,*)

C#### Subroutine: SE2OUTPUT
C###  Description:
C###    SE2OUTPUT updates scale-factors by putting values from
C##     SE2 into SE. Called from UPLINE.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,ne,ner,nl,ns,NLL(12,NEM),
     '  NSB(NKM,NNM,NBFM),VTOE(3,12),VTOE2d(3,4)
      REAL*8 SE(NSM,NBFM,NEM),SE2(18,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nae,nk,nn1,nn2
      LOGICAL PROGRESS

      CALL ENTERS('SE2OUTPUT',*9999)
      SE2(18,ner)=0
      PROGRESS=.TRUE.
      nae=0
      DO WHILE(nae.LE.(NLE(nb)-1).AND.PROGRESS) ! numb of lines for basis nb
        nae=nae+1
        IF(NIT(nb).EQ.3) THEN
          nn1=VTOE(1,nae)
          nn2=VTOE(2,nae)
        ELSE
          nn1=VTOE2d(1,nae)
          nn2=VTOE2d(2,nae)
        ENDIF
        IF(nl.EQ.NLL(nae,ne)) THEN
          DO nk=2,NKT(nn1,nb)
            ns=NSB(nk,nn1,nb)
            SE(ns,nb,ne)=SE2(nk-1,ner)
          ENDDO !nk
          DO nk=2,NKT(nn2,nb)
            ns=NSB(nk,nn2,nb)
            SE(ns,nb,ne)=SE2(6+nk,ner)
          ENDDO !nk
          PROGRESS=.FALSE.
        ENDIF !nl
      ENDDO !nae

      CALL EXITS('SE2OUTPUT')
      RETURN
 9999 CALL ERRORS('SE2OUTPUT',ERROR)
      CALL EXITS('SE2OUTPUT')
      RETURN 1
      END
