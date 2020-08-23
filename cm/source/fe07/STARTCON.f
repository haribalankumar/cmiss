      SUBROUTINE STARTCON(NPNODE,SURF,XP,ZP,ERROR,*)

C#### Subroutine: STARTCON
C###  Description:
C###    STARTCON calculates solution parameters along initial solution
C###    front from given starting condition and stores them in arrays
C###    SURF(1,IX,IS).  This subroutine must be amended whenever ics
C###    are changed.
C###    Calls Functions A_COEFF,B_COEFF,C_COEFF,F_COEFF to calculate
C###    PDE coefficients A,B,C,F at each solution point.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM)
      REAL*8 SURF(100,11,9),XP(NKM,NVM,NJM,NPM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nc,nonode,np
      REAL*8 A_COEFF,B_COEFF,C_COEFF,F_COEFF,QI,RI,TI,UI,XI

      CALL ENTERS('STARTCON',*9999)

      nc=1 !temporary cpb 24/11/94

      DO nonode=1,NPNODE(0,1)
        np=NPNODE(nonode,1)
        XI=XP(1,1,1,np)
        SURF(1,nonode,7)=XI
        TI=0.0D0
        SURF(1,nonode,8)=TI
        RI=1.0D0
        SURF(1,nonode,5)=RI
        QI=1.0D0
        SURF(1,nonode,6)=QI
        UI=ZP(1,1,1,np,nc)
        SURF(1,nonode,9)=UI
        SURF(1,nonode,1)=A_COEFF()
        SURF(1,nonode,2)=B_COEFF(XI)
        SURF(1,nonode,3)=C_COEFF(XI)
        SURF(1,nonode,4)=F_COEFF()
      ENDDO

      CALL EXITS('STARTCON')
      RETURN
 9999 CALL ERRORS('STARTCON',ERROR)
      CALL EXITS('STARTCON')
      RETURN 1
      END


