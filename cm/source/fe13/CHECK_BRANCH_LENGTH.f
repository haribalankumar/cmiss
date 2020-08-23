      SUBROUTINE CHECK_BRANCH_LENGTH(np1,np2,np,XP,ERROR,*)

C#### Subroutine: CHECK_BRANCH_LENGTH
C###  Description:
C###    CHECK_BRANCH_LENGTH calculates the centre of mass of a collection of
C###    random points by averaging their coordinates.
C***  Created by Merryn Howatson Tawhai, February 1997

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INTEGER np1,np2,np
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER nj
      REAL*8 Lbranch1,Lbranch2,Lparent

      CALL ENTERS('CHECK_BRANCH_LENGTH',*9999)

      Lparent=0.d0
      Lbranch1=0.d0
      Lbranch2=0.d0
      DO nj=1,NJT
        Lparent=Lparent+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2
        Lbranch1=Lbranch1+(XP(1,1,nj,np-1)-XP(1,1,nj,np1))**2
        Lbranch2=Lbranch2+(XP(1,1,nj,np)-XP(1,1,nj,np1))**2
      ENDDO
      Lparent=DSQRT(Lparent)
      Lbranch1=DSQRT(Lbranch1)
      Lbranch2=DSQRT(Lbranch2)
      
C.....Check size
      IF(Lbranch1.LT.0.75d0*Lparent)THEN
C        DO nj=1,NJT
C          XP(1,1,nj,np-1)=XP(1,1,nj,np1)+XP(1,2,nj,np-1)*Lparent
C     &      *0.75d0
C        ENDDO
      ELSEIF(Lbranch1.GT.1.5d0*Lparent)THEN
        DO nj=1,NJT
          XP(1,1,nj,np-1)=XP(1,1,nj,np1)+XP(1,2,nj,np-1)*Lparent
     &      *1.5d0
        ENDDO
      ENDIF
      IF(Lbranch2.LT.0.75d0*Lparent)THEN
C        DO nj=1,NJT
C          XP(1,1,nj,np)=XP(1,1,nj,np1)+XP(1,2,nj,np)*Lparent
C     &      *0.75d0
C        ENDDO
      ELSEIF(Lbranch2.GT.1.5d0*Lparent)THEN
        DO nj=1,NJT
          XP(1,1,nj,np)=XP(1,1,nj,np1)+XP(1,2,nj,np)*Lparent
     &      *1.5d0
        ENDDO
      ENDIF

      CALL EXITS('CHECK_BRANCH_LENGTH')
      RETURN
 9999 CALL ERRORS('CHECK_BRANCH_LENGTH',ERROR)
      CALL EXITS('CHECK_BRANCH_LENGTH')
      RETURN 1
      END


