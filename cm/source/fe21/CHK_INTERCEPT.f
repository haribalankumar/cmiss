      SUBROUTINE CHK_INTERCEPT(NPNODE,nr1,nr2,PENALTY,XP,ERROR,*)

C#### Subroutine: CHK_INTERCEPT
C###  Description:
C###    CHK_INTERCEPT checks that the nodes of one region do not cross
C###    over the mesh of a neighbouring region. Finds closest node in
C###    neighbouring region and checks outward normal of outer
C###    region with vector from "inside" node to "outside" node.
C**** Created: J.Crocombe 13/5/97

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM),nr1,nr2
      REAL*8 PENALTY,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj,nonode,nonode2,np,np2,min_np
      REAL*8 CHECK,DISTANCE,OUTNORM(3),min_dist,PROJ_VECT(3)

      CALL ENTERS('CHK_INTERCEPT',*9999)

      PENALTY=0.0d0
      min_dist=0.0d0

      DO nonode=1,NPNODE(0,nr2)
        np=NPNODE(nonode,nr2)
        min_np=0
        DO nonode2=1,NPNODE(0,nr1)
          np2=NPNODE(nonode2,nr1)
          DISTANCE=0.0d0
          DO nj=1,NJT
            DISTANCE=DISTANCE+(XP(1,1,nj,np2)-XP(1,1,nj,np))**2
          ENDDO
          DISTANCE=DSQRT(DISTANCE)
          IF(DISTANCE.LT.min_dist.OR.min_np.EQ.0) THEN
            min_np=np2
            min_dist=DISTANCE
          ENDIF
        ENDDO
        DO nj=1,NJT
          PROJ_VECT(nj)=XP(1,1,nj,min_np)-XP(1,1,nj,np)
        ENDDO !nj
        OUTNORM(1)=XP(2,1,2,min_np)*XP(3,1,3,min_np)-
     '                XP(2,1,3,min_np)*XP(3,1,2,min_np)
        OUTNORM(2)=XP(2,1,3,min_np)*XP(3,1,1,min_np)-
     '                XP(2,1,1,min_np)*XP(3,1,3,min_np)
        OUTNORM(3)=XP(2,1,1,min_np)*XP(3,1,2,min_np)-
     '                XP(2,1,2,min_np)*XP(3,1,1,min_np)
        CHECK=0.0d0
        DO nj=1,NJT
          CHECK=CHECK+OUTNORM(nj)*PROJ_VECT(nj)
        ENDDO !nj
        IF(CHECK.LT.-zero_tol) THEN
            !cross over occurs or too close to other surface - penalise
          PENALTY=PENALTY+min_dist+1
        ELSEIF(min_dist.LT.0.5d0) THEN
          PENALTY=PENALTY+min_dist
        ENDIF
      ENDDO

 9999 CALL ERRORS('CHK_INTERCEPT',ERROR)
      CALL EXITS('CHK_INTERCEPT')
      RETURN 1
      END


