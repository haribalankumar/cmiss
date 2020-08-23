      SUBROUTINE MESH_BRANCH(LD,LD_NUM,LD_TEMP,numzero,nen,np1,np2,
     &  POINT_LIMIT,COFM,FRAC,A_LIM,L_LIM,LENGTH_PARENT,
     &  min_length,XP,XP1,ZD,BRANCH,SS,ERROR,*)

C#### Subroutine: MESH_BRANCH
C###  Description:
C###    MESH_BRANCH is used in the Bifurcating-Distributive method to
C###    create a branch that runs from a defined point to some fraction
C###    along a line towards the centre of mass of a collection of
C###    random points.
C***  Created by Merryn Howatson Tawhai, February 1997

      IMPLICIT NONE
      INCLUDE 'lung00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
      INTEGER LD(NDM),LD_NUM,LD_TEMP(NDM),numzero,nen,np1,np2,
     &  POINT_LIMIT
      REAL*8 COFM(3),FRAC,A_LIM,L_LIM,LENGTH_PARENT,min_length,
     '  XP(NKM,NVM,NJM,NPM),XP1(3),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      LOGICAL BRANCH,SS
      !Local variables
      INTEGER N,NCLOSEST(100),nd,nj,nsp,NUM_CLOSEST,NUM_ND,
     '  number_of_points
      REAL*8 CLOSEST(100),DIST,L_COFM,LENGTH,LU,MIN_DIST,
     '  VECTOR(3),temp

      CALL ENTERS('MESH_BRANCH',*9999)
      
      number_of_points=0
      DO nd=1,NDT
        IF(LD(nd).EQ.nen) number_of_points=number_of_points+1
      ENDDO
      LENGTH=0.0d0
      L_COFM=0.0d0
      DO nj=1,NJT
        XP1(nj)=XP(1,1,nj,np1)+FRAC*(COFM(nj)-XP(1,1,nj,np1))
        VECTOR(nj)=COFM(nj)-XP(1,1,nj,np1)
        LENGTH=LENGTH+(XP1(nj)-XP(1,1,nj,np1))**2.0d0
        L_COFM=L_COFM+VECTOR(nj)**2.D0
      ENDDO !nj
      LENGTH=DSQRT(LENGTH)
      L_COFM=DSQRT(L_COFM)
      CALL NORMALISE2(NJT,VECTOR,temp,ERROR,*9999) !unit vector

      IF(LENGTH.GE.L_LIM)THEN !the branch will not be terminal
        BRANCH=.TRUE.
      ELSE
        BRANCH=.FALSE.
      ENDIF !LENGTH.GE.L_LIM
      IF(.NOT.BRANCH.AND.LENGTH.LT.min_length) THEN !KSB - only used for blood vessels to ensure
 !that the vessels don't go below min level solvable by Navier-Stokes solution technique.
        LENGTH=0.d0
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+(min_length/L_COFM)*
     &      (COFM(nj)-XP(1,1,nj,np1))
          LENGTH=LENGTH+(XP1(nj)-XP(1,1,nj,np1))**2.0d0
        ENDDO
        LENGTH=DSQRT(LENGTH)
      ENDIF
C*** Position branch
      MIN_DIST=1.d6
      NUM_ND=0
      NUM_CLOSEST=1
      NCLOSEST(1)=1
      CLOSEST(1)=1.d6
      DO nd=1,NDT
        nsp=LD(nd) !space # that random point belongs to
        IF(nsp.EQ.nen)THEN !random point belongs to this element space
          DIST=0.d0
          DO nj=1,NJT
            DIST=DIST+(ZD(nj,nd)-XP1(nj))**2.d0
          ENDDO !nj
          DIST=DSQRT(DIST)
          IF(DIST.LT.MIN_DIST)THEN
            MIN_DIST=DIST
          ENDIF
          NUM_ND=NUM_ND+1
          IF(.NOT.SS)THEN
            LD_TEMP(nd)=LD(nd)
            LD(nd)=0 !set to zero to remove from seed points
            numzero=numzero+1
            LD_NUM=LD_NUM-1
          ENDIF
          IF(.NOT.BRANCH)THEN !remove closest data points
            IF(DIST.LT.CLOSEST(NUM_CLOSEST))THEN !store this data point
              IF(NUM_CLOSEST.LT.POINT_LIMIT)THEN 
                NUM_CLOSEST=NUM_CLOSEST+1 !increment number of closest
              ENDIF
              CLOSEST(NUM_CLOSEST)=DIST !store distance
              NCLOSEST(NUM_CLOSEST)=nd !store data point number
            ENDIF !DIST
            CALL RSORT(NUM_CLOSEST,CLOSEST,NCLOSEST) !sort into ascending
          ENDIF !NOT.BRANCH
        ENDIF !nsp
      ENDDO !nd

C*** Check branch angle to parent
c      CALL MESH_ANGLE_CHECK(np1,np2,A_LIM,LENGTH,LU,XP,XP1,ERROR,
c     '  *9999)
      IF(LENGTH.LT.0.5d0*LENGTH_PARENT)THEN
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+VECTOR(nj)*0.5d0*LENGTH_PARENT
        ENDDO
      ENDIF        

      IF(LENGTH.LT.MIN_LENGTH)THEN
        DO nj=1,NJT
          XP1(nj)=XP(1,1,nj,np1)+VECTOR(nj)*MIN_LENGTH
        ENDDO
      ENDIF        

      IF(.NOT.BRANCH)THEN !remove the closest data points
        DO N=1,NUM_CLOSEST
          nd=NCLOSEST(N)
          LD_TEMP(nd)=LD(nd)
          LD(nd)=0
          numzero=numzero+1
          LD_NUM=LD_NUM-1
        ENDDO
      ENDIF !NOT.BRANCH

      CALL EXITS('MESH_BRANCH')
      RETURN
 9999 CALL ERRORS('MESH_BRANCH',ERROR)
      CALL EXITS('MESH_BRANCH')
      RETURN 1
      END


