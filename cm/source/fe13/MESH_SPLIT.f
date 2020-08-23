      SUBROUTINE MESH_SPLIT(LD,nde_change1,nde_change2,
     &  POINT_LIMIT,ne1,ne,np1,np2,np3,COFM,XP,ZD,SS,ERROR,*)

C#### Subroutine: MESH_SPLIT
C###  Description:
C###    MESH_SPLIT divides a set of random points into two sets
C###    depending on which side of the branching plane the point
C###    lies.  Used in the Monte-Carlo method for lung trees.
C***  Created by Merryn Howatson Tawhai, February 1997
C***  Decides which side of a plane a random point is on by calculating
C***  the distance between two parallel planes: one which is defined by
C***  the parent and grandparent branch, and the other which contains a
C***  random point.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INTEGER LD(NDM),nde_change1,nde_change2,
     &  POINT_LIMIT,ne,ne1,np1,np2,np3
      REAL*8 COFM(3),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
      LOGICAL SS(2)
      !Local variables
      INTEGER DAT1,DAT2,nd,ND1_1ST,ND2_1ST,nj,nsp,NPOINTS
      REAL*8 DIST,NORML(4),P(3),Q(3),R(3)
      LOGICAL COLINEAR,COLINEAR2


      CALL ENTERS('MESH_SPLIT',*9999)

      nde_change1=0
      nde_change2=0
      DO nj=1,NJT
        R(nj)=COFM(nj) !split based on cofm and lobar bronchus
      ENDDO
      DO nj=1,NJT
        P(nj)=XP(1,1,nj,np2) !point at start of parent branch
        Q(nj)=XP(1,1,nj,np1) !point at end of parent branch
      ENDDO !nj
      CALL CHECK_COLINEAR(P,Q,R,COLINEAR,ERROR,*9999)
      IF(COLINEAR)THEN
        DO nj=1,NJT
          R(nj)=XP(1,1,nj,np3) !split based on parent and aunt
        ENDDO !nj
        CALL CHECK_COLINEAR(P,Q,R,COLINEAR2,ERROR,*9999)
        IF(COLINEAR2)THEN
          DO nj=1,3
            R(nj)=XP(1,1,nj,np2)+XP(1,2,nj,np2) !splitting plane
         ENDDO !nj
        ENDIF
        CALL CHECK_COLINEAR(P,Q,R,COLINEAR2,ERROR,*9999)
      ENDIF
      CALL PLANE_FROM_3_PTS(NORML,1,P,Q,R,ERROR,*9999) !calculate plane

      NPOINTS=0
      DAT1=0
      DAT2=0
      ND1_1ST=0
      ND2_1ST=0
c      nd_1st=0
      DO nd=1,NDT
        nsp=LD(nd) !space # that random point belongs to
        IF(nsp.EQ.ne1)THEN !random point belongs to this element space
          NPOINTS=NPOINTS+1
c          nd_1st=nd
          DIST=0.0d0
          DO nj=1,NJT
            DIST=DIST+NORML(nj)*ZD(nj,nd) !d in plane containing ZD
          ENDDO !nj
          DIST=-DIST-NORML(4) !distance between two planes
          IF(DIST.GE.0.d0)THEN
            IF(DAT1.EQ.0) ND1_1ST=nd
            DAT1=DAT1+1
            LD(nd)=ne+1
            nde_change1=ne+1
          ELSE IF(DIST.LT.0.d0)THEN
            IF(DAT2.EQ.0) ND2_1ST=nd
            DAT2=DAT2+1
            LD(nd)=ne+2
            nde_change2=ne+2
          ENDIF
        ENDIF
      ENDDO !nd
      
      IF(NPOINTS.LE.POINT_LIMIT)THEN
        SS(1)=.FALSE.
        SS(2)=.FALSE.
c        IF(NPOINTS.GT.0)THEN
c          LD_TEMP(nd_1st)=LD(nd_1st)
c          LD(nd_1st)=0
c          numzero=numzero+1
c        ENDIF
      ELSE
        SS(1)=.TRUE.
        SS(2)=.TRUE.
        
        IF(DAT1.EQ.0)THEN
          IF(ND2_1ST.NE.0)THEN    
c            LD(ND2_1ST)=ne+1
c            DAT1=DAT1+1
c            DAT2=DAT2-1
          ELSE
c            SS(1)=.FALSE.
          ENDIF
        ENDIF
        IF(DAT2.EQ.0)THEN
          IF(ND1_1ST.NE.0)THEN
c            LD(ND1_1ST)=ne+2
c            DAT2=DAT2+1
c            DAT1=DAT1-1
          ELSE
c            SS(2)=.FALSE.
          ENDIF
        ENDIF
        IF(DAT1.LT.POINT_LIMIT) SS(1)=.FALSE.
        IF(DAT2.LT.POINT_LIMIT) SS(2)=.FALSE.
      ENDIF
      
      CALL EXITS('MESH_SPLIT')
      RETURN
 9999 CALL ERRORS('MESH_SPLIT',ERROR)
      CALL EXITS('MESH_SPLIT')
      RETURN 1
      END


