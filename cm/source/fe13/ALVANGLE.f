      SUBROUTINE ALVANGLE(NPLIST,centre,centre2,XP,TRIANGULAR,ERROR,*)

C#### Subroutine: ALVANGLE
C###  Description:
C###    ALVANGLE calculates angles of vector normal to the opening
C###    plane of an alveolus.

C***  Created by Kelly Burrowes, April, 2003.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter variables
      INTEGER NPLIST(0:NPM)
      REAL*8 centre(NJT),centre2(NJT),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL TRIANGULAR
!     Local variables
      INTEGER nj,nonode,np,np1,np2,np3
      REAL*8 A,B,CX,CY,DOT_PROD,N(NJT),NORMAL(NJT),P_21(NJT),
     '  P_31(NJT),SX,SY,theta
      LOGICAL CHANGE_SIGN

      CALL ENTERS('ALVANGLE',*9999)

C... calculating the centre of this entrance
      DO nonode=1,NPLIST(0)
        np=NPLIST(nonode)
        DO nj=1,NJT
          centre2(nj)=centre2(nj)+XP(1,1,nj,np)/NPLIST(0)
        ENDDO
      ENDDO
C... calculating vector from centre2 to centre
      DO nj=1,NJT
        V(nj)=centre(nj)-centre2(nj)
      ENDDO
      CALL NORMALISE(NJT,V,ERROR,*9999) !unit vector
C... finding vector normal to entrance plane
C... should maybe calculate average entrance plane
      np1=NPLIST(1) !NPLIST(NPLIST(0)/3) !3 nodes to calculate plane
      np2=NPLIST(NPLIST(0)/2) !NPLIST((NPLIST(0)/3)*2)
      np3=NPLIST(NPLIST(0))
      
C... find an average entrance plane
C      count1=0
C      count2=0
C      count3=0
C      DO nj=1,NJT
C        XP1(nj)=0.d0
C        XP2(nj)=0.d0
C        XP3(nj)=0.d0
C      ENDDO
C      DO nonode=1,NPLIST(0)
C        np=NPLIST(nonode)
C        IF(nonode.LE.NPLIST(0)/3) THEN
C          count1=count1+1
C          DO nj=1,NJT
C            XP1(nj)=XP1(nj)+XP(1,1,nj,np)
C          ENDDO
C        ELSE IF(nonode.GT.NPLIST(0)/3.AND.nonode.LE.2*NPLIST(0)/3) THEN
C          count2=count2+1
C          DO nj=1,NJT
C            XP2(nj)=XP2(nj)+XP(1,1,nj,np)
C          ENDDO
C        ELSE
C          count3=count3+1
C          DO nj=1,NJT
C            XP3(nj)=XP3(nj)+XP(1,1,nj,np)
C          ENDDO
C        ENDIF
C      ENDDO !nonode
C      DO nj=1,NJT
C        XP1(nj)=XP1(nj)/count1
C        XP2(nj)=XP2(nj)/count2
C        XP3(nj)=XP3(nj)/count3
C      ENDDO
      A=0.d0
      B=0.d0
      DO nj=1,NJT !each direction
C       P_21(nj)=XP2(nj)-XP1(nj) !XP(1,1,nj,np2)-XP(1,1,nj,np1)
        P_21(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
c       P_31(nj)=XP3(nj)-XP1(nj) !XP(1,1,nj,np3)-XP(1,1,nj,np1)
        P_31(nj)=XP(1,1,nj,np3)-XP(1,1,nj,np1)
        A=A+P_21(nj)**2.d0 !length P_21
        B=B+P_31(nj)**2.d0 !length P_31
      ENDDO !nj
      A=DSQRT(A)
      B=DSQRT(B)
      theta=DACOS(DOT_PROD(P_21,P_31)/(A*B))
      CALL CROSS(P_31,P_21,NORMAL) !calculates normal vector (N)
      DO nj=1,NJT
        NORMAL(nj)=NORMAL(nj)/(A*B*DSIN(theta))
      ENDDO
      CHANGE_SIGN=.FALSE.
      IF(DABS(NORMAL(1)).GT.0.1d0) THEN !change direction of N
        IF(V(1).GT.0.d0.AND.NORMAL(1).LT.0.d0.OR.V(1).LT.0.d0.AND.
     '    NORMAL(1).GT.0.d0) CHANGE_SIGN=.TRUE.
      ENDIF
      IF(DABS(NORMAL(2)).GT.0.1d0) THEN
        IF(V(2).GT.0.d0.AND.NORMAL(2).LT.0.d0.OR.V(2).LT.0.d0.AND.
     '    NORMAL(2).GT.0.d0) CHANGE_SIGN=.TRUE.
      ENDIF
      IF(DABS(NORMAL(3)).GT.0.1d0) THEN
        IF(V(3).GT.0.d0.AND.NORMAL(3).LT.0.d0.OR.
     '    V(3).LT.0.d0.AND.NORMAL(3).GT.0.d0) CHANGE_SIGN=
     '    .TRUE.
      ENDIF
      IF(CHANGE_SIGN) THEN !change direction of normal vector
        DO nj=1,NJT
          NORMAL(nj)=-NORMAL(nj)
        ENDDO
      ENDIF
C... to calculate rotation angles for each alveolus
      DO nj=1,NJT
        IF(TRIANGULAR) N(nj)=NORMAL(nj)
        IF(.NOT.TRIANGULAR) N(nj)=V(nj)
      ENDDO
      CX=DSQRT(N(1)**2.d0+N(3)**2.d0)
      SX=N(2)
      ALV_ANGLE(N_ALVEOLI,1)=DATAN2(SX,CX)
      CY=N(3)/DSQRT(N(1)**2.d0+N(3)**2.d0)
      SY=N(1)/DSQRT(N(1)**2.d0+N(3)**2.d0)
      ALV_ANGLE(N_ALVEOLI,2)=DATAN2(SY,CY)
      ALV_ANGLE(N_ALVEOLI,3)=0.d0
      
      CALL EXITS('ALVANGLE')
      RETURN
 9999 CALL ERRORS('ALVANGLE',ERROR)
      CALL EXITS('ALVANGLE')
      RETURN 1
      END
      
      
