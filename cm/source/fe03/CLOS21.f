      SUBROUTINE CLOS21(IBT,IDO,INP,IT,ITMAX,NBJ,SQ,XE,XI,ZD,
     '  INELEM,ERROR,*)

C#### Subroutine: CLOS21
C###  Description:
C###    CLOS21 finds the XI-coordinates at the closest approach of a 2D
C###    element to a data point with coordinates XD using a modified
C###    Newton algorithm.
C###    If INELEM is true then the closest point within the element is
C###    returned; if false then the projection from the element to the
C###    coordinate must be orthogonal.  If such a projection can't be
C###    obtained the xi position returned is outside the element.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 2-D elements in rectangular cartesian coords.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IT,ITMAX,NBJ(NJM)
      REAL*8 SQ,XE(NSM,NJM),XI(3),ZD(NJM)
      LOGICAL INELEM
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER BOUND(2),it2,MI,nb,ni,nifix,nj
      REAL*8 DELTA,DET,D2SQV2,D2SQVW2,D2SQXI(2,2),D2ZXI(3,2,2),DSQXI(2),
     '  DSQXI1,DSQXI2,DSQV,DSQVW,DZ(3),DZXI(3,2),EVMIN,EVMAX,H(2),
     '  MU,PXI,SQLIN,SQDIFF,SQDPRED,TEMP,TEMP1,TEMP2,TOL,
     '  TOL2,V(2),V1,V2,VMAX,W,XILIN(2),Z(3)
      CHARACTER FORMAT*200
      LOGICAL CONVERGED,ENFORCE(2),FREE,NEWTON
      DATA VMAX /1.0d0/

      CALL ENTERS('CLOS21',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(A)') ' >CLOS21 2-D rectangular cartesian'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DELTA=VMAX/4.0d0
      TOL=5.0d0*LOOSE_TOL !must be > sqrt(eps) or SQLIN<=SQ check may not work
      TOL2=TOL**2
      SQ=0.0d0
      DO nj=1,NJT
        nb=NBJ(nj)
        Z(nj)=
     '    PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,XE(1,nj))
        DZ(nj)=Z(nj)-ZD(nj)
        SQ=SQ+DZ(nj)**2
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        FORMAT='(/'' ZD(nj)='',3(E12.6,4X),/'' XI(ni)='',2(E12.6,4X),'
     '       //'/''  Z(nj)='',3(E12.6,4X),/''     SQ='',E12.6)'
        WRITE(OP_STRING,FORMAT) (ZD(nj),nj=1,NJT),(XI(ni),ni=1,2),
     '    (Z(nj),nj=1,NJT),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IT=0
      CONVERGED=.FALSE.
      DO WHILE(.NOT.CONVERGED.AND.IT.LT.ITMAX)
        DSQXI(1)=0.0d0
        DSQXI(2)=0.0d0
C        JJT(1,1)=0.0d0
C        JJT(1,2)=0.0d0
C        JJT(2,2)=0.0d0
        DO nj=1,NJT
          nb=NBJ(nj)
          DZXI(nj,1)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,XE(1,nj))
          DZXI(nj,2)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,XI,XE(1,nj))
          DSQXI(1)=DSQXI(1)+DZXI(nj,1)*DZ(nj)
          DSQXI(2)=DSQXI(2)+DZXI(nj,2)*DZ(nj)
C          JJT(1,1)=JJT(1,1)+DZXI(nj,1)*DZXI(nj,1)
C          JJT(1,2)=JJT(1,2)+DZXI(nj,1)*DZXI(nj,2)
C          JJT(2,2)=JJT(2,2)+DZXI(nj,2)*DZXI(nj,2)
        ENDDO
        DO ni=1,2
          IF(XI(ni).EQ.0.0d0) THEN
            BOUND(ni)=1
            ENFORCE(ni)=DSQXI(ni).GE.0.0d0
          ELSE IF(XI(ni).EQ.1.0d0) THEN
            BOUND(ni)=-1
            ENFORCE(ni)=DSQXI(ni).LE.0.0d0
          ELSE
            BOUND(ni)=0
            ENFORCE(ni)=.FALSE.
          ENDIF
        ENDDO
        IF(ENFORCE(1).AND.ENFORCE(2)) GO TO 9998
        D2SQXI(1,1)=0.0d0
        D2SQXI(1,2)=0.0d0
        D2SQXI(2,2)=0.0d0
        DO nj=1,NJT
          nb=NBJ(nj)
          D2ZXI(nj,1,1)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,3,XI,XE(1,nj))
          D2ZXI(nj,1,2)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,6,XI,XE(1,nj))
          D2ZXI(nj,2,2)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,5,XI,XE(1,nj))
          D2SQXI(1,1)=
     '      D2SQXI(1,1)+DZXI(nj,1)*DZXI(nj,1)+D2ZXI(nj,1,1)*DZ(nj)
          D2SQXI(1,2)=
     '      D2SQXI(1,2)+DZXI(nj,1)*DZXI(nj,2)+D2ZXI(nj,1,2)*DZ(nj)
          D2SQXI(2,2)=
     '      D2SQXI(2,2)+DZXI(nj,2)*DZXI(nj,2)+D2ZXI(nj,2,2)*DZ(nj)
        ENDDO
c cpb 15/9/95 Rolling loops back up as they cause floating exceptions
C from uninitialised variables
C        DSQXI(1)=DZXI(1,1)*DZ(1)+DZXI(2,1)*DZ(2)+DZXI(3,1)*DZ(3)
C        DSQXI(2)=DZXI(1,2)*DZ(1)+DZXI(2,2)*DZ(2)+DZXI(3,2)*DZ(3)
C        D2SQXI(1,1)=DZXI(1,1)*DZXI(1,1)+D2ZXI(1,1,1)*DZ(1)
C     '             +DZXI(2,1)*DZXI(2,1)+D2ZXI(2,1,1)*DZ(2)
C     '             +DZXI(3,1)*DZXI(3,1)+D2ZXI(3,1,1)*DZ(3)
C        D2SQXI(1,2)=DZXI(1,1)*DZXI(1,2)+D2ZXI(1,1,2)*DZ(1)
C     '             +DZXI(2,1)*DZXI(2,2)+D2ZXI(2,1,2)*DZ(2)
C     '             +DZXI(3,1)*DZXI(3,2)+D2ZXI(3,1,2)*DZ(3)
C        D2SQXI(2,2)=DZXI(1,2)*DZXI(1,2)+D2ZXI(1,2,2)*DZ(1)
C     '             +DZXI(2,2)*DZXI(2,2)+D2ZXI(2,2,2)*DZ(2)
C     '             +DZXI(3,2)*DZXI(3,2)+D2ZXI(3,2,2)*DZ(3)

C KAT 20Nov98: Changing to model trust region approach
C        DET=D2SQXI(1,1)*D2SQXI(2,2)-D2SQXI(1,2)**2
C        DSQXI2=DSQXI(1)**2+DSQXI(2)**2
C        NEWTON=.FALSE.
C        IF(DET.GT.0.0d0.AND.D2SQXI(2,2).GT.0.0d0) THEN !try Newton step
C          V(1)=-(D2SQXI(2,2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
C          V(2)=-(DSQXI(2)+D2SQXI(1,2)*V(1))/D2SQXI(2,2)
C          V2=V(1)**2+V(2)**2
C          DSQV=DSQXI(1)*V(1)+DSQXI(2)*V(2)
CC         This checks that numerical errors in the Hessian have not
CC         prevented the step direction being a descent direction.
C          DELTA=DSQXI2*V2*TOL2
C          IF(DSQV**2.GT.DELTA) NEWTON=.TRUE. !have good step direction
C        ENDIF !+ve def.
C        IF(.NOT.NEWTON) THEN
C          FINDSTEP=.TRUE. !look for a search direction
CC KAT 19Nov98: Assume flat surface if Hessian is not +ve definite
CC              This chooses a more optimal direction than simple
CC              steepest descent.
C          DET=JJT(1,1)*JJT(2,2)-JJT(1,2)**2
CC KAT 19Nov98: Steepest descent is now only used if dz/dxi*dz/dxi(trans) is
CC         for some reason singular.  This may happen if there is no
CC         spatial change in some xi direction.
C          IF(DET.GT.0.0d0.AND.JJT(2,2).GT.0.0d0) THEN
C            V(1)=-(JJT(2,2)*DSQXI(1)-JJT(1,2)*DSQXI(2))/DET
C            V(2)=-(DSQXI(2)+JJT(1,2)*V(1))/JJT(2,2)
C            V2=V(1)**2+V(2)**2
C            DSQV=DSQXI(1)*V(1)+DSQXI(2)*V(2)
CC           This checks that numerical errors in the Hessian have not
CC           prevented the step direction being a descent direction.
C            DELTA=DSQXI2*V2*TOL2
CC            IF(DSQV**2.GT.DELTA) FINDSTEP=.FALSE. !have good step direction
C          ENDIF !+ve def.
CCC KAT 18Nov98: Use steepest descent if Hessian is not +ve definite
CCC        IF(DABS(DET).LE.TOL.OR.D2SQXI(2,2).LE.TOL) THEN
C          IF(FINDSTEP) THEN
C            V(1)=-DSQXI(1)
C            V(2)=-DSQXI(2)
C            V2=V(1)**2+V(2)**2
CC KAT 19Nov98: Using maximum step length
C            IF(V2.NE.0.0d0) THEN
C              V1=DSQRT(V2)
C              TEMP=VMAX/V1
C              V(1)=V(1)*TEMP
C              V(2)=V(2)*TEMP
C              V2=VMAX2
C              DSQV=-V1*VMAX
C            ELSE
C              DSQV=0.0d0
C            ENDIF
C          ENDIF
C        ENDIF
CC       Check feasible and limit step size
C        W=1.0d0
C        FREE=.TRUE.
C        IF(INELEM) THEN
C          DO ni=1,2
C            IF(BOUND(ni).NE.0.AND
C     '        .(BOUND(ni).GT.0.EQV.V(ni).LT.0.0d0)) THEN
C              FREE=.FALSE.
C              nifix=ni
C            ENDIF
C          ENDDO
C        ENDIF !INELEM
C        IF(FREE) THEN
C          IF(V2.GT.VMAX2) W=VMAX/DSQRT(V2)
C          IF(.NOT.NEWTON) THEN
CC           don't step further than estimate of minimum along line
C            D2SQV2=V(1)*(V(1)*D2SQXI(1,1)+2.0d0*V(2)*D2SQXI(1,2))
C     '        +V(2)**2*D2SQXI(2,2)
C            IF(D2SQV2.GT.0.0D0) THEN !minimum exists
C              TEMP=-DSQV/D2SQV2 !minimum position
C              IF(TEMP.LT.W) W=TEMP
C            ENDIF
C          ENDIF
C        ELSE
C          IF(ENFORCE(2)) THEN !gradient suggests must use ni=1
C            ni=1
C            V(2)=0.0d0
C          ELSE IF(ENFORCE(2)) THEN !gradient suggests must use ni=2
C            ni=2
C            V(1)=0.0d0
C          ELSE !Gradient points into element.
CC           Fix the direction that prevented the step.
CC           P.d. Hessian guarantees there is only one of these
CC           for this type of gradient.
C            ni=3-nifix
C            V(nifix)=0.0d0
C          ENDIF
C          IF(D2SQXI(ni,ni).GT.0.0D0) THEN !Newton step good
C            V(ni)=-DSQXI(ni)/D2SQXI(ni,ni)
C          ELSE IF(JJT(ni,ni).GT.0.0D0) THEN !Assume flat
C            V(ni)=-DSQXI(ni)/JJT(ni,ni)
C          ELSE !a descent direction
C            V(ni)=-DSIGN(VMAX,DSQXI(ni))
C          ENDIF
C          IF(V(ni).GT.VMAX) W=VMAX/V(ni)
C          DSQV=DSQXI(ni)*V(ni)
C          V2=V(ni)**2
C        ENDIF !free
CC        DSQXIV=DSQXI(1)*V(1)+DSQXI(2)*V(2)    !new AAY 4 Sept 94
C        IF(DOP) THEN
CC$        call mp_setlock()
C          FORMAT='(''    IT1='',I4,10X,''DSQXI(ni)='',2(E12.6,4X),'
C     '      //'/18X,''D2SQXI(MI,ni)='',2(E12.6,4X),/32X,2(E12.6,4X),'
C     '      //'/26X,''V(ni)='',2(E12.6,4X))'
C          WRITE(OP_STRING,FORMAT) it1,(DSQXI(ni),ni=1,2),
C     '      ((D2SQXI(MI,ni),MI=1,2),ni=1,2),(V(ni),ni=1,2)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
C        ENDIF
CC ***   Performs line search: fit a quadratic with the value & deriv at
CC ***   current point IT1 & value at point IT2, then choose minimum point
C        XILIN(1)=XI(1)+V(1)*W
C        XILIN(2)=XI(2)+V(2)*W
C        IF(INELEM) THEN
C          DO ni=1,2
C            IF(XILIN(ni).LT.0.0d0) THEN
C              XILIN(ni)=0.0d0
C              W=XI(ni)/-V(ni)
C              XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
C            ELSE IF(XILIN(ni).GT.1.0d0) THEN
C              XILIN(ni)=1.0d0
C              W=(1.0d0-XI(ni))/V(ni)
C              XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
C            ENDIF
C          ENDDO !ni
C        ENDIF
C        SQLIN=0.d0
C        DO nj=1,NJT
C          nb=NBJ(nj)
C          Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
C     '      XILIN,XE(1,nj))
C          DZ(nj)=Z(nj)-ZD(nj)
C          SQLIN=SQLIN+DZ(nj)**2
C        ENDDO
CC KAT 20Nov98: testing before line search as small SQLIN may always be >
CC              SQ if DSQXI is small (-> V is small).
CC KAT 18Nov98: using v norm as there may be small steps at a bound
C        IF(V2.LT.TOL2.AND.DABS(SQ-SQLIN)/(1.d0+SQ).LE.TOL) THEN
C          SQ=SQLIN
C          XI(1)=XILIN(1)
C          XI(2)=XILIN(2)
C          GO TO 9998
C        ENDIF
C        DO it2=1,ITMAX
C          IF(DOP) THEN
CC$          call mp_setlock()
C            FORMAT='(8X,''IT2='',I4,14X,''W='',E12.6/21X,'
C     '        //'''XILIN(ni)='','//'2(E12.6,4X)/26X,''Z(nj)='','
C     '        //'3(E12.6,4X)/26X,''SQLIN='',E12.6)'
C            WRITE(OP_STRING,FORMAT) it2,W,(XILIN(ni),ni=1,2),
C     '        (Z(nj),nj=1,NJT),SQLIN
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
C          ENDIF
CC KAT 20Nov98: Requesting more from the line search
CC          IF(SQLIN.LE.SQ) GO TO 5
CC KAT 17Nov98: I think this is better as SUM is always -ve here
CC         Added factor of half as this minimizes quadratic
C          DSQVW=DSQV*W
C          IF(SQLIN-SQ.LE.0.25d0*DSQVW) GO TO 5
C          W=0.5d0*W*DSQVW/(DSQVW+SQ-SQLIN)
CC          SUM=DSQXIV*W+SQ-SQLIN
CC          IF(DABS(SUM).GT.1.0d-5) THEN
CC            W=(DSQXIV*W*W)/SUM
CC          ELSE
CC            W=W/2.0d0 !AAY 4 Sept 94
CC          ENDIF
C          XILIN(1)=XI(1)+V(1)*W
C          XILIN(2)=XI(2)+V(2)*W
C          SQLIN=0.d0
C          DO nj=1,NJT
C            nb=NBJ(nj)
C            Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
C     '        XILIN,XE(1,nj))
C            DZ(nj)=Z(nj)-ZD(nj)
C            SQLIN=SQLIN+DZ(nj)**2
C          ENDDO
CC KAT 17Nov98: Uninitialized variables
CC          DZ(1)=Z(1)-ZD(1)
CC          DZ(2)=Z(2)-ZD(2)
CC          DZ(3)=Z(3)-ZD(3)
CC          SQLIN=DZ(1)**2+DZ(2)**2+DZ(3)**2
C        ENDDO

C KAT 20Nov98: A model trust region approach.
C       A Newton step is taken if the condition of the Hessian
C       guarantees that the step will be within the trust region.
C       Otherwise the Hessian is shifted towrds a diagonal matrix to
C       shift the step towards steepest descent.  Usually it is a much
C       better direction than steepest descent.  I think it is close to
C       the best direction in the trust region.
C***    Find the smallest eigen value of the Hessian.
        DSQXI2=DSQXI(1)**2+DSQXI(2)**2
        DSQXI1=DSQRT(DSQXI2)
        TEMP1=(D2SQXI(1,1)+D2SQXI(2,2))/2.0d0
        TEMP2=DSQRT(((D2SQXI(1,1)-D2SQXI(2,2))/2.0d0)**2+D2SQXI(1,2)**2)
        EVMIN=TEMP1-TEMP2
        EVMAX=TEMP1+TEMP2
        IF(DSQXI1.LT.TOL2) GO TO 9998
        DO it2=1,ITMAX
          TEMP=DSQXI1/DELTA
          NEWTON=EVMIN.GE.TEMP
          IF(NEWTON) THEN !Newton is safe
            H(1)=D2SQXI(1,1)
            H(2)=D2SQXI(2,2)
            DET=EVMIN*EVMAX
          ELSE
C***        Shift eigenvalues to restrict step
            MU=TEMP-EVMIN
            H(1)=D2SQXI(1,1)+MU
            H(2)=D2SQXI(2,2)+MU
            DET=TEMP*(EVMAX+MU)
          ENDIF
          V(1)=-(H(2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2)=(D2SQXI(1,2)*DSQXI(1)-H(1)*DSQXI(2))/DET
          V2=V(1)**2+V(2)**2
          DSQV=DSQXI(1)*V(1)+DSQXI(2)*V(2)
C         This checks that numerical errors have not
C         prevented the step direction being a descent direction.

          IF(DSQV**2.LT.DSQXI2*V2*TOL2) THEN !try a smaller trust region
            DELTA=DELTA/10.0d0
          ELSE !step is good
C***        Check feasible and limit step size
            FREE=.TRUE.
            DO ni=1,2
              IF(BOUND(ni).NE.0.AND
     '          .(BOUND(ni).GT.0.EQV.V(ni).LT.0.0d0)) THEN
                FREE=.FALSE.
                nifix=ni
              ENDIF
            ENDDO
            W=1.0d0
            IF(FREE) THEN
              V1=DSQRT(V2) !currently < DELTA
              D2SQV2=V(1)*(V(1)*D2SQXI(1,1)+2.0d0*V(2)*D2SQXI(1,2))
     '          +V(2)**2*D2SQXI(2,2)
              IF(.NOT.NEWTON) THEN
C               Try to step to estimate of minimum along line
                IF(V1.GT.0.0d0) THEN
                  W=DELTA/V1
                  IF(D2SQV2.GT.0.0D0) THEN !minimum exists
                    W=DMIN1(W,-DSQV/D2SQV2) !minimum if within trust region
                  ENDIF
                ENDIF
              ENDIF !newton
            ELSE
              IF(ENFORCE(2)) THEN !gradient suggests must use ni=1
                nifix=2
              ELSE IF(ENFORCE(1)) THEN !gradient suggests must use ni=2
                nifix=1
              !ELSE Gradient points into element.
C               Fix the direction that prevented the step.
C               P.d. Hessian guarantees there is only one of these
C               for this type of gradient.
              ENDIF
              ni=3-nifix
              IF(.NOT.INELEM) THEN
C***            If stepping predominantly out of element then exit
                nifix=3-ni
                IF(DABS(V(nifix)).GT.DABS(DSQXI(ni)/H(ni))) THEN
                  XI(nifix)=XI(nifix)+V(nifix)
                  GO TO 9998
                ENDIF
              ENDIF
              V(nifix)=0.0d0
              IF(D2SQXI(ni,ni).GT.0.0D0) THEN !minimum exists
                V(ni)=-DSQXI(ni)/D2SQXI(ni,ni)
                V1=DABS(V(ni))
                NEWTON=V1.LE.DELTA
              ENDIF
              IF(.NOT.NEWTON) THEN
                V(ni)=-DSIGN(DELTA,DSQXI(ni))
                V1=DELTA
              ENDIF
              V2=V1*V1
              DSQV=DSQXI(ni)*V(ni)
              D2SQV2=V2*D2SQXI(ni,ni)
            ENDIF !free
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              FORMAT='(''    IT1='',I4,10X,''DSQXI(ni)='',2(E12.6,4X),'
     '        //'/18X,''D2SQXI(MI,ni)='',2(E12.6,4X)/48X,E12.6,'
     '        //'/26X,''V(ni)='',2(E12.6,4X))'
              WRITE(OP_STRING,FORMAT) IT,(DSQXI(ni),ni=1,2),
     '          ((D2SQXI(MI,ni),MI=1,ni),ni=1,2),(V(ni),ni=1,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
C***        First half of convergence test.
C           Should be before boundary colllision check
            CONVERGED=V1*W.LT.TOL
C***        Try the step.  (Name: XILIN is historical)
            XILIN(1)=XI(1)+V(1)*W
            XILIN(2)=XI(2)+V(2)*W
C***        Test for boundary collision
            DO ni=1,2
              IF(XILIN(ni).LT.0.0d0) THEN
                XILIN(ni)=0.0d0
                W=XI(ni)/(-V(ni))
                XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
              ELSE IF(XILIN(ni).GT.1.0d0) THEN
                XILIN(ni)=1.0d0
                W=(1.0d0-XI(ni))/V(ni)
                XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
              ENDIF
            ENDDO !ni
C***        Calculate new distance
            SQLIN=0.d0
            DO nj=1,NJT
              nb=NBJ(nj)
              Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XILIN,XE(1,nj))
              DZ(nj)=Z(nj)-ZD(nj)
              SQLIN=SQLIN+DZ(nj)**2
            ENDDO
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              FORMAT='(8X,''IT2='',I4,14X,''W='',E12.6/21X,'
     '          //'''XILIN(ni)='','//'2(E12.6,4X)/26X,''Z(nj)='','
     '          //'3(E12.6,4X)/26X,''SQLIN='',E12.6)'
              WRITE(OP_STRING,FORMAT) it2,W,(XILIN(ni),ni=1,2),
     '          (Z(nj),nj=1,NJT),SQLIN
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
C***        Second half of convergence test.
C KAT 20Nov98: testing before trust region test as small SQLIN may
C              always be > SQ if DSQXI is small (-> V is small).
            CONVERGED=CONVERGED.AND.DABS(SQ-SQLIN)/(1.d0+SQ).LE.TOL
            IF(CONVERGED) GO TO 5
            DSQVW=DSQV*W !<0
            D2SQVW2=0.5d0*D2SQV2*W*W !1/2 for computational efficiency
            SQDIFF=0.5d0*(SQLIN-SQ) !1/2 because derivs are for SQ/2
            SQDPRED=SQDIFF-DSQVW-D2SQVW2
C***        Exit loop if decrease is satisfactory
            IF(SQDIFF.LE.0.25d0*DSQVW) THEN
              IF(NEWTON) THEN
                DELTA=V1 !next step smaller unless this is increased
              ELSE IF(W.GE.1.0d0) THEN
C               If the quadratic model is good increase trust region size
                IF(SQDPRED.LT.-0.1d0*SQDIFF) THEN
                  DELTA=DMIN1(VMAX,DELTA*2.0d0)
                ENDIF
              ENDIF
              GO TO 5
            ENDIF
C***        Calculate new trust region size from an estimate of the
C***        minimum along the step direction using a cubic approximation.
            TEMP=-3.0d0*SQDPRED !<0
            DELTA=
     '        W*V1*(D2SQVW2-DSQRT(D2SQVW2**2+TEMP*DSQVW))/TEMP !>0
C            DELTA=0.5d0*W*V1*DSQVW/(DSQVW-SQDIFF)
          ENDIF !DSQV**2.LT.DSQXI2*V2*TOL2
        ENDDO !it2

 5      SQ=SQLIN
        XI(1)=XILIN(1)
        XI(2)=XILIN(2)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(24X,''Xi(ni)='',2(E12.6,4X))')
     '      (XI(ni),ni=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
C        IF(.NOT.INELEM.AND.(XI(1).LT.-0.1D0.OR.XI(1).GT.1.1D0.OR
C     '    .XI(2).LT.-0.1D0.OR.XI(2).GT.1.1D0)) GO TO 9998
C        DXI(1)=V(1)*W
C        DXI(2)=V(2)*W
C KAT 18Nov98: xi norm is already in unit sized dimensions
C        IF(((DXI(1)**2+DXI(2)**2)/(1.d0+XI(1)**2+XI(2)**2).LE.TOL2)
C     '    .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GO TO 9998
        IT=IT+1
      ENDDO

      IF(.NOT.CONVERGED) THEN
C        TEMP=DABS(SQ-SQLIN)/(1.d0+SQ)
C        IF(TEMP.GT.TOL) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>WARNING!!! Iteration in CLOS21 has not'
     '    //' converged.'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(14X,''Estimate of error magnitude in xi:'',D9.2,''.'')')
     '    W*V1
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
C        ENDIF
      ENDIF

 9998 CALL EXITS('CLOS21')
      RETURN
 9999 CALL ERRORS('CLOS21',ERROR)
      CALL EXITS('CLOS21')
      RETURN 1
      END


