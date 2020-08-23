      SUBROUTINE CLOS31(IBT,IDO,INP,IT,ITMAX,NBJ,SQ,USER_TOL,XE,XI,ZD,
     '  FOUND_E,INELEM,ERROR,*)

C#### Subroutine: CLOS31
C###  Description:
C###    CLOS31 finds the XI-coordinates at the closest approach of a 3D
C###    element to a data point with coordinates XD using a modified
C###    Newton algorithm.
C###    If INELEM is true then the closest point within the element is
C###    returned; if false then the xi coords are only calculated if the
C###    data point lies within the element.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 3-D elements in rectangular cartesian coords.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  IT,ITMAX,NBJ(NJM)
      REAL*8 SQ,USER_TOL,XE(NSM,NJM),XI(3),ZD(NJM)
      LOGICAL FOUND_E,INELEM
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER BOUND(3),it2,MI,nb,ni,ni2,ni_f,NIF(2),nifix,nj
      REAL*8 B,C,D,DELTA,DET,D2SQV2,D2SQVW2,D2SQXI(3,3),D2SQXIF(2,2),
     '  D2ZXI(3,3,3),DSQXI(3),DSQXI1,DSQXI2,DSQXIF(2),DSQXIF1,DSQXIF2,
     '  DSQV,DSQVW,DZ(3),DZXI(3,3),EVMIN,EVFMIN,EVFMAX,H(3,3),HF(2),
     '  HINV(3,3),MAXVIOL,MU,P,Q,PXI,SQLIN,SQDIFF,SQDPRED,TEMP,TEMP1,
     '  TEMP2,TOL,TOL2,V(3),V1,V2,VF(2),VIOL(3),VMAX,W,XILIN(3),Z(3)
      LOGICAL CONVERGED,ENFORCE(3),FREE,NEWTON,TRYSTEP
      DATA VMAX /1.0d0/

      CALL ENTERS('CLOS31',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >CLOS31 3-D rectangular cartesian'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DELTA=VMAX/4.0d0
      TOL=5.0d0*LOOSE_TOL !must be > sqrt(eps) or SQLIN<=SQ check may not work
      TOL2=TOL**2
      SQ=0.0d0
      DO NJ=1,NJT
        NB=NBJ(NJ)
        Z(NJ)=PXI(IBT(1,1,NB),IDO(1,1,0,NB),INP(1,1,NB),NB,1,XI,
     '    XE(1,NJ))
        DZ(NJ)=Z(NJ)-ZD(NJ)
        SQ=SQ+DZ(NJ)**2
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' ZD(NJ)='',3(E12.6,4X),/'' XI(NI)='','
     '    //'3(E12.6,4X),/''  Z(NJ)='',3(E12.6,4X),/''     SQ='','
     '    //'E12.6)') (ZD(NJ),NJ=1,NJT),(XI(NI),NI=1,3),(Z(NJ),
     '    NJ=1,NJT),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      IT=0
      CONVERGED=.FALSE.
      DO WHILE(.NOT.CONVERGED.AND.IT.LT.ITMAX)
        DSQXI(1)=0.0d0
        DSQXI(2)=0.0d0
        DSQXI(3)=0.0d0
        DO nj=1,NJT
          nb=NBJ(nj)
          DZXI(nj,1)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,XI,XE(1,nj))
          DZXI(nj,2)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,XI,XE(1,nj))
          DZXI(nj,3)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,7,XI,XE(1,nj))
          DSQXI(1)=DSQXI(1)+DZXI(nj,1)*DZ(nj)
          DSQXI(2)=DSQXI(2)+DZXI(nj,2)*DZ(nj)
          DSQXI(3)=DSQXI(3)+DZXI(nj,3)*DZ(nj)
        ENDDO
        DO ni=1,3
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
        IF(ENFORCE(1).AND.ENFORCE(2).AND.ENFORCE(3)) GO TO 9998
        D2SQXI(1,1)=0.0d0
        D2SQXI(1,2)=0.0d0
        D2SQXI(1,3)=0.0d0
        D2SQXI(2,2)=0.0d0
        D2SQXI(2,3)=0.0d0
        D2SQXI(3,3)=0.0d0
        DO nj=1,NJT
          nb=NBJ(nj)
          D2ZXI(nj,1,1)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,3,XI,XE(1,nj))
          D2ZXI(nj,1,2)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,6,XI,XE(1,nj))
          D2ZXI(nj,1,3)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,9,XI,XE(1,nj))
          D2ZXI(nj,2,2)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,5,XI,XE(1,nj))
          D2ZXI(nj,2,3)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,10,XI,XE(1,nj))
          D2ZXI(nj,3,3)=
     '      PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,8,XI,XE(1,nj))
          D2SQXI(1,1)=
     '      D2SQXI(1,1)+DZXI(nj,1)*DZXI(nj,1)+D2ZXI(nj,1,1)*DZ(nj)
          D2SQXI(1,2)=
     '      D2SQXI(1,2)+DZXI(nj,1)*DZXI(nj,2)+D2ZXI(nj,1,2)*DZ(nj)
          D2SQXI(1,3)=
     '      D2SQXI(1,3)+DZXI(nj,1)*DZXI(nj,3)+D2ZXI(nj,1,3)*DZ(nj)
          D2SQXI(2,2)=
     '      D2SQXI(2,2)+DZXI(nj,2)*DZXI(nj,2)+D2ZXI(nj,2,2)*DZ(nj)
          D2SQXI(2,3)=
     '      D2SQXI(2,3)+DZXI(nj,2)*DZXI(nj,3)+D2ZXI(nj,2,3)*DZ(nj)
          D2SQXI(3,3)=
     '      D2SQXI(3,3)+DZXI(nj,3)*DZXI(nj,3)+D2ZXI(nj,3,3)*DZ(nj)
        ENDDO !nj
C KAT 15Dec98: A model trust region approach.
C       A Newton step is taken if the condition of the Hessian
C       guarantees that the step will be within the trust region.
C       Otherwise the Hessian is shifted towards a diagonal matrix to
C       shift the step towards steepest descent.  Usually it is a much
C       better direction than steepest descent.  I think it is close to
C       the best direction in the trust region.
        DSQXI2=DSQXI(1)**2+DSQXI(2)**2+DSQXI(3)**2
        DSQXI1=DSQRT(DSQXI2)
        IF(DSQXI1.LT.TOL2) GO TO 9998
C***    Solve a cubic to find the smallest eigenvalue of the Hessian.
        B=-(D2SQXI(1,1)+D2SQXI(2,2)+D2SQXI(3,3))/3.0d0
        C=(D2SQXI(2,2)*D2SQXI(3,3)+D2SQXI(1,1)*D2SQXI(3,3)
     '    +D2SQXI(1,1)*D2SQXI(2,2)
     '    -D2SQXI(2,3)**2-D2SQXI(1,3)**2-D2SQXI(1,2)**2)/3.0d0
        D=D2SQXI(1,1)*D2SQXI(2,3)**2+D2SQXI(2,2)*D2SQXI(1,3)**2
     '    +D2SQXI(3,3)*D2SQXI(1,2)**2
     '    -D2SQXI(1,1)*D2SQXI(2,2)*D2SQXI(3,3)
     '    -(D2SQXI(2,3)*D2SQXI(1,3)*D2SQXI(1,2))*2.0d0
C       Eigenvalues of symmetric matrix are real.
        P=C-B*B !<=0
        IF(P.GE.0.0d0) THEN !include > 0 for numerical error
C         All eigenvalues are the same.
          EVMIN=-B
        ELSE
          Q=D-3.0d0*B*C+2.0d0*B**3
          P=DSQRT(-P)
          TEMP=-Q/(2.0d0*P**3) !-1<=TEMP<=1
          IF(DABS(TEMP).GT.1.0d0) TEMP=DSIGN(1.0d0,TEMP) !for numerical error
C         Smallest root
          EVMIN=2.0d0*P*DCOS((DACOS(TEMP)+2.0d0*PI)/3.0d0)-B
        ENDIF
C       Set up part of an approximate Hessian
        H(1,2)=D2SQXI(1,2)
        H(1,3)=D2SQXI(1,3)
        H(2,3)=D2SQXI(2,3)
        H(2,1)=D2SQXI(1,2)
        H(3,1)=D2SQXI(1,3)
        H(3,2)=D2SQXI(2,3)
        DO it2=1,ITMAX
          TRYSTEP=.TRUE.
          TEMP=DSQXI1/DELTA
          NEWTON=EVMIN.GE.TEMP
          IF(NEWTON) THEN !Newton is safe
            H(1,1)=D2SQXI(1,1)
            H(2,2)=D2SQXI(2,2)
            H(3,3)=D2SQXI(3,3)
          ELSE
C***        Shift eigenvalues to restrict step
            MU=TEMP-EVMIN
            H(1,1)=D2SQXI(1,1)+MU
            H(2,2)=D2SQXI(2,2)+MU
            H(3,3)=D2SQXI(3,3)+MU
          ENDIF
          CALL INVERT(3,H,HINV,DET)
          V(1)=
     '      -(HINV(1,1)*DSQXI(1)+HINV(1,2)*DSQXI(2)+HINV(1,3)*DSQXI(3))
          V(2)=
     '      -(HINV(2,1)*DSQXI(1)+HINV(2,2)*DSQXI(2)+HINV(2,3)*DSQXI(3))
          V(3)=
     '      -(HINV(3,1)*DSQXI(1)+HINV(3,2)*DSQXI(2)+HINV(3,3)*DSQXI(3))
          V2=V(1)**2+V(2)**2+V(3)**2
          DSQV=DSQXI(1)*V(1)+DSQXI(2)*V(2)+DSQXI(3)*V(3)
C         This checks that numerical errors have not
C         prevented the step direction being a descent direction.
          IF(DSQV**2.LT.DSQXI2*V2*TOL2) THEN
            TRYSTEP=.FALSE.
          ELSE !step is good
C***        Check feasible and limit step size
            FREE=.TRUE.
            DO ni=1,3
              IF(BOUND(ni).NE.0) THEN
                VIOL(ni)=-BOUND(ni)*V(ni)
                IF(VIOL(ni).GT.0.0d0) THEN !wants to leave element
                  IF(FREE) THEN
                    FREE=.FALSE.
                    MAXVIOL=VIOL(ni)
                    nifix=ni
                  ELSE IF(VIOL(ni).GT.MAXVIOL) THEN
C                   fix the dirn most strongly suggesting leaving the element
                    MAXVIOL=VIOL(ni)
                    nifix=ni
                  ENDIF
                ENDIF !violation
              ENDIF !bound
            ENDDO !ni
            W=1.0d0
            IF(FREE) THEN
              V1=DSQRT(V2) !currently < DELTA
              D2SQV2=D2SQXI(1,1)*V(1)**2+D2SQXI(2,2)*V(2)**2
     '          +D2SQXI(3,3)*V(3)**2
     '          +2d0*(D2SQXI(1,2)*V(1)*V(2)+D2SQXI(2,3)*V(2)*V(3)
     '          +D2SQXI(1,3)*V(1)*V(3))
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
C             Need to choose a direction to bound. If a gradient
C             suggests leaving the element then fix this direction.
              FREE=.TRUE.
              DO ni=1,3
                IF(ENFORCE(ni)) THEN !gradient suggests leaving element
                  IF(FREE) THEN
                    FREE=.FALSE.
                    MAXVIOL=VIOL(ni)
                    nifix=ni
                  ELSE IF(VIOL(ni).GT.MAXVIOL) THEN
C                   fix the dirn most strongly suggesting leaving the element
                    MAXVIOL=VIOL(ni)
                    nifix=ni
                  ENDIF
                ENDIF !enforce
              ENDDO !ni
C             Directions that are still free make a face
              NIF(1)=1+MOD(nifix,3)
              NIF(2)=1+MOD(nifix+1,3)
              DSQXIF(1)=DSQXI(NIF(1))
              DSQXIF(2)=DSQXI(NIF(2))
              DSQXIF2=DSQXIF(1)**2+DSQXIF(2)**2
              DSQXIF1=DSQRT(DSQXIF2)
              IF(DSQXIF1.LT.TOL2) GO TO 9998
C             The Hessian is now 2x2.
              D2SQXIF(1,1)=D2SQXI(NIF(1),NIF(1))
              D2SQXIF(1,2)=H(NIF(1),NIF(2))
              D2SQXIF(2,2)=D2SQXI(NIF(2),NIF(2))
              IF(.NOT.INELEM) THEN
C***            If stepping predominantly out of element then exit
                TEMP=(H(NIF(1),NIF(1))*DSQXIF(1)
     '            -D2SQXIF(1,2)*DSQXIF(2))**2
     '            +(D2SQXIF(1,2)*DSQXIF(1)
     '            -H(NIF(2),NIF(2))*DSQXIF(2))**2
                IF((V(nifix)*DET)**2.GT.TEMP) THEN
                  XI(nifix)=XI(nifix)+V(nifix)
                  GO TO 9997
                ENDIF
              ENDIF
              V(nifix)=0.0d0
C***          Find the minimum eigenvalue of the new Hessian
              TEMP1=(D2SQXIF(1,1)+D2SQXIF(2,2))/2.0d0
              TEMP2=DSQRT(((D2SQXIF(1,1)-D2SQXIF(2,2))/2.0d0)**2
     '          +D2SQXIF(1,2)**2)
              EVFMIN=TEMP1-TEMP2
              EVFMAX=TEMP1+TEMP2
              TEMP=DSQXIF1/DELTA
              NEWTON=EVMIN.GE.TEMP
              IF(NEWTON) THEN !Newton is safe
                HF(1)=D2SQXIF(1,1)
                HF(2)=D2SQXIF(2,2)
                DET=EVFMIN*EVFMAX
              ELSE
C***            Shift eigenvalues to restrict step
                MU=TEMP-EVFMIN
                HF(1)=D2SQXIF(1,1)+MU
                HF(2)=D2SQXIF(2,2)+MU
                DET=TEMP*(EVFMAX+MU)
              ENDIF
              VF(1)=-(HF(2)*DSQXIF(1)-D2SQXIF(1,2)*DSQXIF(2))/DET
              VF(2)=(D2SQXIF(1,2)*DSQXIF(1)-HF(1)*DSQXIF(2))/DET
              V2=VF(1)**2+VF(2)**2
              DSQV=DSQXIF(1)*VF(1)+DSQXIF(2)*VF(2)
C             This checks that numerical errors have not
C             prevented the step direction being a descent direction.
              IF(DSQV**2.LT.DSQXIF2*V2*TOL2) THEN
                TRYSTEP=.FALSE.
              ELSE !step is good
                V1=DSQRT(V2) !currently < DELTA
C***            Check feasible and limit step size
                FREE=.TRUE.
                DO ni_f=1,2
                  IF(BOUND(NIF(ni_f)).NE.0.AND
     '              .(BOUND(NIF(ni_f)).GT.0.EQV.VF(ni_f).LT.0.0d0)) THEN
                    FREE=.FALSE.
                    nifix=ni_f
                  ENDIF
                ENDDO
                W=1.0d0
                IF(FREE) THEN
                  D2SQV2=VF(1)*(VF(1)*D2SQXIF(1,1)
     '              +2.0d0*VF(2)*D2SQXIF(1,2))+VF(2)**2*D2SQXIF(2,2)
                  V(NIF(1))=VF(1)
                  V(NIF(2))=VF(2)
                  IF(.NOT.NEWTON) THEN
C                   Try to step to estimate of minimum along line
                    IF(V1.GT.0.0d0) THEN
                      W=DELTA/V1
                      IF(D2SQV2.GT.0.0D0) THEN !minimum exists
                        W=DMIN1(W,-DSQV/D2SQV2) !minimum if in trust region
                      ENDIF
                    ENDIF
                  ENDIF !newton
                ELSE
                  IF(ENFORCE(NIF(2))) THEN !gradient => must use ni=1
                    nifix=2
                  ELSE IF(ENFORCE(NIF(1))) THEN !gradient => must use ni=2
                    nifix=1
                  !ELSE Gradient points into element.
C                   Fix the direction that prevented the step.
C                   P.d. Hessian guarantees there is only one of these
C                   for this type of gradient.
                  ENDIF
                  ni_f=3-nifix !free direction in face
                  IF(.NOT.INELEM) THEN
C***                If stepping predominantly out of element then exit
                    IF(DABS(VF(nifix)).GT.
     '                DABS(DSQXIF(ni_f)/HF(ni_f))) THEN
                      XI(NIF(nifix))=XI(NIF(nifix))+VF(nifix)
                      GO TO 9997
                    ENDIF
                  ENDIF
                  ni=NIF(ni_f)
                  nifix=NIF(nifix)
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
              ENDIF !trystep
            ENDIF !free
          ENDIF !trystep
          IF(TRYSTEP) THEN
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI(NI)='','
     '          //'3(E12.6,4X),/18X,''D2SQXI(MI,NI)='',3(E12.6,4X),'
     '          //'/32X,3(E12.6,4X),/32X,3(E12.6,4X),/26X,''V(NI)='','
     '          //'3(E12.6,4X),/26X,''DET,DELTA,DSQV='',3(E12.6,4X))')
     '          IT,(DSQXI(NI),NI=1,3),((D2SQXI(MI,NI),MI=1,3),NI=1,3),
     '          (V(NI),NI=1,3),DET,DELTA,DSQV
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
C***        First half of convergence test.
C           Should be before boundary colllision check
            CONVERGED=V1*W.LT.TOL
C***        Try the step.  (Name: XILIN is historical)
            XILIN(1)=XI(1)+V(1)*W
            XILIN(2)=XI(2)+V(2)*W
            XILIN(3)=XI(3)+V(3)*W
C***        Test for boundary collision
            DO ni=1,3
              IF(XILIN(ni).LT.0.0d0) THEN
                XILIN(ni)=0.0d0
                W=XI(ni)/(-V(ni))
                ni2=1+MOD(ni,3)
                XILIN(ni2)=XI(ni2)+V(ni2)*W
                ni2=1+MOD(ni2,3)
                XILIN(ni2)=XI(ni2)+V(ni2)*W
              ELSE IF(XILIN(ni).GT.1.0d0) THEN
                XILIN(ni)=1.0d0
                W=(1.0d0-XI(ni))/V(ni)
                ni2=1+MOD(ni,3)
                XILIN(ni2)=XI(ni2)+V(ni2)*W
                ni2=1+MOD(ni2,3)
                XILIN(ni2)=XI(ni2)+V(ni2)*W
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
              WRITE(OP_STRING,'(8X,''IT2='',I4,14X,''W='',E12.6/21X,'//
     '          '''XILIN(NI)='',3(E12.6,4X)/26X,''Z(NJ)='',3(E12.6,4X)'
     '          //'/26X,''SQLIN='',E12.6)') IT2,W,(XILIN(NI),NI=1,3),
     '          (Z(NJ),NJ=1,NJT),SQLIN
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
            DELTA=W*V1*(D2SQVW2-DSQRT(D2SQVW2**2+TEMP*DSQVW))/TEMP !>0
          ELSE !not try step
            DELTA=DELTA/10.0d0 !try a smaller trust region
          ENDIF
        ENDDO !it2

 5      SQ=SQLIN
        XI(1)=XILIN(1)
        XI(2)=XILIN(2)
        XI(3)=XILIN(3)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(24X,''Xi(ni)='',3(E12.6,4X))') (XI(NI),
     '      NI=1,3)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IT=IT+1
      ENDDO

      IF(.NOT.CONVERGED) THEN
C        TEMP=DABS(SQ-SQLIN)/(1.d0+SQ)
C        IF(TEMP.GT.TOL) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' >>WARNING!!! Iteration in CLOS31 has not'
     '    //' converged.'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,
     '    '(14X,''Estimate of error magnitude in xi:'',D9.2,''.'')')
     '    W*V1
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
C        ENDIF
        FOUND_E=.FALSE.
      ELSE
        IF(DABS(Z(1)-ZD(1)).LE.USER_TOL.AND.DABS(Z(2)-ZD(2)).LE.USER_TOL
     '    .AND.DABS(Z(3)-ZD(3)).LE.USER_TOL) FOUND_E=.TRUE.
      ENDIF

c      IF(DABS(Z(1)-ZD(1)).LE.USER_TOL.AND.DABS(Z(2)-ZD(2)).LE.USER_TOL
c     '  .AND.DABS(Z(3)-ZD(3)).LE.USER_TOL)THEN
c        FOUND_E=.TRUE.
c      ELSE
c        FOUND_E=.FALSE.
c      ENDIF

      CALL EXITS('CLOS31')
      RETURN
 9997 FOUND_E=.FALSE.
      CALL EXITS('CLOS31')
      RETURN
 9998 IF(DABS(Z(1)-ZD(1)).LE.USER_TOL.AND.DABS(Z(2)-ZD(2)).LE.USER_TOL
     &  .AND.DABS(Z(3)-ZD(3)).LE.USER_TOL) FOUND_E=.TRUE.
      CALL EXITS('CLOS31')
      RETURN
 9999 CALL ERRORS('CLOS31',ERROR)
      CALL EXITS('CLOS31')
      RETURN 1
      END


