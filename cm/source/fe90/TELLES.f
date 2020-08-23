      SUBROUTINE TELLES(IBT,IDO,INP,MAXITS,nb,NLL,NW,D,DL,TOLERANCE,
     '  XE,XIMIN,XPFP,ERROR,*)

C#### Subroutine: TELLES
C###  Description:
C###    TELLES finds the xi location of the point in an element
C###    (specified by XE) that is
C###    the minimum distance from a another point (specfied by XPFP)
C###    for Telles adaptive integration.
C###    It also returns the D values for the integration.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Paramter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),MAXITS,nb,NLL(12),NW
      REAL*8 D(*),DL(3,NLM),TOLERANCE,XE(NSM,NJM),XIMIN(*),XPFP(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ni,NUMITS
      REAL*8 ALPHA,E(3,2),ERR,DELTAXI(2),DETDHDXI,DHDXI(2,2),DIST,
     '  DOTG,DXDXI(3,3,2,2),GRADF(2,2),H(2),NORME,NORMGRAD,PXI,X(3,2),
     '  XI(2,2)
      LOGICAL APEXNODE,CONVERGED,HERMSECTOR

      CALL ENTERS('TELLES',*9999)

C CPB 20/6/95 Using unrolled loops

C*** Find point in the element that is closest to the singularity
      NUMITS=0
      CONVERGED=.FALSE.
      IF(NIT(nb).EQ.1) THEN
C*** 1d bounday element therefore assume 2d

C*** cpb 27/6/95 Use Newtons method to find the minimum point. Note if
C*** a maximum distance is within or close to the element then the
C*** hessian will not be positive definite which is required for
C*** Newton's method. In this case reverse the search direction by
C*** making alpha = -1. This is not strictly correct but given the
C*** minimum point will occur at one end of the element this is a
C*** good approximation.

        XI(1,1)=0.5d0
        X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,1),
     '    XE(1,1))
        X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,1),
     '    XE(1,2))
        E(1,1)=XPFP(1)-X(1,1)
        E(2,1)=XPFP(2)-X(2,1)
        DO WHILE(NUMITS.LT.MAXITS.AND..NOT.CONVERGED)
          DXDXI(1,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '      XI(1,1),XE(1,1))
          DXDXI(2,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '      XI(1,1),XE(1,2))
          DXDXI(1,2,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,3,
     '      XI(1,1),XE(1,1))
          DXDXI(2,2,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,3,
     '      XI(1,1),XE(1,2))
          H(1)=-E(1,1)*DXDXI(1,1,1,1)-E(2,1)*DXDXI(2,1,1,1)
          DHDXI(1,1)=DXDXI(1,1,1,1)*DXDXI(1,1,1,1)+DXDXI(2,1,1,1)*
     '      DXDXI(2,1,1,1)-E(1,1)*DXDXI(1,2,1,1)-E(2,1)*DXDXI(2,2,1,1)
          CALL ASSERT(DABS(DHDXI(1,1)).GT.RDELTA,'>>Zero Hessian?',
     '      ERROR,*9999)
          IF(DHDXI(1,1).LT.0.0d0) THEN
            ALPHA=-1.0d0
          ELSE
            ALPHA=1.0d0
          ENDIF
          DELTAXI(1)=-ALPHA*H(1)/DHDXI(1,1)
          XI(1,2)=XI(1,1)
          XI(1,1)=XI(1,1)+DELTAXI(1)
          IF(XI(1,1).LT.0.0d0) XI(1,1)=0.0d0
          IF(XI(1,1).GT.1.0d0) XI(1,1)=1.0d0
          ERR=DABS(XI(1,1)-XI(1,2))
          CONVERGED=ERR.LE.TOLERANCE
          NUMITS=NUMITS+1
          X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '      XI(1,1),XE(1,1))
          X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '      XI(1,1),XE(1,2))
          E(1,1)=XPFP(1)-X(1,1)
          E(2,1)=XPFP(2)-X(2,1)
          NORME=DSQRT(E(1,1)*E(1,1)+E(2,1)*E(2,1))
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_1)
            WRITE(OP_STRING,'(/'' Iteration '',I2)') NUMITS
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Norm. error : '',D12.4,'
     '        //''', Tolerance : '',D12.4)') ERR,TOLERANCE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' OLD XI   : '',D12.4)') XI(1,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' DELTAXI  : '',D12.4)') DELTAXI(1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' XI       : '',D12.4)') XI(1,1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Distance : '',D12.4)') NORME
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_1)
          ENDIF
        ENDDO
        XIMIN(1)=XI(1,1)
        D(1)=NORME/DL(3,NLL(1))
        D(2)=0.0d0
      ELSE IF(NIT(nb).EQ.2) THEN
C*** 2d bounday element therefore assume 3d

C*** cpb 27/6/95 2d element is more tricky than a 1d element. If a
C*** maximum distance is within, or very close, to the element then
C*** the Hessian will not necessisarily be positive definite which is
C*** required for Newton's method to work. To avoid this problem the
C*** following approach is needed. Firstly use steepest descent to find
C*** the search direction. Next take the maximum step length to go to
C*** the boundary of the element and re-evaluate the gradient at this
C*** point. If these two gradients indicate a minimum point within the
C*** element Newton's method will be used otherwise repeat the above
C*** procedure but restrict the domain to the corresponding line of
C*** the element boundary.

C*** Find gradient at the point xi=(0.5,0.5)
        XI(1,1)=0.5d0
        XI(2,1)=0.5d0
        X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,1),
     '    XE(1,1))
        X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,1),
     '    XE(1,2))
        X(3,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,1),
     '    XE(1,3))
        E(1,1)=XPFP(1)-X(1,1)
        E(2,1)=XPFP(2)-X(2,1)
        E(3,1)=XPFP(3)-X(3,1)
        DIST=DSQRT(E(1,1)*E(1,1)+E(2,1)*E(2,1)+E(3,1)*E(3,1))
        DXDXI(1,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '    XI(1,1),XE(1,1))
        DXDXI(2,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '    XI(1,1),XE(1,2))
        DXDXI(3,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '    XI(1,1),XE(1,3))
        DXDXI(1,1,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,
     '    XI(1,1),XE(1,1))
        DXDXI(2,1,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,
     '    XI(1,1),XE(1,2))
        DXDXI(3,1,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,
     '    XI(1,1),XE(1,3))
        GRADF(1,1)=-(E(1,1)*DXDXI(1,1,1,1)+E(2,1)*DXDXI(2,1,1,1)+
     '    E(3,1)*DXDXI(3,1,1,1))/DIST
        GRADF(2,1)=-(E(1,1)*DXDXI(1,1,2,1)+E(2,1)*DXDXI(2,1,2,1)+
     '    E(3,1)*DXDXI(3,1,2,1))/DIST
        NORMGRAD=DSQRT(GRADF(1,1)*GRADF(1,1)+GRADF(2,1)*GRADF(2,1))
        IF(DABS(NORMGRAD).GT.ZERO_TOL) THEN
          GRADF(1,1)=GRADF(1,1)/NORMGRAD
          GRADF(2,1)=GRADF(2,1)/NORMGRAD
        ENDIF
C*** Step along the opposite of the gradient direction to hit the
C*** element boundary.
        APEXNODE=.FALSE.
        HERMSECTOR=.FALSE.
        IF(NBC(nb).EQ.5) THEN !BEM tensor product basis
          IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
            ni=1
            XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
            IF(GRADF(2,1).LT.0.0d0) THEN
              XI(2,2)=1.0d0
            ELSE
              XI(2,2)=0.0d0
            ENDIF
          ELSE
            ni=2
            XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
            IF(GRADF(1,1).LT.0.0d0) THEN
              XI(1,2)=1.0d0
            ELSE
              XI(1,2)=0.0d0
            ENDIF
          ENDIF
        ELSE IF(NBC(nb).EQ.6) THEN !BEM simplex/sector basis
          IF(IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4) THEN !Herm simplex
            HERMSECTOR=.TRUE.
            IF(NKT(1,nb).EQ.1) THEN !Apex node 1
              IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
                ni=1
                IF(GRADF(2,1).LT.0.0d0) THEN
                  XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
                  XI(2,2)=1.0d0
                ELSE
                  APEXNODE=.TRUE.
                  XI(1,2)=0.0d0
                  XI(2,2)=0.0d0
                ENDIF
              ELSE
                ni=2
                XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
                IF(GRADF(1,1).LT.0.0d0) THEN
                  XI(1,2)=1.0d0
                ELSE
                  XI(1,2)=0.0d0
                ENDIF
              ENDIF
            ELSE IF(NKT(3,nb).EQ.1) THEN !Apex node 3
              IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
                ni=1
                IF(GRADF(2,1).LT.0.0d0) THEN
                  APEXNODE=.TRUE.
                  XI(1,2)=0.0d0
                  XI(2,2)=1.0d0
                ELSE
                  XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
                  XI(2,2)=0.0d0
                ENDIF
              ELSE
                ni=2
                XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
                IF(GRADF(1,1).LT.0.0d0) THEN
                  XI(1,2)=1.0d0
                ELSE
                  XI(1,2)=0.0d0
                ENDIF
              ENDIF
            ELSE
              ERROR='>>Invalid Hermite simplex element'
              GOTO 9999
            ENDIF
          ELSE !Sector element
            IF(IBT(1,1,nb).EQ.5) THEN !Collapsed xi1 at xi2=0
              IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
                ni=1
                IF(GRADF(2,1).LT.0.0d0) THEN
                  XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
                  XI(2,2)=1.0d0
                ELSE
                  APEXNODE=.TRUE.
                  XI(1,2)=0.0d0
                  XI(2,2)=0.0d0
                 ENDIF
              ELSE
                ni=2
                XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
                IF(GRADF(1,1).LT.0.0d0) THEN
                  XI(1,2)=1.0d0
                ELSE
                  XI(1,2)=0.0d0
                ENDIF
              ENDIF
            ELSE IF(IBT(1,1,nb).EQ.6) THEN !Collapsed xi1 at xi2=1
              IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
                ni=1
                IF(GRADF(2,1).LT.0.0d0) THEN
                  APEXNODE=.TRUE.
                  XI(1,2)=0.0d0
                  XI(2,2)=1.0d0
                ELSE
                  XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
                  XI(2,2)=0.0d0
                ENDIF
              ELSE
                ni=2
                XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
                IF(GRADF(1,1).LT.0.0d0) THEN
                  XI(1,2)=1.0d0
                ELSE
                  XI(1,2)=0.0d0
                ENDIF
              ENDIF
            ELSE IF(IBT(1,2,nb).EQ.5) THEN !Collapsed xi2 at xi1=0
              IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
                ni=1
                XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
                IF(GRADF(2,1).LT.0.0d0) THEN
                  XI(2,2)=1.0d0
                ELSE
                  XI(2,2)=0.0d0
                ENDIF
              ELSE
                ni=2
                IF(GRADF(1,1).LT.0.0d0) THEN
                  XI(1,2)=1.0d0
                  XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
                ELSE
                  APEXNODE=.TRUE.
                  XI(1,2)=0.0d0
                  XI(2,2)=0.0d0
                ENDIF
              ENDIF
            ELSE IF(IBT(1,2,nb).EQ.6) THEN !Collapsed xi2 at xi1=1
              IF(DABS(GRADF(2,1)).GT.DABS(GRADF(1,1))) THEN
                ni=1
                XI(1,2)=0.5d0*(1.0d0-GRADF(1,1)/GRADF(2,1))
                IF(GRADF(2,1).LT.0.0d0) THEN
                  XI(2,2)=1.0d0
                ELSE
                  XI(2,2)=0.0d0
                ENDIF
              ELSE
                ni=2
                IF(GRADF(1,1).LT.0.0d0) THEN
                  APEXNODE=.TRUE.
                  XI(1,2)=1.0d0
                  XI(2,2)=0.0d0
                ELSE
                  XI(1,2)=0.0d0
                  XI(2,2)=0.5d0*(1.0d0-GRADF(2,1)/GRADF(1,1))
                ENDIF
              ENDIF
            ELSE
              ERROR='>>Invalid sector element'
              GOTO 9999
            ENDIF
          ENDIF
        ELSE
          ERROR='>>Unknown basis type'
          GOTO 9999
        ENDIF
C*** Find gradient at this boundary point
        X(1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,2),
     '    XE(1,1))
        X(2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,2),
     '    XE(1,2))
        X(3,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI(1,2),
     '    XE(1,3))
        E(1,2)=XPFP(1)-X(1,2)
        E(2,2)=XPFP(2)-X(2,2)
        E(3,2)=XPFP(3)-X(3,2)
        NORME=DSQRT(E(1,2)*E(1,2)+E(2,2)*E(2,2)+E(3,2)*E(3,2))
        DXDXI(1,1,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '    XI(1,2),XE(1,1))
        DXDXI(2,1,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '    XI(1,2),XE(1,2))
        DXDXI(3,1,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '    XI(1,2),XE(1,3))
        DXDXI(1,1,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,
     '    XI(1,2),XE(1,1))
        DXDXI(2,1,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,
     '    XI(1,2),XE(1,2))
        DXDXI(3,1,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,4,
     '    XI(1,2),XE(1,3))
        IF(APEXNODE) THEN
          IF(HERMSECTOR) THEN
            GRADF(1,2)=0.0d0
            IF(XI(2,2).EQ.0.0d0) THEN
              GRADF(2,2)=1.0d0
            ELSE
              GRADF(2,2)=-1.0d0
            ENDIF
          ELSE
            IF(IBT(1,1,nb).EQ.5) THEN
              GRADF(1,2)=0.0d0
              GRADF(2,2)=1.0d0
            ELSE IF(IBT(1,1,nb).EQ.6) THEN
              GRADF(1,2)=0.0d0
              GRADF(2,2)=-1.0d0
            ELSE IF(IBT(1,2,nb).EQ.5) THEN
              GRADF(1,2)=1.0d0
              GRADF(2,2)=0.0d0
            ELSE IF(IBT(1,2,nb).EQ.6) THEN
              GRADF(1,2)=-1.0d0
              GRADF(2,2)=0.0d0
            ENDIF
          ENDIF
        ELSE
          GRADF(1,2)=-(E(1,2)*DXDXI(1,1,1,2)+E(2,2)*DXDXI(2,1,1,2)+
     '      E(3,2)*DXDXI(3,1,1,2))/NORME
          GRADF(2,2)=-(E(1,2)*DXDXI(1,1,2,2)+E(2,2)*DXDXI(2,1,2,2)+
     '      E(3,2)*DXDXI(3,1,2,2))/NORME
        ENDIF
        NORMGRAD=DSQRT(GRADF(1,2)*GRADF(1,2)+GRADF(2,2)*GRADF(2,2))
        IF(DABS(NORMGRAD).GT.ZERO_TOL) THEN
          GRADF(1,2)=GRADF(1,2)/NORMGRAD
          GRADF(2,2)=GRADF(2,2)/NORMGRAD
        ENDIF
        DOTG=GRADF(1,1)*GRADF(1,2)+GRADF(2,1)*GRADF(2,2)
        IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_2)
          WRITE(OP_STRING,'(/'' Interior point'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' XI       :'',2D12.4)') XI(1,1),XI(2,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GRADIENT :'',2D12.4)') GRADF(1,1),
     '      GRADF(2,1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' DISTANCE :'',D12.4)') DIST
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Boundary point'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' XI       :'',2D12.4)') XI(1,2),XI(2,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' GRADIENT :'',2D12.4)') GRADF(1,2),
     '      GRADF(2,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' DISTANCE :'',D12.4)') NORME
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Gradient dot product :'',D12.4)') DOTG
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_2)
        ENDIF
        IF(DOTG.LT.-1.0d-8) THEN
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_3)
            WRITE(OP_STRING,
     '        '(/'' Switching to 2D Newton''''s Method'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_3)
          ENDIF
C*** Minimum is within the element - use Newton's method
          DO WHILE(NUMITS.LT.MAXITS.AND..NOT.CONVERGED)
            DXDXI(1,2,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,3,XI(1,1),XE(1,1))
            DXDXI(2,2,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,3,XI(1,1),XE(1,2))
            DXDXI(3,2,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,3,XI(1,1),XE(1,3))
            DXDXI(1,2,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,5,XI(1,1),XE(1,1))
            DXDXI(2,2,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,5,XI(1,1),XE(1,2))
            DXDXI(3,2,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,5,XI(1,1),XE(1,3))
            DXDXI(1,3,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,6,XI(1,1),XE(1,1))
            DXDXI(2,3,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,6,XI(1,1),XE(1,2))
            DXDXI(3,3,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,6,XI(1,1),XE(1,3))
            H(1)=-E(1,1)*DXDXI(1,1,1,1)-E(2,1)*DXDXI(2,1,1,1)-E(3,1)*
     '        DXDXI(3,1,1,1)
            H(2)=-E(1,1)*DXDXI(1,1,2,1)-E(2,1)*DXDXI(2,1,2,1)-E(3,1)*
     '        DXDXI(3,1,2,1)
            DHDXI(1,1)=DXDXI(1,1,1,1)*DXDXI(1,1,1,1)+DXDXI(2,1,1,1)*
     '        DXDXI(2,1,1,1)+DXDXI(3,1,1,1)*DXDXI(3,1,1,1)-E(1,1)*
     '        DXDXI(1,2,1,1)-E(2,1)*DXDXI(2,2,1,1)-E(3,1)*DXDXI(3,2,1,1)
            DHDXI(1,2)=DXDXI(1,1,1,1)*DXDXI(1,1,2,1)+DXDXI(2,1,1,1)*
     '        DXDXI(2,1,2,1)+DXDXI(3,1,1,1)*DXDXI(3,1,2,1)-E(1,1)*
     '        DXDXI(1,3,1,1)-E(2,1)*DXDXI(2,3,1,1)-E(3,1)*DXDXI(3,3,1,1)
            DHDXI(2,2)=DXDXI(1,1,2,1)*DXDXI(1,1,2,1)+DXDXI(2,1,2,1)*
     '        DXDXI(2,1,2,1)+DXDXI(3,1,2,1)*DXDXI(3,1,2,1)-E(1,1)*
     '        DXDXI(1,2,2,1)-E(2,1)*DXDXI(2,2,2,1)-E(3,1)*DXDXI(3,2,2,1)
            DETDHDXI=DHDXI(1,1)*DHDXI(2,2)-DHDXI(1,2)*DHDXI(1,2)
            IF(DABS(DETDHDXI).GT.ZERO_TOL) THEN
C              IF(DHDXI(1,1).LT.0.0d0.AND.DETDHDXI.LT.0.0d0) THEN
CC$              call mp_setlock()
C                WRITE(OP_STRING,
C     '            '('' >>Warning: Negative definate Hessian'')')
C                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
C              ENDIF
              DELTAXI(1)=-1.0d0/DETDHDXI*(DHDXI(2,2)*H(1)-
     '          DHDXI(1,2)*H(2))
              DELTAXI(2)=-1.0d0/DETDHDXI*(DHDXI(1,1)*H(2)-
     '          DHDXI(1,2)*H(1))
            ELSE !Zero determinant Hessian
C cpb 5/6/97 Hessian could have a zero determinant because there is
C a zero row or zero column caused by one of the xi directions having
C no-bearing on the function to be minimised. This could occur at a
C collapsed node. In this case try and switch to 1d Newtons method
C with the other xi direction.
              IF(DABS(DHDXI(1,1)).GT.ZERO_TOL) THEN !xi 2 not involved
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_4)
                  WRITE(OP_STRING,
     '              '(/'' Xi 2 Hessian zero. Switching to 1D '
     '              //'Newton''''s method in Xi 1'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_4)
                ENDIF
                IF(DHDXI(1,1).LT.0.0d0) THEN
                  ALPHA=-1.0d0
                ELSE
                  ALPHA=1.0d0
                ENDIF
                DELTAXI(1)=-ALPHA*H(1)/DHDXI(1,1)
                DELTAXI(2)=0.0d0
              ELSE IF(DABS(DHDXI(2,2)).GT.ZERO_TOL) THEN !xi 1 not involved
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_5)
                  WRITE(OP_STRING,
     '              '(/'' Xi 1 Hessian zero. Switching to 1D '
     '              //'Newton''''s method in Xi 2'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_5)
                ENDIF
                IF(DHDXI(2,2).LT.0.0d0) THEN
                  ALPHA=-1.0d0
                ELSE
                  ALPHA=1.0d0
                ENDIF
                DELTAXI(1)=0.0d0
                DELTAXI(2)=-ALPHA*H(2)/DHDXI(2,2)
              ELSE !zero hessian
                ERROR='>>Zero Hessian?'
                GOTO 9999
              ENDIF
            ENDIF
            XI(1,2)=XI(1,1)
            XI(2,2)=XI(2,1)
            DOTG=GRADF(1,1)*DELTAXI(1)+GRADF(2,1)*DELTAXI(2)
            IF(DOTG.GT.0.0d0) THEN !Make sure we go against the gradient
              XI(1,1)=XI(1,1)-DELTAXI(1)
              XI(2,1)=XI(2,1)-DELTAXI(2)
            ELSE
              XI(1,1)=XI(1,1)+DELTAXI(1)
              XI(2,1)=XI(2,1)+DELTAXI(2)
            ENDIF
            IF(XI(1,1).LT.0.0d0) XI(1,1)=0.0d0
            IF(XI(1,1).GT.1.0d0) XI(1,1)=1.0d0
            IF(XI(2,1).LT.0.0d0) XI(2,1)=0.0d0
            IF(XI(2,1).GT.1.0d0) XI(2,1)=1.0d0
            ERR=DSQRT((XI(1,1)-XI(1,2))*(XI(1,1)-XI(1,2))+
     '        (XI(2,1)-XI(2,2))*(XI(2,1)-XI(2,2)))
            CONVERGED=ERR.LE.TOLERANCE
            NUMITS=NUMITS+1
            X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI(1,1),XE(1,1))
            X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI(1,1),XE(1,2))
            X(3,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI(1,1),XE(1,3))
            E(1,1)=XPFP(1)-X(1,1)
            E(2,1)=XPFP(2)-X(2,1)
            E(3,1)=XPFP(3)-X(3,1)
            NORME=DSQRT(E(1,1)*E(1,1)+E(2,1)*E(2,1)+E(3,1)*E(3,1))
            IF(NUMITS.LT.MAXITS.AND..NOT.CONVERGED) THEN
              DXDXI(1,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,2,XI(1,1),XE(1,1))
              DXDXI(2,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,2,XI(1,1),XE(1,2))
              DXDXI(3,1,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,2,XI(1,1),XE(1,3))
              DXDXI(1,1,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,4,XI(1,1),XE(1,1))
              DXDXI(2,1,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,4,XI(1,1),XE(1,2))
              DXDXI(3,1,2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,4,XI(1,1),XE(1,3))
              GRADF(1,1)=-(E(1,1)*DXDXI(1,1,1,1)+E(2,1)*DXDXI(2,1,1,1)+
     '          E(3,1)*DXDXI(3,1,1,1))/DIST
              GRADF(2,1)=-(E(1,1)*DXDXI(1,1,2,1)+E(2,1)*DXDXI(2,1,2,1)+
     '          E(3,1)*DXDXI(3,1,2,1))/DIST
              NORMGRAD=DSQRT(GRADF(1,1)*GRADF(1,1)+GRADF(2,1)*
     '          GRADF(2,1))
              IF(DABS(NORMGRAD).GT.ZERO_TOL) THEN
                GRADF(1,1)=GRADF(1,1)/NORMGRAD
                GRADF(2,1)=GRADF(2,1)/NORMGRAD
              ENDIF
            ENDIF
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_6)
              WRITE(OP_STRING,'(/'' Iteration '',I2)') NUMITS
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Norm. error : '',D12.4,'
     '          //''', Tolerance : '',D12.4)') ERR,TOLERANCE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' OLD XI   : '',2D12.4)') XI(1,2),
     '          XI(2,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DELTAXI  : '',2D12.4)') DELTAXI(1),
     '          DELTAXI(2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XI       : '',2D12.4)') XI(1,1),
     '          XI(2,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Distance : '',D12.4)') NORME
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_6)
            ENDIF
          ENDDO
        ELSE
C*** Minimum is along boundary
          X(1,1)=X(1,2)
          X(2,1)=X(2,2)
          X(3,1)=X(3,2)
          DXDXI(1,1,ni,1)=DXDXI(1,1,ni,2)
          DXDXI(2,1,ni,1)=DXDXI(2,1,ni,2)
          DXDXI(3,1,ni,1)=DXDXI(3,1,ni,2)
          E(1,1)=E(1,2)
          E(2,1)=E(2,2)
          E(3,1)=E(3,2)
          XI(1,1)=XI(1,2)
          XI(2,1)=XI(2,2)
          GRADF(1,1)=GRADF(1,2)
          GRADF(2,1)=GRADF(2,2)
          IF(DABS(GRADF(ni,2)).GT.1.0d-6) THEN
C*** Move to end of extremum node of the element along the appropriate
C*** xi direction
            IF(GRADF(ni,2).LT.0.0d0) THEN
              XI(ni,2)=1.0d0
            ELSE
              XI(ni,2)=0.0d0
            ENDIF
C*** Evaluate the gradient at this point
            X(1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI(1,2),XE(1,1))
            X(2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI(1,2),XE(1,2))
            X(3,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XI(1,2),XE(1,3))
            E(1,2)=XPFP(1)-X(1,2)
            E(2,2)=XPFP(2)-X(2,2)
            E(3,2)=XPFP(3)-X(3,2)
            NORME=DSQRT(E(1,2)*E(1,2)+E(2,2)*E(2,2)+E(3,2)*E(3,2))
            DXDXI(1,1,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,2,XI(1,2),XE(1,1))
            DXDXI(2,1,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,2,XI(1,2),XE(1,2))
            DXDXI(3,1,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,2,XI(1,2),XE(1,3))
            DXDXI(1,1,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,4,XI(1,2),XE(1,1))
            DXDXI(2,1,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,4,XI(1,2),XE(1,2))
            DXDXI(3,1,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,4,XI(1,2),XE(1,3))
            GRADF(1,2)=-(E(1,2)*DXDXI(1,1,1,2)+E(2,2)*DXDXI(2,1,1,2)+
     '        E(3,2)*DXDXI(3,1,1,2))/NORME
            GRADF(2,2)=-(E(1,2)*DXDXI(1,1,2,2)+E(2,2)*DXDXI(2,1,2,2)+
     '        E(3,2)*DXDXI(3,1,2,2))/NORME
            NORMGRAD=DSQRT(GRADF(1,2)*GRADF(1,2)+GRADF(2,2)*GRADF(2,2))
            IF(DABS(NORMGRAD).GT.ZERO_TOL) THEN
              GRADF(1,2)=GRADF(1,2)/NORMGRAD
              GRADF(2,2)=GRADF(2,2)/NORMGRAD
            ENDIF
            DOTG=GRADF(1,1)*GRADF(1,2)+GRADF(2,1)*GRADF(2,2)
            IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_7)
              WRITE(OP_STRING,'('' Nodal extremum point'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XI       :'',2D12.4)') XI(1,2),
     '          XI(2,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' GRADIENT :'',2D12.4)') GRADF(1,2),
     '          GRADF(2,2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' DISTANCE :'',D12.4)') NORME
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Gradient dot product :'',D12.4)')
     '          DOTG
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_7)
            ENDIF
            IF(DOTG.LT.-1.0d-8) THEN
C*** Minimum is not at the corner node. Use Newton's method in 1D
C*** to find the mimimum
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_8)
                WRITE(OP_STRING,
     '            '(/'' Switching to 1D Newton''''s Method'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_8)
              ENDIF
              DO WHILE(NUMITS.LT.MAXITS.AND..NOT.CONVERGED)
                DXDXI(1,2,ni,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,2*ni+1,XI(1,1),XE(1,1))
                DXDXI(2,2,ni,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,2*ni+1,XI(1,1),XE(1,2))
                DXDXI(3,2,ni,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '            INP(1,1,nb),nb,2*ni+1,XI(1,1),XE(1,3))
                H(1)=-E(1,1)*DXDXI(1,1,ni,1)-E(2,1)*DXDXI(2,1,ni,1)-
     '            E(3,1)*DXDXI(3,1,ni,1)
                DHDXI(1,1)=DXDXI(1,1,ni,1)*DXDXI(1,1,ni,1)+
     '            DXDXI(2,1,ni,1)*DXDXI(2,1,ni,1)+DXDXI(3,1,ni,1)*
     '            DXDXI(3,1,ni,1)-E(1,1)*DXDXI(1,2,ni,1)-E(2,1)*
     '            DXDXI(2,2,ni,1)-E(3,1)*DXDXI(3,2,ni,1)
                CALL ASSERT(DABS(DHDXI(1,1)).GT.RDELTA,
     '            '>>Zero hessian?',ERROR,*9999)
                IF(DHDXI(1,1).LT.0.0d0) THEN
                  ALPHA=-1.0d0
                ELSE
                  ALPHA=1.0d0
                ENDIF
                DELTAXI(1)=-ALPHA*H(1)/DHDXI(1,1)
                XI(ni,2)=XI(ni,1)
                XI(ni,1)=XI(ni,1)+DELTAXI(1)
                IF(XI(ni,1).LT.0.0d0) XI(ni,1)=0.0d0
                IF(XI(ni,1).GT.1.0d0) XI(ni,1)=1.0d0
                ERR=DABS(XI(ni,1)-XI(ni,2))
                CONVERGED=ERR.LE.TOLERANCE
                NUMITS=NUMITS+1
                X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI(1,1),XE(1,1))
                X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI(1,1),XE(1,2))
                X(3,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI(1,1),XE(1,3))
                E(1,1)=XPFP(1)-X(1,1)
                E(2,1)=XPFP(2)-X(2,1)
                E(3,1)=XPFP(3)-X(3,1)
                NORME=DSQRT(E(1,1)*E(1,1)+E(2,1)*E(2,1)+E(3,1)*E(3,1))
                IF(NUMITS.LT.MAXITS.AND..NOT.CONVERGED) THEN
                  DXDXI(1,1,ni,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,2*ni,XI(1,1),XE(1,1))
                  DXDXI(2,1,ni,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,2*ni,XI(1,1),XE(1,2))
                  DXDXI(3,1,ni,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,2*ni,XI(1,1),XE(1,1))
                ENDIF
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_9)
                  WRITE(OP_STRING,'(/'' Iteration '',I2)') NUMITS
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' Norm. error : '',D12.4,'
     '              //''', Tolerance : '',D12.4)') ERR,TOLERANCE
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' OLD XI   : '',D12.4)') XI(1,2)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' DELTAXI  : '',D12.4)')
     '              DELTAXI(1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' XI       : '',D12.4)') XI(1,1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' Distance : '',D12.4)') NORME
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_9)
                ENDIF
              ENDDO
            ELSE
              XI(1,1)=XI(1,2)
              XI(2,1)=XI(2,2)
            ENDIF
          ENDIF
        ENDIF
        XIMIN(1)=XI(1,1)
        XIMIN(2)=XI(2,1)
        XI(1,2)=1.0d0
        XI(2,2)=XI(2,1)
        X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,1))
        X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,2))
        X(3,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,3))
        XI(1,2)=0.0d0
        XI(2,2)=XI(2,1)
        X(1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,1))
        X(2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,2))
        X(3,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,3))
        DIST=DSQRT((X(1,1)-X(1,2))*(X(1,1)-X(1,2))+(X(2,1)-X(2,2))*
     '    (X(2,1)-X(2,2))+(X(3,1)-X(3,2))*(X(3,1)-X(3,2)))
        IF((DIST.LE.RDELTA).AND.(NW.EQ.14.OR.NW.EQ.15.OR.NW.EQ.16)) THEN
C*** Apex of a simplex or sector element has been found. Set D Large so
C*** that nothing is done in GAUSS11 for this xi direction
          D(1)=10.0d0
        ELSE
          CALL ASSERT(DIST.GE.RDELTA,'>>Zero width element?',
     '      ERROR,*9999)
          D(1)=2.0d0*NORME/DIST
        ENDIF
        XI(1,2)=XI(1,1)
        XI(2,2)=1.0d0
        X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,1))
        X(2,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,2))
        X(3,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,3))
        XI(1,2)=XI(1,1)
        XI(2,2)=0.0d0
        X(1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,1))
        X(2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,2))
        X(3,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI(1,2),XE(1,3))
        DIST=DSQRT((X(1,1)-X(1,2))*(X(1,1)-X(1,2))+(X(2,1)-X(2,2))*
     '    (X(2,1)-X(2,2))+(X(3,1)-X(3,2))*(X(3,1)-X(3,2)))
        IF((DIST.LE.RDELTA).AND.(NW.EQ.14.OR.NW.EQ.15.OR.NW.EQ.16)) THEN
C*** Apex of a simplex or sector element has been found. Set D Large so
C*** that nothing is done in GAUSS11 for this xi direction
          D(2)=10.0d0
        ELSE
          CALL ASSERT(DIST.GE.RDELTA,'>>Zero width element?',
     '      ERROR,*9999)
          D(2)=2.0d0*NORME/DIST
        ENDIF
      ENDIF
      IF(NUMITS.GT.0.AND..NOT.CONVERGED) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_10)
C$      IF(.FALSE.) THEN
        WRITE(OP_STRING,'('' >>Warning: Minimum point iterations have '
     '    //'not converged for Telles rule!'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C$      ENDIF
CC$OMP END CRITICAL(TELLES_10)
      ENDIF
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(TELLES_11)
        WRITE(OP_STRING,'(/'' Telles Rule Parameters:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Rmin='',D12.4)') NORME
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' D(1)='',D12.4,'', D(2)='',D12.4)') D(1),
     '    D(2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(TELLES_11)
      ENDIF

      CALL EXITS('TELLES')
      RETURN
9999  CALL ERRORS('TELLES',ERROR)
      CALL EXITS('TELLES')
      RETURN 1
      END


C!!! cpb 28/6/96 Adding node based outer integration loops hence
C!!! subroutines have been renamed i.e. XEGKGQx -> XEPGKGQx,
C!!! XPGKGQ -> XEGKGQ and a new routine XPGKGQ has been introduced.

