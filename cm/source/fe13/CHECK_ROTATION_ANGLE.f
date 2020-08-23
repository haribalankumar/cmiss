      SUBROUTINE CHECK_ROTATION_ANGLE(ne0,np00,np0,np1,np2,np3,np4,
     &  ROTATION_LIMIT,XP,ERROR,*)

C####  Subroutine: CHECK_ROTATION_ANGLE
C     Calculates using quaternions. For angle ROTATION_ANGLE and unit
C     vector a,b,c , calculate rotation matrix for arbitrary point.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter values
      INTEGER ne0,np00,np0,np1,np2,np3,np4
      REAL*8 ROTATION_LIMIT,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj
      REAL*8 ANGLE0,ANGLE,AXIS(3),DIRECTION(3),NRML(3),NRML_PARENT(3),
     &  ROT_ANGLE,U(3),V(3),Q0,Q1,Q2,Q3,Q(3,3),X(3),temp
      LOGICAL COMPLETE
      INTEGER IT,ITMAX
      REAL*8 CALC_ANGLE,SCALAR,ANGLE_BETWEEN,C,S,T,length,XPOINT(3,5),
     &  norm_1(3),norm_2(3)

      CALL ENTERS('CHECK_ROTATION_ANGLE',*9999)

      ITmAX=10
      
C...parent branching plane, cross-product of parent and grand-parent
      DO nj=1,3
        U(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1) !parent
        V(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np0) !grandparent
      ENDDO !nj
      CALL NORMALISE2(3,U,temp,ERROR,*9999) !unit vector
      CALL NORMALISE2(3,V,temp,ERROR,*9999) !unit vector
      CALL CROSS(U,V,NRML_PARENT) !calculate branching plane
      CALL NORMALISE2(3,NRML_PARENT,temp,ERROR,*9999) !unit vector
      
C...current branching plane
      DO nj=1,3
        U(nj)=XP(1,1,nj,np3)-XP(1,1,nj,np2) !branch
        V(nj)=XP(1,1,nj,np4)-XP(1,1,nj,np2) !sibling
      ENDDO !nj
      CALL NORMALISE2(3,U,temp,ERROR,*9999) !unit vector
      CALL NORMALISE2(3,V,temp,ERROR,*9999) !unit vector
      CALL CROSS(U,V,NRML) !calculate branching plane
      CALL NORMALISE2(3,NRML,temp,ERROR,*9999) !unit vector

C...angle between branching planes
      ANGLE=CALC_ANGLE(NRML,NRML_PARENT)
      ANGLE_BETWEEN=ANGLE
      IF(ANGLE_BETWEEN.GT.PI/2.d0)THEN
        ANGLE_BETWEEN=ANGLE_BETWEEN-PI
        ANGLE=ANGLE_BETWEEN
      ENDIF
c      ANGLE=PI/2.d0-ANGLE
      ANGLE=PI/2.d0+ANGLE

      IF(DABS(ANGLE_BETWEEN).GT.ROTATION_LIMIT.AND.DABS(ANGLE_BETWEEN)
     &  .LT.PI/2.d0-ROTATION_LIMIT)THEN
        IF(ANGLE.LT.0.d0)THEN
          ROT_ANGLE=-(ANGLE+ROTATION_LIMIT)
        ELSE
          ROT_ANGLE=-(ANGLE-ROTATION_LIMIT)
        ENDIF

c        write(*,*) 'angle_between=',angle_between*180.d0/PI,
c     &    ', rot_angle=',rot_angle*180.d0/PI,', angle=',angle*180.d0/PI
        
        ANGLE0=ANGLE

C...if the difference in angles is not correct, rotate branches
        DO nj=1,3
          AXIS(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
        ENDDO !nj
        CALL NORMALISE2(3,AXIS,temp,ERROR,*9999) !unit vector for direction
        
        Q0=DCOS(ROT_ANGLE/2.D0)
        Q1=DSIN(ROT_ANGLE/2.D0)*AXIS(1)
        Q2=DSIN(ROT_ANGLE/2.D0)*AXIS(2)
        Q3=DSIN(ROT_ANGLE/2.D0)*AXIS(3)
        
        Q(1,1)=Q0**2+Q1**2-Q2**2-Q3**2
        Q(1,2)=2*(Q1*Q2-Q0*Q3)
        Q(1,3)=2*(Q1*Q3+Q0*Q2)
        Q(2,1)=2*(Q2*Q1+Q0*Q3)
        Q(2,2)=Q0**2-Q1**2+Q2**2-Q3**2
        Q(2,3)=2*(Q2*Q3-Q0*Q1)
        Q(3,1)=2*(Q3*Q1-Q0*Q2)
        Q(3,2)=2*(Q3*Q2+Q0*Q1)
        Q(3,3)=Q0**2-Q1**2-Q2**2+Q3**2

        DO nj=1,3
          X(nj)=XP(1,2,nj,np3) !unit vector
        ENDDO !nj
        length=DSQRT((XP(1,1,1,np3)-XP(1,1,1,np2))**2+(XP(1,1,2,np3)
     &    -XP(1,1,2,np2))**2+(XP(1,1,3,np3)-XP(1,1,3,np2))**2)
        
        DO nj=1,3
          XP(1,1,nj,np3)=XP(1,1,nj,np2)+length*(Q(nj,1)*X(1)+Q(nj,2)
     &      *X(2)+Q(nj,3)*X(3))
        ENDDO !nj
        DO nj=1,3 !unit vector for branch direction
          DIRECTION(nj)=(XP(1,1,nj,np3)-XP(1,1,nj,np2))
        ENDDO !nj
        CALL NORMALISE2(3,DIRECTION,temp,ERROR,*9999)
        DO nj=1,3 !unit vector for branch direction
          XP(1,2,nj,np3)=DIRECTION(nj)
        ENDDO !nj

        DO nj=1,3
          X(nj)=XP(1,2,nj,np4) !unit vector
        ENDDO !nj
        length=DSQRT((XP(1,1,1,np3)-XP(1,1,1,np2))**2+(XP(1,1,2,np3)
     &    -XP(1,1,2,np2))**2+(XP(1,1,3,np3)-XP(1,1,3,np2))**2)
        DO nj=1,3
          XP(1,1,nj,np4)=XP(1,1,nj,np2)+length*(Q(nj,1)*X(1)+Q(nj,2)
     &      *X(2)+Q(nj,3)*X(3))
        ENDDO !nj
        DO nj=1,3 !unit vector for branch direction
          DIRECTION(nj)=(XP(1,1,nj,np4)-XP(1,1,nj,np2))
        ENDDO !nj
        CALL NORMALISE2(3,DIRECTION,temp,ERROR,*9999)
        DO nj=1,3 !unit vector for branch direction
          XP(1,2,nj,np4)=DIRECTION(nj)
        ENDDO !nj

        DO nj=1,3
          U(nj)=XP(1,2,nj,np3) !direction of a branch
          V(nj)=XP(1,2,nj,np4) !direction of its sibling
        ENDDO !nj
        CALL NORMALISE2(3,U,temp,ERROR,*9999) !unit vector
        CALL NORMALISE2(3,V,temp,ERROR,*9999) !unit vector
        CALL CROSS(U,V,NRML) !calculate branching plane
        CALL NORMALISE2(3,NRML,temp,ERROR,*9999) !unit vector
        
C...angle between branching planes
        ANGLE=CALC_ANGLE(NRML,NRML_PARENT)
c        write(*,*) 'angle after first attempt=',angle*180.d0/PI
        
C should find that 90degrees minus new angle is within the limit range
        IF(DABS(ANGLE_BETWEEN).GT.ROTATION_LIMIT.AND.DABS(ANGLE_BETWEEN)
     &    .LT.PI/2.d0-ROTATION_LIMIT)THEN
          COMPLETE=.TRUE.
        ELSE
          COMPLETE=.FALSE.
        ENDIF
c        IF(DABS(ROTATION_LIMIT-DABS(PI/2.d0-ANGLE)).GT.0.001d0)THEN
c          COMPLETE=.FALSE.
c        ELSE
c          COMPLETE=.TRUE.
c        ENDIF

        DO WHILE(.NOT.COMPLETE)
          IT=IT+1
          ANGLE_BETWEEN=ANGLE
          ANGLE=PI/2.d0-ANGLE
          IF(ANGLE.LT.0.d0)THEN
            ROT_ANGLE=-(ANGLE+ROTATION_LIMIT)
          ELSE
            ROT_ANGLE=-(ANGLE-ROTATION_LIMIT)
          ENDIF
c          write(*,*) 'while loop: rot_angle =',rot_angle*180.d0/PI
          
          Q0=DCOS(ROT_ANGLE/2.D0)
          Q1=DSIN(ROT_ANGLE/2.D0)*AXIS(1)
          Q2=DSIN(ROT_ANGLE/2.D0)*AXIS(2)
          Q3=DSIN(ROT_ANGLE/2.D0)*AXIS(3)
          
          Q(1,1)=Q0**2+Q1**2-Q2**2-Q3**2
          Q(1,2)=2*(Q1*Q2-Q0*Q3)
          Q(1,3)=2*(Q1*Q3+Q0*Q2)
          Q(2,1)=2*(Q2*Q1+Q0*Q3)
          Q(2,2)=Q0**2-Q1**2+Q2**2-Q3**2
          Q(2,3)=2*(Q2*Q3-Q0*Q1)
          Q(3,1)=2*(Q3*Q1-Q0*Q2)
          Q(3,2)=2*(Q3*Q2+Q0*Q1)
          Q(3,3)=Q0**2-Q1**2-Q2**2+Q3**2
          
          DO nj=1,3
            X(nj)=XP(1,2,nj,np3) !unit vector
          ENDDO !nj
          length=DSQRT((XP(1,1,1,np3)-XP(1,1,1,np2))**2+(XP(1,1,2,np3)
     &      -XP(1,1,2,np2))**2+(XP(1,1,3,np3)-XP(1,1,3,np2))**2)
          
          DO nj=1,3
            XP(1,1,nj,np3)=XP(1,1,nj,np2)+length*(Q(nj,1)*X(1)+Q(nj,2)
     &        *X(2)+Q(nj,3)*X(3))
          ENDDO !nj
          DO nj=1,3 !unit vector for branch direction
            DIRECTION(nj)=(XP(1,1,nj,np3)-XP(1,1,nj,np2))
          ENDDO !nj
          CALL NORMALISE2(3,DIRECTION,temp,ERROR,*9999)
          DO nj=1,3 !unit vector for branch direction
            XP(1,2,nj,np3)=DIRECTION(nj)
          ENDDO !nj
          
          DO nj=1,3
            X(nj)=XP(1,2,nj,np4) !unit vector
          ENDDO !nj
          length=DSQRT((XP(1,1,1,np3)-XP(1,1,1,np2))**2+(XP(1,1,2,np3)
     &      -XP(1,1,2,np2))**2+(XP(1,1,3,np3)-XP(1,1,3,np2))**2)
          DO nj=1,3
            XP(1,1,nj,np4)=XP(1,1,nj,np2)+length*(Q(nj,1)*X(1)+Q(nj,2)
     &        *X(2)+Q(nj,3)*X(3))
          ENDDO !nj
          DO nj=1,3 !unit vector for branch direction
            DIRECTION(nj)=(XP(1,1,nj,np4)-XP(1,1,nj,np2))
          ENDDO !nj
          CALL NORMALISE(3,DIRECTION,ERROR,*9999)
          DO nj=1,3 !unit vector for branch direction
            XP(1,2,nj,np4)=DIRECTION(nj)
          ENDDO !nj
          
          DO nj=1,3
            U(nj)=XP(1,2,nj,np3) !direction of a branch
            V(nj)=XP(1,2,nj,np4) !direction of its sibling
          ENDDO !nj
          CALL NORMALISE2(3,U,temp,ERROR,*9999) !unit vector
          CALL NORMALISE2(3,V,temp,ERROR,*9999) !unit vector
          CALL CROSS(U,V,NRML) !calculate branching plane
          CALL NORMALISE2(3,NRML,temp,ERROR,*9999) !unit vector
          
C...angle between branching planes
          ANGLE=CALC_ANGLE(NRML,NRML_PARENT)

c          write(*,*) 'testing:',(ROTATION_LIMIT-DABS(PI/2.d0-ANGLE))
c     &      *180.d0/PI
          
          IF(DABS(ROTATION_LIMIT-DABS(PI/2.d0-ANGLE)).LE.0.001d0)THEN
            COMPLETE=.TRUE.
          ENDIF

          IF(IT.GT.ITMAX)THEN
            WRITE(*,*) 'WARNING!!!! rotation angle = ',ANGLE*180.d0/PI
          ENDIF

        ENDDO !do while not found
c        pause
C.......Alternate calculation for the rotation angle
        CALL ROTATION_ANGLE(np1,np2,np00,np3,np4,ANGLE,XP,ERROR,*9999)
      ELSE
c        write(*,*) 'Not',ANGLE_BETWEEN*180.d0/PI
      ENDIF
        
      CALL EXITS('CHECK_ROTATION_ANGLE')
      RETURN
 9999 CALL ERRORS('CHECK_ROTATION_ANGLE',ERROR)
      CALL EXITS('CHECK_ROTATION_ANGLE')
      RETURN 1
      END

      

