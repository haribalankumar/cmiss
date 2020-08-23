      SUBROUTINE ANGLE(TN1,TN2,THETA,XN1,XN2,ERROR,*)

C#### Subroutine: ANGLE
C###  Description:
C###    ANGLE calculates the angle between 2 planes whose normals are
C###    XN1 and XN2 and whose tangents are TN1 and TN2 respectively.

C**** For each normal two tangents are given - one in the +/- xi1
C**** direction of the plane and the other in the +/- xi2 direction
C**** of the plane.  Both are needed to allow for the case when the
C**** first tangent is end on after transforming to the (x_hat,y_hat)
C**** space.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      REAL*8 TN1(3,*),TN2(3,*),THETA,XN1(*),XN2(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,nj
      REAL*8 A_MULT,B,B_MULT,CROSS_LENGTH,DOT,DOT2,DOT3,
     '     S(3),S_LEN,T_STAR(3),TN2_HAT(2,3),
     '     XN2_HAT(3),XN3(3)

      CALL ENTERS('ANGLE',*9999)
      DOT=0.0D0
      DO nj=1,NJ_LOC(NJL_GEOM,0,1) !1 is temporary RGB
        DOT=XN1(nj)*XN2(nj)+DOT
      ENDDO
      IF(DABS(DOT).LT.RDELTA) DOT=0.0D0
      THETA=DACOS(DOT) !Returns angle between 0 and PI
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ANGLE_1)
        WRITE(OP_STRING,'('' Dot product of xn1 and xn2='',D12.4)') DOT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Inverse cos of this angle='',D12.4)')
     '    THETA
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(ANGLE_1)
      ENDIF
      IF(NJ_LOC(NJL_GEOM,0,1).EQ.2) THEN !2d space. IN this case t_star=tn1
        DO nj=1,NJ_LOC(NJL_GEOM,0,1)
          T_STAR(nj)=TN1(nj,1)
        ENDDO
      ELSE !3d space.  Need to rotate axes.
        !Find xn1 cross xn2
        XN3(1)=XN1(2)*XN2(3)-XN1(3)*XN2(2)
        XN3(2)=XN1(3)*XN2(1)-XN1(1)*XN2(3)
        XN3(3)=XN1(1)*XN2(2)-XN1(2)*XN2(1)
        CROSS_LENGTH=DSQRT(XN3(1)**2+XN3(2)**2+XN3(3)**2)
        DO nj=1,3
          XN3(nj)=XN3(nj)/CROSS_LENGTH
        ENDDO
        !Find T_STAR - the tangent vector in the plane containing xn1
        !and xn2.
        DOT2=0.0D0
        DOT3=0.0D0
        DO nj=1,NJ_LOC(NJL_GEOM,0,1)
          DOT2=DOT2+XN3(nj)*TN1(nj,1) !xn3.tn1_1
          DOT3=DOT3+XN3(nj)*TN1(nj,2) !xn3.tn1_2
        ENDDO
        IF(DABS(DOT2).LE.RDELTA) THEN
          DO nj=1,NJ_LOC(NJL_GEOM,0,1)
            T_STAR(nj)=TN1(nj,1)
          ENDDO
        ELSE
          B=DOT3/DOT2
          S_LEN=0.0D0
          DO nj=1,NJ_LOC(NJL_GEOM,0,1)
            S(nj)=TN1(nj,2)-B*TN1(nj,1)
            S_LEN=S_LEN+S(nj)*S(nj)
          ENDDO
          S_LEN=DSQRT(S_LEN)
          IF(S_LEN.LE.RDELTA) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ANGLE_2)
            WRITE(OP_STRING,'('' >>S=0 in ANGLE subroutine ????'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(ANGLE_2)
            GOTO 9999
          ENDIF
          B_MULT=1.0D0/S_LEN
          A_MULT=-B*B_MULT
          IF(A_MULT.LT.0.0D0)A_MULT=-A_MULT !Keep A_MULT and B_MULT > 0
          DO nj=1,NJ_LOC(NJL_GEOM,0,1)
            T_STAR(nj)=A_MULT*TN1(nj,1)+B_MULT*TN1(nj,2)
          ENDDO
        ENDIF
      ENDIF !End of 2d/3d choice.
      !Transform to a new coordinate system (x_hat,y_hat,z_hat) with the
      !x_hat axis lying along XN1 and the y_hat axis given by t_star.
      !To see if an external or internal angle is required we check
      !whether XN2 or TN2 lie in the first quadrant of (x-hat,y_hat)
      !space.
      XN2_HAT(1)=DOT
      XN2_HAT(2)=0.0D0
      DO i=1,2
        DO j=1,2
          TN2_HAT(i,j)=0.0D0
        ENDDO
      ENDDO
      DO nj=1,NJ_LOC(NJL_GEOM,0,1)
        XN2_HAT(2)=XN2_HAT(2)+T_STAR(nj)*XN2(nj)
        DO i=1,NJ_LOC(NJL_GEOM,0,1)-1
          TN2_HAT(1,i)=TN2_HAT(1,i)+XN1(nj)*TN2(nj,i)
          TN2_HAT(2,i)=TN2_HAT(2,i)+T_STAR(nj)*TN2(nj,i)
        ENDDO
      ENDDO
      !If either XN2_HAT or TN2_HAT lie in the first quadrant require
      !external angle
      IF((XN2_HAT(1).GT.0.0D0.AND.XN2_HAT(2).GT.0.0D0).OR.
     ' (TN2_HAT(1,1).GT.0.0D0.AND.TN2_HAT(2,1).GT.0.0D0).OR.
     ' (TN2_HAT(1,2).GT.0.0D0.AND.TN2_HAT(2,2).GT.0.0D0)) THEN
        THETA=PI+THETA !External angle
      ELSE
        THETA=PI-THETA !Internal angle
      ENDIF
      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(ANGLE_3)
        WRITE(OP_STRING,'('' Required angle='',D12.4)') THETA
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(ANGLE_3)
      ENDIF
      CALL EXITS('ANGLE')
      RETURN
9999  CALL ERRORS('ANGLE',ERROR)
      CALL EXITS('ANGLE')
      RETURN 1
      END


