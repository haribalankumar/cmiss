      SUBROUTINE BRANCH3(nq,nr,NXQ,NYNQ,CQ,XQ,YQ,
     '  HALF_TIME_STEP,ERROR,*)


C#### Subroutine: BRANCH3
C###  Description:
C###    BRANCH3 calculates the pressure of the grid points where
C###    transition occurs from an artery to a vein. This grid point
C###    is treated as a reservoir and fluid accumulated together with
C###    elastic properties will detemine the pressure. The velocity
C###    of such a point is equal to zero. Inflow/outflow is determined
C###    by the upstream,downstream pressures and the pressure of the
C###    reservoir itself. Since the radii of these points are undefined
C###    spatial derivatives of radius are computed using one-sided,
C###    second order difference scheme.
C**** Written by Kumar Mithraratne, Nov. 2001.

      IMPLICIT NONE
      INCLUDE 'b12.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER  nq,nr,NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNQ(NHM,NQM,0:NRCM)
      REAL*8  CQ(NMM,NQM),XQ(NJM,NQM),YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
      LOGICAL HALF_TIME_STEP

!     Local Variables
      INTEGER  nj,njtot,nq_up1,nq_up2,nq_up3,nq_dw1,nq_dw2,nq_dw3,
     '  ny_p1,ny_p2,ny_p3,ny_p4,ny_p5,
     '  ny_r1,ny_r2,ny_r3,ny_r4,ny_r5,ny_r6,ny_r7,
     '  ny_v2,ny_v3,ny_v4,ny_v5,
     '  yp_index,yp_value
      REAL*8 Cart,Cvein,DELTA_X,DELX_MEAN,dPdt0,dPdt1,P1,P2,P3,P4,
     '  P5,R1,R2,R3,R5,R6,R7,Rart,Rvein,TERM1,TERM2,TERM3,
     '  TERM4,TERM5,V2,V3,V5
C      INTEGER ny_p6,ny_p7,ny_v1,ny_v6,ny_v7
C      REAL*8 P6,P7,V1,V4,V6,V7

      CALL ENTERS('BRANCH3',*9999)

      IF(HALF_TIME_STEP) THEN
        yp_index=1 !index for previuos time step
        yp_value=5 !index for current half-time step
      ELSE   ! full-time step
        yp_index=5 !index for previous half-time step
        yp_value=1 !index for previous full-time step for RHS variables
                   !          current full time step for LHS variables
      ENDIF

      njtot=NJ_LOC(NJL_GEOM,0,nr)

      IF(HALF_TIME_STEP) THEN

        nq_up1=NXQ(-1,1,nq,1)
        nq_up2=NXQ(-1,1,nq_up1,1)
        nq_up3=NXQ(-1,1,nq_up2,1)
        nq_dw1=NXQ(1,1,nq,1)
        nq_dw2=NXQ(1,1,nq_dw1,1)
        nq_dw3=NXQ(1,1,nq_dw2,1)

        ny_p1=NYNQ(1,nq_up3,0)
        ny_p2=NYNQ(1,nq_up2,0)
        ny_p3=NYNQ(1,nq_up1,0)
        ny_p4=NYNQ(1,nq,0)
        ny_p5=NYNQ(1,nq_dw1,0)
C       ny_p6=NYNQ(1,nq_dw2,0)
C       ny_p7=NYNQ(1,nq_dw3,0)
        ny_r1=NYNQ(2,nq_up3,0)
        ny_r2=NYNQ(2,nq_up2,0)
        ny_r3=NYNQ(2,nq_up1,0)
        ny_r4=NYNQ(2,nq,0)
        ny_r5=NYNQ(2,nq_dw1,0)
        ny_r6=NYNQ(2,nq_dw2,0)
        ny_r7=NYNQ(2,nq_dw3,0)
C       ny_v1=NYNQ(3,nq_up3,0)
        ny_v2=NYNQ(3,nq_up2,0)
        ny_v3=NYNQ(3,nq_up1,0)
        ny_v4=NYNQ(3,nq,0)
        ny_v5=NYNQ(3,nq_dw1,0)
C       ny_v6=NYNQ(3,nq_dw2,0)
C       ny_v7=NYNQ(3,nq_dw3,0)

        P1=YQ(ny_p1,yp_index)
        P2=YQ(ny_p2,yp_index)
        P3=YQ(ny_p3,yp_index)
        P4=YQ(ny_p4,yp_index)
        P5=YQ(ny_p5,yp_index)
C       P6=YQ(ny_p6,yp_index)
C       P7=YQ(ny_p7,yp_index)
        R1=YQ(ny_r1,yp_index)
        R2=YQ(ny_r2,yp_index)
        R3=YQ(ny_r3,yp_index)
C       R4=YQ(ny_r4,yp_index)
        R5=YQ(ny_r5,yp_index)
        R6=YQ(ny_r6,yp_index)
        R7=YQ(ny_r7,yp_index)
C       V1=YQ(ny_v1,yp_index)
        V2=YQ(ny_v2,yp_index)
        V3=YQ(ny_v3,yp_index)
C       V4=YQ(ny_v4,yp_index)
        V5=YQ(ny_v5,yp_index)
C       V6=YQ(ny_v6,yp_index)
C       V7=YQ(ny_v7,yp_index)

C half-time step values at half point upstream of terminal grid point:

        DELTA_X=0.0d0
        DO nj=1,njtot
          DELTA_X=DELTA_X+((XQ(nj,nq)-XQ(nj,nq_up1))**2.0d0)
        ENDDO
        DELTA_X=DELTA_X**0.5d0

        DELX_MEAN=0.0d0
        DO nj=1,njtot
          DELX_MEAN=DELX_MEAN+(XQ(nj,NQ_START(nr))-
     '      XQ(nj,NXQ(1,1,NQ_START(nr),1)))**2.0d0
        ENDDO
        DELX_MEAN=DELX_MEAN**0.5d0

C This is not required if the element (1D arc) is discretised with
C the uniform arc-length spacing rather than equal xi spacing.
        DELTA_X=MAX(DELTA_X,0.75d0*DELX_MEAN)

        TERM1=R3
        TERM2=(-0.125d0)*(TINCR/DELTA_X)*(3.0d0*R3-4.0d0*R2+R1)*V3
        TERM3=0.25d0*(TINCR/DELTA_X)*R3*V3
        YQ(ny_r3,yp_value)=TERM1+TERM2+TERM3                ! Radius

        TERM1=0.5d0*V3
        TERM2=0.25d0*(2.0d0*CQ(3,nq_up1)-1.0d0)*
     '    (TINCR/DELTA_X)*(V3**2.0d0)
        TERM3=(-0.125d0)*(CQ(3,nq_up1)-1.0d0)*(TINCR/DELTA_X)*
     '    (3.0d0*R3-4.0d0*R2+R1)*(V2**2.0d0)/R3
        TERM4=(-0.5d0)*(1.0d0/CQ(1,nq_up1))*(TINCR/DELTA_X)*(P4-P3)
        TERM5=(-0.5d0*CQ(2,nq_up1)*CQ(3,nq_up1))/
     '    (CQ(3,nq_up1)-1.0d0)*TINCR*V3/(R3**2.0d0)
        YQ(ny_v3,yp_value)=TERM1+TERM2+TERM3+TERM4+TERM5    ! Velocity

        TERM1=CQ(5,nq_up1)
        TERM2=(YQ(ny_r3,yp_value)/CQ(4,nq_up1))**
     '    CQ(6,nq_up1)-1.00d0
          YQ(ny_p3,yp_value)=TERM1*TERM2                    ! Pressure

C half-time step values at half point downstream of terminal grid point:

        DELTA_X=0.0d0
        DO nj=1,njtot
          DELTA_X=DELTA_X+((XQ(nj,nq_dw1)-XQ(nj,nq))**2.0d0)
        ENDDO
        DELTA_X=DELTA_X**0.5d0

        DELX_MEAN=0.0d0
        DO nj=1,njtot
          DELX_MEAN=DELX_MEAN+(XQ(nj,NQ_START(nr))-
     '      XQ(nj,NXQ(1,1,NQ_START(nr),1)))**2.0d0
        ENDDO
        DELX_MEAN=DELX_MEAN**0.5d0

C This is not required if the element (1D arc) is discretised with
C the uniform arc-length spacing rather than equal xi spacing.
        DELTA_X=MAX(DELTA_X,0.75d0*DELX_MEAN)

        TERM1=R5
        TERM2=(-0.125d0)*(TINCR/DELTA_X)*(-3.0d0*R5+4.0d0*R6-R7)*V5
        TERM3=(-0.25d0)*(TINCR/DELTA_X)*R5*V5
        YQ(ny_r4,yp_value)=TERM1+TERM2+TERM3                ! Radius

        TERM1=0.5d0*V5
        TERM2=(-0.25d0)*(2.0d0*CQ(3,nq_dw1)-1.0d0)*
     '    (TINCR/DELTA_X)*(V5**2.0d0)
        TERM3=(-0.125d0)*(CQ(3,nq_dw1)-1.0d0)*(TINCR/DELTA_X)*
     '    (-3.0d0*R5+4.0d0*R6-R7)*(V5**2.0d0)/R5
        TERM4=(-0.5d0)*(1.0d0/CQ(1,nq_dw1))*(TINCR/DELTA_X)*(P5-P4)
        TERM5=(-0.5d0*CQ(2,nq_dw1)*CQ(3,nq_dw1))/
     '    (CQ(3,nq_dw1)-1.0d0)*TINCR*V5/(R5**2.0d0)
        YQ(ny_v4,yp_value)=TERM1+TERM2+TERM3+TERM4+TERM5    ! Velocity

        TERM1=CQ(5,nq_dw1)
        TERM2=(YQ(ny_r4,yp_value)/CQ(4,nq_dw1))**
     '    CQ(6,nq_dw1)-1.00d0
        YQ(ny_p4,yp_value)=TERM1*TERM2                      ! Pressure

      ELSE  ! full time step

        nq_up1=NXQ(-1,1,nq,1)
        nq_dw1=NXQ(1,1,nq,1)

        ny_p1=NYNQ(1,nq_up1,0)
        ny_p2=NYNQ(1,nq,0)
        ny_p3=NYNQ(1,nq_dw1,0)
        ny_r2=NYNQ(2,nq,0)
        ny_v2=NYNQ(3,nq,0)

        P1=YQ(ny_p1,yp_value)
        P2=YQ(ny_p2,yp_value)
        P3=YQ(ny_p3,yp_value)

C The values Rart and Rvein depend on the flow and Cart and Cvein depend
C on the pressure. In order to calculate these parameters flow/velocity and
C pressure from the previous time step could be used. Currently, however
C following values are used.

        Rart=1000.0d0
        Cart=100.0d0
        Rvein=1.0d0
        Cvein=100.0d0

        TERM1=(P1-P2)/(Rart*CQ(1,nq)*Cart)
        TERM2=(-1.0d0)*(P2-P3)/(Rvein*CQ(1,nq)*Cvein)
        dPdt0=TERM1+TERM2
        P4=P2+dPdt0*TINCR
        TERM1=(P1-P4)/(Rart*CQ(1,nq)*Cart)
        TERM2=(-1.0d0)*(P4-P3)/(Rvein*CQ(1,nq)*Cvein)
        dPdt1=TERM1+TERM2
        YQ(ny_p2,yp_value)=P2+0.5d0*(dPdt0+dPdt1)*TINCR

        YQ(ny_r2,yp_value)=0.000d0  ! radius at a terminal point is undefined.
        YQ(ny_v2,yp_value)=0.000d0  ! velocity at a terminal point is zero.

      ENDIF ! half-time/full-time


      CALL EXITS('BRANCH3')
      RETURN
 9999 CALL ERRORS('BRANCH3',ERROR)
      CALL EXITS('BRANCH3')
      RETURN 1
      END
      

