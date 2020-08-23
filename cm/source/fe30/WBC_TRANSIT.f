      SUBROUTINE WBC_TRANSIT(flow,length,P1,P2,P_cr,R_WBC,R_vessel,
     '  t,te,tp,ERROR,*)

C#### Subroutine: WBC_TRANSIT
C###  Description:
C###    WBC_TRANSIT calculates the white blood cell (WBC) transit time
C###    through the pulmonary capillary network. Total transit time
C###    is the sum of entrance time (te) and passing time (tp).
C###    Model based on neutrophils which make up about 70% of WBCs.

C***   Created by KSB, 21st June 2002
C***   Based on model by Huang et al., 2001.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b00.cmn'

!     Parameter list
      REAL*8 flow,length,P1,P2,P_cr,R_WBC,R_vessel,t,te,tp
      CHARACTER ERROR*(*)
!     Local variables
      REAL*8 DP,I,m,R1,R2,R_WBC_in,R_min,Ro_star,te_star,
     '  u_cyto

      CALL ENTERS('WBC_TRANSIT',*9999)

      DP=P2-P1 !(Pa) pressure drop over element
      m=6.d0 !from Huang:2001
      u_cyto=135.d0 !Pa.s cytoplasmic viscos. (value from Needham:1990)
      Ro_star=R_WBC/R_vessel !R*=Ro* at t*=0
      IF(Ro_star.GE.1.d0.AND.Ro_star.LT.1.2d0) THEN
        te_star=m*(1.d0/3.d0+(2.d0/3.d0)*(Ro_star**3.d0)-Ro_star**2.d0+
     '    (2.d0/3.d0)*(Ro_star**2.d0-1.d0)**(3.d0/2.d0)+
     '    (Ro_star**2.d0-1.d0)**(1.d0/2.d0)-(Ro_star**2.d0*
     '    DSIN(ACOS(1.d0/Ro_star)))+DLOG((Ro_star*DSIN(ACOS(1/Ro_star))+
     '    Ro_star)/(Ro_star+(Ro_star**2.d0-1.d0)**(1.d0/2.d0))))!eq (41)
      ELSE IF(Ro_star.LT.1.0d0) THEN
        te_star=0.d0 !if R_WBC < R_vessel then entrance time =0.d0
      ENDIF
      IF(Ro_star.LT.1.2d0) THEN !eqn (44), Huang (2001)
        te=u_cyto/(DP-P_cr)*te_star !entrance time (s)
      ELSE IF(Ro_star.GE.1.2d0) THEN
        te=0.107d0*m*u_cyto/(DP-33.d0)*10.d0**(3.68d0*(Ro_star-1.2d0))
      ENDIF
      IF(te.LT.0.d0) te=0.d0 !if te negative, make te=0
C... to consider capillary with different R at inlet & outlet
C... currently radius defined over element therefore constant
C      R_x=R1+((R2-R1)/length*x
c mht R_x set but not used 
c     R_x=R_vessel !currently R is constant over element
      R1=R_vessel
      R2=R_vessel

      IF(R1.EQ.R2) THEN !vessel radius constant
        IF(R_WBC.LE.R1) THEN
          R_WBC_in=R_WBC !radius of WBC in vessel
        ELSE IF(R_WBC.GT.R1) THEN
          R_WBC_in=R1
        ENDIF
        tp=(PI*R1**2.d0/FLOW)*length/2.d0*(1.d0+(R_WBC_in/R1)**2.42d0)
      ELSE
        R_min=MIN(R1,R2) !smallest radius of vessel
        IF(R_WBC.LE.R_min) THEN
          I=1.d0/0.58d0*R_WBC**2.42d0*(R2**0.58d0-R1**0.58d0)
        ELSE IF(R_WBC.GT.R_min.AND.R1.GT.R2) THEN
          I=1.d0/0.58d0*R_WBC**2.42d0*(R2**0.58d0-R1**0.58d0)+
     '      1.d0/3.d0*(R2**3.d0-R_WBC**3.d0)
        ELSE IF(R_WBC.GT.R_min.AND.R1.LE.R2) THEN
          I=1.d0/0.58d0*R1**2.42d0*(R2**0.58d0-R1**0.58d0)
        ENDIF
        tp=PI/(2.d0*FLOW)*length/(R2-R1)*(1.d0/3.d0*
     '    (R2**3.d0-R1**3.d0)+I) !tp=passage time (s)
      ENDIF
C... when WBC deforms to pass through vessel it retains smaller diameter
C... for rest of pathway through capillaries
      R_WBC=R_WBC_in
      t=te+tp !total transit time for WBC through segment (ne)

      CALL EXITS('WBC_TRANSIT')
      RETURN

 9999 CALL ERRORS('WBC_TRANSIT',ERROR)
      CALL EXITS('WBC_TRANSIT')
      RETURN 1
      END


