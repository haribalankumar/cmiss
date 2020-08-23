      SUBROUTINE AM_INTERPOLATE(K,MAX_ORDER,MAX_RND,
     '  NUMBER_EQN,NUM_VARS,DYOUT,PHI,PSI,T,TOUT,Y,YOUT)

C#### Subroutine: AM_INTERPOLATE
C###  Description:
C###    Interpolates the solution at TOUT using the Kth order Adams
C###    polynomial calculated near T in AM_STEP. On output YOUT
C###    contains the solution at TOUT and DYOUT contains the
C###    derivative of the solution at TOUT.

      IMPLICIT NONE
!     Parameter list
      INTEGER K,MAX_ORDER,MAX_RND,NUMBER_EQN,NUM_VARS
      REAL*8 PHI(NUMBER_EQN,MAX_ORDER+2+MAX_RND),PSI(MAX_ORDER),T,
     '  TOUT,Y(NUM_VARS),YOUT(NUM_VARS),DYOUT(NUM_VARS)
!     Local Variables
      INTEGER i,j,l
      REAL*8 ETA,G(13),GAMMA,H_I,RHO(13),TERM,W(13)

C KAT 2003-05-30: Not thread-safe as it places the array in global storage.
C       DATA G(1) /1.0d0/
C       DATA RHO(1) /1.0d0/
      G(1)=1d0
      RHO(1)=1d0

C     The formula for interpolation is
C       y_out=y_(n+1)+h_I\sum{i=1}{K+1}[g_i,1^I.Phi_i(n+1)]
C       dy_out=\sum{i=1}{K+1}[rho_i^I.Phi_i(n+1)]
C     where
C       h_I=t_out-t_(n+1)
C       g_i,q^I=1/q                                       i=1
C       g_i,q^I=GAMMA_(I-1)(1).g_(i-1),q^I-
C         eta_(i-1).g_(i-1),(q+1)^I                       i>=2
C       eta_i=h_I/Psi_i(n+1)
C       GAMMA_i(s)=s.h_I/Psi_1(n+1)                       i=1
C       GAMMA_i(s)=(s.h_I+Psi_(i-1)(n+1))/Psi_i(n+1)      i>=2
C       rho_1=1.0                                         i=1
C       rho_i=rho_(i-1).GAMMA_(i-1)(1)                    i=2,..,K+1

      H_I=TOUT-T

C     Initialise W for computing G
      DO i=1,K+1
        W(i)=1.0d0/DBLE(i)
      ENDDO !i
      TERM=0.0d0

C     Compute G
      DO j=2,K+1
        GAMMA=(H_I+TERM)/PSI(j-1)
        ETA=H_I/PSI(j-1)
        DO i=1,K+2-j
          W(i)=GAMMA*W(i)-ETA*W(i+1)
        ENDDO !i
        G(j)=W(1)
        RHO(j)=GAMMA*RHO(j-1)
        TERM=PSI(j-1)
      ENDDO !j

C     Interpolate
      DO l=1,NUMBER_EQN
        YOUT(l)=0.0d0
        DYOUT(l)=0.0d0
      ENDDO !l
      DO j=1,K+1
        i=K+2-j
        DO l=1,NUMBER_EQN
          YOUT(l)=YOUT(l)+G(i)*PHI(l,i)
          DYOUT(l)=DYOUT(l)+RHO(i)*PHI(l,i)
        ENDDO !l
      ENDDO !j
      DO l=1,NUMBER_EQN
        YOUT(l)=Y(l)+H_I*YOUT(l)
      ENDDO !l

      RETURN
      END


