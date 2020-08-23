      REAL*8 FUNCTION PF1(na,nu,XI)

C#### Function: PF1
C###  Type: REAL*8
C###  Description:
C###    PF1 evaluates Fourier basis function number na at Xi.  nu
C###    specifies time derivative. OMEGA is angular frequency.
C###    If nu=-1 PF1 returns integral terms (from xi=0 so integral of
C###    sin term has integration constant added).

      IMPLICIT NONE
      INCLUDE 'four00.cmn'
!     Parameter List
      INTEGER na,nu
      REAL*8 XI

      IF(na.EQ.1) THEN              !constant term
        IF(nu.EQ.1) THEN            !0th-order time derivative
          PF1=1.0D0
        ELSE IF(nu.EQ.2) THEN     !1st-order time derivative
          PF1=0.0D0
        ELSE IF(nu.EQ.3) THEN     !2nd-order time derivative
          PF1=0.0D0
        ELSE IF(nu.EQ.-1) THEN     !time integral
          PF1=XI
        ENDIF

      ELSE IF(MOD(na,2).EQ.0) THEN  !cosine terms
        IF(nu.EQ.1) THEN            !0th-order time derivative
          PF1=DCOS((DBLE(na)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.2) THEN     !1st-order time derivative
          PF1=-(DBLE(na)/2.0D0)*OMEGA*DSIN((DBLE(na)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.3) THEN     !2nd-order time derivative
          PF1=-(DBLE(na)/2.0D0)*OMEGA*(DBLE(na)/2.0D0)*OMEGA*
     '      DCOS((DBLE(na)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.-1) THEN     !time integral
          PF1=DSIN((DBLE(na)/2.0D0)*OMEGA*XI)/(DBLE(na)/2.0D0*OMEGA)
        ENDIF

      ELSE                         !sine terms
        IF(nu.EQ.1) THEN            !0th-order time derivative
          PF1=DSIN((DBLE(na-1)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.2) THEN     !1st-order time derivative
          PF1=(DBLE(na-1)/2.0D0)*OMEGA*DCOS((DBLE(na-1)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.3) THEN     !2nd-order time derivative
          PF1=-(DBLE(na-1)/2.0D0)*OMEGA*(DBLE(na-1)/2.0D0)*OMEGA*
     '      DSIN((DBLE(na-1)/2.0D0)*OMEGA*XI)
        ELSE IF(nu.EQ.-1) THEN     !time integral
          PF1=(1.0D0-DCOS(DBLE(na-1)/2.0D0*OMEGA*XI))/
     '      (DBLE(na-1)/2.0D0*OMEGA)
        ENDIF
      ENDIF

      RETURN
      END


