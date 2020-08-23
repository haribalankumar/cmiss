      SUBROUTINE WBC_BLOCK(nb,NEELEM,NPNE,NYNP,CE,YP,ERROR,*)

C#### Subroutine: WBC_BLOCK
C###  Description:
C###    WBC_BLOCK calculates all possible capillary segments which
C###    could be blocked with white blood cells (WBCs) then uses
C###    a probability to allocate a certain amount as blocked.

C*** NB/ still under development - should remove blocked segments
C*** from solution rather than just increase resistance. (24/10/02)

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'pulm00.cmn'

!     Parameter list
      INTEGER nb,NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM)
      REAL*8 CE(NMM,NEM),YP(NYM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER ISEED,ne,noelem,np1,np2,TRAPS(0:NE_R_M)
      REAL*8 CM_RANDOM_NUMBER,DP,m,MAX_WBC,MIN_WBC,NUM_BLOCKED,p,P_cr,
     '  Ro_star,R_vessel,R_WBC,tau,te,te_star,u_cyto


      CALL ENTERS('WBC_BLOCK',*9999)

      NUM_BLOCKED=BLOCKED*NEELEM(0) !total # segments blocked
      ISEED=1
      TRAPS(0)=0
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
C... Determine which segments are potentially blocked with WBC's &
C... increase the resistance effectively to infinity (so zero flow)
        MAX_WBC=6.8d-3/2.d0 !mm
        MIN_WBC=8.3d-3/2.d0 !mm radius of white blood cell
        R_WBC=CM_RANDOM_NUMBER(ISEED)*(MAX_WBC-MIN_WBC)+MIN_WBC
        tau=0.035d0*1.d-3 !(N/m) avge tension in cell cortex
        R_vessel=DSQRT(CE(nm_a,ne)*CE(nm_b,ne)) !mm
        P_cr=(2.d0*tau/(R_WBC/1000.d0))*((R_WBC/R_vessel)-1.d0) !Pa
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        DP=DABS(YP(NYNP(1,1,1,np1,0,1))-
     '    YP(NYNP(1,1,1,np2,0,1))) !pressure difference
        IF(DP.GT.P_cr) THEN !calculate entrance time (te)
          m=6.d0
          u_cyto=135.d0 !cytoplasmic viscosity of WBC
          Ro_star=R_WBC/R_vessel
          IF(Ro_star.GE.1.d0.AND.Ro_star.LT.1.2d0) THEN
            te_star=m*(1.d0/3.d0+(2.d0/3.d0)*(Ro_star**3.d0)-
     '        Ro_star**2.d0+(2.d0/3.d0)*(Ro_star**2.d0-1.d0)**
     '        (3.d0/2.d0)+(Ro_star**2.d0-1.d0)**(1.d0/2.d0)-
     '        (Ro_star**2.d0*DSIN(ACOS(1.d0/Ro_star)))+
     '        DLOG((Ro_star*DSIN(ACOS(1/Ro_star))+Ro_star)/
     '        (Ro_star+(Ro_star**2.d0-1.d0)**(1.d0/2.d0)))) !eq (41)
          ELSE IF(Ro_star.LT.1.0d0) THEN
            te_star=0.d0 !if R_WBC < R_vessel then entrance time =0.d0
          ENDIF
          IF(Ro_star.LT.1.2d0) THEN !eqn (44), Huang (2001)
            te=u_cyto/(DP-P_cr)*te_star !entrance time (s)
          ELSE IF(Ro_star.GE.1.2d0) THEN
            te=0.107d0*m*u_cyto/(DP-33.d0)*10.d0**(3.68d0*
     '        (Ro_star-1.2d0))
          ENDIF
          IF(te.GE.1.d0) THEN
            TRAPS(0)=TRAPS(0)+1 !stores all the possible element traps
            TRAPS(TRAPS(0))=ne
          ENDIF
        ELSE
          TRAPS(0)=TRAPS(0)+1
          TRAPS(TRAPS(0))=ne
        ENDIF
      ENDDO !noelem
      NE_BLOCK(0)=0
      IF(TRAPS(0).GT.0) THEN
        p=NUM_BLOCKED/TRAPS(0) !probability of blockage
      ELSE
        p=0.d0
      ENDIF
      DO WHILE(NE_BLOCK(0).LT.NUM_BLOCKED)
        noelem=0
        DO WHILE(NE_BLOCK(0).LT.NUM_BLOCKED.AND.noelem.LT.TRAPS(0))
          noelem=noelem+1
          ne=TRAPS(noelem)
          IF(CM_RANDOM_NUMBER(ISEED).LE.p) THEN
            NE_BLOCK(0)=NE_BLOCK(0)+1 !allocates only certain percentage
            NE_BLOCK(NE_BLOCK(0))=ne !as blocked at this time
            CE(nm_Rseg,ne)=1.d10 !effectively infinite resistance
          ENDIF
        ENDDO !WHILE
      ENDDO !WHILE
C... NB/ should remove blocked segments rather than just increasing
C... resistance. Also the firat pass through here no solution has
C... yet been carried out therefore all DP = 0.d0


      CALL EXITS('WBC_BLOCK')
      RETURN

 9999 CALL ERRORS('WBC_BLOCK',ERROR)
      CALL EXITS('WBC_BLOCK')
      RETURN 1
      END


