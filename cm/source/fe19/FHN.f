      SUBROUTINE FHN(TIME,Y,DY,CONTROL,MODEL,SIZES,VARIANT,DERIVED,
     '  PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR)

C#### Subroutine: FHN
C###  Description:
C###    FHN returns the right hand side vector for the Fitzhugh-Nagumo
C###    cellular model.
C**** (re)Written by Martin Buist, 7 May 1999

C***  SIZES(1) is 3 equations (Vm,Recov,Cai).
C***  SIZES(2) is not used.
C***  SIZES(3) is not used.
C***  SIZES(4) is not used.
C***  SIZES(5) is 11 parameters.
C***  PARAM(1) is the membrane capacitance Cm (corrected by 10^-6).
C***  PARAM(2) is the surface to volume ration Am.
C***  PARAM(3) is resting potential in mV.
C***  PARAM(4) is plateau potential in mV.
C***  PARAM(5) is threshold potential in mV.
C***  PARAM(6) is excitation rate.
C***  PARAM(7) is excitation decay.
C***  PARAM(8) is recovery rate.
C***  PARAM(9) is recovery decay.
C***  PARAM(10) is Panfilov FHN mu2.
C***  PARAM(11) is Ca2+ time constant.
C***  SIZES(6) is 4 protocols.
C***  PROTOCOL(1) is pseudo stimulus current.
C***  PROTOCOL(4) is applied stimulus current.
C***  SIZES(7) is not used.
C***  SIZES(8) is not used.
C***  SIZES(9) is not used.
C***  SIZES(10) is not used.
C***  SIZES(11) is not used.

      IMPLICIT NONE

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cell_fhn.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR,SIZES(11)
      REAL*8 TIME(*),Y(*),DY(*),DERIVED(*),PARAM(*),PROTOCOL(*),ARI(*),
     '  ARO(*)
C      INTEGER SIZES(10)
C      INTEGER CONTROL(SIZES(NUM_CONTROL)),MODEL(SIZES(NUM_MODEL)),
C     '  VARIANT,AII(SIZES(NUM_AII)),AIO(SIZES(NUM_AIO)),ERR
C      REAL*8 T,Y(SIZES(NUM_EQN)),DY(SIZES(NUM_EQN)),
C     '  DERIVED(SIZES(NUM_DERIVED)),PARAM(SIZES(NUM_PARAM)),
C     '  PROTOCOL(SIZES(NUM_PROTOCOL)),ARI(SIZES(NUM_ARI)),
C     '  ARO(SIZES(NUM_ARO))
!     Local Variables
      REAL*8 DIFF,EPS,FHN1,NORM_PHI,PHI,RATE,THRESHOLD,VSTAR
      CHARACTER ERROR_DUMMY*20

      IF(DOP) THEN
        WRITE(OP_STRING,'('' REST='',D12.4)') PARAM(Vmrest)
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
        WRITE(OP_STRING,'('' PLATEAU='',D12.4)') PARAM(Vmplat)
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
        WRITE(OP_STRING,'('' THRESHOLD='',D12.4)') PARAM(Vmthres)
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
        WRITE(OP_STRING,'('' RATE='',D12.4)') PARAM(C1)
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
        WRITE(OP_STRING,'('' DECAY='',D12.4)') PARAM(C2)
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
      ENDIF

C     Normalise wrt resting and plateau potentials
      DIFF=PARAM(Vmplat)-PARAM(Vmrest)
      PHI=(Y(Vm)-PARAM(Vmrest))/DIFF
      THRESHOLD=(PARAM(Vmthres)-PARAM(Vmrest))/DIFF

      FHN1=PARAM(C1)*PHI*(PHI-THRESHOLD)*(PHI-1.d0)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' DIFF='',D12.4)') DIFF
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
        WRITE(OP_STRING,'('' THRESHOLD='',D12.4)') THRESHOLD
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
        WRITE(OP_STRING,'('' FHN='',D12.4)') FHN1
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
      ENDIF

      IF(KTYP33.EQ.1) THEN !Standard FHN
        FHN1=FHN1+PARAM(C2)*Y(Recov)
      ELSE !Roger's FHN or Panfilov FHN
        FHN1=FHN1+PARAM(C2)*PHI*Y(Recov)
      ENDIF

      DY(Vm)=FHN1
      IF(DOP) THEN
        WRITE(OP_STRING,'('' FHN_ION='',D12.4)') FHN1
        CALL WRITES(IODI,OP_STRING,ERROR_DUMMY,*9999)
      ENDIF

! recovery
C      PHI=Y(Vm)+DY(Vm)*DT
      PHI=Y(Vm)
      DIFF=PARAM(Vmplat)-PARAM(Vmrest)
      NORM_PHI=(PHI-PARAM(Vmrest))/DIFF
      IF(KTYP33.LT.3) THEN !Standard or Roger's FHN
        DY(Recov)=PARAM(B)*(NORM_PHI-PARAM(D)*Y(Recov))
      ELSE !Panfilov FHN
        EPS=PARAM(B)+PARAM(D)*Y(Recov)/(NORM_PHI+
     '    PARAM(Mu2))
        RATE=PARAM(C1)/DIFF/DIFF
        VSTAR=PHI+PARAM(Vmrest)-PARAM(Vmplat)-
     '    PARAM(Vmthres)
        DY(Recov)=EPS*(-Y(Recov)-RATE*VSTAR*(PHI-PARAM(Vmrest)))
      ENDIF

! calcium
      IF(DABS(PARAM(CaTime)).GT.ZERO_TOL) THEN
        IF(NORM_PHI.LT.ZERO_TOL) NORM_PHI=0.0d0
        DY(Cai)=(NORM_PHI-Y(Cai))/PARAM(CaTime)
      ELSE
        Y(Cai)=0.0d0
        DY(Cai)=0.0d0
      ENDIF

      DY(Vm)=-(DY(Vm)/PARAM(Cm))+((PROTOCOL(PseudoIs)+
     '  PROTOCOL(Is1current))/(PARAM(Cm)*PARAM(Am)))

 9999 RETURN
      END


