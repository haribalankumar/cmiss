! cell_hh.inc
!  Array elements for the Hodgkin-Huxley ionic current model.

C *** Y,DY arrays
      INTEGER h,m,n
c      PARAMETER(Vm=1) !y,dy given in cell_reserved.inc
      PARAMETER(h=2) !y,dy
      PARAMETER(m=3) !y,dy
      PARAMETER(n=4) !y,dy

C *** DERIVED
      INTEGER DIStim,DINa,DIK,DIl
      PARAMETER(DIStim=1) !derived
      PARAMETER(DINa=2) !derived
      PARAMETER(DIK=3) !derived
      PARAMETER(DIl=4) !derived

C *** Parameter array
      INTEGER Er,gNa_max,gK_max,gl_max
      !PARAMETER(Cm=1) - from cell_reserved.inc
      !PARAMETER(Am=2) - from cell_reserved.inc
      PARAMETER(Er=3) ! Resting potential (mV)
      PARAMETER(gNa_max=4) ! peak Na conductance (mS.mm^-2)
      PARAMETER(gK_max=5) ! peak K conductance (mS.mm^-2)
      PARAMETER(gl_max=6) ! peak leak conductance (mS.mm^-2)

C *** Control variables
      INTEGER RETURN_CURRENTS
      !PARAMETER(ODE=1) - from cell_reserved.inc
      PARAMETER(RETURN_CURRENTS=2)
