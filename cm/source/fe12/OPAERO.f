      SUBROUTINE OPAERO(NPL,ERROR,*)

C#### Subroutine: OPAERO
C###  Description:
C###    OPAERO outputs aerofoil parameters.
C**** Xn(nn,nj)      are nodal coordinates
C**** dXndXi(nn,nj)  are nodal coordinate derivatives
C**** Zn(nn)         are nodal values of veloc potential Phi
C**** dZndXi(nn)     are nodal Phi derivatives
C**** X(nj)          are coordinates at Xi point
C**** dXdXi(nj,ni)   are derivs of Xj wrt Xi(ni)
C**** dSdXi(ni)      are arclength derivs wrt Xi(ni)
C**** Z              is Phi at Xi(ni)
C**** dZdXi(ni)      is deriv of Phi wrt Xi(ni) at Xi(ni)
C**** NORMAL(nj)     is unit normal vector to aerofoil

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NPL(5,0:3,NLM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,ne1,ne2,ng,nl,no_aero,no_entry,no_exit,no_wake,
     '  NP1,NP2

      CALL ENTERS('OPAERO',*9999)
      CALL ASSERT(CALL_AERO,
     '  '>>Aerofoil parameters not defined',ERROR,*9999)

! Aerofoil calculations
      WRITE(OP_STRING,'(/'' Aerofoil surface lines: '')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      DO no_aero=1,NL_AERO(0,1)+NL_AERO(0,2)
        IF(no_aero.LE.NL_AERO(0,1)) THEN !upper surface
          nl=NL_AERO(no_aero,1)
          ne=NE_AERO(no_aero,1)
        ELSE                             !lower surface
          nl=NL_AERO(no_aero-NL_AERO(0,1),2)
          ne=NE_AERO(no_aero-NL_AERO(0,1),2)
        ENDIF

        NP1=NPL(2,1,nl) !is node # at Xi=0 end of line
        NP2=NPL(3,1,nl) !is node # at Xi=1 end of line

        IF(no_aero.LE.NL_AERO(0,1)) THEN !upper surface
          WRITE(OP_STRING,'('' Line '',I5,'' Element '',I5,'
     '      //'''   Upper surface Nodes: '',2I6)') nl,ne,NP1,NP2
        ELSE                             !lower surface
          WRITE(OP_STRING,'('' Line '',I5,'' Element '',I5,'
     '      //'''   Lower surface Nodes: '',2I6)') nl,ne,NP1,NP2
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ENDDO

      WRITE(OP_STRING,'(/'' Leading  edge node = '',I5)') NP_aero_LE
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'( '' Trailing edge node = '',I5)') NP_aero_TE1
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

! Upstream inflow face nodes
      WRITE(OP_STRING,'(/'' Upstream entry face nodes: '',20I6)')
     '  (NP_ENTRY(no_entry),no_entry=1,NP_ENTRY(0))
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

! Downstream outflow face parameters
      WRITE(OP_STRING,'(/'' Downstream exit face lines: '')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      DO no_exit=1,NL_EXIT(0,1)
        nl=NL_EXIT(no_exit,1)

        NP1=NPL(2,1,nl) !is node # at Xi=0 end of line
        NP2=NPL(3,1,nl) !is node # at Xi=1 end of line

        WRITE(OP_STRING,'('' Line '',I5,'' Downstream exit Nodes: '','
     '    //'2I6)') nl,NP1,NP2
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO !no_exit

! Wake calculations
      WRITE(OP_STRING,'(/'' Wake surface lines: '')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      DO no_wake=1,NL_WAKE(0,1)
        nl=NL_WAKE(no_wake,1)
        ne1=NE_WAKE(no_wake,1)
        ne2=NE_WAKE(no_wake,2)
        NP1=NPL(2,1,nl) !is node # at Xi=0 end of line
        NP2=NPL(3,1,nl) !is node # at Xi=1 end of line
        WRITE(OP_STRING,'('' Line '',I5,'' Elements: '',2I6,'
     '    //'''   Wake nodes: '',2I6)') nl,ne1,ne2,NP1,NP2
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO !no_wake

      WRITE(OP_STRING,'(/'' Nodal d(PHI)_residuals '
     '  //'& Pressure diffs (at Xi1=1) across wake:'')')
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO no_wake=1,NL_WAKE(0,1)
        NP1=NPL(3,1,NL_WAKE(no_wake,1))
        NP2=NPL(3,1,NL_WAKE(no_wake,2))
        WRITE(OP_STRING,'('' Nodes '',2I5,'' dPHI_residual ='',E12.5,'
     '    //''' Press_diff. ='',E12.5)')
     '    NP1,NP2,dPHI_resid(no_wake),PRESS_DIFF_WAKE(no_wake)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO !no_wake

      WRITE(OP_STRING,'(''  Trailing edge pressure diff. = '',E12.5)')
     '  PRESS_DIFF_WAKE(NL_WAKE(0,1)+1)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      WRITE(OP_STRING,'(/'' Reference velocity = '',E12.5)')
     '  REF_VELOC
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      Fluid density = '',E12.5)')
     '  FLUID_DENSITY
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      Circumference = '',E12.5)') AERO_PERIM
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''        Circulation = '',E12.5)')
     '  CIRCULATION
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      Aerofoil lift = '',E12.5)') TOT_LIFT
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''      Aerofoil drag = '',E12.5)') TOT_DRAG
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'(''   Lift coefficient = '',E12.5)') LIFT_COEFF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(DABS(TOT_LIFT).GT.1.D-6) THEN
        WRITE(OP_STRING,'(''    Ratio drag/lift = '',E12.5)')
     '    TOT_DRAG/TOT_LIFT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      WRITE(OP_STRING,'(''  TE velocity diff. = '',E12.5)')
     '  TE_VELOC_DIFF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(CALL_COUP) THEN !coupling problem setup
        WRITE(OP_STRING,'(/'' Pressure diff.s across aerofoil:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO no_aero=1,NL_AERO(0,1)
          WRITE(OP_STRING,'('' Aerofoil line '',I2,'
     '      //''' Gauss point pressure diff.s = '',10E12.5)') no_aero,
     '      (PRESS_DIFF_AERO(ng,no_aero),ng=1,NGT(NB_AERO_PRESS))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF !call_coup

      CALL EXITS('OPAERO')
      RETURN
 9999 CALL ERRORS('OPAERO',ERROR)
      CALL EXITS('OPAERO')
      RETURN 1
      END


