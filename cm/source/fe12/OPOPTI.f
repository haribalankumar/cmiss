      SUBROUTINE OPOPTI(LDR,NLNO,NMNO,NP1OPT,NP2OPT,NP3OPT,
     '  NYNO,PAOPTY,CM,CONTR,PAOPTI,PMIN,PMAX,RESID,RESIDM,
     '  SUMMARY,ERROR,*)

C#### Subroutine: OPOPTI
C###  Description:
C###    OPOPTI outputs optimisation parameters.

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'opti00.cmn'
      INCLUDE 'trsf00.cmn'
!     Parameter List
      INTEGER LDR(0:NDM),NLNO(NOPM),NMNO(1:2,0:NOPM),NP1OPT(NOPM),
     '  NP2OPT(NOPM),NP3OPT(NOPM),NYNO(0:NYOM,NOOPM,NRCM),PAOPTY(NOPM)
      REAL*8 CM(*),CONTR(*),PAOPTI(*),PMIN(*),PMAX(*),
     '  RESID(*),RESIDM(*)
      CHARACTER ERROR*(*)
      LOGICAL SUMMARY
!     Local Variables
      INTEGER nd,nl,nocont,nojoin,noopti,nores,NP1,NP2,NPR,ny
      REAL*8 DATARESID,SUMSQUARE,RMSERR
      CHARACTER TITLE1(2)*20,TITLE2(9)*41,TITLE3(13)*41,TITLE4(0:4)*13,
     '  TITLE5(2)*5,TITLE6(3)*52

      DATA TITLE1/'material parameters',
     '            'geometric parameters'/

      DATA TITLE2/'square of max princ. stress diff.',!ktyp27=1,ktyp26=1
     '            'sum of squared reaction diff.s',   !ktyp27=2    "
     '            'zero flux differences',            !ktyp27=3    "
     '            'hydrostatic pressure condition',   !ktyp27=4    "
     '            'geometric least squares',          !ktyp27=5    "
     '            ' ',                                !ktyp27=6    "
     '            ' ',                                !ktyp27=7    "
     '            ' ',                                !ktyp27=8    "
     '            ' '/                                !ktyp27=9    "

      DATA TITLE3/'squared area of trapezoids',       !ktyp27=1,ktyp26=2
     '            'curvature of stripes',             !ktyp27=2    "
     '            'zero flux differences',            !ktyp27=3    "
     '            'hydrostatic pressure condition',   !ktyp27=4    "
     '            'data fitting',                     !ktyp27=5    "
     '            'fluid interface condition',        !ktyp27=6    "
     '            'aerofoil wake pressure difference',!ktyp27=7    "
     '            'aerofoil lift and wake',           !ktyp27=8    "
     '            'boundary layer thickness condition',!ktyp27=9   "
     '            'customiasation measurments',       !ktyp27=10   "
     '            'difference in second moments',     !ktyp27=11   "
     '            'activation times',                 !ktyp27=12   "
     '            'dipole source'/                    !ktyp27=13   "
      DATA TITLE4/'none         ', !CONSTRAINT_TYPE=0
     '            'circumference', !CONSTRAINT_TYPE=1
     '            '             ', !CONSTRAINT_TYPE=2
     '            '             ', !CONSTRAINT_TYPE=3
     '            '             '/ !CONSTRAINT_TYPE=4

      DATA TITLE5/'NPSOL', !KTYP29=1
     '            'MINOS'/ !KTYP29=2
      DATA TITLE6/'s are components of data projections',      !KTYP29B=1
     '  ' is the Euclidean norm of data projection magnitudes',!KTYP29B=2
     '  's are squares of data projections magnitudes'/        !KTYP29B=3

      CALL ENTERS('OPOPTI',*9999)

      WRITE(OP_STRING,'(''  Optimisation package is '',A)')
     ' TITLE5(KTYP29)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      IF(KTYP29.EQ.2) THEN
        IF(SPARSEJAC) THEN
          WRITE(OP_STRING,'(''  Problem is sparse'')')
        ELSE
          WRITE(OP_STRING,'(''  Problem is not sparse'')')
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(KTYP26.EQ.1) THEN   !optimising material parameters
        WRITE(OP_STRING,'(''  Optimise '',A,'' by minimizing '',A)')
     '    TITLE1(KTYP26),TITLE2(KTYP27)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(KTYP27.EQ.1) THEN      !Square of max princ stress diff

        ELSE IF(KTYP27.EQ.2) THEN !Squared reaction diffs
          WRITE(OP_STRING,'(''  Material params in optimisation '','
     '      //''': '',12I3)') (NMNO(1,noopti),noopti=1,NTOPTI)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            WRITE(OP_STRING,'(''  Parameter #'',I2,'': '',3E13.5)')
     '        NMNO(1,noopti),PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'(''  Number of sets of apriori '','
     '      //'''measurements = '',I2)') KTYP28
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(KTYP27.EQ.3) THEN

        ELSE IF(KTYP27.EQ.4) THEN !Hydrostatic press condition
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            NP1=NP1OPT(noopti)
            WRITE(OP_STRING,
     '        '(''  Parameter #'',I2,'' (node '',I3,''): '',3E13.5)')
     '        noopti,NP1,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(KTYP27.EQ.5) THEN !Geometric least squares
        ENDIF

      ELSE IF(KTYP26.EQ.2) THEN !optimising geometric parameters
        WRITE(OP_STRING,'(''  Optimise '',A,'' by minimizing '',A)')
     '    TITLE1(KTYP26),TITLE3(KTYP27)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(KTYP27.EQ.1) THEN !Squared area of trapezoids
          WRITE(OP_STRING,'(''  Material params in optimisation '
     '      //'are: '',12I3)') (NMNO(1,noopti),noopti=1,NTOPTI)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            WRITE(OP_STRING,'(''  Parameter #'',I5,'': '',3E13.5)')
     '        NMNO(1,noopti),PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(KTYP27.EQ.2) THEN !Curvature of stripes
          WRITE(OP_STRING,'(''  Material params in optimisation '','
     '      //'''are: '',12I3)') (NMNO(1,noopti),noopti=1,NTOPTI)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            WRITE(OP_STRING,'(''  Parameter #'',I2,'': '',3E13.5)')
     '        NMNO(1,noopti),PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'(''  Number of sets of apriori '','
     '      //'''measurements = '',I2)') KTYP28
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(KTYP27.EQ.3) THEN !Zero flux diffs
          WRITE(OP_STRING,'('' Node num.s in optimisation are: '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            NP1=NPJOIN(noopti,1)
            NP2=NPRAD(-1,noopti,1)
            WRITE(OP_STRING,'('' Interface node '',I3)')NP1
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Corresponding fixed surface node '','
     '        //'I3)') NP2
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Nodes on same radial line  '')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)
     '        (NPRAD(NPR,noopti,1),NPR=1,NPRAD(0,noopti,1))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(KTYP27.EQ.4) THEN !Hydrostatic pressure condition
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            NP1=NP1OPT(noopti)
            WRITE(OP_STRING,
     '        '(''  Parameter #'',I2,'' (node '',I3,''): '',3E13.5)')
     '        noopti,NP1,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(KTYP27.EQ.5) THEN !Data fitting
          WRITE(OP_STRING,'(''  Residual'',A)') TITLE6(KTYP29B)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(.NOT.SUMMARY) THEN
            WRITE(OP_STRING,
     '        '(''  Current, minimum & maximum parameter values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noopti=1,NTOPTI
              IF(PAOPTY(noopti).EQ.1) THEN !Parameter is geometric dof
                ny=NYNO(1,noopti,2)
                WRITE(OP_STRING,
     '            '(''  Parameter #'',I5,'' (ny = '',I4,''): '','
     '      //'3D13.5)')
     '            noopti,ny,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE IF(PAOPTY(noopti).EQ.2) THEN !Parameter is a line
                nl=NLNO(noopti)
          WRITE(OP_STRING,
     '      '(''  Parameter #'',I5,'' (nl = '',I4,''): '','
     '      //'3D13.5)')
     '            noopti,nl,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            WRITE(OP_STRING,
     '        '(''  Constraint values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nocont=1,NTCNTR
              IF(KTYP29.EQ.1) THEN
                WRITE(OP_STRING,
     '      '(''  Constraint number '',I5,'' : '',E13.5)') nocont,
     '      CONTR(nocont)
              ELSE IF(KTYP29.EQ.2) THEN
                WRITE(OP_STRING,
     '      '(''  Constraint number '',I5,'' : '',E13.5)') nocont,
     '      CM(nocont)
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
          IF(KTYP1B.EQ.1) THEN !unit length components
            WRITE(OP_STRING,
     '        '(''  Unit length derivatives are '
     '        //'optimised as components'')')
          ELSE IF(KTYP1B.EQ.2) THEN !derivative as angle
            WRITE(OP_STRING,
     '        '(''  Derivatives are optimised as angles'')')
          ELSE IF(KTYP1B.EQ.3) THEN !non-unit length
            WRITE(OP_STRING,
     '        '(''  Non-unit length derivatives are '
     '        //'optimised as components'')')
          ELSE !invalid
            WRITE(OP_STRING,
     '        '(''  Invalid KTYP1B'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DATARESID=0.0d0
          SUMSQUARE=0.0d0
          IF(KTYP29.EQ.1) THEN
            DO nd=1,NT_RES-KTYP12
              DATARESID=DATARESID+RESID(nd)**2
              SUMSQUARE=SUMSQUARE+RESID(nd)
            ENDDO
          ELSE IF(KTYP29.EQ.2) THEN
            DO nd=1,NT_RES-KTYP12
              DATARESID=DATARESID+RESIDM(nd)**2
              SUMSQUARE=SUMSQUARE+RESIDM(nd)
            ENDDO
          ENDIF
          IF(KTYP29B.NE.3.AND.KTYP29.EQ.1) SUMSQUARE=DATARESID
          RMSERR=DSQRT(SUMSQUARE/DBLE(LDR(0)))
          WRITE(OP_STRING,'(''  Data residual       = '',D12.5)')
     '      DATARESID
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  RMS error           = '',D12.5)')
     '      RMSERR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP12.EQ.1) THEN !Sobolev smoothing
            IF(KTYP29.EQ.1) THEN
              WRITE(OP_STRING,'(''  Sobolev value       = '',D12.5)')
     '    RESID(NT_RES)**2
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''  Total residual      = '',D12.5)')
     '          DATARESID+RESID(NT_RES)**2
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
               IF(RESID(NT_RES).ne.0.0d0) THEN
          WRITE(OP_STRING,'(''  Smoothing ratio     = '',D12.5)')
     '      (RESID(NT_RES)**2)/DATARESID
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''  Log Smoothing ratio = '',D12.5)')
     '            DLOG10((RESID(NT_RES)**2)/DATARESID)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
           ELSE IF(KTYP29.EQ.2) THEN
              WRITE(OP_STRING,'(''  Sobolev value       = '',D12.5)')
     '    RESIDM(NT_RES)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''  Total residual      = '',D12.5)')
     '          DATARESID+RESIDM(NT_RES)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(RESIDM(NT_RES).ne.0.0d0) THEN
          WRITE(OP_STRING,'(''  Smoothing ratio     = '',D12.5)')
     '      RESIDM(NT_RES)/DATARESID
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''  Log Smoothing ratio = '',D12.5)')
     '            DLOG10((RESIDM(NT_RES)**2)/DATARESID)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
          ELSE
            WRITE(OP_STRING,'(''  Total residual      = '',D12.5)')
     '        DATARESID
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          WRITE(OP_STRING,'(''  Number of optimisation variables = '','
     '      //'I4)') NTOPTI
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Number of constraints            = '','
     '      //'I4)') NTCNTR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Number of residuals              = '','
     '      //'I4)') NT_RES
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(KTYP27.EQ.6) THEN !Fluid interface condition
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nojoin=1,NPJOIN(0,1)
            NP1=NPJOIN(nojoin,1)
            noopti=2*nojoin-1
            WRITE(OP_STRING,
     '        '(''  Parameter #'',I2,'' (node '',I3,''): '',3E13.5)')
     '        noopti,NP1,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Parameter #'',I2,11X,'': '',3E13.5)')
     '        noopti+1,PAOPTI(noopti+1),PMIN(noopti+1),
     '        PMAX(noopti+1)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(KTYP27.EQ.7) THEN !Aerofoil wake pressure difference
!         & (if sail stress eqtns defined) sail stress residuals
          WRITE(OP_STRING,'(''  #wake parameters = '',I2)') N_OPTI(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  #sail parameters = '',I2)') N_OPTI(2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(N_OPTI(1).GT.0) THEN !sail optimisation defined
           WRITE(OP_STRING,'('' Wake parameters (nodes then deriv):'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noopti=1,N_OPTI(1)
              WRITE(OP_STRING,'(''  Parameter #'',I2,'': '','
     '          //'''(Nodes '',2I5,'')'',3E13.5)')
     '          noopti,NP1OPT(noopti),NP2OPT(noopti),
     '          PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

          IF(N_OPTI(2).GT.0) THEN !sail optimisation defined
            WRITE(OP_STRING,'('' Sail stress parameters:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noopti=N_OPTI(1)+1,NTOPTI
              WRITE(OP_STRING,'(''  Parameter #'',I2,'': '','
     '          //'''(Node '',I5,'' var '',I1,'' deriv '',I1,'')'','
     '          //'3E13.5)')
     '          noopti,NP1OPT(noopti),NP2OPT(noopti),
     '          NP3OPT(noopti),PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

        ELSE IF(KTYP27.EQ.8) THEN !Aerofoil lift and wake
          WRITE(OP_STRING,
     '      '(''  Initial, minimum & maximum parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(N_OPTI(1).GE.NL_WAKE(0,1)) THEN
            WRITE(OP_STRING,'('' Wake params (nodes then deriv):'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noopti=1,NL_WAKE(0,1)
              WRITE(OP_STRING,'(''  Parameter #'',I2,'': '','
     '          //'''(Nodes '',2I5,'')'',3E13.5)')
     '          noopti,NP1OPT(noopti),NP2OPT(noopti),
     '          PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF

          WRITE(OP_STRING,'('' Aero parameters:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=N_OPTI(1)+1,NTOPTI
            WRITE(OP_STRING,'(''  Parameter #'',I2,'': '','
     '        //'''(Nodes '',2I5,'')'',3E13.5)')
     '        noopti,NP1OPT(noopti),NP2OPT(noopti),
     '        PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
          WRITE(OP_STRING,'('' Constraint is: '',A,'', value='',E12.5)')
     '      TITLE4(CONSTRAINT_TYPE),PMIN(NTOPTI+1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING,'('' Lift upper bound = '',E12.5)')
     '      UPPER_BOUND
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(KTYP27.EQ.10) THEN !torso measurements
          WRITE(OP_STRING,'('' Number of mesh parameters = '',I2)')
     '      NTOPTI
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '(''  Current, minimum & maximum mesh parameter values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noopti=1,NTOPTI
            WRITE(OP_STRING,'(''  Mesh parameter #'',I3,'
     '        //''' : '',3D13.5)')
     '        noopti,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !noopti
          IF(CUSTOMISATION_TYPE.EQ.1) THEN
          WRITE(OP_STRING,'('' Using circumferential measurements '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(CUSTOMISATION_TYPE.EQ.2) THEN
          WRITE(OP_STRING,'('' Using volume differences '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(CUSTOMISATION_TYPE.EQ.3) THEN
          WRITE(OP_STRING,'('' Using ellipse approx. to circumference '
     '       //''')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(CUSTOMISATION_TYPE.EQ.4) THEN
          WRITE(OP_STRING,'('' Using 2D measurements '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' Number of measurements = '',I2)') NT_RES
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(CUSTOMISATION_TYPE.EQ.1.OR.CUSTOMISATION_TYPE.EQ.3) THEN
          WRITE(OP_STRING,
     '      '(''  Circumferential measurement & measurement height:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nores=1,NT_RES
            WRITE(OP_STRING,'(''  Measurement #'',I3,'' : '',2D13.5)')
     '        nores,CIRMEASURE(1,nores),CIRMEASURE(2,nores)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nores
          ELSEIF(CUSTOMISATION_TYPE.EQ.4)THEN
          WRITE(OP_STRING,
     '      '(''  2-d in x and y with height of measurements:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
           DO nores=1,NT_RES/2
            WRITE(OP_STRING,'(''  Measurement #'',I3,'' : '',3D13.5)')
     '        nores,WDMEASURE(1,nores),WDMEASURE(2,nores),
     '       WDMEASURE(3,nores)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
           ENDDO
          ENDIF
          WRITE(OP_STRING,
     '      '(''  Residuals:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nores=1,NT_RES
            WRITE(OP_STRING,'(''  Residual #'',I5,'' : '',1D13.5)')
     '        nores,RESID(nores)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nores
        ELSE IF(KTYP27.EQ.11) THEN !second moments
          WRITE(OP_STRING,'(''  Current theta value:'',D12.4)')
     '      PAOPTI(1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP27.EQ.12) THEN !Activation inverse
          IF(ACTN_OPTI_WAVE_PARAM) THEN
            WRITE(OP_STRING,'(''  Wavefront parameters are '
     '        //'optimised'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Current transmembrane jump '
     '        //'potential: '',D12.3)') TRSF_ACTN_WAVE_JUMP
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Current wavefront width: '',D12.3)')
     '        TRSF_ACTN_WAVE_WIDTH
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'(''  Wavefront parameters are not '
     '        //'optimised'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'(''  Sum-of-squares objective weight:'''
     '      //',D12.4)') SS_OBJ_WEIGHT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Correlation coefficient objective '
     '      //'weight:'',D12.4)') CC_OBJ_WEIGHT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ACTN_IREGULARISE.EQ.1) THEN ! Additional constraints
            WRITE(OP_STRING,'(''  Additional contraint is : none'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(ACTN_IREGULARISE.EQ.2) THEN
            WRITE(OP_STRING,'(''  Additional constraint is :'','
     '        //''' a surface Laplacian'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            ERROR='>>Unknown additional constraints'
            GOTO 9999
          ENDIF

          IF(ACTN_IREGULARISE.EQ.2) THEN
            WRITE(OP_STRING,'(''  The regularisation parameter : '',
     '        F11.4)') ACTN_REG_PARAM_LAPLACE
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          IF(ACTN_IRESFUN.EQ.1) THEN
            WRITE(OP_STRING,
     '        '(''  Using inverse routines : normal'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSEIF(ACTN_IRESFUN.EQ.2) THEN
            WRITE(OP_STRING,
     '        '(''  Using inverse routines : optimised'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            ERROR='>>Code needs updating - Unknown code type'
            GOTO 9999
          ENDIF

        ENDIF

        IF(KTYP27.NE.5) THEN !not already output
          IF(.NOT.SUMMARY) THEN
            WRITE(OP_STRING,
     '        '(''  Current, minimum & maximum parameter values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noopti=1,NTOPTI
              IF(PAOPTY(noopti).EQ.1) THEN !Parameter is geometric dof
                ny=NYNO(1,noopti,2)
                WRITE(OP_STRING,
     '            '(''  Parameter #'',I5,'' (ny = '',I4,''): '','
     '      //'3D13.5)')
     '            noopti,ny,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE IF(KTYP26.EQ.2.AND.KTYP27.EQ.12) THEN
C               Wave parameters for activation inverse optimisation
                WRITE(OP_STRING,
     '            '(''  Parameter #'',I5,''            : '','
     '      //'3D13.5)')
     '            noopti,PAOPTI(noopti),PMIN(noopti),PMAX(noopti)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
            WRITE(OP_STRING,
     '        '(''  Constraint values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nocont=1,NTCNTR
              WRITE(OP_STRING,
     '          '(''  Constraint number '',I3,'' : '',E13.5)') nocont,
     '          CONTR(nocont)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            DATARESID=0.0d0
            SUMSQUARE=0.0d0
            DO nd=1,NT_RES
              WRITE(OP_STRING,'(''  Residual '',I5,'' = '',D12.5)')
     '          nd,RESID(nd)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DATARESID=DATARESID+RESID(nd)**2
              SUMSQUARE=SUMSQUARE+RESID(nd)
            ENDDO
          ELSE ! SUMMARY - only calculate residual
            DATARESID=0.0d0
            SUMSQUARE=0.0d0
            DO nd=1,NT_RES
              DATARESID=DATARESID+RESID(nd)**2
              SUMSQUARE=SUMSQUARE+RESID(nd)
            ENDDO
          ENDIF
          RMSERR=DSQRT(SUMSQUARE/DBLE(NT_RES))
          WRITE(OP_STRING,'(''  Total residual       = '',D12.5)')
     '      DATARESID
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  RMS error            = '',D12.5)')
     '      RMSERR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING,'(''  Number of optimisation variables = '','
     '      //'I4)') NTOPTI
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Number of constraints            = '','
     '      //'I4)') NTCNTR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Number of residuals              = '','
     '      //'I4)') NT_RES
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !ktyp27.ne.5
      ENDIF

      IF(   (KTYP26.EQ.1.AND.KTYP27.EQ.2)       !mat pars sqd reac diffs
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.5)       !Data Fitting
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.7)       !aerofoil wake & stress
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.8)       !aerofoil lift
     '  .OR.(KTYP26.EQ.2.AND.KTYP27.EQ.12)) THEN ! activation timing

C***    optimisations that call E04UPF

        WRITE(OP_STRING,'(''  Print level = '',I5)') IPPLEV
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP29.EQ.2) THEN !MINOS
          IF(SUMMFILE) THEN
            WRITE(OP_STRING,'(''  Problem has a summary file'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        WRITE(OP_STRING,'(''  Optimality tolerance = '',E10.3)') OPTTOL
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  Linesearch tolerance = '',E10.3)') LNSTOL
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  Step limit = '',E10.3)') STEPLIM
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(''  Derivative level = '',I2)') IPDLEV
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)


C LKC 23-MAR-1999 Adding maximum major iteration output
        IF(MAX_MAJOR_ITER.LT.0) THEN
          WRITE(OP_STRING,
     '      '(''  Using default number of maximum major iterations '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,
     '      '(''  Maximum number of major iterations is '',I3)')
     '      MAX_MAJOR_ITER
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF


C LKC 21-NOV-1999 Adding maximum minor iteration output
        IF(MAX_MINOR_ITER.LT.0) THEN
          WRITE(OP_STRING,
     '      '(''  Using default number of maximum minor iterations '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,
     '      '(''  Maximum number of minor iterations is '',I3)')
     '      MAX_MINOR_ITER
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF


        IF(IPDLEV.EQ.0) THEN
           IF(DIFFER.GT.0.D0) THEN
            WRITE(OP_STRING,'(''  Finite difference interval = '','
     '        //'E10.3)') DIFFER
          ELSE
            WRITE(OP_STRING,'(''  Finite difference interval to be '
     '        //'computed'')')
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'(''  Verify Level = '',I2)') IPVLEV
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(IPVLEV.EQ.1.OR.IPVLEV.EQ.3) THEN
            WRITE(OP_STRING,'(''  Starting objective check var. = '','
     '        //'I3)') ISROCV
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Stopping objective check var. = '','
     '        //'I3)') ISPOCV
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(NTCNTR.GT.0) THEN
          WRITE(OP_STRING,'(''  Nonlinear feasibility tolerance = '','
     '      //'E10.3)') NLFTOL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(IPVLEV.EQ.2.OR.IPVLEV.EQ.3) THEN
            WRITE(OP_STRING,'(''  Starting constraint check var. = '','
     '        //'I3)') ISRCCV
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  Stopping constraint check var. = '','
     '        //'I3)') ISPCCV
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        IF(KTYP29.EQ.2) THEN !MINOS
          WRITE(OP_STRING,'(''  Debug level = '',I3)') DBGLEV
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Iterations limit = '',I5)') ITERLM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('OPOPTI')
      RETURN
 9999 CALL ERRORS('OPOPTI',ERROR)
      CALL EXITS('OPOPTI')
      RETURN 1
      END


