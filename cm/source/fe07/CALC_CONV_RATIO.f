      SUBROUTINE CALC_CONV_RATIO(ITER1,IOUNIT,NBH,ne,NEELEM,
     &  NONY,NPNY,nr,nr_solve,nx,NYNO,NYNR,
     &  CONY,ERRMAX,GRR,RATIO,YG,YP,CONTACT,
     &  OUTPUT,
     &  ERROR,*)

C 01/08/07 XSL
C#### Subroutine: CALC_CONV_RATIO
C###  Description:
C###   Calculates and outputs the convergence criteria for nonlinear problems

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'nonl00.cmn'

!     Parameter List
      INTEGER ITER1,IOUNIT,NBH(NHM,NCM,NEM),ne,
     &  NEELEM(0:NE_R_M,0:NRM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     &  NPNY(0:6,NYM,0:NRCM),nr,nr_solve,nx,
     &  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     &  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),ERRMAX,GRR(NOM),
     &  RATIO(0:6),YG(NIYGM,NGM,NEM),YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL CONTACT,OUTPUT

!     Local Variables
      INTEGER ng,nh,nk,no,noelem,no_nynr,no_resid,
     &  noy,ny1,ny2,nyo
      REAL*8 co,RSUM_CONSTRAINED(0:6),
     &  RSUM_UNCONSTRAIN(0:6),ZERO_TOL

      DATA ZERO_TOL /1.0d-8/

      CALL ENTERS('CALC_CONV_RATIO',*9999)

C news VJ. To sum all the residuals is dimensionally incorrect as each residual has different set of units
C So, now distinguishing between force residuals (index 1) 
C                                moment residuals (index 2)
C                                2nd moment residuals (index 3)        
C                                3rd moment residuals (index 4)
C                                hydrostatic pressure residuals (index 5)
C                                element based residuals (index 6)
C reinitialise residual sums and ratios and sum of soln increments
      RATIO(0)=0.0d0
C 21/07/08 XSL Added new option energy norm for convergence check
      RSUM_CONSTRAINED(0)=0.0d0
      RSUM_UNCONSTRAIN(0)=0.0d0
      DO no_resid=1,6
        RSUM_CONSTRAINED(no_resid)=0.0D0
        RSUM_UNCONSTRAIN(no_resid)=0.0D0
        RATIO(no_resid)=0.0d0
      ENDDO
C         Calc sum of constrained and unconstrained residuals
C NEWS 28/06/05 JHC only compute GRR in NONLIN for 1st Newton step, otherwise they
C are calculated through SOLVE5        
      IF(ITER1.EQ.0) THEN
        DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows
C         using GRR as temp storage for constraint reduced eqn resids
          GRR(no)=0.0d0
        ENDDO !no
      ENDIF
C NEWE JHC          
      DO no_nynr=1,NYNR(0,1,1,nr_solve) !loop over rows
        ny1=NYNR(no_nynr,1,1,nr_solve) !is row number         

C NEWS 28/06/05 JHC JHC only compute GRR in NONLIN for 1st Newton step, otherwise they
C are calculated through SOLVE5            
        IF(NONY(0,ny1,1,nr_solve).GT.0) THEN !free dependent variable
          IF(ITER1.EQ.0) THEN
C NEW 01/07/05 no need to compute ny2 as YP(ny1,4) has already been subtracted off from solve5
C                ny2=GETNYR(2,NPNY,nr_solve,0,1,ny1,NYNE,NYNP) !is RHS var #
            DO noy=1,NONY(0,ny1,1,nr_solve) !loop on rows assoc with ny1
              no=NONY(noy,ny1,1,nr_solve) !is row number for ny1
              co=CONY(noy,ny1,1,nr_solve) !is coupling coeff for ny1
C NEW 01/07/05 YP(ny1,4) has already been subtracted off from solve5
C                  GRR(no)=GRR(no)+(YP(ny1,4)-YP(ny2,1))*co
              GRR(no)=GRR(no)+YP(ny1,4)*co
            ENDDO !noy
          ENDIF
C NEWE JHC            
C 21/07/08 XSL Added new option energy norm for convergence check
        ELSE !bdry cond applied to dependent variable
          IF((KTYP007.EQ.1).OR.(KTYP007.EQ.2)) THEN !differentiated ratios and max gauss
C news VJ: differentiating between residuals of different orders of magnitude
C This is to non-dimensionalise each order before arriving at a ratio of unconstrained to constrained
C It is not dimensionally correct to sum all the residuals as they are of different units - force/moment/2nd moment etc
            IF(NPNY(0,ny1,0).EQ.1) THEN !node based
              nk=NPNY(1,ny1,0)
              nh=NPNY(3,ny1,0)
              IF(nh.LE.3) THEN!dealing with nodal attributes
                IF(nk.EQ.1) THEN!then dealing with nodal coordinates
                  RSUM_CONSTRAINED(1)=RSUM_CONSTRAINED(1)+
     &              DABS(YP(ny1,4))
                ELSE IF((nk.EQ.2).OR.(nk.EQ.3).OR.(nk.EQ.5)) THEN !dealing with 1st derivs
                  RSUM_CONSTRAINED(2)=RSUM_CONSTRAINED(2)+
     &              DABS(YP(ny1,4))
                ELSE IF((nk.EQ.4).OR.(nk.EQ.6).OR.(nk.EQ.7)) THEN ! 2nd derivs
                  RSUM_CONSTRAINED(3)=RSUM_CONSTRAINED(3)+
     &              DABS(YP(ny1,4))
                ELSE IF(nk.EQ.8) THEN ! 3rd derivs
                  RSUM_CONSTRAINED(4)=RSUM_CONSTRAINED(4)+
     &              DABS(YP(ny1,4))
                ENDIF ! nk, derivative type
              ELSE !dealing with hydrostatic pressures...
                RSUM_CONSTRAINED(5)=RSUM_CONSTRAINED(5)+
     &            DABS(YP(ny1,4))
              ENDIF ! nh, dependent variable type
            ELSE !dealing with element based params
              RSUM_CONSTRAINED(6)=RSUM_CONSTRAINED(6)+
     &          DABS(YP(ny1,4))
            ENDIF ! node/elem based parameter
C 21/07/08 XSL Added new option energy norm for convergence check
          ELSE IF(KTYP007.EQ.3) THEN !sum of unconstr R/sum of constr R
            RSUM_CONSTRAINED(0)=RSUM_CONSTRAINED(0)+DABS(YP(ny1,4))
          ENDIF !ktyp007
        ENDIF ! free or fixed variable
      ENDDO !no_nynr
      DO no=1,NOT(1,1,nr_solve,nx) !loop over global soln rows             
        DO nyo=1,NYNO(0,no,1,nr_solve)
          ny1=NYNO(nyo,no,1,nr_solve) !is row number
C 21/07/08 XSL Added new option energy norm for convergence check
          ny2=NYNO(nyo,no,2,nr_solve) !is global variable number
          IF ((KTYP007.EQ.1).OR.(KTYP007.EQ.2)) THEN !differentiated ratios and max gauss
C NEW JHC 04/07/05 not necessary to calculate ny2 as FIX(ny2,1) is no longer required
C              ny2=GETNYR(2,NPNY,nr_solve,0,1,ny1,NYNE,NYNP) !is RHS var #
C             Only add in residual if no force bc is applied
C NEWS 28/06/05 JHC force residual component should also be added to unconstrained residual
C             IF(.NOT.FIX(ny2,1)) THEN
            IF(NPNY(0,ny1,0).EQ.1) THEN !node based
              nk=NPNY(1,ny1,0)
              nh=NPNY(3,ny1,0)
              IF(nh.LE.3) THEN!dealing with nodal attributes
                IF(nk.EQ.1) THEN!then dealing with nodal coordinates
                  RSUM_UNCONSTRAIN(1)=RSUM_UNCONSTRAIN(1)+
     &              DABS(GRR(no))
                ELSE IF((nk.EQ.2).OR.(nk.EQ.3).OR.(nk.EQ.5)) THEN !dealing with 1st derivs
                  RSUM_UNCONSTRAIN(2)=RSUM_UNCONSTRAIN(2)+
     &              DABS(GRR(no))
                ELSE IF((nk.EQ.4).OR.(nk.EQ.6).OR.(nk.EQ.7)) THEN
                  RSUM_UNCONSTRAIN(3)=RSUM_UNCONSTRAIN(3)+
     &              DABS(GRR(no))
                ELSE IF(nk.EQ.8) THEN
                  RSUM_UNCONSTRAIN(4)=RSUM_UNCONSTRAIN(4)+
     &              DABS(GRR(no))
                ENDIF
              ELSE !dealing with hydrostatic pressures...
                RSUM_UNCONSTRAIN(5)=RSUM_UNCONSTRAIN(5)+
     &            DABS(GRR(no))
              ENDIF 
            ELSE !dealing with element based params
              RSUM_UNCONSTRAIN(6)=RSUM_UNCONSTRAIN(6)+
     &          DABS(GRR(no))
            ENDIF !NPNY
C             ENDIF ! .NOT.FIX
C 21/07/08 XSL NEWS Added new option energy norm for convergence check
          ELSE IF(KTYP007.EQ.3) THEN !sum of unconstr/sum of constr
            RSUM_UNCONSTRAIN(0)=RSUM_UNCONSTRAIN(0)+DABS(GRR(no))
          ELSE IF(KTYP007.EQ.4) THEN ! scalar product R.delta
            RATIO(0)=RATIO(0)+YP(ny1,4)*YP(ny2,5)
          ENDIF
C NEWE
C NEWE JHC
        ENDDO !nyo
      ENDDO !no

C 21/07/08 XSL NEWS Added new option energy norm for convergence check
      IF(KTYP007.EQ.4) THEN
        RATIO(0)=DABS(RATIO(0)) ! use the magnitude of the scalar product
        IF((.NOT.CONTACT.AND.ITER1.EQ.0).OR. ! 1st iteration, not contact
     &    (CONTACT.AND.CONV_IT.EQ.1)) THEN ! 1st iteration, contact
          RATIO(0)=2.0d0*ERRMAX ! for printing to terminal and loop can keep on going
        ELSE IF(INITIAL_ENERGY.LT.ZERO_TOL) THEN
C 04/03/09 XSL NEWS Check if INITIAL_ENERGY is zero
C if it is, then there is nothing to solve for, terminate the process
C if not add this additional check, R.deltau/INITIAL_ENERGY will be a non-zero number
C because the denominator is almost zero
C However, this has to be checked after first Newton step
C when YP5 has been calculated, so YP1 will be altered by a almost zero value
C which hopefully shouldnt matter
C Otherwise the convergence check for enery norm has to be done after SOLVE5,
C which will change how convergence check has been handled in NONLIN
          RATIO(0)=0.5d0*ERRMAX !set this to be half of the ERRMAX^2 so that it will terminate
        ELSE
          RATIO(0)=RATIO(0)/INITIAL_ENERGY ! normalise with respect to initial values
        ENDIF !check if it's the 1st iteration
      ENDIF !ktyp007
C NEWE

      IF(KTYP007.EQ.1) THEN !sum of ratios
        DO no_resid=1,6
          IF(RSUM_CONSTRAINED(no_resid).GT.ZERO_TOL) THEN !only include the ratio of a particular residual type if some d.o.f are constrained
            RATIO(no_resid)=
     &        RSUM_UNCONSTRAIN(no_resid)/RSUM_CONSTRAINED(no_resid)
          ELSE IF(RSUM_UNCONSTRAIN(no_resid).LT.ZERO_TOL) THEN !If both unconstrained and constrained are practically zero
            RATIO(no_resid)=0.0d0 ! needs to be smaller than zero_tol as they get added up
          ELSE !If unconstrained resid is high and constrained is zero
C NEWS 01/07/05 JHC for some cases, the unconstrained is slightly bigger than zero_tol, so it will have a very small number divided by zero.
C It is not necessary to make it very large (although it is mathematically more valid) as long as the convergence criteria is not met and 
C performs another Newton step. So make the ratio as small as possible but still prevents from converging
C                RATIO(no_resid)=1.0d0/ZERO_TOL
            RATIO(no_resid)=2.0d0*ERRMAX
C NEWE JHC
          ENDIF
          RATIO(0)=RATIO(0)+RATIO(no_resid)
        ENDDO

      ELSE IF(KTYP007.EQ.2) THEN ! Gauss point value check
C ***       DPN 09 May 2001 - allow a different measure of convergence by
C           comparing the unconstrained residuals to the some maximum value
C           from YG - i.e. compare the generated forces to the maximum active
C           tension rather than the generated, constrained boundary forces -
C           which are going to be zero for isotonic contraction...
C new       Apr15 VJ: ratio calc made dimensionally consistent by dividing sum of unconstr. force residual by max gauss point tension val 
        RSUM_CONSTRAINED(1)=0.0d0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ng=1,NGT(NBH(NH_LOC(1,nx),1,ne))
            RSUM_CONSTRAINED(1)=MAX(RSUM_CONSTRAINED(1),
     &        DABS(YG(KTYP007a,ng,ne)))
          ENDDO !ng
        ENDDO !noelem (ne)
  
        IF(RSUM_CONSTRAINED(1).GT.ZERO_TOL) THEN
          RATIO(0)=RSUM_UNCONSTRAIN(1)/RSUM_CONSTRAINED(1)
        ELSE IF(RSUM_UNCONSTRAIN(1).LT.ZERO_TOL) THEN !If both unconstrained and constrained are practically zero
          RATIO(0)=0.0d0
        ELSE !If unconstrained resid is high and constrained is zero
          RATIO(0)=1.0d0/ZERO_TOL
        ENDIF

C 21/07/08 XSL NEWS Added new option energy norm for convergence check
      ELSE IF(KTYP007.EQ.3) THEN !ratio of sums
        IF(RSUM_CONSTRAINED(0).GT.ZERO_TOL) THEN !only include the ratio of a particular residual type if some d.o.f are constrained
            RATIO(0)=RSUM_UNCONSTRAIN(0)/RSUM_CONSTRAINED(0)
        ELSE IF(RSUM_UNCONSTRAIN(0).LT.ZERO_TOL) THEN !If both unconstrained and constrained are practically zero
            RATIO(0)=0.0d0 ! needs to be smaller than zero_tol as they get added up
        ELSE !If unconstrained resid is high and constrained is zero
C NEWS 01/07/05 JHC for some cases, the unconstrained is slightly bigger than zero_tol, so it will have a very small number divided by zero.
C It is not necessary to make it very large (although it is mathematically more valid) as long as the convergence criteria is not met and 
C performs another Newton step. So make the ratio as small as possible but still prevents from converging
C                RATIO(no_resid)=1.0d0/ZERO_TOL
            RATIO(0)=2.0d0*ERRMAX
C NEWE JHC
        ENDIF !check RSUM_CONSTRAIN
C NEWE
      ENDIF !KTYP007
  
      IF(OUTPUT.OR.IWRIT4(nr_solve,nx).GT.1) THEN !output residuals
        IF((IWRIT4(nr_solve,nx).GT.1).AND.KTYP007.EQ.1) THEN !output the different components of residuals, ratios, and solution inc.
          WRITE(OP_STRING,
     &      '(/'' Sum of unconstrained value residuals    '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of constrained value residuals      '',
     &        ''           ='',D11.4,'
     &      //'/'' Their ratio (unconstr/constr)           '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of unconstrained 1st deriv residuals'',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of constrained 1st deriv residuals  '',
     &        ''           ='',D11.4,'
     &      //'/'' Their ratio (unconstr/constr)           '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of unconstrained 2nd deriv residuals'',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of constrained 2nd deriv residuals  '',
     &        ''           ='',D11.4,'
     &      //'/'' Their ratio (unconstr/constr)           '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of unconstrained 3rd deriv residuals'',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of constrained 3rd deriv residuals'',
     &        ''             ='',D11.4,'
     &      //'/'' Their ratio (unconstr/constr)           '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of unconstrained hydrostatic pressure '',
     &        ''residuals='',D11.4,'
     &      //'/'' Sum of constrained hydrostatic pressure '',
     &        ''residuals  ='',D11.4,'
     &      //'/'' Their ratio (unconstr/constr)           '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of unconstrained element based '',
     &        ''residuals       ='',D11.4,'
     &      //'/'' Sum of constrained element based '',
     &        ''residuals         ='',D11.4,'
     &      //'/'' Their ratio (unconstr/constr)           '',
     &        ''           ='',D11.4,'
     &      //'/'' Sum of ratios (unconstr/constr)         '',
     &        ''           ='',D11.4)')
     &      RSUM_UNCONSTRAIN(1),RSUM_CONSTRAINED(1),RATIO(1),
     &      RSUM_UNCONSTRAIN(2),RSUM_CONSTRAINED(2),RATIO(2),
     &      RSUM_UNCONSTRAIN(3),RSUM_CONSTRAINED(3),RATIO(3),
     &      RSUM_UNCONSTRAIN(4),RSUM_CONSTRAINED(4),RATIO(4),
     &      RSUM_UNCONSTRAIN(5),RSUM_CONSTRAINED(5),RATIO(5),
     &      RSUM_UNCONSTRAIN(6),RSUM_CONSTRAINED(6),RATIO(6),
     &      RATIO(0)
C newe VJ
        ELSE IF(OUTPUT.AND.(KTYP007.EQ.1)) THEN !only output sum of ratios and sum of soln inc.
          WRITE(OP_STRING,
     &      '('' Sum of (unconstr/constr) ratios''
     &        '' of the degrees of freedom='',D11.4)')
     &      RATIO(0)   
        ELSE IF(KTYP007.EQ.2) THEN ! Gauss point value check
          WRITE(OP_STRING,
     &      '(/'' Sum of unconstrained residual forces''
     &        '' ='',D11.4,'
     &      //'/'' Maximum Gauss point active tension ''
     &        ''  ='',D11.4,'
     &      //'/'' Their ratio''
     &        ''                          ='',D11.4)')
     &      RSUM_UNCONSTRAIN(1),RSUM_CONSTRAINED(1),
     &      RATIO(0)
C 21/07/08 XSL NEWS Added new option energy norm for convergence check
        ELSE IF(OUTPUT.AND.(KTYP007.EQ.3)) THEN !sum of unconstrained R/ sum of constrained R
           WRITE(OP_STRING,
     &      '('' Sum of unconstrained residuals over sum of ''
     &        ''constrained residuals='',D11.4)')
     &      RATIO(0)   
         ELSE IF(OUTPUT.AND.(KTYP007.EQ.4)) THEN !scalar product R.delta
          WRITE(OP_STRING,
     &      '('' Normalised scalar product of residuals ''
     &        ''and solution increments='',D11.4)')
     &      RATIO(0)   
        ENDIF
C NEWE
        CALL WRITES(IOUNIT,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('CALC_CONV_RATIO')
      RETURN
 9999 CALL ERRORS('CALC_CONV_RATIO',ERROR)
      CALL EXITS('CALC_CONV_RATIO')
      RETURN 1
      END
