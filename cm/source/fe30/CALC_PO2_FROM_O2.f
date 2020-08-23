      SUBROUTINE CALC_PO2_FROM_O2(content,Hb,pH,PCO2,PO2,SHbO2,Temp,
     &  ERROR,*)

C#### Subroutine: CALC_PO2_FROM_O2
C###  Description:
C###     Calculates the partial pressure from the oxygen content.
C###     Variable stepping iterative procedure is used for inversion. 
C###  Inputs:
C###     content -> O2 content ml O2 / ml blood
C###     Hb -> Hemoglobin content molar (standard value = 2.33E-3 = 0.15g/ml)
C###     pH -> pH of plasma
C###     PCO2 -> CO2 partial pressure mmHg (standard value = 40)
C###     Temp -> temperature degrees celcius
C###  Outputs:
C###     PO2 -> O2 partial pressure mmHg
C###     SO2 -> fractional oxyhemoglobin saturation
C**** Created by AJS, Nov 2009

     
      IMPLICIT NONE

!     Parameter List
      REAL*8 content,Hb,pH,PCO2,PO2,SHbO2,Temp
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,max_iterations
      REAL*8 diff_new,diff_old,inc,O2_old,O2_new,O2_next,O2step,PO2_old,
     &  PO2_new,PO2_next,SHbO2_old,SHbO2_new,SHbO2_next,tolerance
      LOGICAL CONVERGED

      CALL ENTERS('CALC_PO2_FROM_O2',*9999)
!       write(*,*) 'For content=',content

      CALL ASSERT(content.GE.0.0d0,'>>Error: O2 content<0',ERROR,*9999)

      IF(content.EQ.0.0d0)THEN
        PO2=0.0d0
        SHbO2=0.0d0
      ELSE

C     Parameters
!       Temp=37.0d0  !degrees C
!       pH=7.40d0    !pH of plasma (pHrbc=7.24)
!       Hb=2.33d-3   !hemoglobin content in blood, molar (=15 g/100ml)
!       PCO2=40.0d0  !CO2 partial pressure mmHg    
      max_iterations=100 !terminate when reached
      CONVERGED=.FALSE.
      i=1
      tolerance=1.0d-3 ! convergence criteria
      inc=10.0d0  ! increment for PO2, mmHg 

C     Choose start point
      PO2_new=50.0d0            ! mmHg
      CALL CALC_O2_KELMAN(O2_new,Hb,pH,PCO2,PO2_new,SHbO2_new,Temp,
     &  ERROR,*9999)

C     Check convergence
      IF(ABS((O2_new-content)/content).LT.(tolerance*content))
     &  CONVERGED=.TRUE.

C     Loop to find PO2 value
      DO WHILE (.NOT.CONVERGED.AND.(i.LT.max_iterations))
    
C       Modify increment size
        IF(O2_new.GT.content)THEN
          inc=-1.0d0*ABS(inc)
        ELSEIF(O2_new.LT.content)THEN
          inc=ABS(inc)
        ENDIF
        IF(i.GT.1)THEN 
          diff_new=O2_new-content
          diff_old=O2_old-content
          O2step=ABS(O2_new-O2_old)
          IF((diff_old.GT.0.0d0.AND.diff_new.LT.0.0d0).OR.
     &      (diff_old.LT.0.0d0.AND.diff_new.GT.0.0d0))THEN !last 2 steps straddle point
            inc=inc/2.0d0
          ELSEIF(ABS(diff_new).GT.O2step)THEN
            inc=inc*2.0d0
          ENDIF
        ENDIF

C       Increment to find new PO2, SO2, O2
        PO2_next=PO2_new+inc
        CALL CALC_O2_KELMAN(O2_next,Hb,pH,PCO2,PO2_next,SHbO2_next,Temp,
     &    ERROR,*9999)
!       write(*,*) 'i=',i,' inc=',inc,' PO2=',PO2_next,' SO2=',
!      &  SHbO2_next,' O2=',O2_next

C       Check convergence
        IF(ABS((O2_next-content)/content).LT.(tolerance*content))
     &    CONVERGED=.TRUE.
        i=i+1
        PO2_old=PO2_new
        O2_old=O2_new
        SHbO2_old=SHbO2_new
        PO2_new=PO2_next
        O2_new=O2_next
        SHbO2_new=SHbO2_next

      ENDDO !while
      
      CALL ASSERT(CONVERGED,'>>Error: PO2 value not found',ERROR,*9999)
      PO2=PO2_new
      SHbO2=SHbO2_new

      ENDIF !content=0

      CALL EXITS('CALC_PO2_FROM_O2')
      RETURN
 9999 CALL ERRORS('CALC_PO2_FROM_O2',ERROR)
      CALL EXITS('CALC_PO2_FROM_O2')
      RETURN 1
      END

