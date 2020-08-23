      SUBROUTINE SET_COLOUR_LUT(LUT,ERROR,*)

C#### Subroutine: SET_COLOUR_LUT
C###  Description:
C###    SET_COLOUR_LUT defines colour lookup table.

      IMPLICIT NONE
      INCLUDE 'colo00.cmn'
!     Parameter List
      REAL LUT(3,0:*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j
      REAL THETA

      CALL ENTERS('SET_COLOUR_LUT',*9999)

C **  Set up colours for colour indices 0(grey), 1(black) and 2(white)
      DO j=1,3
        LUT(j,0)=0.65 !grey
        LUT(j,1)=0.0  !black
        LUT(j,2)=1.0  !white
      ENDDO

C **  Set up colours for colour indices 3 (red) to MAXCOLOURS (blue)
      DO i=3,MAXCOLOURS
        THETA=REAL(i-3)*6.0/REAL(MAXCOLOURS-3)    !THETA range
        IF(THETA.LT.2.0) THEN                     !0 - 2
          LUT(1,i)=1.0
          LUT(3,i)=0.0
          IF(THETA.LT.1.0) THEN                   !  0 - 1
            LUT(2,i)=THETA*0.75
          ELSE                                    !  1 - 2
            LUT(2,i)=0.75+(THETA-1.0)/4.0
          ENDIF
        ELSE IF(THETA.LT.4.0) THEN                !2 - 4
          LUT(1,i)=(4.0-THETA)/2.0
          LUT(2,i)=1.0
          LUT(3,i)=(THETA-2.0)/2.0
        ELSE                                      !4 - 6
          LUT(1,i)=0.0
          LUT(3,i)=1.0
          IF(THETA.LT.5.0) THEN                   !  4 - 5
            LUT(2,i)=1.0-(THETA-4.0)/4.0
          ELSE                                    !  5 - 6
            LUT(2,i)=0.75-(THETA-5.0)*0.75
          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('SET_COLOUR_LUT')
      RETURN
 9999 CALL ERRORS('SET_COLOUR_LUT',ERROR)
      CALL EXITS('SET_COLOUR_LUT')
      RETURN 1
      END


