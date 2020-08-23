      SUBROUTINE PREODES(N,ERROR,*)

C#### Subroutine: PREODES. Sets parameter values before calling RADAU5 subroutine.


      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'odes00.cmn'
      
!     Parameter List
      INTEGER N
      CHARACTER ERROR*(*)
!----------------------------------------------------------------
!                      Local Variables
!---------------------------------------------------------------- 
      INTEGER K        
!     LOGICAL 
!     CHARACTER 
!-----------------------------------------------------------------
                        !ENTER PREODES
!----------------------------------------------------------------
      CALL ENTERS('PREODES',*9999)             
      
      !SWITCH VARIABLES
      IMAS=0
      ITOL=0
      IJAC=0
      MLJAC=N
      MLMAS=N
      MUMAS=1
      IOUT=0
      IDID=0
      H=0.1D-5
      !TOLERANCES AND STEP SIZE
      RTOL=0.1D-10!0.1D-10
      ATOL=0.1D-10!0.1D-10
      LWORK=1000
      LIWORK=1000   
      DO K=1,LIWORK
        IWORK(K)=0
        WORK(K)=0.0D1
      ENDDO   
      !WORK(1)=1.0D-18
!-----------------------------------------------------------------  

      CALL EXITS('PREODES')
      RETURN
 9999 CALL ERRORS('PREODES',ERROR)
      CALL EXITS('PREODES')
      RETURN 1
      END


