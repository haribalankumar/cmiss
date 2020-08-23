      SUBROUTINE ONE_2_TWO(nq_1d,nq,nx,YQ_1D,YQS_1D,YQ,YQS,ERROR,*)

C#### Subroutine: ONE_2_TWO
C###  Description:
C###    Used to copy 1-D strips into 2-D sheets

      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nq_1d,nq,nx
      REAL*8 YQ_1D(NYQM,NIQM),YQS_1D(NIQSM,NQM),
     '  YQ(NYQM,NIQM,NAM,NXM),YQS(NIQSM,NQM)
      CHARACTER ERROR*(*)
      !Local variables
      INTEGER niqs

      CALL ENTERS('ONE_2_TWO',*9999)
      DO niqs=1,NIQM
        YQ(nq,niqs,1,nx)=YQ_1D(nq_1d,niqs)
      ENDDO
      DO niqs=1,NIQSM
        YQS(niqs,nq)=YQS_1D(niqs,nq_1d)
      ENDDO
      CALL EXITS('ONE_2_TWO')
      RETURN
 9999 CALL ERRORS('ONE_2_TWO',ERROR)
      CALL EXITS('ONE_2_TWO')
      RETURN 1
      END


