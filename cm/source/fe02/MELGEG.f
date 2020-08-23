      SUBROUTINE MELGEG(LGE,nq,NHST,NHQ,nr,NYNQ,ERROR,*)

C#### Subroutine: MELGEG
C###  Description:
C###    MELGEG calculates the row numbers (LGE(*,1)) and column numbers
C###    (LGE(*,2)) in the matrix grid variables nhs.
C###    It also returns the total number of element variables NHST(nrc).
C###  See-Also: MELGE

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER LGE(NHM*NSM,NRCM),nq,NHST(NRCM),
     '  NHQ(NRM),nr,NYNQ(NHM,NQM,0:NRCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nh,nrc

      CALL ENTERS('MELGEG',*9999)

      DO nrc=1,2
        NHST(nrc)=0
        DO nh=1,NHQ(nr)
          NHST(nrc)=NHST(nrc)+1
          LGE(NHST(nrc),nrc)=NYNQ(nh,nq,nrc)
        ENDDO !nh
      ENDDO !nrc

      CALL EXITS('MELGEG')
      RETURN
 9999 CALL ERRORS('MELGEG',ERROR)
      CALL EXITS('MELGEG')
      RETURN 1
      END


