      SUBROUTINE EVALDW(En,EG,n,nr,ne,ng,NQNE,dWdE,ERROR,*)
      IMPLICIT NONE

C#### Subroutine: EVALDW
C###  Description:
C###    EVALDW evaluates the derivative of the strain energy function
C###      for a single component of the Green strain tensor while holding all
C###      other components constant. This is used for integrating the derivative
C###      of the strain energy function to obtain the strain energy.

      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='EVALDW')

      INCLUDE 'error0.inc'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp50.cmn'

C     Parameter List
      INTEGER n,nr,ne,ng,NQNE(NEQM,NQEM)
      REAL*8 En, !strain value of component n of Green's strain
     & EG(3,3),DW(6),dWdE
      CHARACTER ERROR*(*)

C     Local Variables
      INTEGER IDXI(6),IDXJ(6),i,j

C     Map n (index of DW) to i and j (indices of EG).
C     This avoids a chain of IF ELSE.
      DATA IDXI /1,2,3,1,1,2/
      DATA IDXJ /1,2,3,2,3,3/

      CALL ENTERS(ROUTINENAME,*9999)

      dWdE=0.0d0

      EG(IDXI(n),IDXJ(n))=En
      EG(IDXJ(n),IDXI(n))=En

!       PRINT *,'EVALDW EG=',
!      &((EG(i,j),i=1,3),j=1,3)

      IF((KTYP54(nr).EQ.3).AND.(KTYP3B.EQ.2)) THEN
      ! calculate DW from EG using cellml material law
        CALL EGDW_CELLML(EG,%VAL(ICQS_SPATIAL_PTR),
     &      %VAL(IRCQS_SPATIAL_PTR),ne,ng,
     &      %VAL(NQNE_PTR),%VAL(RCQS_SPATIAL_PTR),
     &      %VAL(YQS_PTR),DW,ERROR,*9999)

      ELSEIF((KTYP54(nr).EQ.1).AND.(KTYP55(nr).EQ.3)) THEN
        CALL ENERGY(nr,%VAL(CG_PTR),DW,EG(1,1),EG(2,2),EG(3,3),
     &    EG(1,2),EG(1,3),EG(2,3),%VAL(YG_PTR),ERROR,*9999)
      ELSE
        CALL ASSERT(.FALSE.,
     &    '>>Strain energy can only be calculated if the material'
     &    //'law is in terms of the fibre & transverse strains.',
     &    ERROR,*9999)
C       GR If you want to implement other cases then you will need to
C       calculate the invariants from the strain tensors. See ZGTG53
C       for examples.
      ENDIF
! 
!       PRINT *,'EVALDW DW=',
!      &(DW(i),i=1,6)

      dWdE=DW(n)

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END
