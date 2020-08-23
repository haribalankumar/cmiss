      REAL*8 FUNCTION PSI20(IDO,INP,nb,nu,nk,nn,XI)

C#### Function: PSI20
C###  Type: REAL*8
C###  Description:
C###    PSI20 evaluates tensor product Lagrange and Hermite basis
C###    functions at XI.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nk,nn,nu
      REAL*8 XI(*)
!     Local Variables
      INTEGER IPD(9,3),ni
      REAL*8 PHT3
      DATA IPD/4,3,2,1,3,2,1,1,1, 1,2,3,4,1,1,1,3,2, 1,1,1,1,2,3,4,2,3/

C     CALL ENTERS('PSI20',*9999)
      PSI20=1.0d0
      DO ni=1,NIT(nb)
        PSI20=PSI20*PHT3(INP(nn,ni),IDO(nk,nn,ni),IPD(nu,ni),XI(ni))
      ENDDO

C     CALL EXITS('PSI20')
      RETURN
      END

C FE50 Functions
C ==============

