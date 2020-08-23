      REAL*8 FUNCTION RAN0(idum)
C#### Function: RAN0
C###  Type: REAL*8
C###  Description:
C###    "Minimal" Random number generator of Park and Miller.
C###    Returns a uniform random deviate between 0.0 and 1.0. Set or
C###    reset idum to and integer value (except the unlikely value
C###    of MASK) to initialise the sequence.
C###
C###    Note : The seed - idum must not be altered
C###    between calls for successive divates in a sequence.
C***    Referenced from Numerical Recepies in Fortran, Second Edition

      IMPLICIT NONE
!     Parameter List
      INTEGER idum
!     Local Variables
      INTEGER iA,IM,IQ,IR,k,MASK
      REAL*8 AM

      PARAMETER(IA=16807,IM=2147483647,AM=1.D0/2147483647.0d0,
     '  IQ=127773, IR=2836,MASK=123459876)

      idum=ieor(idum,MASK) !XOR
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      IF(idum.LT.0.D0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      RETURN
      END

