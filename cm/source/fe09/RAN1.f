      REAL*8 FUNCTION RAN1(idum)
C#### Function: RAN1
C###  Type: REAL*8
C###  Description:
C###    "Minimal" Random number generator of Park and Miller with
C###    Bays-Durham shuffle and added saeguards. Returns a uniform
C###    random deviate between 0.0 and 1.0 (exclusive of the endpoint
C###    values). RNMX should approximate the largest floating value
C###    that is less than 1.
C###
C###    Note : Call with idum a negative integer to initialise;
C###    thereafter, do not alter idum between successive deviates in a
C###    sequence.
C***  Referenced from Numerical Recepies in Fortran, Second Edition

      IMPLICIT NONE
!     Paramter
      INTEGER idum
!     Local Variables
      REAL*8 AM,EPS,RNMX
      INTEGER IA,IM,IQ,IR,NTAB,NDIV
      PARAMETER(IA=16807,IM=2147483647,IQ=127773,
     '  IR=2836,NTAB=32,EPS=1.2D-08,
     '  RNMX=1.D0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/

      AM=1.D0/DBLE(IM)
      NDIV=1+(IM-1)/NTAB
      IF(idum.LE.0.OR.iy.EQ.0) THEN    !initialise
        idum=max(-idum,1)              !prevent idum=0
        DO j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          IF(idum.LT.0) idum=idum+IM
          IF(j.LE.NTAB) iv(j)=idum
        ENDDO
        iy=iv(1)
      ENDIF
      k=idum/IQ                        !start here when not init.
      idum=IA*(idum-k*IQ)-IR*k
      IF(idum.LT.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)             !don't want endpoints
      ran1=min(AM*iy,RNMX)
      RETURN
      END


