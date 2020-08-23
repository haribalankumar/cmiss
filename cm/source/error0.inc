!  error0.cmn
C     Suggested length of the ERROR string.
C     Note that the interpreter does not use the ERROR string supplied
C     to it  but instead uses a string of its own length.  It would be
C     sensible to make this value the same as the length of the
C     interpreter's error length.
      INTEGER ERRSTRLEN
      PARAMETER (ERRSTRLEN=1024)
