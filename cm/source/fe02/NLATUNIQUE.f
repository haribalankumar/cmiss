      INTEGER FUNCTION NLATUNIQUE(L,start,step,end,NQNLAT)

C#### Function: NLATUNIQUE
C###  Type: INTEGER
C###  Description:
C###    NLATUNIQUE takes the start and end lattice positions
C###    of a lattice based grid in an xi direction and then
C###    counts how many grid points are below and above L,
C###    where L is the most extreme lattice point mapping
C###    to grid point nq. This allows the xi location of
C###    nq to be calculated.      
      
      IMPLICIT NONE

      INCLUDE 'grid00.cmn'       
      INCLUDE 'geom00.cmn'
!     Parameter List      
      INTEGER end,extreme,L,NQNLAT(NEQM*NQEM),start,step

!     Local Variables
      INTEGER i

      extreme=L
      
      DO i=start,end,step
        IF(NQNLAT(i).NE.NQNLAT(i+step)) THEN
          extreme=extreme+1
        ENDIF
      ENDDO

      NLATUNIQUE=extreme
      
      RETURN
      END


