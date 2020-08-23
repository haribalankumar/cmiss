      SUBROUTINE ILISTRMDUP(N,IDATA,ERROR,*)

C#### Subroutine: ILISTRMDUP
C###  Description:
C###    ILISTRMDUP sorts N integer IDATA values into a non-decreasing
C###    sequence using IHEAPSORT (N>50) or ISHELLSORT (N>50) and then 
C###    removes all duplicates from the list. On exit N contains the 
C###    number of unique elements in the list.
C***  Rewritten, Martin Buist, 2004

      IMPLICIT NONE

!     Parameter List
      INTEGER IDATA(*),N
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER index,nolist

      CALL ENTERS('ILISTRMDUP',*9999)

      CALL ISORT(N,IDATA)

      index=0
      DO nolist=2,N
        IF(IDATA(nolist).EQ.IDATA(nolist-1)) THEN
          index=index+1
        ELSE
          IDATA(nolist-index)=IDATA(nolist)
        ENDIF
      ENDDO !nolist

      N=N-index

      CALL EXITS('ILISTRMDUP')
      RETURN
 9999 CALL ERRORS('ILISTRMDUP',ERROR)
      CALL EXITS('ILISTRMDUP')
      RETURN 1
      END

