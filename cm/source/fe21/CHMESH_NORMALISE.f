      SUBROUTINE CHMESH_NORMALISE(XP,nk,np,ERROR,*)

C#### Subroutine: CHMESH_NORMALISE
C###  Description:
C###    CHMESH_NORMALISE is J.Crocombes normalise subroutine.  It
C###    normalises the derivative vector passed to it.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nk,np
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nj
      REAL*8 TEMP

      CALL ENTERS('CHMESH_NORMALISE',*9999)

      TEMP = 0.0d0
      DO nj = 1,NJT
        TEMP = TEMP+(XP(nk,1,nj,np)*XP(nk,1,nj,np))
      ENDDO !nj
      IF(TEMP.GT.0.1D-6)THEN
C GMH 8/1/97 Update cmgui link
        CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
        DO nj = 1,NJT
          XP(nk,1,nj,np) = XP(NK,1,NJ,NP)/DSQRT(TEMP)
        ENDDO !nj
      ENDIF

      CALL EXITS('CHMESH_NORMALISE')
      RETURN
 9999 CALL ERRORS('CHMESH_NORMALISE',ERROR)
      CALL EXITS('CHMESH_NORMALISE')
      RETURN 1
      END


