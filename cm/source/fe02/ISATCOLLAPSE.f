      LOGICAL FUNCTION ISATCOLLAPSE(IBT,INP,nb,nn)

C#### Function: ISATCOLLAPSE
C###  Type: LOGICAL
C###  Description:
C###    ISATCOLLAPSE returns .TRUE. if the local node nn is a collapsed
C###    node in a sector, .FALSE. if not.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      ! Parameter list
      INTEGER IBT(3,NIM),INP(NNM,NIM),nb,nn
      ! Local variables
      INTEGER ni
      LOGICAL COLLAPSED

      COLLAPSED=.FALSE.
      DO ni=1,NIT(nb)
        IF(IBT(1,ni).EQ.3.AND.IBT(2,ni).EQ.4) THEN
          IF(NKT(nn,nb).EQ.1) COLLAPSED=.TRUE.
        ELSE IF(IBT(1,ni).EQ.5) THEN
          IF(INP(nn,IBT(3,ni)).EQ.1) COLLAPSED=.TRUE.
        ELSE IF(IBT(1,ni).EQ.6) THEN
          IF(IBT(2,ni).EQ.4) THEN
            IF(INP(nn,IBT(3,ni)).EQ.2) COLLAPSED=.TRUE.
          ELSE
            IF(INP(nn,IBT(3,ni)).EQ.IBT(2,ni)+1) COLLAPSED=.TRUE.
          ENDIF
        ENDIF
      ENDDO

      ISATCOLLAPSE=COLLAPSED

      RETURN
      END


