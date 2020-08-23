      SUBROUTINE DATXY(AXIS,DEPTH,WIDTH,Z,ZD,ERROR,*)

C#### Subroutine: DATXY
C###  Description:
C###    DATXY calculates the width and depth of a data set.
C    Created: JPC 20/5/97

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      REAL*8 DEPTH,WIDTH,Z,ZD(NJM,NDM)
      INTEGER AXIS
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER d,D1nd,D2nd,nd,w,W1nd,W2nd
      REAL*8 DEPTH1,DEPTH2,WIDTH1,WIDTH2,ZDAT

      CALL ENTERS('DATXY',*9999)
        WIDTH1=1000.0d0
        WIDTH2=1000.0d0
        DEPTH1=1000.0d0
        DEPTH2=1000.0d0
        IF(AXIS.EQ.1) THEN
          w=2
          d=3
        ELSEIF(AXIS.EQ.2) THEN
          w=1
          d=3
        ELSE
          w=1
          d=2
        ENDIF
        DO nd=1,NDT
          IF(NJT.EQ.2) THEN
            IF(AXIS.EQ.3) THEN
              ZDAT=0.0d0
            ELSE
              ZDAT=ZD(AXIS,nd)
            ENDIF
          ELSEIF(NJT.EQ.3) THEN
            ZDAT=ZD(AXIS,nd)
          ENDIF
          IF(NJT.EQ.3.OR.AXIS.EQ.3) THEN
            IF(ZD(w,nd).GE.zero_tol) THEN
              IF(((ZDAT-Z)**2+(ZD(d,nd)-0.0d0)**2).LT.DEPTH1) THEN
                DEPTH1=(ZDAT-Z)**2+(ZD(d,nd)-0.0d0)**2
                D1nd=nd
              ENDIF
            ELSE
              IF(((ZDAT-Z)**2+(ZD(d,nd)-0.0d0)**2).LT.DEPTH2) THEN
                DEPTH2=(ZDAT-Z)**2+(ZD(d,nd)-0.0d0)**2
                D2nd=nd
              ENDIF
            ENDIF
            IF(ZD(d,nd).GE.zero_tol) THEN
              IF(((ZDAT-Z)**2+(ZD(w,nd)-0.0d0)**2).LT.WIDTH1) THEN
                WIDTH1=(ZDAT-Z)**2+(ZD(w,nd)-0.0d0)**2
                W1nd=nd
              ENDIF
            ELSE
              IF(((ZDAT-Z)**2+(ZD(w,nd)-0.0d0)**2).LT.WIDTH2) THEN
                WIDTH2=(ZDAT-Z)**2+(ZD(w,nd)-0.0d0)**2
                W2nd=nd
              ENDIF
            ENDIF
          ELSE
            IF(ZD(w,nd).GE.zero_tol) THEN
              IF((ZDAT-Z)**2.LT.WIDTH1) THEN
                WIDTH1=(ZDAT-Z)**2
                W1nd=nd
              ENDIF
            ELSE
              IF((ZDAT-Z)**2.LT.WIDTH2) THEN
                WIDTH2=(ZDAT-Z)**2
                W2nd=nd
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        IF(NJT.EQ.3.OR.AXIS.EQ.3) DEPTH=(ZD(d,W1nd)-ZD(d,W2nd))
        WIDTH=(ZD(w,D1nd)-ZD(w,D2nd))

      CALL EXITS('DATXY')
      RETURN
 9999 CALL ERRORS('DATXY',ERROR)
      CALL EXITS('DATXY')
      RETURN 1
      END


