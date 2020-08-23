      SUBROUTINE CALCULATE_FSPACE(ELEM,IREGION,TETRA,FSPACE,NODE,ERROR,
     '  *)

C#### Subroutine: CALCULATE_FSPACE
C###  Description:
C###    Calulates the f-space. Genmesh routine.
CC JMB 20-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),IREGION(0:N_GM),TETRA(0:LDTETRA,4)
      REAL*8 FSPACE(0:N_GM),NODE(LDNODE,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,REF
      REAL*8 SPACE

      CALL ENTERS('CALCULATE_FSPACE',*9999)

      ! Initialize
      DO I=1,IREGION(0)
        FSPACE(I)=9.9d+20
      ENDDO !i
      REF=TETRA(0,HEAD)
      DO WHILE(REF.NE.0)
        DO I=1,NPTS
          IF(IREGION(ELEM(REF,I)).GT.1)THEN
            DO J=I+1,NPTS
              IF(IREGION(ELEM(REF,J)).GT.1)THEN
                SPACE=ZERO
                DO K=1,NJT
                  SPACE=SPACE+
     '              (NODE(ELEM(REF,I),K)-NODE(ELEM(REF,J),K))**2.d0
                ENDDO !k
                SPACE=SPACE*0.6123d0**2.d0
                IF(FSPACE(ELEM(REF,I)).GT.SPACE)THEN
                  FSPACE(ELEM(REF,I))=SPACE
                ENDIF !fspace
                IF(FSPACE(ELEM(REF,J)).GT.SPACE)THEN
                  FSPACE(ELEM(REF,J))=SPACE
                ENDIF !fspace
              ENDIF !iregion
            ENDDO !j
          ENDIF !iregion
        ENDDO !i
        ! Increment
        REF=TETRA(REF,NEXT)
      ENDDO !while ref

      CALL EXITS('CALCULATE_FSPACE')
      RETURN
 9999 CALL ERRORS('CALCULATE_FSPACE',ERROR)
      CALL EXITS('CALCULATE_FSPACE')
      RETURN 1
      END


