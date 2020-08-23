      SUBROUTINE CHECK_INTEGRITY(ELEM,FACE,TETRA,CIRC,NODE,RSQ,ERROR,*)

C#### Subroutine: CHECK_INTEGRITY
C###  Description:
C###    Checks the Delaunay mesh integrity.
CC JMB 27-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),FACE(LDFACE,NFACES),TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ(N_GM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,ISUM,J,K,NEIB,REF,VERT
      REAL*8 SUM
!     Local Functions
      REAL*8 DLAMCH

      CALL ENTERS('CHECK_INTEGRITY',*9999)

      REF=TETRA(0,HEAD)
      DO WHILE(REF.NE.0)
        DO I=1,NFACES
          NEIB=FACE(REF,I)
          IF(NEIB.NE.0) THEN
            ISUM=0
            DO J=1,NPTS
              IF(I.NE.J)THEN
                ISUM=ISUM+ELEM(REF,J)
              ENDIF !i
            ENDDO !j
            J=1
            DO WHILE((FACE(NEIB,J).NE.REF).AND.(J.LT.NFACES))
              J=J+1
            ENDDO !while
            DO K=1,NPTS
              IF(J.NE.K)THEN
                ISUM=ISUM-ELEM(NEIB,K)
              ENDIF !j
            ENDDO !k
            VERT=ELEM(NEIB,J)
            CALL ASSERT_VORO(ISUM.GT.0,
     '        'Face point convention irregularity',ERROR,*9999)
            SUM=ZERO
            DO J=1,NJT
              SUM=SUM+(CIRC(REF,J)-NODE(VERT,J))**2.d0
            ENDDO
            SUM=RSQ(REF)-SUM
            CALL ASSERT_VORO(SUM.LT.DLAMCH('S'), 'Non DELAUNAY face',
     '        ERROR,*9999)
          ENDIF !neib
        ENDDO !i
        ! Increment
        REF=TETRA(REF,NEXT)
      ENDDO !while ref

      CALL EXITS('CHECK_INTEGRITY')
      RETURN
 9999 CALL ERRORS('CHECK_INTEGRITY',ERROR)
      CALL EXITS('CHECK_INTEGRITY')
      RETURN 1
      END


