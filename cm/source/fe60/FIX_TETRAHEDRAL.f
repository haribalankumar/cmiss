      SUBROUTINE FIX_TETRAHEDRAL(ELEM,FACE,TETRA,CIRC,NODE,RSQ,
     '  PLANE,ERROR,*)

C#### Subroutine: FIX_TETRAHEDRAL
C###  Description:
C###    Genmesh routine.
CC JMB 18-11-01

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),FACE(LDFACE,NFACES),TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ(N_GM)
      LOGICAL PLANE
      CHARACTER ERROR*(*)
!     Local Varaiables
      INTEGER COUNT,I,AXIS(2),ITETRA(3),J,K,L,M,NEIB,INODE(4)
      LOGICAL BOOL
!     Local Functions
      REAL*8 DLAMCH

      CALL ENTERS('FIX_TETRAHEDRAL',*9999)

      COUNT=0
      I=TETRA(0,HEAD)
      DO WHILE(I.NE.0)
        IF(DABS(RSQ(I)).LT.SQRT(DLAMCH('E')))THEN
          COUNT=COUNT+1 ! Increment counter
          ! Get the 4 nodes at the tip of the neighbouring tetrahedral
          DO J=1,NPTS
            NEIB=FACE(I,J)
            IF(NEIB.EQ.0)THEN
              INODE(J)=0
            ELSE
              DO K=1,NPTS
                BOOL=.TRUE.
                DO L=1,NPTS
                  BOOL=BOOL.AND.(ELEM(NEIB,K).NE.ELEM(I,L))
                ENDDO !l
                IF(BOOL)THEN
                  INODE(J)=ELEM(NEIB,K)
                ENDIF !bool
              ENDDO !k
            ENDIF !neib
          ENDDO !j
          BOOL=.FALSE.
          DO J=1,NPTS
            BOOL=BOOL.OR.(INODE(J).EQ.0)
          ENDDO !j
          IF(BOOL)THEN
            CALL PEEL_TETRAHEDRAL(FACE,I,TETRA,ERROR,*9999)
          ELSE ! Find any two matching nodes
            DO J=1,NJT
              DO K=J+1,NJT
                IF(INODE(J).EQ.INODE(K))THEN
                  ITETRA(1)=I
                  ITETRA(2)=J
                  ITETRA(3)=K
                  M=1
                  DO L=1,NFACES
                    IF((L.NE.J).AND.(L.NE.K))THEN
                      AXIS(M)=L
                      M=M+1
                    ENDIF !l.ne.j
                  ENDDO !l
                  CALL FLIP_EDGE_TO_FACE(AXIS,ELEM,FACE,ITETRA,
     '              TETRA,CIRC,NODE,RSQ,PLANE,ERROR,*9999)
                ENDIF !inode
              ENDDO !k
            ENDDO !j
          ENDIF !bool
        ENDIF !dabs(rsq(i))
        I=TETRA(I,NEXT) ! Increment counter
      ENDDO !while(i.ne.0)

      CALL EXITS('FIX_TETRAHEDRAL')
      RETURN
 9999 CALL ERRORS('FIX_TETRAHEDRAL',ERROR)
      CALL EXITS('FIX_TETRAHEDRAL')
      RETURN 1
      END


