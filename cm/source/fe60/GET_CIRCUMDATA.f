      SUBROUTINE GET_CIRCUMDATA(ELEM,CIRC,NODE,RSQ,PLANE,ERROR,*)

C#### Subroutine: GET_CIRCUMDATA
C###  Description:
C###    Determines the tetrahedral circum-data. Genmesh routine.
CC JMB 15-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ
      LOGICAL PLANE
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 A(3,3),B(3),IPIV(3),SUM,SQR,VOLUME
      INTEGER I,INFO,J,K
!     Local Functions
      REAL*8 DETERMINANT,DLAMCH

      CALL ENTERS('GET_CIRCUMDATA',*9999)

      ! Determine the coordinates of the tetrahedral
      VOLUME=DABS(DETERMINANT(ELEM(1,2),ELEM(1,3),ELEM(1,4),NODE,LDNODE)
     '  -DETERMINANT(ELEM(1,1),ELEM(1,3),ELEM(1,4),NODE,LDNODE)
     '  +DETERMINANT(ELEM(1,1),ELEM(1,2),ELEM(1,4),NODE,LDNODE)
     '  -DETERMINANT(ELEM(1,1),ELEM(1,2),ELEM(1,3),NODE,LDNODE))
      IF(VOLUME.GT.SQRT(DLAMCH('E')))THEN
        I=1
        DO K=2,NPTS
          DO J=1,NJT
            A(I,J)=NODE(ELEM(1,K),J)-NODE(ELEM(1,1),J)
          ENDDO !j
          ! Increment counter
          I=I+1
        ENDDO !k
        SQR=ZERO
        DO I=1,NJT
          SQR=SQR+NODE(ELEM(1,1),I)**2.d0
        ENDDO
        J=2
        DO I=1,NJT
          B(I)=ZERO
          DO K=1,NJT
            B(I)=B(I)+NODE(ELEM(1,J),K)**2.d0
          ENDDO !k
          B(I)=0.5d0*(B(I)-SQR)
          J=J+1
        ENDDO !i
        ! Solve the set of linear equations
        CALL DGESV(NJT,1,A,NJT,IPIV,B,NJT,INFO)
        CALL ASSERT_VORO(INFO.EQ.0,
     '    'DGESV failed in routine GET_CIRCUMDATA',ERROR,*9999)
        RSQ=ZERO
        DO I=1,NJT
          RSQ=RSQ+(B(I)-NODE(ELEM(1,1),I))**2.d0
        ENDDO !i
        ! Tetrahedral circumdata
        CALL DCOPY(NJT,B,1,CIRC(1,1),LDCIRC)
      ELSE
        ! The tetrahedra has zero volume
        PLANE=.TRUE.
        DO J=1,NJT
          SUM=ZERO
          DO I=1,NFACES
            SUM=SUM+NODE(ELEM(1,I),J)
          ENDDO !i
          CIRC(1,J)=0.25d0*SUM
        ENDDO !j
        RSQ=ZERO
      ENDIF !volume

      CALL EXITS('GET_CIRCUMDATA')
      RETURN
 9999 CALL ERRORS('GET_CIRCUMDATA',ERROR)
      CALL EXITS('GET_CIRCUMDATA')
      RETURN 1
      END


