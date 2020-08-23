
      SUBROUTINE AMG_COPY_MATRIX(N,NDA,NDU,NZA,
     &     ISC_A,ISR_A,A,AMG_A,AMG_IA,AMG_JA,AMG_BDY,
     &     SPARSE_A,ERROR,*)

C#### Subroutine: AMG_COPY_MATRIX
C###  Description:
C###    AMG_COPY_MATRIX copies the matrix in CMISS format into the
C###    matrix style used by AMG.  The arrays used for storage for
C###    amg contain space for all grids so we place it in first part
C###    of array.
C###  Written by Travis Austin 07/12/04

      IMPLICIT NONE
C     Parameter List
      INTEGER N,NDA,NDU,NZA,SPARSE_A
      INTEGER ISC_A(*),ISR_A(*)
      REAL*8  A(*)
      INTEGER AMG_IA(*),AMG_JA(*),AMG_BDY(*)
      REAL*8  AMG_A(*)
      CHARACTER ERROR*(*)
C     Local Variables
      INTEGER I,ITMP,J,SGN
      REAL*8  RTMP

      CALL ENTERS('AMG_COPY_MATRIX',*9999)

C ==============================================
C First zero out coarse level data (necessary?)
C ==============================================

C$OMP PARALLEL DO
      DO I=N+2,NDU
         AMG_IA(I) = 0
      ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
      DO I=NZA+1,NDA
        AMG_A(I)  = 0.D0
        AMG_JA(I) = 0
      ENDDO
C$OMP END PARALLEL DO
      
C ============================================================
C Copy Matrix into AMG equivalents (depends on storage format)
C ============================================================


      IF( SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4 ) THEN

C$OMP PARALLEL DO
        DO I=1,N+1
          AMG_IA(I) = ISR_A(I)
        ENDDO
C$OMP END PARALLEL DO
              
C$OMP PARALLEL DO
        DO I=1,N
          AMG_BDY(I) = 0
        ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I,ITMP,J,RTMP,SGN)
C$OMP& SHARED(A,AMG_A,AMG_BDY,
C$OMP&        AMG_IA,AMG_JA,ISC_A,N)
        DO I=1,N

           SGN = 1
           DO J=AMG_IA(I),AMG_IA(I+1)-1
              AMG_JA(J) = ISC_A(J)
              AMG_A(J)  = A(J)
              IF( AMG_JA(J).EQ.I .AND. AMG_A(J).LE.0.D0 ) SGN = -1
           ENDDO

           ! Move diagonal to first entry in row
           DO J=AMG_IA(I),AMG_IA(I+1)-1
              AMG_A(J) = SGN*AMG_A(J)              
              IF( AMG_JA(J).EQ.I ) THEN
                 RTMP = AMG_A(AMG_IA(I))
                 ITMP = AMG_JA(AMG_IA(I))
                 AMG_A(AMG_IA(I))  = AMG_A(J)
                 AMG_JA(AMG_IA(I)) = AMG_JA(J)
                 AMG_A(J)  = RTMP 
                 AMG_JA(J) = ITMP 
              ENDIF
           ENDDO

           DO J=AMG_IA(I)+1,AMG_IA(I+1)-1
              IF( AMG_A(J).NE.0.D0 ) AMG_BDY(I) = SGN
           ENDDO
           
        ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO
C$OMP& PRIVATE(I,J)
C$OMP& SHARED(AMG_A,AMG_BDY,AMG_IA,AMG_JA)
        DO I=1,N
           DO J=AMG_IA(I)+1,AMG_IA(I+1)-1
              IF( AMG_BDY(AMG_JA(J)).EQ.0 ) AMG_A(J) = 0.D0
           ENDDO
        ENDDO
C$OMP END PARALLEL DO

      ENDIF

      CALL EXITS('AMG_COPY_MATRIX')
      RETURN

 9999 CALL ERRORS('AMG_COPY_MATRIX',ERROR)
      CALL EXITS('AMG_COPY_MATRIX')
      RETURN 1
      END
