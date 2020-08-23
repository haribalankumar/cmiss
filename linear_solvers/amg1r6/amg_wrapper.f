
      SUBROUTINE AMG_WRAPPER(AA,IA,JA,U,F,IG,
     &                       CMISS_AA,CMISS_IA,CMISS_JA,
     &                       CMISS_X,CMISS_B, 
     &                       AMG_BDY,AMG_NDIMS,AMG_IPARMS,
     &                       AMG_RPARMS,ERROR,*)

C#### Subroutine: AMG_WRAPPER
C###  Description:
C###    AMG_WRAPPER makes the call to AMG1R6 look simpler and hides some
C###    details from the user. 
C###  Written by Travis Austin 02/12/04

      IMPLICIT NONE
!     Parameter List
      INTEGER IA(*),JA(*),IG(*)
      REAL*8  AA(*),U(*),F(*)
      INTEGER CMISS_IA(*),CMISS_JA(*)
      REAL*8  CMISS_AA(*),CMISS_X(*),CMISS_B(*)
      INTEGER AMG_BDY(*),AMG_NDIMS(7), AMG_IPARMS(14)
      REAL*8  AMG_RPARMS(4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,I1,IERR,ISWTCH,J,J1,J2,JBEG,JEND,JBEG1,JEND1,
     &        JBEG2,JEND2,N,NDU,IRS,SPARSE_A
      INTEGER NDIFF
      REAL*8  ROWSUM
      
      ISWTCH = AMG_IPARMS(2)

      NDU = AMG_NDIMS(4)
      N   = AMG_NDIMS(7)

      SPARSE_A = 1

      IF( ISWTCH.NE.1 .AND. ISWTCH.NE.5 ) THEN

C$OMP PARALLEL DO 
         DO I=1,N
            U(I) = CMISS_X(I)
         ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO 
         DO I=N+1,NDU
            U(I) = 0.D0
            F(I) = 0.D0
         ENDDO
C$OMP END PARALLEL DO

C$OMP PARALLEL DO 
C$OMP&  PRIVATE(I,J,JBEG,JEND)
C$OMP&  SHARED(CMISS_IA,CMISS_JA,CMISS_AA,
C$OMP&         CMISS_B,AMG_BDY,F,N)
         DO I=1,N
            JBEG = CMISS_IA(I)
            JEND = CMISS_IA(I+1)-1
            F(I) = CMISS_B(I)
            DO J=JBEG,JEND
               IF( AMG_BDY(CMISS_JA(J)).EQ.0 .AND. ! UNKNOWN IS DIRICHLET
     &             I.NE.CMISS_JA(J) ) THEN         ! NOT DIAGONAL
                  F(I) = F(I) - CMISS_AA(J)*CMISS_B(CMISS_JA(J))
               ENDIF
            ENDDO
            IF( AMG_BDY(I).LT.0 ) THEN
               F(I) = -F(I)
            ENDIF
         ENDDO
C$OMP END PARALLEL DO

      ENDIF

      IRS = 0
      IF( ISWTCH.GE.4 ) THEN

C ----------------------------------------
C        CHECK FOR ZERO ROW SUM
C ----------------------------------------

C$OMP PARALLEL DO 
C$OMP& PRIVATE(I,J,JBEG,JEND,ROWSUM)
C$OMP& SHARED(IRS,AA,N)
         DO I=1,N
            ROWSUM=0.D0
            JBEG = IA(I)
            JEND = IA(I+1)-1
            DO J=JBEG,JEND
               ROWSUM = ROWSUM + AA(J)
            ENDDO
            IF( IRS.EQ.0 .AND. ABS(ROWSUM).GE.1E-10 ) IRS = 1
         ENDDO
C$OMP END PARALLEL DO

C ----------------------------------------
C        CHECK FOR SYMMETRY
C ----------------------------------------

         NDIFF = 0
C$OMP PARALLEL DO REDUCTION(+:NDIFF)
CC$OMP& DEFAULT(NONE)
C$OMP& PRIVATE(I1,JBEG1,JEND1,JBEG2,JEND2,J1,J2)
C$OMP& SHARED(AA,IA,JA,N)
         DO I1=1,N
            JBEG1 = IA(I1)+1
            JEND1 = IA(I1+1)-1
            DO J1=JBEG1,JEND1
               IF(JA(J1).GT.I1) THEN
                  JBEG2 = IA(JA(J1))+1
                  JEND2 = IA(JA(J1))-1
                  DO J2 = JBEG2,JEND2
                     IF( JA(J2).EQ.I1 .AND. AA(J2).NE.AA(J1) ) THEN
                        NDIFF = NDIFF + 1
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
C$OMP END PARALLEL DO

         IF( IRS.EQ.0 ) THEN
            IF( NDIFF.GT.0 ) THEN
               ! Zero Row Sum & Nonsymmetric
               AMG_IPARMS(1) = 21
            ELSEIF( NDIFF.EQ.0 ) THEN
               ! Zero Row Sum & Symmetric
               AMG_IPARMS(1) = 11
            ENDIF
         ELSEIF( IRS.EQ.1 ) THEN
            IF( NDIFF.GT.0 ) THEN
               ! Nonzero Row Sum & Nonsymmetric
               AMG_IPARMS(1) = 22
            ELSEIF( NDIFF.EQ.0 ) THEN
               ! Nonzero Row Sum & Symmetric
               AMG_IPARMS(1) = 12
            ENDIF
         ENDIF

      ENDIF

      CALL AMG1R6(AA,IA,JA,U,F,IG, 
     &            AMG_NDIMS(1),AMG_NDIMS(2),AMG_NDIMS(3),AMG_NDIMS(4), 
     &            AMG_NDIMS(5),AMG_NDIMS(6),AMG_NDIMS(7),
     &            AMG_IPARMS(1),AMG_IPARMS(2),AMG_IPARMS(3), 
     &            AMG_IPARMS(4),AMG_IPARMS(5),AMG_IPARMS(6), 
     &            AMG_IPARMS(7),AMG_IPARMS(8),AMG_RPARMS(1), 
     &            AMG_IPARMS(9),AMG_IPARMS(10),AMG_IPARMS(11), 
     &            AMG_IPARMS(12),AMG_RPARMS(2),AMG_RPARMS(3), 
     &            AMG_RPARMS(4),AMG_IPARMS(13),AMG_IPARMS(14),
     &            IERR)       



      IF( IERR.GT.0 ) GOTO 9999 

      IF( ISWTCH.NE.1 .AND. ISWTCH.NE.5 ) THEN 

C$OMP PARALLEL DO 
         DO I=1,N 
            CMISS_X(I) = U(I)  
         ENDDO 
C$OMP END PARALLEL DO 

      ENDIF 


      CALL EXITS('AMG_WRAPPER') 
      RETURN 
      
9999  CALL ERRORS('AMG_WRAPPER',ERROR) 
      RETURN 1
      END
