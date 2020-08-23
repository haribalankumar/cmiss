      SUBROUTINE BMG2_PerSymStd_SETUP_cg_LU(
     &                SO, II, JJ, NStncl, ABD, NABD1, NABD2, JPN 
     &                )

C ===========================================================================
C
C***BEGIN PROLOGUE  BMG2_PerSymStd_SETUP_cg_LU
C***SUBSIDIARY
C***PURPOSE  BMG2_PerSymStd_SETUP_cg_LU copies the matrix on the 
C            coarsest grid in an LAPACK array, and then uses an
C            LAPACK routine to form the l-u decomposition.
C
C***LIBRARY   SLATEC
C***AUTHOR  DENDY, J. E. JR.
C             LOS ALAMOS NATIONAL LABORATORY
C           VOYTKO, M. H.
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C
C          Copyright, 1988. The Regents of the University of California.
C          This software was produced under a U.S. Government contract
C          (W-7405-ENG-36) by Los Alamso National Laboratory, which is
C          operated by the University of California for the U.S. Depart-
C          ment of Energy. The U.S. Government is licensed to use,
C          reproduce, and distribute this software. Permission is
C          granted to the public to copy and use this software without
C          charge, provided that this Notice and any statement of
C          authorship are reproduced on all copies. Neither the
C          Government nor the University makes any warranty, expressed
C          or implied, or assumes any liability or responsibility for
C          the use of this software.
C
C***PARAMETERS
C***INPUT
C   SO        Refer to BOXMG.
C   II        Number of grid points in x direction, including
C             two fictitious points.
C   JJ        Number of grid points in y direction, including
C             two fictitious points.
C   ABD       Refer to BOXMG.
C   IPN       Refer to BOXMG.
C
C***ROUTINES CALLED  SPBFA, SPOFA
C
C***REVISION HISTORY  (YYYYMMDD)
C   1983/09/25  * DATE WRITTEN
C   1990/06/27  - modified to conform to the 4/10/90 "Guide to the SLATEC
C                 Common Mathematical Library" by Victor A. Bandy
C   2000/01/10  - updated declarations, mostly cosmetic 
C               - added INCLUDE to initialize constants (J. David Moulton) 
C
C***END PROLOGUE  MGSETP
C
C ===========================================================================

      IMPLICIT NONE

C ----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C ---------------------------
C     Arguments Declarations
C
      INTEGER JPN, II, JJ, NABD1, NABD2, NStncl
      REAL*8  ABD(NABD1,NABD2), SO(II,JJ,NStncl)

C ---------------------------
C     Local Declarations
C
      INTEGER I, INFO, IPN, I1, I2, J, J1, J2, KK, N

C ==========================================================================

C -------------------------------------------------------
C     Copy the operator on the coarsest grid into ABD 
C -------------------------------------------------------

      IPN=IABS(JPN)
      I1=II-1
      J1=JJ-1
      I2=I1-1
      N=I2*(J1-1)
      KK=0

      IF ( NStncl.EQ.5 ) THEN

         IF ( IPN.EQ.4 ) THEN
            !
            ! Nonperiodic, but possibly indefinite
            !
            DO J=2, J1
               DO I=2, I1
                  KK=KK+1
                  ABD(II,KK) =  SO(I,J,KO)
                  ABD(I1,KK) = -SO(I,J,KW)
                  ABD(3,KK)  = -SO(I+1,J,KNW)
                  ABD(2,KK)  = -SO(I,J,KS)
                  ABD(1,KK)  = -SO(I,J,KSW)
               ENDDO
            ENDDO
            !
            ! Indefinite ...
            !
            IF ( JPN.LT.0 ) THEN
               ABD(II,KK)=ABD(II,KK)+SO(I1,J1,KO)
            ENDIF
            !
            ! Factor using the LAPACK routine DPBTRF
            !
            CALL DPBTRF('U',N,I1,ABD,NABD1,INFO)
            !
         ELSEIF( (JPN.EQ.(-3)).AND.(II.EQ.4).AND.(JJ.EQ.4)) THEN
            !
            ! Periodic in both x and y, but only 4 computational nodes
            !
            ABD(1,1) = SO(2,2,KO)
            ABD(2,2) = SO(3,2,KO)
            ABD(1,2) = -SO(2,2,KW) - SO(3,2,KW)
            ABD(3,3) = SO(2,3,KO)
            ABD(1,3) = -SO(2,2,KS) - SO(2,3,KS)
            ABD(2,3) = - SO(3,3,KNW) - SO(4,2,KNW)
     &                 - SO(4,3,KSW) - SO(3,2,KSW)
            ABD(4,4) = SO(3,3,KO)
            ABD(1,4) = - SO(3,3,KSW) - SO(2,2,KSW)
     &                 - SO(2,3,KNW) - SO(3,2,KNW)
            ABD(2,4) = -SO(3,2,KS) - SO(3,3,KS)
            ABD(3,4) = -SO(2,3,KW) - SO(3,3,KW)
            !
            ! Indefinite ...
            !
            ABD(4,4) = ABD(4,4) + SO(4,4,KO)
            !
            ! Factor using the LAPACK routine DPOTRF
            !
            CALL DPOTRF('U',N,ABD,NABD1,INFO)
            !
         ELSEIF( IPN.EQ.1.OR.IPN.EQ.2.OR.IPN.EQ.3) THEN
            !
            ! Periodic in y, x or both
            !
            KK=1
            ABD(1,1) = SO(2,2,KO)
            DO I=3,I1
               KK=KK+1
               ABD(KK,KK) = SO(I,2,KO)
               ABD(KK-1,KK) = -SO(I,2,KW)
            ENDDO
            
            IF ( IPN.EQ.2.OR.IPN.EQ.3 ) THEN
               ABD(KK-I2+1,KK) = -SO(II,2,KW)
            ENDIF

            DO J=3,J1
               
               IF ( IPN.EQ.2.OR.IPN.EQ.3 ) THEN
                  ABD(KK,KK+1) = -SO(2,J,KSW)
               ENDIF
               
               DO I=2,I1
                  KK=KK+1
                  ABD(KK,KK) = SO(I,J,KO)
                  IF (I.NE.2) THEN 
                     ABD(KK-1,KK) = -SO(I,J,KW)
                     ABD(KK-I2-1,KK) = -SO(I,J,KSW)
                  ENDIF
                  ABD(KK-I2+1,KK) = -SO(I+1,J,KNW)
                  ABD(KK-I2,KK) = -SO(I,J,KS)
               ENDDO

               IF( IPN.EQ.2 .OR. IPN.EQ.3 ) THEN
                  ABD(KK-I2+1,KK) = -SO(II,J,KW)
                  ABD(KK-2*I2+1,KK) = -SO(II,J,KNW)
               ENDIF

            ENDDO
            
            IF ( IPN.EQ.1 .OR. IPN.EQ.3 ) THEN 
               
               KK=KK-I2
               J2=(J1-2)*I2
               KK=KK+1
               ABD(KK-J2,KK) = -SO(2,JJ,KS)
               ABD(KK-J2+1,KK) = -SO(3,JJ,KSW)

               IF ( IPN.EQ.3 ) THEN
                  ABD(I2,KK)=-SO(2,JJ,KNW)
               ENDIF

               DO I=3,I1
                  KK=KK+1
                  ABD(KK-J2,KK)   = -SO(I,JJ,KS)
                  ABD(KK-J2-1,KK) = -SO(I,JJ,KNW)
                  ABD(KK-J2+1,KK) = -SO(I+1,JJ,KSW)
               ENDDO
               
               ABD(KK-J2+1,KK) = RZERO
               
               IF( IPN.EQ.3 ) THEN
                  ABD(1,KK) = -SO(II,JJ,KSW)
               ENDIF
               
               IF( IPN.EQ.2.OR.IPN.EQ.3 ) THEN
                  ABD(KK-2*I2+1,KK) = -SO(II,J1,KNW)
               ENDIF
               
            ENDIF
            
            IF ( JPN.LT.0 ) THEN
               ABD(KK,KK) = ABD(KK,KK) + SO(I1,J1,KO)
            ENDIF

            !
            ! Factor using the LAPACK routine DPOTRF
            !
            CALL DPOTRF('U',N,ABD,NABD1,INFO)
            
         ENDIF

      ELSEIF ( NStncl.EQ.3 ) THEN

         IF ( IPN.EQ.4 ) THEN
           !
           ! Nonperiodic, but possibly indefinite
           !
            DO J=2, J1
               DO I=2, I1
                  KK=KK+1
                  ABD(II,KK) =  SO(I,J,KO)
                  ABD(I1,KK) = -SO(I,J,KW)
                  ABD(3,KK)  =  rZERO
                  ABD(2,KK)  = -SO(I,J,KS)
                  ABD(1,KK)  =  rZERO
               ENDDO
            ENDDO
            !
            ! Indefinite ...
            !
            IF ( JPN.LT.0 ) THEN
               ABD(II,KK) = ABD(II,KK) + SO(I1,J1,KO)
            ENDIF
            !
            ! Factor using the LAPACK routine DPBTRF
            !
            CALL DPBTRF('U',N,I1,ABD,NABD1,INFO)
            !
         ELSEIF( IPN.EQ.1.OR.IPN.EQ.2.OR.IPN.EQ.3) THEN
            !
            ! Periodic in y, x or both
            !      
            KK=1
            ABD(1,1) = SO(2,2,KO)
            DO I=3,I1
               KK=KK+1
               ABD(KK,KK) = SO(I,2,KO)
               ABD(KK-1,KK) = -SO(I,2,KW)
            ENDDO
            
            IF ( IPN.EQ.2.OR.IPN.EQ.3 ) THEN
               ABD(KK-I2+1,KK) = -SO(II,2,KW)
            ENDIF

            DO J=3,J1
               
               IF ( IPN.EQ.2.OR.IPN.EQ.3 ) THEN
                  ABD(KK,KK+1) = rZERO
               ENDIF
               
               DO I=2,I1
                  KK=KK+1
                  ABD(KK,KK) = SO(I,J,KO)
                  IF ( I.NE.2 ) THEN
                     ABD(KK-1,KK) = -SO(I,J,KW)
                     ABD(KK-I2-1,KK) = rZERO
                  ENDIF
                  ABD(KK-I2+1,KK) = rZERO
                  ABD(KK-I2,KK) = -SO(I,J,KS)
               ENDDO

               IF( IPN.EQ.2 .OR. IPN.EQ.3 ) THEN
                  ABD(KK-I2+1,KK) = -SO(II,J,KW)
                  ABD(KK-2*I2+1,KK) = rZERO
               ENDIF

            ENDDO
            
            IF ( IPN.EQ.1 .OR. IPN.EQ.3 ) THEN 
               
               KK=KK-I2
               J2=(J1-2)*I2
               KK=KK+1
               ABD(KK-J2,KK) = -SO(2,JJ,KS)
               ABD(KK-J2+1,KK) = rZERO

               IF ( IPN.EQ.3 ) THEN
                  ABD(I2,KK) = rZERO
               ENDIF

               DO I=3,I1
                  KK=KK+1
                  ABD(KK-J2,KK) = -SO(I,JJ,KS)
                  ABD(KK-J2-1,KK) = rZERO
                  ABD(KK-J2+1,KK) = rZERO
               ENDDO

               ABD(KK-J2+1,KK) = rZERO
               
               IF( IPN.EQ.3 ) THEN
                  ABD(1,KK) = rZERO
               ENDIF
               
               IF( IPN.EQ.2.OR.IPN.EQ.3 ) THEN
                  ABD(KK-2*I2+1,KK) = rZERO
               ENDIF
               
            ENDIF
            
            IF ( JPN.LT.0 ) THEN
               ABD(KK,KK) = ABD(KK,KK) + SO(I1,J1,KO)
            ENDIF

            !
            ! Factor using the LAPACK routine DPOTRF
            !
            CALL DPOTRF('U',N,ABD,NABD1,INFO)
            
         ENDIF
         !
      ENDIF





      RETURN
      END
      
      
