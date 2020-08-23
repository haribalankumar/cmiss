      SUBROUTINE ZGTG55(nr,AZL,CG,EG,TG,ERROR,*)

C#### Subroutine: ZGTG55
C###  Description:
C###    ZGTG55 evaluates components of 2nd Piola-Kirchhoff stress
C###    tensor TG at current Gauss point for string problems
C###    (KTYP51(nr)=5).
C**** Stress components must always be with respect to
C**** (undeformed) Nu (body/fibre) coords (KTYP53(nr)>1).
C**** TGP are principal second Piola Kirchhoff stresses

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nr
      REAL*8 AZL(3,3),CG(NMM),EG(3,3),TG(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 STRAIN

      CALL ENTERS('ZGTG55',*9999)

      IF(KTYP53(nr).EQ.1) THEN      !stresses in theta coords
        ERROR='>>Stresses must be wrt body/fibre coords'
        GO TO 9999
      ELSE IF(KTYP53(nr).GT.1) THEN !stresses in body/fibre Nu coords
        STRAIN=0.50d0*(AZL(1,1)-1.0d0)
        IF(DABS(STRAIN).LT.1.0d+02) THEN
          EG(1,1)=STRAIN
        ELSE
          EG(1,1)=100.0d0
C          ERROR='>>>Strain cpt out of bounds'
C          GO TO 9999
        ENDIF
      ENDIF

      IF(KTYP55(nr).EQ.1) THEN !TG is fn of the principal invariants
        ERROR='>>Stresses must be fn of fibre strain'
        GO TO 9999
      ELSE IF(KTYP55(nr).EQ.2) THEN
        !TG is a fn of the principal stretches
        ERROR='>>Stresses must be fn of fibre strain'
        GO TO 9999
      ELSE IF(KTYP55(nr).EQ.3) THEN
        !TG is func of the fibre transv strains
        IF(KTYP52(nr).EQ.1) THEN !compressible
          IF(EG(1,1).GT.0.0d0) THEN !stretch
            TG(1,1)=CG(1)*EG(1,1)
          ELSE                     !compression
            TG(1,1)=0.0d0
          ENDIF
        ELSE                 !incompressible
          ERROR='>>Incompressible case not implemented'
          GO TO 9999
        ENDIF
      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,
     '    '('' AZL(1,1)='',D12.3,'' Strain=EG(1,1)='',D12.3,'
     '    //''' Stress=TG(1,1)='',D12.3)') AZL(1,1),EG(1,1),TG(1,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZGTG55')
      RETURN
 9999 CALL ERRORS('ZGTG55',ERROR)
      CALL EXITS('ZGTG55')
      RETURN 1
      END


