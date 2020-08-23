      SUBROUTINE OPNOIS(ERROR,*)

C#### Subroutine: OPNOIS
C###  Description:
C###    OPNOIS outputs noise parameters

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'nois00.cmn'
C LKC 27-OCT-97 not required
C      INCLUDE 'cmiss$reference:geom00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables


      CALL ENTERS('OPNOIS',*9999)

C***  Write noise types
      CALL ASSERT(2.GE.DEFNOIS_TYPE,
     '  '>>Code Needs Updating',ERROR,*9999)
      IF(DEFNOIS_TYPE.EQ.1) THEN
        WRITE(OP_STRING,'(/''  The type of noise is ABSOLUTE'')')
      ELSE
        WRITE(OP_STRING,'(/''  The type of noise is RELATIVE'')')
      ENDIF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL ASSERT(3.GE.NOIS_DIST,'>>Code Needs Updating',ERROR,*9999)
C***  Write distribution types
      IF(NOIS_DIST.EQ.1) THEN
        WRITE(OP_STRING,'(''  with GAUSSIAN distribution '')')
      ELSE IF(NOIS_DIST.EQ.2) THEN
        WRITE(OP_STRING,'(''  with UNIFORM distribution (ran0)'')')
      ELSE IF(NOIS_DIST.EQ.3) THEN
        WRITE(OP_STRING,'(''  with UNIFORM distribution (ran1)'')')
       ENDIF
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

C***  Write Noise Level
      WRITE(OP_STRING,'(''  and a noise level of '',D12.4)')
     '  NOIS_LEV
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPNOIS')
      RETURN
 9999 CALL ERRORS('OPNOIS',ERROR)
      CALL EXITS('OPNOIS')
      RETURN 1
      END


