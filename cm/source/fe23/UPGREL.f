      SUBROUTINE UPGREL(ne,NE_NEW,ERROR,*)

C#### Subroutine: UPGREL
C###  Description:
C###    UPGREL updates group element lists when mesh has been refined.

      IMPLICIT NONE
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER ne,NE_NEW
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER noelem,nogrel
      INTEGER*4 ILIST_PTR
!     Functions
      LOGICAL INLIST
      INTEGER*4 ILISTLOC

      CALL ENTERS('UPGREL',*9999)
      DO nogrel=1,NTGREL !Loop over groups
        IF(INLIST(ne,%VAL(LIGREL_PTR(nogrel)),NLIGREL(nogrel),noelem))
     '    THEN
          ILIST_PTR=0
          CALL ALLOCATE_MEMORY(NLIGREL(nogrel)+1,0,INTTYPE,ILIST_PTR,
     '      MEM_INIT,ERROR,*9999)
          CALL ILIST_COPY(NLIGREL(nogrel),
     '      %VAL(LIGREL_PTR(nogrel)),%VAL(ILIST_PTR))
          NLIGREL(nogrel)=NLIGREL(nogrel)+1
          CALL ILIST_COPY(1,
     '      NE_NEW,%VAL(ILISTLOC(%VAL(ILIST_PTR),NLIGREL(nogrel))))
          CALL FREE_MEMORY(LIGREL_PTR(nogrel),ERROR,*9999)
          LIGREL_PTR(nogrel)=ILIST_PTR
        ENDIF
C!!! This needs updating
C        DO noelem=1,0!LIGREL(0,nogrel) !Loop over elements in each group
C          IF(ne.EQ.LIGREL(noelem,nogrel))THEN
C!           !Element ne is in group nogrel
C            LIGREL(0,nogrel)=LIGREL(0,nogrel)+1
C            LIGREL(LIGREL(0,nogrel),nogrel)=NE_NEW
C            GOTO 100
C          ENDIF
C        ENDDO
C 100    CONTINUE
      ENDDO

      CALL EXITS('UPGREL')
      RETURN
 9999 CALL ERRORS('UPGREL',ERROR)
      CALL EXITS('UPGREL')
      RETURN 1
      END


