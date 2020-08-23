      SUBROUTINE GRELEM_SUB(NELIST,STRING2,AS,ERROR,*)

C#### Subroutine: GRELEM_SUB
C###  Description:
C###    GRELEM_SUB groups elements without necessarily going through
C###    the command line.  This routine is called from GRELEM to group
C###    elements in NELIST, or it can be called from another subroutine
C###    to set up a group of elements, eg. during mesh generation of
C###    specific structures.

C****
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NELIST(0:NEM)
      LOGICAL AS
      CHARACTER ERROR*(*),STRING2*(255)
!     Local Variables
      INTEGER IBEG,IBEG2,IEND,IEND2,nogrel,N1GREL
      CHARACTER CHAR*30,CHAR2*2,LABEL*30

      CALL ENTERS('GRELEM_SUB',*9999)

      IF(AS)THEN
C     Check whether group name already exists
        CALL STRING_TRIM(STRING2,IBEG,IEND)
        CALL CUPPER(STRING2(IBEG:IEND),CHAR)
        N1GREL=0
        DO nogrel=1,NTGREL
          CALL CUPPER(LAGREL(nogrel),LABEL)
          CALL STRING_TRIM(LABEL,IBEG2,IEND2)
          IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
            N1GREL=nogrel !is existing group label ID
            GO TO 100
          ENDIF
        ENDDO
 100    IF(N1GREL.EQ.0) THEN !need new group label
          NTGREL=NTGREL+1 !is new total #groups
          CALL ASSERT(NTGREL.LE.GREL_MAXGRP,
     '      '>>Increase array sizes in grou00.cmn',ERROR,*9999)
          N1GREL=NTGREL
          LIGREL_PTR(N1GREL)=0 !hasn't already been allocated
c          CALL STRING_TRIM(STRING(N3CO+1),IBEG,IEND)
          LAGREL(N1GREL)=STRING2(IBEG:IEND) !new elem group label
        ENDIF
      ELSE
        NTGREL=NTGREL+1 !is new total
        CALL ASSERT(NTGREL.LE.GREL_MAXGRP,
     '    '>>Increase array sizes in grou00.cmn',ERROR,*9999)
        N1GREL=NTGREL
        LIGREL_PTR(N1GREL)=0 !hasn't already been allocated
C          CHAR2=CFROMI(N1GREL,'(I2)')
        WRITE(CHAR2,'(I2)') N1GREL
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        LAGREL(N1GREL)='element_'//CHAR2(IBEG2:IEND2) !new label
      ENDIF

      NLIGREL(N1GREL)=NELIST(0)
      CALL ALLOCATE_MEMORY(NELIST(0),0,INTTYPE,LIGREL_PTR(N1GREL),
     '  MEM_INIT,ERROR,*9999)
C     Copy NELIST into LIGREL
      CALL ILIST_COPY(NELIST(0),NELIST(1),%VAL(LIGREL_PTR(N1GREL)))

      CALL EXITS('GRELEM_SUB')
      RETURN
 9999 CALL ERRORS('GRELEM_SUB',ERROR)
      CALL EXITS('GRELEM_SUB')
      RETURN 1
      END


