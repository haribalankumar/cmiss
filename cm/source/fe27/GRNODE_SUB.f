      SUBROUTINE GRNODE_SUB(NPLIST,STRING2,AS,ERROR,*)

C#### Subroutine: GRNODE_SUB
C###  Description:
C###    GRNODE_SUB groups nodes without necessarily going through
C###    the command line.  This routine is called from GRNODE to group
C###    nodes in NPLIST, or it can be called from another subroutine
C###    to set up a group of nodes, eg. during mesh generation of
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
      INTEGER NPLIST(0:NPM)
      CHARACTER ERROR*(*),STRING2*(255)
      LOGICAL AS
!     Local Variables
      INTEGER IBEG,IBEG2,IEND,IEND2,nogrno,N1GRNO
      CHARACTER CHAR*30,CHAR2*2,LABEL*30


      CALL ENTERS('GRNODE_SUB',*9999)

      IF(AS)THEN
 !Check whether group name already exists
        CALL STRING_TRIM(STRING2,IBEG,IEND)
C          CHAR=CUPPER(CO(N3CO+1)(IBEG:IEND))
        CALL CUPPER(STRING2(IBEG:IEND),CHAR)
        N1GRNO=0
        DO nogrno=1,NTGRNO
C            LABEL=CUPPER(LAGRNO(nogrno))
          CALL CUPPER(LAGRNO(nogrno),LABEL)
          CALL STRING_TRIM(LABEL,IBEG2,IEND2)
          IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
            N1GRNO=nogrno !is existing group label ID
            GO TO 100
          ENDIF
        ENDDO
 100    IF(N1GRNO.EQ.0) THEN !need new group label
          NTGRNO=NTGRNO+1 !is new total #groups
          CALL ASSERT(NTGRNO.LE.GRNO_MAXGRP,
     '      '>>Increase array sizes in grou00.cmn',ERROR,*9999)
          N1GRNO=NTGRNO
C          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          LAGRNO(N1GRNO)=STRING2(IBEG:IEND) !new node group label
        ENDIF

      ELSE !use the total #groups as default label
        NTGRNO=NTGRNO+1 !is new total
        CALL ASSERT(NTGRNO.LE.GRNO_MAXGRP,
     '    '>>Increase array sizes in grou00.cmn',ERROR,*9999)
        N1GRNO=NTGRNO
C          CHAR2=CFROMI(N1GRNO,'(I2)')
        WRITE(CHAR2,'(I2)') N1GRNO
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
        LAGRNO(N1GRNO)='node_'//CHAR2(IBEG2:IEND2) !new label
      ENDIF

C TVK 02July1999: dynamic groups
C        LIGRNO(0,N1GRNO)=NPLIST(0)
C        DO nolist=1,NPLIST(0)
C          LIGRNO(nolist,N1GRNO)=NPLIST(nolist)
C        ENDDO

      NLIGRNO(N1GRNO)=NPLIST(0)
      CALL ALLOCATE_MEMORY(NPLIST(0),0,INTTYPE,LIGRNO_PTR(N1GRNO),
     '  MEM_INIT,ERROR,*9999)
      CALL ILIST_COPY(NPLIST(0),NPLIST(1),%VAL(LIGRNO_PTR(N1GRNO)))

      CALL EXITS('GRNODE_SUB')
      RETURN
 9999 CALL ERRORS('GRNODE_SUB',ERROR)
      CALL EXITS('GRNODE_SUB')
      RETURN 1
      END


