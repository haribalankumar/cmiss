      SUBROUTINE CALC_GLOBAL_NODES(NPNODE,ERROR,*)

C#### Subroutine: CALC_GLOBAL_NODES
C###  Description:
C###    CALC_GLOBAL_NODES calculates athe global nodes based on the
C###    nodes from each region (1..NRT).  The list is checked for
C###    duplicates.

C***  GMH 22/7/96 Why is npnode not dimensioned to NPM - the global
C***  list could contain more than NP_R_M nodes.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nr,nonode

      CALL ENTERS('CALC_GLOBAL_NODES',*9999)
C     Initialise the list
      NPNODE(0,0)=0
C     Add nodes from each region
      DO nr=1,NRT
        DO nonode=1,NPNODE(0,nr)
          NPNODE(0,0)=NPNODE(0,0)+1
          IF(NPNODE(0,0).LE.NP_R_M) THEN
            NPNODE(NPNODE(0,0),0)=NPNODE(nonode,nr)
          ENDIF
        ENDDO !nonode
      ENDDO !nr
      CALL ASSERT(NPNODE(0,0).LE.NP_R_M,'>> Increase NP_R_M',ERROR,
     '  *9999)
C     Check the list for duplicates
      CALL ILISTRMDUP(NPNODE(0,0),NPNODE(1,0),ERROR,*9999)
C     Warn if we did not find any nodes
      IF(NPNODE(0,0).EQ.0) THEN
        WRITE(OP_STRING,'('' >>WARNING!!! No global nodes found'')')
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('CALC_GLOBAL_NODES')
      RETURN
 9999 CALL ERRORS('CALC_GLOBAL_NODES',ERROR)
      CALL EXITS('CALC_GLOBAL_NODES')
      RETURN 1
      END

      
