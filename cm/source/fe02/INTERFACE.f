      SUBROUTINE INTERFACE(NP_INTERFACE,NPNODE,ERROR,*)

C#### Subroutine: INTERFACE
C###  Description:
C###    INTERFACE finds regions that a node belongs to.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NP_INTERFACE(0:NPM,0:3),NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n1node,n2node,nonode,np,nr,nr1,nr2,nr3
      LOGICAL FOUND

      CALL ENTERS('INTERFACE',*9999)
      DO nr=1,NRT
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          NP_INTERFACE(np,0)=1 !each node in at least one region
          NP_INTERFACE(np,1)=0 !initialize to no interface
        ENDDO !nonode (np)
      ENDDO !nr

      DO nr1=1,NRT
        DO n1node=1,NPNODE(0,nr1)
          np=NPNODE(n1node,nr1)
          DO nr2=1,nr1-1
            DO n2node=1,NPNODE(0,nr2)
              IF(np.EQ.NPNODE(n2node,nr2)) THEN !np shared between nr1/2
                FOUND=.FALSE.
                DO nr3=1,NP_INTERFACE(np,0) !check not already listed
                  IF(nr1.EQ.NP_INTERFACE(np,nr3)) FOUND=.TRUE.
                ENDDO
                IF(.NOT.FOUND) THEN
                  NP_INTERFACE(np,0)=NP_INTERFACE(np,0)+1
                  IF(NP_INTERFACE(np,0).LE.3)
     '              NP_INTERFACE(np,NP_INTERFACE(np,0))=nr1
                ENDIF
              ENDIF
            ENDDO !n2node
          ENDDO !nr2
          IF(NP_INTERFACE(np,1).EQ.0) NP_INTERFACE(np,1)=nr1
          CALL ASSERT(NP_INTERFACE(np,0).LE.3,'>>np is defined in more '
     '      //'than three regions: Increase 2nd dim of NP_INTERFACE',
     '      ERROR,*9999)
        ENDDO !n1node
      ENDDO !nr1

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nr=1,NRT
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            WRITE(OP_STRING,*)' NP_INTERFACE(',np,',0)=',
     '        NP_INTERFACE(np,0)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            DO nr1=1,NP_INTERFACE(np,0)
              WRITE(OP_STRING,*)' NP_INTERFACE(np,',nr1,')=',
     '          NP_INTERFACE(np,nr1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nr1
          ENDDO !nonode (np)
        ENDDO !nr
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('INTERFACE')
      RETURN
 9999 CALL ERRORS('INTERFACE',ERROR)
      CALL EXITS('INTERFACE')
      RETURN 1
      END


