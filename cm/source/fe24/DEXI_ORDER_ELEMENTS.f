      SUBROUTINE DEXI_ORDER_ELEMENTS(NEELEM,NELIST2,NENP,NPNODE,X,XP,
     '  ERROR,*)

C#### Subroutine: DEXI_ORDER_ELEMENTS
C###  Description:
C###  DEXI_ORDER_ELEMENTS sets up a search order for the host elements,
C###  by first calculating the closest node in host region to the
C###  current point. All elements containing this node are set to be the
C###  first elements that will be searched, followed by all remaining
C###  elements in the region, in their defaul listed order.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NEELEM(0:NE_R_M),NELIST2(0:1000),
     '  NENP(NPM,0:NEPM),NPNODE(0:NP_R_M)
      REAL*8 X(3),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ne,ne2,nj,noelem,noelem2,noelem_count,nonode,np,np_closest
      REAL*8  DIST,MIN_DIST
      LOGICAL DONE_ELEMENT

      CALL ENTERS('DEXI_ORDER_ELEMENTS',*9999)

      MIN_DIST=1.d6
      DO nonode=1,NPNODE(0) !for each node in host
        np=NPNODE(nonode)
        DIST=0.d0
        DO nj=1,NJT
          DIST=DIST+(X(nj)-XP(1,1,nj,np))**2.d0
        ENDDO !nj
        DIST=DSQRT(DIST)
        IF(DIST.LT.MIN_DIST)THEN
          MIN_DIST=DIST
          np_closest=np
        ENDIF !DIST
      ENDDO !nonode (np)
      DO noelem=1,NENP(np_closest,0)
        ne=NENP(np_closest,noelem)
        NELIST2(noelem)=ne
      ENDDO !noelem (ne)
      NELIST2(0)=NENP(np,0)
      noelem_count=NELIST2(0)
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        DONE_ELEMENT=.FALSE.
        DO noelem2=1,NELIST2(0)
          ne2=NELIST2(noelem2)
          IF(ne.EQ.ne2) DONE_ELEMENT=.TRUE. !already entered into NELIST
        ENDDO !noelem2
        IF(.NOT.DONE_ELEMENT)THEN
          noelem_count=noelem_count+1
          NELIST2(noelem_count)=ne
        ENDIF
      ENDDO !noelem (ne)
      NELIST2(0)=noelem_count

      CALL EXITS('DEXI_ORDER_ELEMENTS')
      RETURN
 9999 CALL ERRORS('DEXI_ORDER_ELEMENTS',ERROR)
      CALL EXITS('DEXI_ORDER_ELEMENTS')
      RETURN 1
      END



