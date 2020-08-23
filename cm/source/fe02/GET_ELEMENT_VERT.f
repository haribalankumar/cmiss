      SUBROUTINE GET_ELEMENT_VERT(NBJ,NNB,NPNE,
     '  NVERTNE,ERROR,*)
      
C#### Subroutine: GET_ELEMENT_VERT
C###  Description:
C###    GET_ELEMENT_VERT finds the vertices of element ne and
C###    returns them in the vector NVERTNE.   

C#### Variable: NVERTNE
C###  Type: INTEGER      
C###  Set_up: GET_ELEMENT_VERT      
C###  Description:
C###    NVERTNE holds the vertices of an element
C###    calculated in the subroutine GET_ELEMENT_VERT
C###    The 0th position
C###    of NVERTNE stores the dimensionality of the element.
C###    1, 2 and 3d elements are all given as having 8
C###    vertices, with the vertices repeated for 1 and 2d
C###    elements. Column 1 of NVERTNE is the local node numers
C###    and column 2 is the global node numbers.
      
      IMPLICIT NONE

      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM),NNB(4,4,4,NBFM),NPNE(NNM,NBFM),NVERTNE(0:8,2)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,INPI(3),INPIT(3),j,nb,ni,NITB
      LOGICAL NOTFOUND      

      CALL ENTERS('GET_ELEMENT_VERT',*9999)

      nb=NBJ(1)
      NITB=NIT(nb) !assume <= 3

C     Find the number of nodes in each direction.
C     Assumes that the largest number of nodes in each
C     direction can be found by searching the first line in NNB.
      DO ni=1,3
        INPI(ni)=1
      ENDDO !ni
      DO ni=1,3
        INPI(ni)=2
        NOTFOUND=.TRUE.
        DO WHILE(INPI(ni).LE.4.AND.NOTFOUND)
          IF(NNB(INPI(1),INPI(2),INPI(3),nb).EQ.0) THEN
            NOTFOUND=.FALSE.
          ELSE
            INPI(ni)=INPI(ni)+1
          ENDIF
        ENDDO !INPI(ni)
        INPIT(ni)=INPI(ni)-1
        INPI(ni)=1
      ENDDO !ni     

C     Write to NVERTNE

      NVERTNE(0,1)=NITB ! 0th indice stores dimension of element
C     Write to each of the eight vertex positions in NVERTNE(nn,1)
C     the appropriate local node number taken from NNB, utilizing
C     the vertex indices calculated above.
      NVERTNE(1,1)=NNB(1,1,1,nb)
      NVERTNE(2,1)=NNB(INPIT(1),1,1,nb) 
      NVERTNE(3,1)=NNB(1,INPIT(2),1,nb)
      NVERTNE(4,1)=NNB(INPIT(1),INPIT(2),1,nb)
      NVERTNE(5,1)=NNB(1,1,INPIT(3),nb)
      NVERTNE(6,1)=NNB(INPIT(1),1,INPIT(3),nb)
      NVERTNE(7,1)=NNB(1,INPIT(2),INPIT(3),nb) 
      NVERTNE(8,1)=NNB(INPIT(1),INPIT(2),INPIT(3),nb)      
      
C     Global vertex numbers
C     Fill in the global node numbers for the vertices,
C     using NPNE.
      DO i=1,8
        NVERTNE(i,2)=NPNE(NVERTNE(i,1),nb)
      ENDDO

C     Now loop over the local and global vertex numbers,
C     finding the local vertex repitition. This loops
C     over the eight vertices global numbers
C     and compares them to the unchecked vertices global
C     numbers . If global numbers are repeated then the
C     local vertex number is written over to reflect
C     this repetition.
      
      DO i=1,8
        DO j=i+1,8
          IF(NVERTNE(i,2).EQ.NVERTNE(j,2)) THEN
            NVERTNE(j,1)=NVERTNE(i,1)
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('GET_ELEMENT_VERT')
      RETURN
 9999 CALL ERRORS('GET_ELEMENT_VERT',ERROR)
      CALL EXITS('GET_ELEMENT_VERT')
      RETURN 1
      END

      
