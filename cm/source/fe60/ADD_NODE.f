      SUBROUTINE ADD_NODE(ELEM,FACE,FOUND,INODE,TETRA,CIRC,NODE,
     '  RSQ,PLANE,ERROR,*)

C#### Subroutine: ADD_NODE
C###  Desription:
C###    Adds a node to the geometry. Genmesh routine.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),FACE(LDFACE,NFACES),FOUND(0:N_GM),
     '  INODE, TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ(N_GM)
      LOGICAL PLANE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,L,M,MATCH,NBOUR,REF,WORK(3)
      INTEGER FLIST(N_GM,3),IPTR(N_GM),NEIGH(N_GM),PTR(N_GM)
      LOGICAL ACTIVE(N_GM),BOOL
      CHARACTER CHAR*1
!     External functions
      INTEGER IDIGITS

      CALL ENTERS('ADD_NODE',*9999)

      ! Initialization
      PLANE=.FALSE.
      ! Construct facelist including duplicates
      M=1
      DO I=1,FOUND(0)
        REF=FOUND(I)
        ! Circulate through the faces of the tetrahedral
        DO J=NFACES,1,-1
          L=1
          DO K=1,NFACES
            IF(J.NE.K)THEN
              WORK(L)=ELEM(REF,K)
              L=L+1
            ELSE
              ! Neighbouring face
              NBOUR=FACE(REF,K)
            ENDIF !j.ne.k
          ENDDO !k
          ! Sort the global face node numbers
          CALL INSERTION(NJT,WORK,ERROR,*9999)
          DO K=1,NJT
            FLIST(M,K)=WORK(K)
          ENDDO !k
          NEIGH(M)=NBOUR
          ACTIVE(M)=.TRUE.
          IF(NBOUR.NE.0)THEN
            K=1
            CALL ASSERT(NBOUR.LE.LDFACE,
     '        '>>Increase LDFACE in genmesh.cmn',ERROR,*9999)
            CALL ASSERT(K.LE.NFACES,
     '        '>>Increase NFACES in genmesh.cmn',ERROR,*9999)
            DO WHILE((FACE(NBOUR,K).NE.REF).AND.(K.LT.NFACES))
              K=K+1
            ENDDO !while
            PTR(M)=NBOUR
            IPTR(M)=K
          ENDIF !nbour
          M=M+1! Increment face list
        ENDDO !j
      ENDDO !i
      ! Eliminate duplicate faces
      DO I=1,NFACES*FOUND(0)
        IF(ACTIVE(I))THEN
          DO J=1,NFACES*FOUND(0)
            IF(I.NE.J)THEN
              BOOL=.TRUE.
              DO K=1,NJT
                BOOL=BOOL.AND.(FLIST(I,K).EQ.FLIST(J,K))
              ENDDO !k
              IF(BOOL)THEN ! Duplicate face
                ACTIVE(I)=.FALSE.
                ACTIVE(J)=.FALSE.
              ENDIF !bool
            ENDIF !i.ne.j
          ENDDO !j
        ENDIF !active
      ENDDO !i
      ! Delete the insertion polygon
      DO I=1,FOUND(0)
        REF=FOUND(I)
        IF(REF.EQ.TETRA(0,NEXT))THEN
          TETRA(0,NEXT)=TETRA(REF,NEXT)
        ENDIF !ref
        ! Delete the tetrahedral
        CALL DELETE_TETRAHEDRAL(REF,TETRA,ERROR,*9999)
      ENDDO !i
      ! Replace existing triangles in insertion polygon with new ones
      DO I=1,NFACES*FOUND(0)
        IF(ACTIVE(I))THEN
          ! A new tetrahedral is to be added with this face
          CALL ADD_TETRAHEDRAL(TETRA,ERROR,*9999)
          DO J=1,(NPTS-1)
            ELEM(IELEM,J)=FLIST(I,J)
          ENDDO !j
          IF(IELEM.GT.LDELEM) THEN
            WRITE(CHAR,'(I1)') IDIGITS(IELEM)
            WRITE(ERROR,'(''>>Increase LDELEM to '',I'//CHAR//')') IELEM
            GO TO 9999
          ENDIF
c          CALL ASSERT(IELEM.LE.LDELEM,'>>Increase LDELEM',ERROR,*9999)
          ELEM(IELEM,NPTS)=INODE
          
          CALL ASSERT(IELEM.LE.LDFACE,'>>Increase LDFACE',ERROR,*9999)
          FACE(IELEM,NFACES)=NEIGH(I)
          IF(NEIGH(I).NE.0)THEN
            ! Assign the curent tetrehedral to the pointer face
            FACE(PTR(I),IPTR(I))=IELEM
          ENDIF !neigh
          NEIGH(I)=IELEM
          ! Determine the circumdata
          CALL ASSERT(IELEM.LE.LDCIRC,'>>Increase LDCIRC',ERROR,*9999)
          CALL ASSERT(IELEM.LE.N_GM,'>>Increase N_GM',ERROR,*9999)
          CALL GET_CIRCUMDATA(ELEM(IELEM,1),CIRC(IELEM,1),NODE,
     '      RSQ(IELEM),PLANE,ERROR,*9999)
        ENDIF !active
      ENDDO !i
      ! Complete adjacency structure matching
      DO I=1,NFACES*FOUND(0)
        IF(ACTIVE(I))THEN
          DO J=1,(NFACES)*FOUND(0)
            IF(ACTIVE(J).AND.(I.NE.J))THEN
              MATCH=0
              DO K=1,(NPTS-1)
                WORK(K)=0
                DO L=1,(NPTS-1)
                  IF(ELEM(NEIGH(I),K).EQ.ELEM(NEIGH(J),L))THEN
                    WORK(K)=1
                    MATCH=MATCH+1 ! Increment match
                  ENDIF !elem
                ENDDO !l
              ENDDO !k
              IF(MATCH.EQ.2)THEN
                K=1
                DO WHILE(WORK(K).NE.0)
                  K=K+1
                ENDDO !while
                FACE(NEIGH(I),K)=NEIGH(J)
              ENDIF !match
            ENDIF !active
          ENDDO !j
        ENDIF !active
      ENDDO !i

      CALL EXITS('ADD_NODE')
      RETURN
 9999 CALL ERRORS('ADD_NODE',ERROR)
      CALL EXITS('ADD_NODE')
      RETURN 1
      END


