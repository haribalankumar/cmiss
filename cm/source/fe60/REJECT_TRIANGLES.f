      SUBROUTINE REJECT_TRIANGLES(BOUNDARY2,NPI,TRIANGLES,ONBOUND,
     '  ERROR,*)

C#### Subroutine: REJECT_TRIANGLES
C###  Description:
C###    REJECT_TRIANGLES removes boundary triangles (those involving
C###    3 B nodes), and replaces triangles that involve 1 or 2 B nodes
C###    with a structure that is appropriate for defining the boundary.

C***  The B and IB nodes define the object (in this case a sphere)
C***  boundary by spanning the boundary, lying on the object surface,
C***  and being orthogonal to one another.  To exploit the B and IB
C***  node boundary definition, the final triangulation requires B
C***  nodes to only be in triangles with B-B-IB or B-IB-IB structures.
C***  See diagram below.

C***             IN
C***             /\
C***            /  \
C***           /    \
C***          /      \
C***         /        \
C***        /          \
C***     IB ------------ IB
C***       |\           |
C***       |  \         |
C***       |    \       |
C***       |      \     |
C***       |        \   |
C***       |          \ |
C***        ------------
C***       B            B

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER BOUNDARY2(N_BOUNDARY*2,2),NPI(N_VERT_EST),
     '  TRIANGLES(3*MAX_TRIANGLES)
      LOGICAL ONBOUND(MAX_TRIANGLES)
      CHARACTER ERROR*(*)
      INTEGER BNODE(2),i,IBNODE(2),index,INODE,j,k,nbound,NFOUND,
     '  NIBOUND,nn,np,ntri,ntri2,TRI_A(3),TRI_ADJ(2),TRI_B(3)
      LOGICAL IB,trib(3),NFOUND1,NFOUND2

      CALL ENTERS('REJECT_TRIANGLES',*9999)

      ntri=0
      DO WHILE(ntri.LT.N_TRIANGLES) !replace boundary triangles
        ntri=ntri+1
        TRI_A(1)=NPI(TRIANGLES(3*ntri-2))
        TRI_A(2)=NPI(TRIANGLES(3*ntri-1))
        TRI_A(3)=NPI(TRIANGLES(3*ntri))

        ONBOUND(ntri)=.FALSE.
        NBOUND=0
        NIBOUND=0
        DO j=1,3 !for each triangle vertex
          TRIB(J)=.FALSE.
          DO i=1,N_BOUND !for each boundary node
            IF(TRI_A(j).EQ.BOUNDARY2(I,2))THEN
              NBOUND=NBOUND+1
              TRIB(J)=.TRUE.
            ENDIF
          ENDDO
          DO i=1,N_IBOUND !for each internal boundary node
            IF(TRI_A(j).EQ.BOUNDARY2(I,1))THEN
              NIBOUND=NIBOUND+1
            ENDIF
          ENDDO !i
        ENDDO !j

        IF(NBOUND.EQ.2)THEN
          ONBOUND(ntri)=.TRUE.
C       Check whether other node is IN or IB
          IB=.FALSE.
          DO j=1,3
            DO i=1,N_IBOUND
              IF(TRI_A(j).EQ.BOUNDARY2(i,1)) IB=.TRUE.
            ENDDO !i
          ENDDO !j
C         If not IB node, then replace triangle with composite triangle
          IF(.NOT.IB)THEN
            i=0
            DO j=1,3 !sort out which nodes are B and IN
              IF(TRIB(j))THEN
                i=i+1
                BNODE(i)=TRI_A(j) !store boundary nodes
              ELSE
                INODE=TRI_A(j)
              ENDIF
            ENDDO !j
 !find the adjacent triangles
            NFOUND1=.FALSE.
            NFOUND2=.FALSE.
            NFOUND=0
            ntri2=0
            IBNODE(1)=0
            IBNODE(2)=0
            DO WHILE((ntri2.LT.N_TRIANGLES).AND.((.NOT.NFOUND1).OR.
     '        (.NOT.NFOUND2)))
C Probably a better way to structure the following, more compact.
              ntri2=ntri2+1
              IF(ntri2.NE.ntri) THEN !can't use current triangle
                TRI_B(1)=NPI(TRIANGLES(3*ntri2-2))
                TRI_B(2)=NPI(TRIANGLES(3*ntri2-1))
                TRI_B(3)=NPI(TRIANGLES(3*ntri2))
                IF(TRI_B(1).EQ.BNODE(1))THEN
                  IF(TRI_B(2).EQ.INODE)THEN
                    IF(.NOT.NFOUND1) THEN
                      NFOUND1=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(1)=TRI_B(3) !other (IB) node
                    ENDIF
                  ELSE IF(TRI_B(3).EQ.INODE)THEN
                    IF(.NOT.NFOUND1) THEN
                      NFOUND1=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(1)=TRI_B(2) !other (IB) node
                    ENDIF
                  ENDIF
                ELSE IF(TRI_B(2).EQ.BNODE(1))THEN
                  IF(TRI_B(1).EQ.INODE)THEN
                    IF(.NOT.NFOUND1) THEN
                      NFOUND1=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(1)=TRI_B(3) !other (IB) node
                    ENDIF
                  ELSE IF(TRI_B(3).EQ.INODE)THEN
                    IF(.NOT.NFOUND1) THEN
                      NFOUND1=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(1)=TRI_B(1) !other (IB) node
                    ENDIF
                  ENDIF
                ELSE IF(TRI_B(3).EQ.BNODE(1))THEN
                  IF(TRI_B(1).EQ.INODE)THEN
                    IF(.NOT.NFOUND1) THEN
                      NFOUND1=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(1)=TRI_B(2) !other (IB) node
                    ENDIF
                  ELSE IF(TRI_B(2).EQ.INODE)THEN
                    IF(.NOT.NFOUND1) THEN
                      NFOUND1=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(1)=TRI_B(1) !other (IB) node
                    ENDIF
                  ENDIF
                ELSE IF(TRI_B(1).EQ.BNODE(2))THEN
                  IF(TRI_B(2).EQ.INODE)THEN
                    IF(.NOT.NFOUND2) THEN
                      NFOUND2=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(2)=TRI_B(3) !other (IB) node
                    ENDIF
                  ELSE IF(TRI_B(3).EQ.INODE)THEN
                    IF(.NOT.NFOUND2) THEN
                      NFOUND2=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(2)=TRI_B(2) !other (IB) node
                    ENDIF
                  ENDIF
                ELSE IF(TRI_B(2).EQ.BNODE(2))THEN
                  IF(TRI_B(1).EQ.INODE)THEN
                    IF(.NOT.NFOUND2) THEN
                      NFOUND2=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(2)=TRI_B(3) !other (IB) node
                    ENDIF
                  ELSE IF(TRI_B(3).EQ.INODE)THEN
                    IF(.NOT.NFOUND2) THEN
                      NFOUND2=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(2)=TRI_B(1) !other (IB) node
                    ENDIF
                  ENDIF
                ELSE IF(TRI_B(3).EQ.BNODE(2))THEN
                  IF(TRI_B(1).EQ.INODE)THEN
                    IF(.NOT.NFOUND2) THEN
                      NFOUND2=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(2)=TRI_B(2) !other (IB) node
                    ENDIF
                  ELSE IF(TRI_B(2).EQ.INODE)THEN
                    IF(.NOT.NFOUND2) THEN
                      NFOUND2=.TRUE.
                      NFOUND=NFOUND+1
                      TRI_ADJ(NFOUND)=ntri2
                      IBNODE(2)=TRI_B(1) !other (IB) node
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF !ntri2.NE.ntri
            ENDDO !WHILE

C... replace triangles with composites
            DO nn=1,2
              k=0
              np=0
              DO WHILE(np.NE.IBNODE(nn))
                k=k+1
                np=NPI(k)
              ENDDO !WHILE
C           Replace 2B + IN triangle
              TRIANGLES(3*ntri-(3-nn))=k
C           Replace second adjacent triangle
              TRIANGLES(3*TRI_ADJ(2)-(3-nn))=k
            ENDDO !nn
            ONBOUND(ntri)=.TRUE.
            ONBOUND(TRI_ADJ(2))=.TRUE.
            DO nn=1,2
              k=0
              np=0
              DO WHILE(np.NE.BNODE(nn))
                k=k+1
                np=NPI(k)
              ENDDO !WHILE
C           Replace first adjacent triangle
              TRIANGLES(3*TRI_ADJ(1)-(3-nn))=k
            ENDDO !nn
            ONBOUND(TRI_ADJ(1))=.TRUE.
            k=0
            np=0
            DO WHILE(np.NE.INODE)
              k=k+1
              np=NPI(k)
            ENDDO !WHILE
            TRIANGLES(3*ntri)=k
            k=0
            np=0
            DO WHILE(np.NE.IBNODE(1))
              k=k+1
              np=NPI(k)
            ENDDO !WHILE
            TRIANGLES(3*TRI_ADJ(1))=k
            k=0
            np=0
            DO WHILE(np.NE.BNODE(2))
              k=k+1
              np=NPI(k)
            ENDDO !WHILE
            TRIANGLES(3*TRI_ADJ(2))=k
C.. old stuff below
C           Replace 2B + IN triangle
C            TRIANGLES(3*ntri-2)=IBNODE(1)
C            TRIANGLES(3*ntri-1)=IBNODE(2)
C            TRIANGLES(3*ntri)=INODE
C            ONBOUND(ntri)=.TRUE. !note, this needed only sometimes
C           Replace first adjacent triangle
C            ntri2=TRI_ADJ(1)
C            TRIANGLES(3*ntri2-2)=BNODE(1)
C            TRIANGLES(3*ntri2-1)=BNODE(2)
C            TRIANGLES(3*ntri2)=IBNODE(1)
C            ONBOUND(ntri2)=.TRUE.
C           Replace second adjacent triangle
C            ntri2=TRI_ADJ(2)
C            TRIANGLES(3*ntri2-2)=IBNODE(1)
C            TRIANGLES(3*ntri2-1)=IBNODE(2)
C            TRIANGLES(3*ntri2)=BNODE(2)
C            ONBOUND(ntri2)=.TRUE.
            ntri=ntri-1 !forces the triangle to be re-checked
          ENDIF
        ENDIF !NBOUND
      ENDDO !WHILE

      ntri=0
      DO WHILE(ntri.LT.N_TRIANGLES) !replace boundary triangles
        ntri=ntri+1
        TRI_A(1)=NPI(TRIANGLES(3*ntri-2))
        TRI_A(2)=NPI(TRIANGLES(3*ntri-1))
        TRI_A(3)=NPI(TRIANGLES(3*ntri))

        NBOUND=0
        NIBOUND=0
        DO j=1,3 !for each triangle vertex
          TRIB(J)=.FALSE.
          DO i=1,N_BOUND !for each boundary node
            IF(TRI_A(j).EQ.BOUNDARY2(I,2))THEN
              NBOUND=NBOUND+1
              TRIB(J)=.TRUE.
            ENDIF
          ENDDO
          DO i=1,N_IBOUND !for each internal boundary node
            IF(TRI_A(j).EQ.BOUNDARY2(I,1))THEN
              NIBOUND=NIBOUND+1
            ENDIF
          ENDDO !i
        ENDDO !j

        IF(NBOUND.EQ.1.AND.NIBOUND.EQ.1)THEN
          index=0
          DO j=1,3
            IF(TRIB(j))THEN
              BNODE(1)=TRI_A(j) !boundary node
              index=index+j
            ENDIF
            DO i=1,N_IBOUND
              IF(TRI_A(j).EQ.BOUNDARY2(i,1))THEN
                IBNODE(1)=TRI_A(j)
                index=index+j
              ENDIF
            ENDDO
          ENDDO !j
          INODE=TRI_A(6-index)
 !Find triangle adjacent to B-IN edge
          ntri2=0
          NFOUND=0
          DO WHILE(NFOUND.LT.2.AND.ntri2.LT.N_TRIANGLES)
            ntri2=ntri2+1
            IF(ntri2.NE.ntri)THEN
              TRI_B(1)=NPI(TRIANGLES(3*ntri2-2))
              TRI_B(2)=NPI(TRIANGLES(3*ntri2-1))
              TRI_B(3)=NPI(TRIANGLES(3*ntri2))
              IF(TRI_B(1).EQ.BNODE(1))THEN
                IF(TRI_B(2).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  IBNODE(2)=TRI_B(3) !other (IB) node
                ELSE IF(TRI_B(3).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  IBNODE(2)=TRI_B(2) !other (IB) node
                ENDIF
              ELSE IF(TRI_B(2).EQ.BNODE(1))THEN
                IF(TRI_B(1).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  IBNODE(2)=TRI_B(3) !other (IB) node
                ELSE IF(TRI_B(3).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  IBNODE(2)=TRI_B(1) !other (IB) node
                ENDIF
              ELSE IF(TRI_B(3).EQ.BNODE(1))THEN
                IF(TRI_B(1).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  IBNODE(2)=TRI_B(2) !other (IB) node
                ELSE IF(TRI_B(2).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  IBNODE(2)=TRI_B(1) !other (IB) node
                ENDIF
              ELSE IF(TRI_B(1).EQ.IBNODE(1))THEN
                IF(TRI_B(2).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  BNODE(2)=TRI_B(3) !other (IB) node
                ELSE IF(TRI_B(3).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  BNODE(2)=TRI_B(2) !other (IB) node
                ENDIF
              ELSE IF(TRI_B(2).EQ.IBNODE(1))THEN
                IF(TRI_B(1).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  BNODE(2)=TRI_B(3) !other (IB) node
                ELSE IF(TRI_B(3).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  BNODE(2)=TRI_B(1) !other (IB) node
                ENDIF
              ELSE IF(TRI_B(3).EQ.IBNODE(1))THEN
                IF(TRI_B(1).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  BNODE(2)=TRI_B(2) !other (IB) node
                ELSE IF(TRI_B(2).EQ.INODE)THEN
                  NFOUND=NFOUND+1
                  TRI_ADJ(NFOUND)=ntri2
                  BNODE(2)=TRI_B(1) !other (IB) node
                ENDIF
              ENDIF
            ENDIF
          ENDDO !WHILE
C... fixed so correct ref # used in TRIANGLES array, KSB June 2002
          DO nn=1,2
            k=0
            np=0
            DO WHILE(np.NE.IBNODE(nn))
              k=k+1
              np=NPI(k)
            ENDDO !WHILE
C           Replace B + IN + IB triangle
            TRIANGLES(3*ntri-(3-nn))=k
C           Replace first adjacent triangle
            TRIANGLES(3*TRI_ADJ(1)-(3-nn))=k
          ENDDO !nn
          ONBOUND(ntri)=.TRUE.
          ONBOUND(TRI_ADJ(1))=.TRUE.
          k=0
          np=0
          DO WHILE(np.NE.IBNODE(1))
            k=k+1
            np=NPI(k)
          ENDDO !WHILE
C           Replace second adjacent triangle
          TRIANGLES(3*TRI_ADJ(2)-2)=k
          k=0
          np=0
          DO WHILE(np.NE.INODE)
            k=k+1
            np=NPI(k)
          ENDDO !WHILE
          TRIANGLES(3*TRI_ADJ(2)-1)=k
          ONBOUND(TRI_ADJ(2))=.TRUE.
          k=0
          np=0
          DO WHILE(np.NE.INODE)
            k=k+1
            np=NPI(k)
          ENDDO !WHILE
          TRIANGLES(3*ntri)=k
          k=0
          np=0
          DO WHILE(np.NE.BNODE(1))
            k=k+1
            np=NPI(k)
          ENDDO !WHILE
          TRIANGLES(3*TRI_ADJ(1))=k
          k=0
          np=0
          DO WHILE(np.NE.BNODE(2))
            k=k+1
            np=NPI(k)
          ENDDO !WHILE
          TRIANGLES(3*TRI_ADJ(2))=k
C           Replace B + IN + IB triangle
C          TRIANGLES(3*ntri-2)=IBNODE(1)
C          TRIANGLES(3*ntri-1)=IBNODE(2)
C          TRIANGLES(3*ntri)=INODE
C          ONBOUND(ntri)=.TRUE.
C           Replace first adjacent triangle
C          ntri2=TRI_ADJ(1)
C          TRIANGLES(3*ntri2-2)=IBNODE(1)
C          TRIANGLES(3*ntri2-1)=IBNODE(2)
C          TRIANGLES(3*ntri2)=BNODE(1)
C          ONBOUND(ntri2)=.TRUE.
C           Replace second adjacent triangle
C          ntri2=TRI_ADJ(2)
C          TRIANGLES(3*ntri2-2)=IBNODE(1)
C          TRIANGLES(3*ntri2-1)=INODE
C          TRIANGLES(3*ntri2)=BNODE(2)
C          ONBOUND(ntri2)=.TRUE.
          ntri=ntri-1 !forces the triangle to be re-checked

 !Find triangle adjacent to B-IB edge
        ENDIF
      ENDDO !WHILE

      ntri=0
      DO WHILE(ntri.LT.N_TRIANGLES) !delete 3 B noded triangles
        ntri=ntri+1
        TRI_A(1)=NPI(TRIANGLES(3*ntri-2))
        TRI_A(2)=NPI(TRIANGLES(3*ntri-1))
        TRI_A(3)=NPI(TRIANGLES(3*ntri))

        NBOUND=0
        NIBOUND=0
        DO j=1,3 !for each triangle vertex
          TRIB(J)=.FALSE.
          DO i=1,N_BOUND !for each boundary node
            IF(TRI_A(j).EQ.BOUNDARY2(I,2))THEN
              NBOUND=NBOUND+1
              TRIB(J)=.TRUE.
            ENDIF
          ENDDO
          DO i=1,N_IBOUND
            IF(TRI_A(j).EQ.BOUNDARY2(I,1))THEN !IB node
              NIBOUND=NIBOUND+1
            ENDIF
          ENDDO !i
        ENDDO !j

        IF(NIBOUND.EQ.2.AND.NBOUND.EQ.1) THEN
          ONBOUND(ntri)=.TRUE.
        ELSE IF(NBOUND.EQ.3)THEN !delete triangle
          DO ntri2=ntri+1,N_TRIANGLES
            TRIANGLES(3*(ntri2-1)-2)=TRIANGLES(3*ntri2-2)
            TRIANGLES(3*(ntri2-1)-1)=TRIANGLES(3*ntri2-1)
            TRIANGLES(3*(ntri2-1))=TRIANGLES(3*ntri2)
            ONBOUND(ntri2-1)=ONBOUND(ntri2)
          ENDDO !ntri2
          N_TRIANGLES=N_TRIANGLES-1
          ntri=ntri-1
        ENDIF
      ENDDO !WHILE

      CALL EXITS('REJECT_TRIANGLES')
      RETURN
 9999 CALL ERRORS('REJECT_TRIANGLES',ERROR)
      CALL EXITS('REJECT_TRIANGLES')
      RETURN 1
      END


