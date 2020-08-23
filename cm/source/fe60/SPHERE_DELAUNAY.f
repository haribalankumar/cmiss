      SUBROUTINE SPHERE_DELAUNAY(VERTICES_XYZ,TRIANGLES,RWORKING,
     '  IERROR,ERROR,*)

C#### Subroutine: SPHERE_DELAUNAY
C###  Description:
C###    Calculates the Delaunay triangulation of the
C###    vertices on a unit sphere.

C***  Created by: David Bullivant, January 2002
C***  Last modified: 10 April 2002

C***  N_VERTICES is the number of vertices to be triangulated
C***  VERTICES_XYZ is an array of length 3*N_VERTICES, containing the
C***    x,y,z coordinates of the vertices (X1,Y1,Z1, X2,Y2,Z2, ...).
C***    After successful completion, the vertices will be projected onto
C***    the unit sphere
C***  MAX_TRIANGLES is the number of triangles that can be held in
C***    TRIANGLES
C***  TRIANGLES.  After successful completion, TRIANGLES will contain
C***    the vertex numbers for the N_TRIANGLES.  It needs to be at least
C***    3*MAX_TRIANGLES long
C***  N_TRIANGLES.  After successful completion, this is the number of
C***    triangles in the triangulation
C***  RWORKING a real working array that is at least 4*MAX_TRIANGLES
C***    long
C***  IERROR not zero indicates an error
C***    IERROR=1 more than MAX_TRIANGLES
C***    IERROR=2 error in SPHERE_PROJECT
C***    IERROR=3 invalid arguments
C***    IERROR=4 could not find initial tetrahedron

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER IERROR,TRIANGLES(3*MAX_TRIANGLES)
      REAL*8 RWORKING(4*MAX_TRIANGLES),VERTICES_XYZ(3*N_VERT_EST)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,INITIAL_TETRAHEDRON(4),ITEMP,J,K,L,M,N,N_ADJACENT_LINES,
     '  N1,N2
      REAL*8 LENGTH,PI,RTEMP,XYZ(3),X1,X2,X3,X4,Y1,Y2,Y3,Y4,Z1,Z2,Z3,Z4
!     Constants
      REAL*8 TOLERANCE
      PARAMETER (TOLERANCE=1.d-8)

      CALL ENTERS('SPHERE_DELAUNAY',*9999)

      PI=4.d0*ATAN(1.d0)
C     normalize vertices
      I=1
      DO J=1,N_VERTICES
        CALL SPHERE_PROJECT(VERTICES_XYZ(I),XYZ,IERROR,ERROR,*9999)
        VERTICES_XYZ(I)=XYZ(1)
        I=I+1
        VERTICES_XYZ(I)=XYZ(2)
        I=I+1
        VERTICES_XYZ(I)=XYZ(3)
        I=I+1
      ENDDO !J
      IERROR=0
      IF (4.LE.N_VERTICES) THEN
        IF (4.LE.MAX_TRIANGLES) THEN
C         find a tetrahedron for the initial triangulation (4 triangles)
C         first corner
          INITIAL_TETRAHEDRON(1)=1
          I=1
C         second corner must be different from first corner and not
C           opposite to first corner
   10     CONTINUE
            I=I+1
            CALL SPHERE_CALC_DIST(VERTICES_XYZ(1),
     '        VERTICES_XYZ(3*I-2),LENGTH,IERROR,ERROR,*9999)
          IF ((0.EQ.IERROR).AND.(I.LT.N_VERTICES-2).AND.
     '      ((LENGTH.LE.TOLERANCE).OR.(LENGTH.GE.PI-TOLERANCE))) GOTO 10
          IF ((0.EQ.IERROR).AND.(TOLERANCE.LT.LENGTH).AND.
     '      (LENGTH.LT.PI-TOLERANCE)) THEN
            INITIAL_TETRAHEDRON(2)=I
C           third corner can't be on the same line as first two corners
            X1=VERTICES_XYZ(1)
            Y1=VERTICES_XYZ(2)
            Z1=VERTICES_XYZ(3)
            J=3*I-2
            X2=VERTICES_XYZ(J)
            Y2=VERTICES_XYZ(J+1)
            Z2=VERTICES_XYZ(J+2)
            X3=Y1*Z2-Z1*Y2
            Y3=Z1*X2-X1*Z2
            Z3=X1*Y2-Y1*X2
   20       CONTINUE
              I=I+1
              J=3*I-2
              X2=VERTICES_XYZ(J)
              Y2=VERTICES_XYZ(J+1)
              Z2=VERTICES_XYZ(J+2)
              X4=Y1*Z2-Z1*Y2
              Y4=Z1*X2-X1*Z2
              Z4=X1*Y2-Y1*X2
              LENGTH=SQRT(X4*X4+Y4*Y4+Z4*Z4)
              RTEMP=SQRT((X3-X4)*(X3-X4)+(Y3-Y4)*(Y3-Y4)+
     '          (Z3-Z4)*(Z3-Z4))
            IF ((I.LT.N_VERTICES-1).AND.
     '        ((LENGTH.LE.TOLERANCE).OR.(RTEMP.LE.TOLERANCE))) GOTO 20
            IF (((LENGTH.GT.TOLERANCE).AND.(RTEMP.GT.TOLERANCE))) THEN
              INITIAL_TETRAHEDRON(3)=I
C             fourth corner can't be in same plane as other 3
              CALL SPHERE_CALC_CIRCUM(
     '          VERTICES_XYZ(3*INITIAL_TETRAHEDRON(1)-2),
     '          VERTICES_XYZ(3*INITIAL_TETRAHEDRON(2)-2),
     '          VERTICES_XYZ(3*INITIAL_TETRAHEDRON(3)-2),XYZ,IERROR,
     '          ERROR,*9999)
              IF (0.EQ.IERROR) THEN
                X1=XYZ(1)
                Y1=XYZ(2)
                Z1=XYZ(3)
   30           CONTINUE
                  I=I+1
C                 make sure that first face is the smaller triangle
                  CALL SPHERE_CALC_CIRCUM(
     '              VERTICES_XYZ(3*INITIAL_TETRAHEDRON(1)-2),
     '              VERTICES_XYZ(3*INITIAL_TETRAHEDRON(2)-2),
     '              VERTICES_XYZ(3*I-2),XYZ,IERROR,ERROR,*9999)
                  IF (0.EQ.IERROR) THEN
                    X2=XYZ(1)
                    Y2=XYZ(2)
                    Z2=XYZ(3)
                    X4=Y1*Z2-Z1*Y2
                    Y4=Z1*X2-X1*Z2
                    Z4=X1*Y2-Y1*X2
                    LENGTH=SQRT(X4*X4+Y4*Y4+Z4*Z4)
                  ENDIF
                IF ((0.EQ.IERROR).AND.(I.LT.N_VERTICES).AND.
     '            (LENGTH.LE.TOLERANCE)) GOTO 30
                IF ((0.EQ.IERROR).AND.(LENGTH.GT.TOLERANCE)) THEN
                  INITIAL_TETRAHEDRON(4)=I
C                 determine which side of the first face the fourth
C                   corner is
                  J=3*I-2
                  K=3*INITIAL_TETRAHEDRON(1)-2
                  X2=VERTICES_XYZ(J)-VERTICES_XYZ(K)
                  Y2=VERTICES_XYZ(J+1)-VERTICES_XYZ(K+1)
                  Z2=VERTICES_XYZ(J+2)-VERTICES_XYZ(K+2)
                  LENGTH=X1*X2+Y1*Y2+Z1*Z2
                  IF (LENGTH.GT.0.D0) THEN
                    J=INITIAL_TETRAHEDRON(2)
                    INITIAL_TETRAHEDRON(2)=INITIAL_TETRAHEDRON(3)
                    INITIAL_TETRAHEDRON(3)=J
                  ENDIF
                ELSE
                  IERROR=4
                ENDIF
              ENDIF
            ELSE
              IERROR=4
            ENDIF
          ELSE
            IF (0.EQ.IERROR) THEN
              IERROR=4
            ENDIF
          ENDIF
          IF (0.EQ.IERROR) THEN
            N_TRIANGLES=4
C           initial triangulation covers the sphere
C           triangle 1
            TRIANGLES(1)=INITIAL_TETRAHEDRON(1)
            TRIANGLES(2)=INITIAL_TETRAHEDRON(2)
            TRIANGLES(3)=INITIAL_TETRAHEDRON(3)
            CALL SPHERE_CALC_CIRCUM(VERTICES_XYZ(3*TRIANGLES(1)-2),
     '        VERTICES_XYZ(3*TRIANGLES(2)-2),
     '        VERTICES_XYZ(3*TRIANGLES(3)-2),RWORKING(1),IERROR,ERROR,
     '        *9999)
            IF (0.EQ.IERROR) THEN
C             Calculate radius of triangle
              CALL SPHERE_CALC_DIST(RWORKING(1),
     '          VERTICES_XYZ(3*TRIANGLES(1)-2),RWORKING(4),IERROR,ERROR,
     '          *9999)
              IF (0.EQ.IERROR) THEN
C               triangle 2
                TRIANGLES(4)=INITIAL_TETRAHEDRON(1)
                TRIANGLES(5)=INITIAL_TETRAHEDRON(4)
                TRIANGLES(6)=INITIAL_TETRAHEDRON(2)
                CALL SPHERE_CALC_CIRCUM(VERTICES_XYZ(3*TRIANGLES(4)-2),
     '            VERTICES_XYZ(3*TRIANGLES(5)-2),
     '            VERTICES_XYZ(3*TRIANGLES(6)-2),RWORKING(5),IERROR,
     '            ERROR,*9999)
                IF (0.EQ.IERROR) THEN
C                 Calculate radius of triangle
                  CALL SPHERE_CALC_DIST(RWORKING(5),
     '              VERTICES_XYZ(3*TRIANGLES(4)-2),RWORKING(8),IERROR,
     '              ERROR,*9999)
                  IF (0.EQ.IERROR) THEN
C                   triangle 3
                    TRIANGLES(7)=INITIAL_TETRAHEDRON(2)
                    TRIANGLES(8)=INITIAL_TETRAHEDRON(4)
                    TRIANGLES(9)=INITIAL_TETRAHEDRON(3)
                    CALL SPHERE_CALC_CIRCUM(
     '                VERTICES_XYZ(3*TRIANGLES(7)-2),
     '                VERTICES_XYZ(3*TRIANGLES(8)-2),
     '                VERTICES_XYZ(3*TRIANGLES(9)-2),RWORKING(9),IERROR,
     '                ERROR,*9999)
                    IF (0.EQ.IERROR) THEN
C                     Calculate radius of triangle
                      CALL SPHERE_CALC_DIST(RWORKING(9),
     '                  VERTICES_XYZ(3*TRIANGLES(7)-2),RWORKING(12),
     '                  IERROR,ERROR,*9999)
                      IF (0.EQ.IERROR) THEN
C                       triangle 4
                        TRIANGLES(10)=INITIAL_TETRAHEDRON(3)
                        TRIANGLES(11)=INITIAL_TETRAHEDRON(4)
                        TRIANGLES(12)=INITIAL_TETRAHEDRON(1)
                        CALL SPHERE_CALC_CIRCUM(
     '                    VERTICES_XYZ(3*TRIANGLES(10)-2),
     '                    VERTICES_XYZ(3*TRIANGLES(11)-2),
     '                    VERTICES_XYZ(3*TRIANGLES(12)-2),RWORKING(13),
     '                    IERROR,ERROR,*9999)
                        IF (0.EQ.IERROR) THEN
C                         Calculate radius of triangle
                          CALL SPHERE_CALC_DIST(RWORKING(13),
     '                      VERTICES_XYZ(3*TRIANGLES(10)-2),
     '                      RWORKING(16),IERROR,ERROR,*9999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            INITIAL_TETRAHEDRON(1)= 3*INITIAL_TETRAHEDRON(1)-2
            INITIAL_TETRAHEDRON(2)= 3*INITIAL_TETRAHEDRON(2)-2
            INITIAL_TETRAHEDRON(3)= 3*INITIAL_TETRAHEDRON(3)-2
            INITIAL_TETRAHEDRON(4)= 3*INITIAL_TETRAHEDRON(4)-2
            I=1
            DO WHILE ((0.EQ.IERROR).AND.(I.LE.3*N_VERTICES))
              IF ((I.NE.INITIAL_TETRAHEDRON(1)).AND.
     '          (I.NE.INITIAL_TETRAHEDRON(2)).AND.
     '          (I.NE.INITIAL_TETRAHEDRON(3)).AND.
     '          (I.NE.INITIAL_TETRAHEDRON(4))) THEN
C               delete the triangles that no longer satisfy the
C                 in-circle criterion
                N_ADJACENT_LINES=0
                J=1
                K=0
                DO WHILE ((0.EQ.IERROR).AND.(J.LE.N_TRIANGLES))
C                 Calculate distance from vertex point to circumcentre
C                   of triangle J
                  CALL SPHERE_CALC_DIST(RWORKING(4*J-3),VERTICES_XYZ(I),
     '              LENGTH,IERROR,ERROR,*9999)
C                 If length.LE.radius, then vertex is within the circle
C                 Triangle will be discarded, replaced by new ones
                  IF (0.EQ.IERROR) THEN
                    IF (LENGTH.LE.RWORKING(4*J)+TOLERANCE) THEN
C                     add to the list of adjacent vertices
C                     use the end of the triangles array as working
C                       storage for the adjacent lines
                      N2=TRIANGLES(3*J)
                      N=1
                      DO WHILE ((0.EQ.IERROR).AND.(N.LE.3))
                        N1=N2
                        N2=TRIANGLES(3*J-3+N)
C                       check for duplicates
                        L=1
                        M=3*MAX_TRIANGLES
                        DO WHILE ((L.LE.N_ADJACENT_LINES).AND.
     '                    .NOT.(((N1.EQ.TRIANGLES(M)).AND.
     '                    (N2.EQ.TRIANGLES(M-1))).OR.
     '                    ((N2.EQ.TRIANGLES(M)).AND.
     '                    (N1.EQ.TRIANGLES(M-1)))))
                          L=L+1
                          M=M-2
                        ENDDO
                        IF (L.LE.N_ADJACENT_LINES) THEN
                          IF (L.LT.N_ADJACENT_LINES) THEN
                            L=3*MAX_TRIANGLES-2*N_ADJACENT_LINES+2
                            TRIANGLES(M)=TRIANGLES(L)
                            TRIANGLES(M-1)=TRIANGLES(L-1)
                          ENDIF
                          N_ADJACENT_LINES=N_ADJACENT_LINES-1
                        ELSE
                          IF (3*N_TRIANGLES+2*(N_ADJACENT_LINES+1).LE.
     '                      3*MAX_TRIANGLES) THEN
C                           add to the list of adjacent vertices
                            L=3*MAX_TRIANGLES-2*N_ADJACENT_LINES
                            TRIANGLES(L)=N2
                            TRIANGLES(L-1)=N1
                            N_ADJACENT_LINES=N_ADJACENT_LINES+1
                          ELSE
                            IERROR=1
                          ENDIF
                        ENDIF
                        N=N+1
                      ENDDO
                      IF (0.EQ.IERROR) THEN
C                       remove triangle
                        IF (J.LT.N_TRIANGLES) THEN
                          L=3*J-2
                          M=3*N_TRIANGLES-2
                          DO N=1,3
                            ITEMP=TRIANGLES(L)
                            TRIANGLES(L)=TRIANGLES(M)
                            TRIANGLES(M)=ITEMP
                            L=L+1
                            M=M+1
                          ENDDO
                          L=4*J-3
                          M=4*N_TRIANGLES-3
                          DO N=1,4
                            RTEMP=RWORKING(L)
                            RWORKING(L)=RWORKING(M)
                            RWORKING(M)=RTEMP
                            L=L+1
                            M=M+1
                          ENDDO
                        ENDIF
                        K=K+1
                        N_TRIANGLES=N_TRIANGLES-1
                      ENDIF
                    ELSE
                      J=J+1
                    ENDIF
                  ENDIF
                ENDDO
                IF (0.EQ.IERROR) THEN
C                 determine new triangles
                  IF (N_TRIANGLES+N_ADJACENT_LINES.LE.MAX_TRIANGLES)
     '              THEN
                    L=3*MAX_TRIANGLES-2*N_ADJACENT_LINES
                    M=3*N_TRIANGLES+1
                    K=1
                    DO WHILE ((0.EQ.IERROR).AND.(K.LE.N_ADJACENT_LINES))
                      TRIANGLES(M)=I/3+1
                      M=M+1
                      L=L+1
                      TRIANGLES(M)=TRIANGLES(L)
                      M=M+1
                      L=L+1
                      TRIANGLES(M)=TRIANGLES(L)
                      M=M+1
                      N_TRIANGLES=N_TRIANGLES+1
                      CALL SPHERE_CALC_CIRCUM(
     '                  VERTICES_XYZ(3*TRIANGLES(3*N_TRIANGLES-2)-2),
     '                  VERTICES_XYZ(3*TRIANGLES(3*N_TRIANGLES-1)-2),
     '                  VERTICES_XYZ(3*TRIANGLES(3*N_TRIANGLES)-2),
     '                  RWORKING(4*N_TRIANGLES-3),IERROR,ERROR,*9999)
                      IF (0.EQ.IERROR) THEN
C                       Calculate radius of triangle
                        CALL SPHERE_CALC_DIST(
     '                    RWORKING(4*N_TRIANGLES-3),VERTICES_XYZ(
     '                    3*TRIANGLES(3*N_TRIANGLES-2)-2),
     '                    RWORKING(4*N_TRIANGLES),IERROR,ERROR,*9999)
                      ENDIF
                      K=K+1
                    ENDDO
                  ELSE
                    IERROR=1
                  ENDIF !N_TRIANGLES
                ENDIF !((0.EQ.IERROR)..)
              ENDIF
              I=I+3
            ENDDO !WHILE((0.EQ.IERROR)...)
          ENDIF
        ELSE
          IERROR=1
        ENDIF !(2.LT.MAX..)
      ELSE
        IF (3.EQ.N_VERTICES) THEN
          IF (2.LE.MAX_TRIANGLES) THEN
C           triangulation covers the sphere
            N_TRIANGLES=2
C           triangle 1
            TRIANGLES(1)=1
            TRIANGLES(2)=2
            TRIANGLES(3)=3
C           triangle 2
            TRIANGLES(4)=1
            TRIANGLES(5)=3
            TRIANGLES(6)=2
          ELSE
            IERROR=1
          ENDIF
        ELSE
          IF (0.LT.N_VERTICES) THEN
            N_TRIANGLES=0
          ELSE
            IERROR=3
          ENDIF
        ENDIF
      ENDIF !(2.LT.NUM...)

      CALL EXITS('SPHERE_DELAUNAY')
      RETURN
 9999 CALL ERRORS('SPHERE_DELAUNAY',ERROR)
      CALL EXITS('SPHERE_DELAUNAY')
      RETURN 1
      END


C      SUBROUTINE SPHERE_DELAUNAY_DYNAM(IERROR,IWORKING,TRIANGLES,
C     '  RWORKING,VERTICES_XYZ,ERROR,*)
C
CC#### Subroutine: SPHERE_DELAUNAY_DYNAM
CC###  Description:
CC###    Calculates the Delaunay triangulation of the
CC###    vertices on a unit sphere.
C
CC***  Created by: David Bullivant, January 2002
CC***  Last modified: 3 February 2002
C
CC***  N_VERTICES is the number of vertices to be triangulated
CC***  VERTICES_XYZ is an array of length 3*N_VERTICES,
CC***    containing the x,y,z coordinates of the vertices (X1,Y1,Z1,
CC***    X2,Y2,Z2, ...). After successful completion, the vertices will
CC***    be projected onto the unit sphere
CC***  MAX_TRIANGLES is the number of triangles that can be
CC***    held in TRIANGLES
CC***  TRIANGLES.  After successful completion, TRIANGLES will contain
CC***    the vertex numbers for the N_TRIANGLES
CC***  N_TRIANGLES.  After successful completion, this is the
CC***    number of triangles in the triangulation
CC***  IWORKING an integer working array that is at least
CC***    N_VERTICES long
CC***  RWORKING a real working array that is at least N_VERTICES
CC***    long
CC***  IERROR not zero indicates an error.  IERROR=1 more than
CC***    MAX_TRIANGLES
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:genmesh.cmn'
C!     Parameter List
C      INTEGER IERROR,IWORKING(10*N_VERTICES),TRIANGLES(MAX_TRIANGLES)
C      REAL*8 RWORKING(10*N_VERTICES),VERTICES_XYZ(3*N_VERTICES)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER I,J,K,L,M,NN,N_ADJACENT_VERTICES,VERTEX_NUMBER
C      LOGICAL ADD_VERTEX
C      REAL*8 ANGLE,LENGTH,RADIUS,XYZ(3),X1,X2,X3,X_AXIS_1,X_AXIS_2,Y1,
C     '  Y2,Y3,Y_AXIS_1,Y_AXIS_2,Z1,Z2,Z3,Z_AXIS_1,Z_AXIS_2
C
C      CALL ENTERS('SPHERE_DELAUNAY_DYNAM',*9999)
C
CC     Initialise working arrays
C      DO I=1,10*N_VERTICES
C        RWORKING(I)=0.d0
C        IWORKING(I)=0
C      ENDDO !I
C
CC     normalize vertices
C      I=1
C      DO J=1,N_VERTICES
C        CALL SPHERE_PROJECT(VERTICES_XYZ(I),XYZ,IERROR,ERROR,*9999)
C        VERTICES_XYZ(I)=XYZ(1)
C        I=I+1
C        VERTICES_XYZ(I)=XYZ(2)
C        I=I+1
C        VERTICES_XYZ(I)=XYZ(3)
C        I=I+1
C      ENDDO !I
C      IERROR=0
C      IF (2.LT.N_VERTICES) THEN
C        IF (2.LT.MAX_TRIANGLES) THEN
C          N_TRIANGLES=2
CC         initial triangulation covers the sphere
CC         triangle 1
C          TRIANGLES(1)=1
C          TRIANGLES(2)=2
C          TRIANGLES(3)=3
CC         triangle 2
C          TRIANGLES(4)=1
C          TRIANGLES(5)=3
C          TRIANGLES(6)=2
C          I=10
C          DO WHILE ((0.EQ.IERROR).AND.(I.LE.3*N_VERTICES))
CC           delete the triangles that no longer satisfy the in-circle
CC             criterion
C            J=1
C            K=0
C            N_ADJACENT_VERTICES=0
C            X1=VERTICES_XYZ(I)
C            Y1=VERTICES_XYZ(I+1)
C            Z1=VERTICES_XYZ(I+2)
CC           calculate a coordinate vectors for the plane perpendicular
CC             to (X1,Y1,Z1)
C            IF (X1.LT.0.) THEN
C              X2= -X1
C            ELSE
C              X2=X1
C            ENDIF
C            IF (Y1.LT.0.) THEN
C              Y2= -Y1
C            ELSE
C              Y2=Y1
C            ENDIF
C            IF (Z1.LT.0.) THEN
C              Z2= -Z1
C            ELSE
C              Z2=Z1
C            ENDIF
C            IF (X2.GT.Y2) THEN
C              IF (X2.GT.Z2) THEN
C                IF (Y2.GT.Z2) THEN
C                  X_AXIS_1=Y1
C                  Y_AXIS_1= -X1
C                  Z_AXIS_1=0
C                ELSE
C                  X_AXIS_1=Z1
C                  Y_AXIS_1=0
C                  Z_AXIS_1= -X1
C                ENDIF
C              ELSE
C                X_AXIS_1= -Z1
C                Y_AXIS_1=0
C                Z_AXIS_1=X1
C              ENDIF !(X2.GT.Z2)
C            ELSE
C              IF (Y2.GT.Z2) THEN
C                IF (X2.GT.Z2) THEN
C                  X_AXIS_1= -Y1
C                  Y_AXIS_1=X1
C                  Z_AXIS_1=0
C                ELSE
C                  X_AXIS_1=0
C                  Y_AXIS_1=Z1
C                  Z_AXIS_1= -Y1
C                ENDIF
C              ELSE
C                X_AXIS_1=0
C                Y_AXIS_1= -Z1
C                Z_AXIS_1=Y1
C              ENDIF
C            ENDIF !(X2.GT.Y2)
C            LENGTH=SQRT(X_AXIS_1*X_AXIS_1+Y_AXIS_1*Y_AXIS_1+
C     '        Z_AXIS_1*Z_AXIS_1)
C            X_AXIS_1=X_AXIS_1/LENGTH
C            Y_AXIS_1=Y_AXIS_1/LENGTH
C            Z_AXIS_1=Z_AXIS_1/LENGTH
C            X_AXIS_2=Y1*Z_AXIS_1-Z1*Y_AXIS_1
C            Y_AXIS_2=Z1*X_AXIS_1-X1*Z_AXIS_1
C            Z_AXIS_2=X1*Y_AXIS_1-Y1*X_AXIS_1
C            DO WHILE ((0.EQ.IERROR).AND.(J.LE.N_TRIANGLES))
CC             Get circumcentre of triangle J, put into XYZ
C              CALL SPHERE_CALC_CIRCUM(VERTICES_XYZ(3*TRIANGLES(3*J-2)
C     '          -2),VERTICES_XYZ(3*TRIANGLES(3*J-1)-2),VERTICES_XYZ(3*
C     '          TRIANGLES(3*J)-2),XYZ,IERROR,ERROR,*9999)
C              IF (0.EQ.IERROR) THEN
CC               Calculate radius of triangle
C                CALL SPHERE_CALC_DIST(XYZ,VERTICES_XYZ(3*
C     '            TRIANGLES(3*J-2)-2),RADIUS,IERROR,ERROR,*9999)
C                IF (0.EQ.IERROR) THEN
CC                 Calculate distance from vertex point to circumcentre
C                  CALL SPHERE_CALC_DIST(XYZ,VERTICES_XYZ(I),
C     '              LENGTH,IERROR,ERROR,*9999)
C                ENDIF
C              ENDIF !(0.EQ.IERROR)
CC             If length.LE.radius, then vertex is within the circle
CC             Triangle will be discarded, replaced by new ones
C              IF ((0.EQ.IERROR).AND.(LENGTH.LE.RADIUS)) THEN
C                K=K+1
CC               add the vertices to the list of adjacent vertices,
CC                 ordering the adjacent vertices so that they go
CC                 anti-clockwise around the new vertex
C                DO L=1,3
C                  VERTEX_NUMBER=TRIANGLES(3*J-3+L)
C                  M=3*VERTEX_NUMBER-2
C                  X2=VERTICES_XYZ(M)
C                  Y2=VERTICES_XYZ(M+1)
C                  Z2=VERTICES_XYZ(M+2)
CC                 cross-product
C                  X3=Y1*Z2-Z1*Y2
C                  Y3=Z1*X2-X1*Z2
C                  Z3=X1*Y2-Y1*X2
C                  ANGLE=ATAN2(X3*X_AXIS_1+Y3*Y_AXIS_1+Z3*Z_AXIS_1,
C     '              X3*X_AXIS_2+Y3*Y_AXIS_2+Z3*Z_AXIS_2)
C                  M=1
C                  DO WHILE ((M.LE.N_ADJACENT_VERTICES).AND.
C     '              (VERTEX_NUMBER.NE.IWORKING(M)).AND.
C     '              (ANGLE.LT.RWORKING(M)))
C                    M=M+1
C                  ENDDO !while((m..))
C                  IF (M.GT.N_ADJACENT_VERTICES) THEN
C                    ADD_VERTEX=.TRUE.
C                  ELSE
C                    IF (VERTEX_NUMBER.NE.IWORKING(M)) THEN
C                      ADD_VERTEX=.TRUE.
C                    ELSE
C                      ADD_VERTEX=.FALSE.
C                    ENDIF
C                  ENDIF !(m.gt.num...)
C                  IF (ADD_VERTEX) THEN
C                    DO NN=N_ADJACENT_VERTICES,M,-1
C                      IWORKING(NN+1)=IWORKING(NN)
C                      RWORKING(NN+1)=RWORKING(NN)
C                    ENDDO !NN
C                    IWORKING(M)=VERTEX_NUMBER
C                    RWORKING(M)=ANGLE
C                    N_ADJACENT_VERTICES=
C     '                N_ADJACENT_VERTICES+1
C                  ENDIF !ADD_VERTEX
C                ENDDO !L
C              ELSE
C                IF (K.GT.0) THEN
C                  L=3*J-2
C                  M=L-3*K
C                  TRIANGLES(M)=TRIANGLES(L)
C                  TRIANGLES(M+1)=TRIANGLES(L+1)
C                  TRIANGLES(M+2)=TRIANGLES(L+2)
C                ENDIF !K.GT.0
C              ENDIF !((0.EQ.IERROR)...
C              J=J+1
C            ENDDO !WHILE((0.EQ.IERROR)...
C            IF ((0.EQ.IERROR).AND.(K.GT.0)) THEN
C              N_TRIANGLES=N_TRIANGLES-K
CC             determine new triangles
C              IF (N_TRIANGLES+N_ADJACENT_VERTICES.LE.
C     '          MAX_TRIANGLES) THEN
C                K=1
C                DO WHILE (K.LT.N_ADJACENT_VERTICES)
C                  L=3*N_TRIANGLES+1
C                  TRIANGLES(L)=I/3+1
C                  TRIANGLES(L+1)=IWORKING(K)
C                  TRIANGLES(L+2)=IWORKING(K+1)
C                  N_TRIANGLES=N_TRIANGLES+1
C                  K=K+1
C                ENDDO !WHILE(K.LT...)
C                L=3*N_TRIANGLES+1
C                TRIANGLES(L)=I/3+1
C                TRIANGLES(L+1)=IWORKING(N_ADJACENT_VERTICES)
C                TRIANGLES(L+2)=IWORKING(1)
C                N_TRIANGLES=N_TRIANGLES+1
C              ELSE
C                IERROR=1
C              ENDIF !N_TRIANGLES
C            ENDIF !((0.EQ.IERROR)..)
C            I=I+3
C          ENDDO !WHILE((0.EQ.IERROR)...)
C        ELSE
C          IERROR=1
C        ENDIF !(2.LT.MAX..)
C      ELSE
C        N_TRIANGLES=0
C      ENDIF !(2.LT.NUM...)
C
C      CALL EXITS('SPHERE_DELAUNAY_DYNAM')
C      RETURN
C 9999 CALL ERRORS('SPHERE_DELAUNAY_DYNAM',ERROR)
C      CALL EXITS('SPHERE_DELAUNAY_DYNAM')
C      RETURN 1
C      END
C
C
