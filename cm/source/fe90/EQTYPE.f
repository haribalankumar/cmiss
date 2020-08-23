      SUBROUTINE EQTYPE(IBT,NBH,NEELEM,nr,NW,nx,ERROR,*)

C#### Subroutine: EQTYPE
C###  Description:
C###    EQTYPE decides which equation (and hence Green's funtion) is
C###    required for each boundary element region.  This subroutine
C###    returns a Green's function identifier for each Boundary element
C###    domain, IGREN.

C**** COMPLEX is true if the Green's function is complex valued.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'eqt000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),nr,NW(NEM,3),nx
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,ne,NE_TYPE,noelem

      CALL ENTERS('EQTYPE',*9999)
      COMPLEX=.FALSE.
      HERMITE=.FALSE.
      noelem=NEELEM(0,nr)
      ne=NEELEM(noelem,nr)! This identifies an element in the BE region.
      !Code below assumes that we are solving the same equation in
      !each element
      IF(ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4) THEN  !static or quasi
        GOTO(10,20,30,40,50)ITYP2(nr,nx) !DECIDING WHAT EQ IS TO BE SOLVED.
C ***           ITYP2(nr,nx)  1 : LINEAR ELASTICITY
C ***                      2 :
C ***                      3 : LAPLACE'S EQUATION
C ***                      4 : HELMHOLTZ EQUATION
C ***                      5 : POISSON'S EQUATION

10      CONTINUE
          GOTO(12,14,16)NJ_LOC(NJL_GEOM,0,nr)
12        CONTINUE
            GOTO 1000
14        CONTINUE
            IF(NW(1,1).EQ.11) THEN !Plane stress element
              IGREN(nr)=14
            ELSEIF(NW(1,1).EQ.12) THEN !Plane strain
              IGREN(nr)=13
            ENDIF
            GOTO 1000
16        CONTINUE
            IGREN(nr)=15
            GOTO 1000

20      CONTINUE
          GOTO(22,24,26)NJ_LOC(NJL_GEOM,0,nr)
22        CONTINUE
            GOTO 1000
24        CONTINUE
            GOTO 1000
26        CONTINUE
            GOTO 1000

30      CONTINUE
          GOTO(32,34,36)NJ_LOC(NJL_GEOM,0,nr)
32        CONTINUE
            GOTO 1000
34        CONTINUE
            IF(ITYP3(nr,nx).EQ.1) THEN !Standard Laplace
              IGREN(nr)=1
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN !div(sigma*grad(phi))
              IGREN(nr)=7
            ENDIF
            GOTO 1000
36        CONTINUE
            IF(ITYP3(nr,nx).EQ.1) THEN !Standard Laplace
              IGREN(nr)=2
            ELSE IF(ITYP3(nr,nx).EQ.2) THEN !div(sigma*grad(phi))
              IGREN(nr)=8
            ENDIF
            GOTO 1000

40      CONTINUE
          GOTO (41,45)ITYP3(nr,nx)
C ***            ITYP3(nr,nx) 1 : STANDARD HELMHOLTZ
C ***                      2 : MODIFIED HELMHOLTZ (YUKAWA)

41        CONTINUE
            GOTO(42,43,44)NJ_LOC(NJL_GEOM,0,nr)
42          CONTINUE
              GOTO 1000
43          CONTINUE
              IF(JTYP4.EQ.2) THEN
                !3d,Cylindrically symmetric about x (or r) axis
    !(i.e 2d overall)
                IGREN(nr)=11
                COMPLEX=.TRUE.
C ***           Get weights and Gauss pts required for GREENC and
C ***           DGREENC
                ng=8
C old MLB 18/3/97 CALL GAUSSPWB(0.0d0,1.0d0,ng,-1,W1,D1,10,ERROR,*9999)
                CALL GAUSSLEG(0.0d0,1.0d0,D1,W1,ng,ERROR,*9999)
              ELSE
                IGREN(nr)=3
              ENDIF
              GOTO 1000
44          CONTINUE
              IF(JTYP4.EQ.1) THEN !3d,Unsymmetric
                IGREN(nr)=4
                COMPLEX=.TRUE.
              ENDIF
              GOTO 1000
45        CONTINUE
            GOTO(46,47,48)NJ_LOC(NJL_GEOM,0,nr)
46          CONTINUE
              GOTO 1000
47          CONTINUE
              IF(JTYP4.EQ.3) THEN
                !3d,Cylindrically symmetric about y (or z) axis
    !(i.e 2d overall)
                IGREN(nr)=12
C ***           Get weights and Gauss pts required for GREENA and
C ***           DGREENA
                ng=8
C old MLB 18/3/97 CALL GAUSSPWB(0.0d0,1.0d0,ng,-1,W1,D1,10,ERROR,*9999)
                CALL GAUSSLEG(0.0d0,1.0d0,D1,W1,ng,ERROR,*9999)
              ELSE
                IGREN(nr)=5
              ENDIF
              GOTO 1000
48          CONTINUE
              IF(JTYP4.EQ.1) THEN !3d,Unsymmetric
                IGREN(nr)=6
              ENDIF
              GOTO 1000
50      CONTINUE
          IF(ITYP3(nr,nx).EQ.2) THEN
            GOTO(52,54,56)NJ_LOC(NJL_GEOM,0,nr)
52            CONTINUE
              GOTO 1000
54          CONTINUE
              IGREN(nr)=9
              GOTO 1000
56          CONTINUE
              IGREN(nr)=10
              GOTO 1000
          ELSE
            ERROR='>>Invalid set of ITYPEs'
            GOTO 9999
          ENDIF
      ENDIF

1000  CONTINUE
      !Define NW for each element
      !AJP 2-3-94
      !NW NEEDS TO BE THE SAME FOR EACH ELEMENT FOR THE BE ROUTINES.
      !THEREFORE DON'T NEED TO DO IT FOR EACH ELEMENT

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        ! This identifies an element in the BE region.
        nb=NBASEF(NBH(NH_LOC(1,nx),1,ne),1) ! First dependent variable
        IF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
C ***   Two-dimensional or 3-d and cylindrically sym i.e 1D integrations
          IF(IBT(1,1,nb).EQ.1) THEN !Lagrange
            IF(IBT(2,1,nb).EQ.0) THEN !1D Constant
              NW(ne,2)=1
            ELSE IF(IBT(2,1,nb).EQ.1) THEN !1D Linear
              NW(ne,2)=2
            ELSE IF(IBT(2,1,nb).EQ.2) THEN !1D Quadratic
              NW(ne,2)=3
            ELSE IF(IBT(2,1,nb).EQ.3) THEN !1D Cubic
              NW(ne,2)=4 !1D Cubic Lagrange
            ENDIF
          ELSE IF(IBT(1,1,nb).EQ.2) THEN !Hermite
            NW(ne,2)=5 !1D Cubic Hermite
            HERMITE=.TRUE.
          ENDIF
        ELSE IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
C ***   Three-dimensional i.e. 2D integrations
          IF(IBT(1,1,nb).EQ.1.AND.
     '       IBT(1,2,nb).EQ.1) THEN
            !Lagrange in first and second directions
            IF(IBT(2,1,nb).EQ.0.AND.
     '         IBT(2,2,nb).EQ.0) THEN
              NW(ne,2)=1 !2D constant
            ELSE IF(IBT(2,1,nb).EQ.1) THEN !Linear in first dir
              IF(IBT(2,2,nb).EQ.1) THEN !Linear in 2nd dir.
                NW(ne,2)=6 !Bilinear
              ELSE IF(IBT(2,2,nb).EQ.2) THEN !Quad in 2nd dir.
                NW(ne,2)=10 !Linear-quadratic
              ELSE IF(IBT(2,2,nb).EQ.3) THEN !Cubic in 2nd dir.
C                  NW(ne,2)= !Linear - Cubic Lagrange
              ENDIF
            ELSE IF(IBT(2,1,nb).EQ.2) THEN !Quad in first dir.
              IF(IBT(2,2,nb).EQ.1) THEN !Linear
C                NW(ne,2)= !Quadratic-linear
              ELSE IF(IBT(2,2,nb).EQ.2) THEN !Quad in second dir.
                NW(ne,2)=7 !Biquadratic
              ELSE IF(IBT(2,2,nb).EQ.3) THEN !Cubic in 2nd dir.
C                NW(ne,2)= !Quadratic - Cubic Lagrange
              ENDIF
            ELSE IF(IBT(2,1,nb).EQ.3) THEN !Cubic in first dir.
              IF(IBT(2,2,nb).EQ.1) THEN !Linear in 2nd dir
C                NW(ne,2)= !Cubic Lagrange - linear
              ELSE IF(IBT(2,2,nb).EQ.2) THEN !Quad in 2nd dir
C                NW(ne,2)= !Cubic Lagrange - Quadratic
              ELSE IF(IBT(2,2,nb).EQ.3) THEN !Cubic in 2nd dir
                NW(ne,2)=8 !Bicubic Lagrange
              ENDIF
            ENDIF !End of bi-Lagrange
          ELSE IF(IBT(1,1,nb).EQ.2.AND.
     '           IBT(1,2,nb).EQ.2) THEN
            !Hermite in each direction
            NW(ne,2)=9 !Bicubic Hermite
            HERMITE=.TRUE.
          ELSE IF(IBT(1,1,nb).EQ.1.AND.
     '           IBT(1,2,nb).EQ.2) THEN
            !Lagrange in first direction and Hermite in the second
            !direction
            IF(IBT(2,1,nb).EQ.1) THEN !Linear in 1st dir
              NW(ne,2)=11 !linear-cubic hermite
              HERMITE=.TRUE.
            ELSE
              ERROR='>>This basis function not implemented'
              GOTO 9999
            ENDIF
          ELSE IF(IBT(1,1,nb).EQ.2.AND.
     '           IBT(1,2,nb).EQ.1) THEN
            !Hermite in first direction and Lagrange in the second
            !direction
            IF(IBT(2,2,nb).EQ.1) THEN !Linear in 2nd dir
              NW(ne,2)=13 !cubic hermite-linear
              HERMITE=.TRUE.
            ELSE
              ERROR='>>This basis function not implemented'
              GOTO 9999
            ENDIF
          ELSE IF(IBT(1,1,nb).EQ.3) THEN !Simplex
            IF(IBT(2,1,nb).EQ.4) THEN !Special hermite simplex
              IF(NKT(1,nb).EQ.1) THEN !Apex at node 1
                NW(ne,2)=14
                HERMITE=.TRUE.
              ELSE IF(NKT(3,nb).EQ.1) THEN !Apex at node 3
                NW(ne,2)=15
                HERMITE=.TRUE.
              ELSE
                ERROR='>>Error in Hermite simplex'
                GOTO 9999
              ENDIF
            ELSE
              ERROR='>>This simplex element not implemented'
            ENDIF
          ELSE IF(IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6.OR.
     '        IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6) THEN !3 Node Sector
            NW(ne,2)=16
            IF(IBT(2,1,nb).EQ.4.OR.IBT(2,2,nb).EQ.4) THEN
              HERMITE=.TRUE.
            ENDIF
          ELSE !A mixture of Hermite and Lagrange
            ERROR='>>This basis function not implemented'
            GOTO 9999
          ENDIF
        ENDIF
      ENDDO

!     !Set up DLIM - the limits of the various integration schemes
      DO nb=1,12
        DO ne_type=1,15
          DLIM(nb,ne_type)=RMAX
        ENDDO
      ENDDO
      DLIM(1,1)=1.0d-2 !Constant elements
      DLIM(2,1)=1.0d0
c cpb 23/6/95 Changing DLIM for linear elements
C      DLIM(1,2)=1.0d-2 !Linear elements
C      DLIM(2,2)=1.0d0
c cpb 6/11/96 adding log quadrature
C      DLIM(3,2)=1.0d-2 !Linear elements
C      DLIM(4,2)=5.0d0 !for the lack of a better number cpb.
      DLIM(5,2)=1.0d-2 !Linear elements
      DLIM(6,2)=5.0d0 !for the lack of a better number cpb.
      DO nb=1,5
        DLIM(nb,3)=1.0d-2 !Quadratic elements
      ENDDO
      DLIM(6,3)=1.0d0
      DO nb=1,3
        DLIM(nb,4)=1.0d-2 !Cubic Lagrange elements
      ENDDO
      DLIM(4,4)=1.0d0
c cpb 23/6/95 Changing DLIM for cubic Hermite elements
C      DLIM(1,5)=1.0d-2 !Cubic hermite elements
C      DLIM(2,5)=1.0d0
c cpb 6/11/96 adding log quadrature
C      DLIM(3,5)=1.0d-2 !Cubic hermite elements
C      DLIM(4,5)=5.0d0 !cpb 23/6/95 For the lack of a better number
      DLIM(5,5)=1.0d-2 !Cubic hermite elements
      DLIM(6,5)=5.0d0 !cpb 23/6/95 For the lack of a better number
      DO nb=1,9
        DLIM(nb,6)=1.0d0-2 !2D bilinear
      ENDDO
      DLIM(10,6)=1.0d0
      DLIM(1,7)=1.0d-2 !2D biquadratic
      DLIM(1,8)=1.0d-2 !2D bicubic Lagrange
      DO nb=1,9
        DLIM(nb,9)=1.0d-2 !2D bicubic hermite
      ENDDO
      DLIM(10,9)=1.0d0
      DLIM(1,10)=1.0d-2 !2D Linear-quadratic
      DO nb=1,9
        DLIM(nb,11)=1.0d-2 !2D linear-cubic hermite
      ENDDO
      DLIM(10,11)=1.0d0
      DLIM(1,12)=1.0d-2 !2D quadratic-cubic hermite
      DO nb=1,9
        DLIM(nb,13)=1.0d-2 !2D cubic hermite-linear
      ENDDO
      DLIM(10,13)=1.0d0
      DO nb=1,7
        DLIM(nb,14)=1.0d-2 !Hermite simplex (Apex at node 1)
      ENDDO
      DLIM(8,14)=1.0d0
      DO nb=1,7
        DLIM(nb,15)=1.0d-2 !Hermite simple (Apex at node 3)
      ENDDO
      DLIM(8,15)=1.0d0
      DO nb=1,7
        DLIM(nb,16)=1.0d-2 !2D 3 node sector element
      ENDDO
      DLIM(8,16)=1.0d0

      CALL EXITS('EQTYPE')
      RETURN
9999  CALL ERRORS('EQTYPE',ERROR)
      CALL EXITS('EQTYPE')
      RETURN 1
      END


