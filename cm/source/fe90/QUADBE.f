      SUBROUTINE QUADBE(intscheme,nb,nbbem,ERROR,*)

C#### Subroutine: QUADBE
C###  Description:
C###    QUADBE identifies the appropriate basis function class to use
C###    based on the integration scheme used. Returns local child
C###    number nbbem for the basis families used to interpolate
C###    geometry, dependent and normal derivatives [same child number
C###    for each family].

C**** DLIM stores the values of the bounds on DMIN between which
C**** different quadrature schemes will be used. These are defined in
C**** EQTYPE and are in the BEM000 common block.
C**** The programmer can experiment with these values to try and
C**** obtain optimum efficiency (initialised in EQTYPE).
C**** DLIM currently dimensioned as DLIM(12,15) i.e. a maximum of
C**** 12 different basis functions for any particular element type
C**** and a maximum of 15 different element types.
C**** nbbt=nbasef(nb,0) - the number of basis functions in family nb.
C**** DLIM(i,1),i=1,nbbt are the limits for 1D or 2D constant geometric
C****                    interpolation
C****
C**** DLIM(i,2),i=1,nbbt are the limits for 1D linear geometric inter.
C****
C**** DLIM(i,3),i=1,nbbt are the limits for 1D quadratic geometric
C**** interpolation
C**** (note that for quadratic basis functions three of the values in
C****  DLIM are equal.  This allows the following code to step over
C****  those quadrature schemes for the 0 to 1/2 and 1/2 to 1
C****  integrals - stored as the 2nd and 3rd basis functions.
C****  The same is true for the cubic basis functions.)
C****
C**** DLIM(i,4),i=1,nbbt are the limits for 1D cubic Lagrange geometric
C****                    interpolation
C****
C**** DLIM(i,5),i=1,nbbt are the limits for 1D cubic Hermite geometric
C****                    interpolation
C****
C**** DLIM(i,6),i=1,nbbt are the limits for Bilinear geometric inter
C****
C**** DLIM(i,7),i=1,nbbt are the limits for Biquadratic geometric inter
C****
C**** DLIM(i,8),i=1,nbbt are the limits for Bicubic Lagrange
C****                    geometric interpolation
C****
C**** DLIM(i,9),i=1,nbbt are the limits for Bicubic Hermite geometric
C****                    interpolation
C****
C**** DLIM(i,10),i=1,nbbt are the limits for 2d linear-quadratic
C****                     geometric interpolation
C****
C**** DLIM(i,11),i=1,nbbt are the limits for 2d linear-cubic Hermite
C****                     geometric interpolation
C****
C**** DLIM(i,12),i=1,nbbt are the limits for 2d quadratic-cubic Hermite
C****                     geometric interpolation
C**** DLIM(i,13)   2D cubic Hermite-linear interpolation
C**** DLIM(i,14)   2D cubic Hermite simplex (apex at node 1)
C**** DLIM(i,15)   2D cubic Hermite simplex (apex at node 3)
C**** The user is prompted for 2 quadrature schemes at the problem
C**** setup stage and NBBT basis functions are defined with quadrature
C**** schemes varying between the low and high order. This routine
C**** identifies the appropriate
C**** scheme to use based on the distance from the singularity to the
C**** element under consideration.
C**** The subroutine return the value of nbbem, from which
C**** nbuh, nbqh and nb1j are found.
C**** These are the appropriate basis function numbers for the dependent
C**** variable, the normal deriviative and the geometric coordinate
C**** respectively. The same type of element splitting must be used
C**** for each basis function.
C**** Code assumes basis function numbers are the same for each
C**** dependent variable nh=1,nhp(np)
C**** If adaptive is .TRUE. then an adaptive basis function is set up
C**** based on Telles's rule
C**** Assumes the interpolation is the same in each direction.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER intscheme,nb,nbbem
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('QUADBE',*9999)

      IF(intscheme.EQ.1) THEN !element splitting
        nbbem=1
      ELSE IF(intscheme.EQ.2) THEN !adaptive
        nbbem=NBASEF(nb,0)
      ELSE IF(intscheme.GE.3) THEN !high/medium/low gauss scheme
        nbbem=NBASEF(nb,0)-6+intscheme
      ENDIF

      CALL EXITS('QUADBE')
      RETURN
9999  CALL ERRORS('QUADBE',ERROR)
      CALL EXITS('QUADBE')
      RETURN 1
      END


