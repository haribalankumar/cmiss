eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

use strict;
BEGIN {
  die "Environment variable CMISS_ROOT not defined\n"
	unless exists $ENV{CMISS_ROOT};
}
use lib "$ENV{CMISS_ROOT}/cmiss_perl/lib";
use CmUtils::IpFileFix qw(fixfiles);

fixfiles(
  filter => sub {
    my $file = shift;
    my $krylov = 0;


    # Reorder solvers, and add solver names:
    # First fix broken existing examples ...
    $file->replace( target =>  "*(2) SVD Factorisation",
                    replace => " (2) SVD Factorisation" );

    $file->replace(
      target => <<EOF
   (4) Generalised minimum residual
   (5) Unused
EOF
      , replace => <<EOF
   (4) Generalised minimum residual
   (5) Least Squares
EOF
    );
    $file->replace(
      target => <<EOF
   (5) Least Squares
   (6) Unused
EOF
      , replace => <<EOF
   (5) Least Squares
   (6) Conjugate gradient
EOF
    );
    $file->replace(
      target => <<EOF
   (6) Conjugate gradient
   (7) Unused
EOF
      , replace => <<EOF
   (6) Conjugate gradient
   (7) Cholesky Decomposition
EOF
    );

    $file->replace(
      target => <<EOF
 Specify type of linear solution procedure [1]:
   (1) LU Factorisation
   (2) SVD Factorisation
   (3) Biconjugate gradient
   (4) Generalised minimum residual
   (5) Least Squares
   (6) Conjugate gradient
   (7) Cholesky Decomposition
EOF
      , replace => <<EOF
 Specify type of linear solution procedure [1]:
   (1)  LU Decomposition
   (2)  Single Value Decomposition
   (3)  Least Squares
   (4)  Cholesky Decomposition
   (5)  Jacobi Iteration
   (6)  Succesive Over Relaxation
   (7)  Incomplete LU Decomposition(0)
   (8)  Incomplete LU Decomposition(1)
   (9)  Conjugate Gradient
   (10) Biconjugate Gradient Stabilised
   (11) Generalised Minimum Residual
EOF
    );

    # Get the solver type -- they have been renumbers, so we have to 
    # kludge a wee bit
    $file->replace( target =>  "   (11) Generalised Minimum Residual\n    ",
                    replace => "   (11) Generalised Minimum Residual\n    SOLNAME" );

    # Least Squares solver
    if ($file->match( pattern => "SOLNAME5" )) {
       $file->replace( target =>  "SOLNAME5",
                       replace => "3" );
    }
    # Cholesky
    elsif ($file->match( pattern => "SOLNAME7" )) {
      $file->replace( target =>  "SOLNAME7",
                      replace => "4" );
    }
    # Krylov space methods
    elsif ($file->match( pattern => "SOLNAME6" )) {
      $file->replace( target =>  "SOLNAME6",
                      replace => "9" );
    }
    elsif ($file->match( pattern => "SOLNAME3" )) {
      $file->replace( target =>  " SOLNAME3",
                      replace => "10" );
    }
    elsif ($file->match( pattern => "SOLNAME4" )) {
      $file->replace( target =>  " SOLNAME4",
                      replace => "11" );
    }
    else {
      $file->replace( target =>  "SOLNAME",
                      replace => "" );
    }

    $file->replace(
      target => <<EOF
   (3) Incomplete LU
   (4) Row scale and swap
EOF
    , replace => <<EOF
   (3) Symmetric succesive relaxation
   (4) Row scale and swap
EOF
    );
    $file->replace(
      target => <<EOF
   (4) Row scale and swap
   (5) Unused
EOF
    , replace => <<EOF
   (4) Row scale and swap
  *(5) Incomplete LU
EOF
    );
    $file->replace(
      target => <<EOF
   (1) None
   (2) Diagonal
   (3) Symmetric succesive relaxation
   (4) Row scale and swap
  *(5) Incomplete LU
EOF
      , replace => <<EOF
   (1) None
   (2) Point Jacobi
   (3) Jacobi Iteration
   (4) Symmetric Succesive Over Relaxation
   (5) Incomplete LU Decomposition(0)
   (6) Incomplete LU Decomposition(1)
   (7) Row scale
EOF
    );
		$file->replace( target => "Row scale and swap", replace => "Row scale" );

    if( $file->match( pattern => qr/(?m)^ Enter the tolerance/ ) ) {
      $file->replace( target => qr/(?m)^ The number of orthogonalisations before restarting.*\n/,
		      replace => "" );
      $file->replace( target => qr/(?m)^ Enter the tolerance.*$/,
		      replace => " Enter the solver tolerance [1.0d-6]: 1.0d-6\n The number of orthogonalisations before restarting is [10]: 10");
    }

    if ($file->match( pattern => " Specify type of preconditioning" )) {
      $file->add_before( target => " Specify option for linear solution [0]:",
                         insert => " The number of preconditioner iterations per loop [2]: 2\n",
                         unless => " The number of preconditioner iterations per loop" );
    }

    # Get rid of the Harwell sparse LU solver, and add options
    # for SuperLU and Umfpack:
    if ($file->match( pattern => "Specify the LU solver [1]:")) {

      $file->replace( target => <<EOF
 Specify the LU solver [1]:
   (1) Harwell MA28
   (2) SuperLU
   (3) Umfpack
EOF
        , replace => <<EOF
 Specify the LU solver [2]:
   (1) SuperLU
   (2) Umfpack
EOF
      );
      $file->replace( target => <<EOF
 Specify the LU solver [1]:
   (1) Harwell MA28
   (2) SuperLU
EOF
        , replace => <<EOF
 Specify the LU solver [2]:
   (1) SuperLU
   (2) Umfpack
EOF
      );

      $file->replace( target =>  "(2) Umfpack\n    ",
                      replace => "(2) Umfpack\n    LUSOLVER" );
      # If we have the SuperLU solver
      if ($file->match( pattern => "LUSOLVER2")) {
        $file->replace( target =>  "LUSOLVER2",
                        replace => "1" );
        $file->replace( target => <<EOF
   (1) Natural ordering
   (2) Minimum degree ordering from A^T.A
   (3) Minimum degree ordering from A+A^T
EOF
         , replace => <<EOF
   (1) Natural ordering
   (2) Minimum degree ordering from A^T.A
   (3) Minimum degree ordering from A+A^T
   (4) Column Minimum degree ordering
EOF
        );
      }
      # Otherwise we have Umfpack
      else {
        $file->replace( target =>  "LUSOLVER1", replace => "LUSOLVER3" );
        $file->replace( target =>  "LUSOLVER3",
                        replace => "2\n Specify the pivot threshold [0.1D0]: 0.1" );
      }
    }

    # Add an error check to the solution process
    $file->replace(
      target =>  "   (1) Timing only",
      replace => "   (1) Timing and error check"
    );

  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
