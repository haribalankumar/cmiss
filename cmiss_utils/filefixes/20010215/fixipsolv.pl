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
    # David Nickerson, 07 September 2001
    #   Adds question asking for the convergence criteria to use for the
    #   solution of nonlinear equations
    $file->add_before(
      target => " Do you want to use parallel element stiffness computations",
      insert => <<EOF,
   Specify the convergence criteria [1]:
     (1) Ratio of unconstrained to constrained residuals
     (2) Ratio of unconstrained residuals to maximum Gauss point value
      1
EOF
      unless => "Specify the convergence criteria"
    );
    # Richard Boyes, 3rd October 2001
    #   Add question for LSODA
    my $new;
    if ($file->match(pattern => qr/\(2\) Automatic stepping\n\s*2/)) {
      $new = "  *(6) LSODA (adaptive step, stiff switching)\n";
    } else {
      $new = "   (6) LSODA (adaptive step, stiff switching)\n";
    }
    $file->add_after(
      target => "(5) Adams-Moulton (variable order, adaptive time step)\n",
      insert => $new,
      unless => "(6) LSODA (adaptive step, stiff switching)"
    );
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
