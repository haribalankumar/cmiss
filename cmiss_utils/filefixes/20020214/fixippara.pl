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

    # Created: Leo Cheng 11-MAR-2002
    $file->add_before(
      target => " USE_MAPS      (0 or 1)",
      insert => " USE_MAGNETIC  (0 or 1)[0]: 0\n",
      unless => " USE_MAGNETIC  (0 or 1)"
    );

    # Added: Leo Cheng 11-MAR-2002
    $file->add_before(
      target => " Max# grid degrees of freedom      (NYQM)",
      insert => " Max# signal sets                  (NSSM)[1]:         1\n",
      unless => " Max# signal sets                  (NSSM)[1]:"
    );

    # Added: Carey Stevens 23-June-2002
    $file->add_before(
      target => " Size iter. solver array (NZ_ITERATIVE_M)[1]:",
      insert => " Max# mesh dofs map to 1 mesh dof  (NYYM)[1]:         1\n",
      unless => " Max# mesh dofs map to 1 mesh dof  (NYYM)[1]:"
    );

    # S Norris: Removed the unwanted ITERATIVE parameter values.
    $file->delete( string => qr/(?m)^\ Size\ iter\.\ solver\ array\ \(NZ_ITERATIVE_M\).*\n/ );
    $file->delete( string => qr/(?m)^\ 2nd\ dimen\.\ iter\.\ solver\ \ \(N_ITERATIVE_M\).*\n/ );
    $file->delete( string => qr/(?m)^\ USE_ITERATIVE\ \(0\ or\ 1\)\[1\]\:.*\n/ );
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);









