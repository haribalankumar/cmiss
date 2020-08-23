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

    #
    # Created: Travis Austin 31-May-2005
    #   Adds option to solver list of using Algebraic Multigrid
    #
 
    {
        $file->add_after
        (
          target =>  qr/(?m)^   \(11\) Generalised Minimum Residual/,
          insert =>  "\n   (14) Algebraic Multigrid (AMG1R6)",
          unless =>  qr/(?m)^   \(14\) Algebraic Multigrid \(AMG1R6\)/
        );
    }
  }
)

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
