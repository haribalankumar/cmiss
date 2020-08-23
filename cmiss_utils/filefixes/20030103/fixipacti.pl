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

    {
      # EWR: Added regional variation of calcium level parameter
      $file->add_before(
        target => " Enter initial calcium level",
        insert => " Specify whether the initial calcium level [Ca]i is [1]:\n  (1) Constant spatially\n  (2) Piecewise constant (defined by elements)\n  (3) Defined by Gauss points\n   1\n",
        unless => " Specify whether the initial calcium level [Ca]i is [1]:"
      );
    }
    1;
  }
)

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
