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
    # Created: Carey Stevens
    $file->add_after(
      target => "   (3) Micro-structure parameters\n",
      insert => "   (4) Holmes constitutive law constants\n",
      unless => "   (4) Holmes constitutive law constants"
    );
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
