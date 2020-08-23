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

    # Created: David Hall 3-MAR-2003
    # Question added asking if lattice grid points are to be used.

    $file->add_before(
      target => qr/(?m)^ Enter element #s.*\n/,
      insert => " Use lattice based grid points [N]? N\n",
      unless => qr/Use lattice/
    );
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
