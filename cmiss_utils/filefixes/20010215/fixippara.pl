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
    # Added: Merryn Tawhai 9 Oct 2001
    $file->add_before(
      target => " Max# global face segments          (NFM)",
      insert => " Max# adjacent elements in Xi      (NEIM)[1]:         4\n",
      unless => " Max# adjacent elements in Xi      (NEIM)"
    );
    # Created: David Nickerson
    $file->add_before(
      target => " USE_MINOS     (0 or 1)",
      insert => " USE_MAPS      (0 or 1)[1]: 0\n",
      unless => " USE_MAPS      (0 or 1)"
    );
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
