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
    # Added: Merryn Tawhai, 23 Oct 2001
    #   Changes option from 'Humidity in lung' to 'Pulmonary transport'
    $file->replace(
      target => " (11) Humidity transport in lung\n",
      replace => " (11) Pulmonary transport\n"
    );
    # Created: David Nickerson, 07 May 2001
    #   Adds an option for user defined, cellular based, mechanics models
    $file->add_after(
      target => "   (4) Infarct\n",
      insert => "   (5) User defined\n",
      unless => "   (4) Infarct\n   (5) User defined\n"
    );
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
