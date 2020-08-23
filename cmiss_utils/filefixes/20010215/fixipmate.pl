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
    # Created: Carey Stevens 24 May 2001
    #     Adds question for inputting pole-zero law parameters
    my $new_string = <<EOF;
 The shear terms of the pole-zero are [1]:
   (1) Determined using the fibre distribution model
   (2) Input
    1
EOF
    $file->add_before(
      target => " Enter the number of constitutive law parameters",
      unless => "The shear terms of the pole-zero are",
      insert => $new_string
    );

    {
    # Created: Ben Wright 8 Nov 2001
    # Change in file format for specifying element names and #'s individually
    # 	
      my $pattern = qr/The value in element (\d.*) is .*?:\s*(.*)\n/;
      while (my ($elem, $value) = $file->match(pattern => $pattern)) {
        my $old_estring = qr/(?m)^ The value in element.*\n/;
        my $new_estring = <<EOF;
 Enter element #s/name [EXIT]: $elem
 Enter value for element(s) [0.00]: $value
EOF
        $file->replace(target => $old_estring, replace => $new_estring);
      }
      # Place an element=0 for exit crteria, after each block.
      $pattern = qr/( Enter value for element.*)\n\s*\n/;
      while (my ($line) = $file->match(pattern => $pattern)) {
        $file->replace(target => $pattern, replace => <<EOF);
$line
 Enter element #s/name [EXIT]: 0

EOF
      }
    }
    1;
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
