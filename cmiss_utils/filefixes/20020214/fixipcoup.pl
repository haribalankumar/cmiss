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
       # Duane Malcolm, 20 August 2002
       # Updates option (9) from "Unused" to "Coupled tube flow and mechanics"
       # and adds option "(10) Cellular processes" which is required for 
       # coupling advection-diffusion problems.
       my $one = qr/Specify the problem setup  \[1\]:\n/;
       my $two = qr/\(0\) Reset\n/;
       my $three = qr/\(1\) Coupled saturated-unsaturated\n/;
       my $four = qr/\(2\) Coupled Laplace\n/;
       my $five = qr/\(3\) Coupled tree growth problem\n/;

       if( $file->match(pattern=>$one) && $file->match(pattern=>$two) && $file->match(pattern=>$three)  && $file->match(pattern=>$four) && $file->match(pattern=>$five) ){
         print " .. Match found for ipcoup problem setup";
         $file->replace(target => "(9) Unused", replace => "(9) Coupled tube flow and mechanics");
         $file->add_after(target => " (9) Coupled tube flow and mechanics\n", insert => "  (10) Cellular processes\n", unless => "(10) Cellular processes");
       }else{
         print ".. No Match for ADV/DIFF";
       }
    }
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);


