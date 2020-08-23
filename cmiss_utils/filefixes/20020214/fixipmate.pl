eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

use strict;
use lib "/product/cmiss/cmiss_perl/lib";
use CmUtils::IpFileFix qw(fixfiles);

fixfiles(
  filter => sub {
    my $file = shift;

    # D Hall: Extend the exponential form question for skin funcitons.
    #         Should only be used if membrane theory is used.
    $file->add_after( 
      target => "   (2) Dr J.W. Holmes distributed fibre formulation\n",
      insert => "   (3) Tong & Fung skin function\n",
      unless => "   (3) Tong & Fung skin function"
    );


    if( $file->match( pattern => " Specify whether the strain energy W" )) {
      $file->{text} .= "\n"; 
      $file->{text} .= " Specify whether the density is [1]:\n"; 
      $file->{text} .= "  (1) Constant spatially\n"; 
      $file->{text} .= "  (2) Piecewise constant (defined by elements)\n"; 
      $file->{text} .= "  (3) Piecewise linear   (defined by nodes)\n"; 
      $file->{text} .= "  (4) Defined by Gauss points\n"; 
      $file->{text} .= "  (5) Defined by Grid points\n"; 
      $file->{text} .= "    1\n"; 
      $file->{text} .= " The value is [0.000D+00]: 0\n"; 
    }
  }
)

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
