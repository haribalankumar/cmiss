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
      # Justin Fernandez, 08 March 2002
      #   Adds a question for contact mechanics if the problem
      #   involves 3D static finite elasticity with Galerkin FE.
      #   Will only update once.
  
      my $one = qr/Specify whether \[1\]:\n(\s*\(\d\)(\s|\w|\d|\-)*\n){6}\s*1/;
      my $two = qr/Specify equation \[1\]:\n(\s*\(\d\)(\s|\w|\d|\-|\(|\)|\.|=)*\n){9}\s*2/;
      my $three = qr/Specify problem type \[1\]:\n([\s\*]*\(\d\)(\s|\w|\d|\-)*\n){6}\s*3/;
      my $four = qr/Specify whether solution is by \[1\]:\n(\s*\(\d\)(\s|\w|\d|\-)*\n){6}\s*1/;
      my $new = " Does this problem involve contact mechanics [N]? N\n";
  
      if($file->match(pattern =>$new)) {
        print "Already Modified";
      }else{
        if(($file->match(pattern=>$one))&&($file->match(pattern=>$two))&&($file->match(pattern=>$three))&&($file->match(pattern=>$four))){
         $file->add_after(
            target => qr/Do you want the global matrices stored as sparse matrices \[Y\]\? [Nn]\n/,
            insert => $new,
            unless => qr/Do you want the global matrices stored as sparse matrices \[Y\]\? [Yy]\n/
          );
          $file->add_after(
            target => qr/Do you want to calculate the sparsity pattern for the global matrices \[Y\]\? [Yy]\n/,
            insert => $new,
            unless => qr//
          );     
        }else{
          print " .. No Match for Contact";
        }
      }
    }
     {
       # Duane Malcolm 20/08/02
       my $one = qr/Specify whether \[1\]:\n(\s*\(\d+\).*\n){6}\s*2\s*\n/;
       my $two = qr/Specify equation \[1\]:\n(\s*\*?\(\d+\).*\n){12}\s*3\s*\n/;
 
       if( $file->match(pattern=>$one) && $file->match(pattern=>$two) ){
         print " .. Match found for ADV/DIFF";
         $file->replace(
           target => " Specify the number of dependent variables (0 for cell) [1]: 1\n",
           replace => " Specify the number of dependent variables [1]: 1\n"
         );
       }else{
         print ".. No Match for ADV/DIFF";
       }
     }
     
     # Justin Fernandez, 29 Oct 2002
     #   Adds a question for contact mechanics if the problem
     #   involves 3D static linear elasticity with Galerkin FE and zero thermal effects
     #   Will only update once.

     my $one = qr/Specify whether \[1\]:\n(\s*\(\d\)(\s|\w|\d|\-)*\n){6}\s*1/;
     my $two = qr/Specify equation \[1\]:\n(\s*\(\d\)(\s|\w|\d|\-|\(|\)|\.|=)*\n){9}\s*1/;
     my $three = qr/Element type is \[9\]:*\s*9/;    
     my $four = qr/Specify whether the thermal strain effects are \[0\]:\n(\s*\(\d\)(\s|\w|\d|\.)*\n){4}\s*0/;
     my $five = qr/Specify whether solution is by \[1\]:\n(\s*\(\d\)(\s|\w|\d|\-)*\n){6}\s*1/;
     my $new = " Does this problem involve contact mechanics [N]? N\n";

     if($file->match(pattern =>$new)) {
       print "Already Modified";
     }else{
      if(($file->match(pattern=>$one))&&($file->match(pattern=>$two))&&($file->match(pattern=>$three))&&($file->match(pattern=>$four))&&($file->match(pattern=>$five))){
         print "Match";
        $file->add_after(
           target => qr/Do you want the global matrices stored as sparse matrices \[Y\]\? [Nn]\n/,
           insert => $new,
           unless => qr/Do you want the global matrices stored as sparse matrices \[Y\]\? [Yy]\n/
         );
         $file->add_after(
           target => qr/Do you want to calculate the sparsity pattern for the global matrices \[Y\]\? [Yy]\n/,
           insert => $new,
           unless => qr//
         );     
       }else{
         print "No Match for Contact";
       }
     }   
     
     
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);


