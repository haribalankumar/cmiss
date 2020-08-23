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
    
    # Created: Mark Trew 17-APRIL-2003
    #   Adds grid-based finite volume option to grid solution methods 
    #   Also adds grid-based finite element option if it does not
    #   exist.

    # The option must be able to be added at various places through 
    # an irequa file.

    # Regular expression place identifiers
    my $Option1 = qr/ {3}\(6\) Grid-based Finite Element *\n/;
    my $Option2 = qr/(?s) {3}\(4\) Collocation.*\(5\) *\n/;
    my $Option3 = qr/(?s) {3}\(4\) Collocation.*\(5\) Finite volume technique *\n/;
    
    # Strings to use for line deletion and replacement
    my $Colline      = "   (4) Collocation\n";
    my $FVline       = "   (5) Finite volume technique\n";
    my $BlankFVline  = "   (5)\n";
    my $GFEline      = "   (6) Grid-based Finite Element\n";
    my $GFVline      = "   (7) Grid-based Finite Volume\n";

    # Look to see if Option1 exists
    if ($file->match(pattern => $Option1)) {

      # Globally delete any occurrences of the string to be added.
      while ($file->match(pattern => $GFVline)) {
        $file->delete(string => $GFVline);
      }

      # Add new line
      $file->{text} =~ s/($Option1)/$GFEline$GFVline/g;
    }

    # Look to see if Option2 or Option3 exists
    elsif ($file->match(pattern => $Option2) | ($file->match(pattern => $Option3))) {

      # Globally delete any occurrences of the string to be added.
      my $DeleteLine = $GFEline.$GFVline;
      print "$DeleteLine";
      while ($file->match(pattern => $DeleteLine)) {
        $file->delete(string => $DeleteLine);
      }

      # Add new line
      if ($file->match(pattern => $Option3)) { 
        $file->{text} =~ s/($Option3)/$Colline$FVline$GFEline$GFVline/g;
      }
      else {
        $file->{text} =~ s/($Option2)/$Colline$BlankFVline$GFEline$GFVline/g;
      }      
    }    
    else {
      print " Could not locate position to insert. ";
    }

    # EWR: Added option for using pressure bcs with compressive material
    $file->add_after
      (
	target =>  qr/(?gm)^   \(6\) Compressible \+ fluid for lung/,
	insert =>  "\n   (7) Compressible + face pressure bcs",
	unless =>  qr/(?m)^   \(7\) Compressible \+ face pressure bcs/
       );

  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);

