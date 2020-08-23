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
      # Added: Greg Sands, 30 Nov 2001
      #   Whole lot of changes to ipinve files
      my $orig =<<EOF;
 Specify the imaging approach  [1]
   (1) Zero-crossing
   (2) Unused
   (3) Unused
    1
EOF
      my $new =<<EOF;

 Enter the regularisation scheme [1]:
   (1) None
   (2) Surface Laplacian
   (3) Unused
    2

 Specify regularisation parameter choice [1]:
   (1) Constant
   (2) Ratio
   (3) L-curve zero
   (4) Minimum Residual Norm
   (5) Unused
    3

 Enter the initial regularisation parameter [0.5]: 1.0d-9
EOF
      $file->add_after(
        target => $orig,
        insert => $new,
        unless => "Minimum Residual Norm"
      );
      $file->replace(
        target  => "imaging approach  [",
        replace => "imaging approach ["
      );
      $file->replace(
        target  => "\n Specify inverse approach [",
        replace => " Specify inverse approach ["
      );
      $file->{text} =~ s/: \[/ [/g;
      $file->{text} =~ s/\]\n/]:\n/g;
    }
    1;
  }
);

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
