package CmUtils::RotationTransform;

use strict;
use warnings;

use Carp;

=head1 CmUtils::RotationTransform

Converts a 6 element translation/axis-rotation vector (as provided by the
'rotate' program) into a 16 element transformation matrix.  This module
requires the C<PDL> and C<Math::Trig> modules to be installed.  The 
transformation matrix is returned as a 1D 16-element list.

=head2 B<angle2trans>

Expects an array/list/ptr-to-array of 6 numbers:
  x-trans y-trans z-trans rot/x rot/y rot/z
where the rotations about the axes are given in radians.

=head2 B<degangle2trans>

Expects an array/list/ptr-to-array of 6 numbers:
  x-trans y-trans z-trans rot/x rot/y rot/z
where the rotations about the axes are given in degrees.

=head2 B<readDat>

Expects a file name.  File contains 6 numbers:
  x-trans y-trans z-trans rot/x rot/y rot/z
as output from the 'rotate' program.

=head2 VARIABLES

=item $CmUtils::RotationTransform::Transpose

Set this to '0' for returned matrix
not to be transposed.  By default the matrix is returned by 
column (as required for CMGUI).

=cut

BEGIN {
  use CmUtils::Exporter ();
  our ( $VERSION, @ISA, @EXPORT, @EXPORT_OK );

=head1 VERSION

1.02 (28 March 2001)

=head1 CHANGES

1.02 - now using PDL::Lite for speed

=cut

  $VERSION     = 1.02;
  @ISA         = qw( CmUtils::Exporter );
  @EXPORT      = qw( degangle2trans angle2trans );
  @EXPORT_OK   = qw( readDat );
}

# Set $Transpose to "0" for no transposition
our $Transpose = 1;

sub readDat {
  my $file = shift;
  open DAT, "<$file" or croak "File $file not found";
  my $line = <DAT>;
  close DAT;
  return split ' ', $line;
}

sub degangle2trans (@) {
  use Math::Trig qw( deg2rad );
  my @transform;

  if ($_[0] =~ /ARRAY/) {
    @transform=@{$_[0]};
  } else {
    @transform = @_;
  }

  do { $_= deg2rad($_) } foreach @transform[3..5];
  angle2trans(@transform);
}

sub angle2trans (@) {
  use PDL::Lite;
  
  my ( @transform );
  my ( $m, $s, $c );
  
  if ($_[0] =~ /ARRAY/) {
    @transform=@{$_[0]};
  } else {
    @transform = @_;
  }
  unless (@transform == 6) {return undef};
  
  my $trans = pdl( 
    [1,0,0,$transform[0]],
    [0,1,0,$transform[1]],
    [0,0,1,$transform[2]],
    [0,0,0,1]);
  
  $s = sin ($transform[3]);
  $c = cos ($transform[3]);
  $m = pdl(
    [1, 0,  0,0],
    [0,$c,-$s,0],
    [0,$s, $c,0],
    [0, 0,  0,1]);
  $trans = $m x $trans;
  
  $s = sin ($transform[4]);
  $c = cos ($transform[4]);
  $m = pdl(
    [ $c,0,$s,0],
    [  0,1, 0,0],
    [-$s,0,$c,0],
    [0, 0,  0,1]);  
  $trans = $m x $trans;
  
  $s = sin ($transform[5]);
  $c = cos ($transform[5]);  
  $m = pdl(
    [$c,-$s,0,0],
    [$s, $c,0,0],
    [ 0,  0,1,0],
    [ 0,  0,0,1]);
  $trans = $m x $trans;
  
  # transpose for CMGUI (read down columns)
  return $Transpose ? $trans->transpose()->list()
                    : $trans->list();
}

1;
__END__
