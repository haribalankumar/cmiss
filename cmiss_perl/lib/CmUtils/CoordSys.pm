package CmUtils::CoordSys;
require 5.006;

=head1 CmUtils::CoordSys

Routines to convert from various coordinate systems to RC.  Angles are 
expected in degrees and are converted to radians.  Either two or three 
coordinates may be specified, as either an array, a list, or an
array reference.

The following four calls return the same answer.

  use CmUtils::CoordSys qw/Prolate_To_RC/;
  $focus = 1;

  @rc = Prolate_To_RC($focus,5,30,90);

  @pcoord = (5, 30, 90);
  @rc = Prolate_To_RC($focus,@pcoord);

  $pref = \@pcoord;
  @rc = Prolate_To_RC($focus,$pref);

  $pref = [5, 30, 90];
  @rc = Prolate_To_RC($focus,$pref);

The wrapper routines each call the base routine B<Convert_To_RC()>.  The
C<Math::Trig> module must be available.

=head1 TO DO

Conversion to other coordinate systems (not just rc).

=cut

use strict;
use warnings;

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK);

=head1 VERSION

1.01 (28 February 2001)

=cut

  $VERSION     = 1.01;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw();
  @EXPORT_OK   = qw(
    &Convert_To_RC 
    &Cylindrical_To_RC &Spherical_To_RC &Prolate_To_RC &Oblate_To_RC
  );
}

# GBS
sub Convert_To_RC ($$@) {

=head2 B<Convert_To_RC>

Called as

  @RCpt = Convert_To_RC($CoordSys, $Focus, @Coords);

CoordSys = 1 for rectangular cartesian coordinates
           2 for cylindrical polar coordinates
           3 for spherical polar coordinates
           4 for prolate spheriodal coordinates
           5 for oblate spheroidal coordinates

Note that angles are expected in degrees and are converted to radians.  Either
two or three coordinates may be specified, as either an array, a list, or an
array reference.

=cut

  use Math::Trig;

  my($CoordSys) = shift;
  my($Focus) = shift;
  my(@Pt,@nPt);
  if ($_[0] =~ /ARRAY/) {
    @Pt=@{$_[0]};
  } else {
    @Pt = @_;
  }

  if ($CoordSys == 1) {
    @nPt=@Pt;
  } 
  elsif ($CoordSys == 2) {
    if (@Pt == 1) {
      warn ('>>Error: Requires at least 2 dimensions');
    } 
    elsif (@Pt == 2) {
      my ($r,$theta) = ($Pt[0], deg2rad($Pt[1]));
      $nPt[0] = $r * cos($theta);
      $nPt[1] = $r * sin($theta);
    } 
    elsif (@Pt == 3) {
      my($r,$theta,$z) = ($Pt[0], deg2rad($Pt[1]), $Pt[2]);
      $nPt[0] = $r * cos($theta);
      $nPt[1] = $r * sin($theta);
      $nPt[2] = $z;
    } 
    else {
      warn ('>>Error: Too many coordinates');
    }
  } 
  elsif ($CoordSys == 3) {
    if (@Pt == 1) {
      warn ('>>Error: Requires at least 2 dimensions');
    }
    elsif (@Pt == 2) {
      my ($r,$theta) = ($Pt[0], deg2rad($Pt[1]));
      $nPt[0] = $r * cos($theta);
      $nPt[1] = $r * sin($theta);
    } 
    elsif (@Pt == 3) {
      my($r,$theta,$phi) = ($Pt[0], deg2rad($Pt[1]), deg2rad($Pt[2]));
      $nPt[0] = $r * cos($theta) * cos($phi); 
      $nPt[1] = $r * sin($theta) * cos($phi);
      $nPt[2] = $r * sin($phi);
    } 
    else {
      warn ('>>Error: Too many coordinates');
    }
  } 
  elsif ($CoordSys == 4) {
    if (@Pt == 1) {
      warn ('>>Error: Requires at least 2 dimensions');
    } 
    elsif (@Pt == 2) {
      my($lambda,$mu) = ($Pt[0], deg2rad($Pt[1]));
      $nPt[0] = $Focus * cosh($lambda) * cos($mu);
      $nPt[1] = $Focus * sinh($lambda) * sin($mu);
    } 
    elsif (@Pt == 3) {
      my($lambda,$mu,$theta) = ($Pt[0], deg2rad($Pt[1]), deg2rad($Pt[2]));
      $nPt[0] = $Focus * cosh($lambda) * cos($mu);
      $nPt[1] = $Focus * sinh($lambda) * sin($mu) * cos($theta);
      $nPt[2] = $Focus * sinh($lambda) * sin($mu) * sin($theta);
    } 
    else {
      warn ('>>Error: Too many coordinates');
    }
  } 
  elsif ($CoordSys == 5) {
    if (@Pt == 1) {
      warn ('>>Error: Requires at least 2 dimensions');
    } 
    elsif (@Pt == 2) {
      my($lambda,$mu) = ($Pt[0], deg2rad($Pt[1]));
      $nPt[0] = $Focus * cosh($lambda) * cos($mu);
      $nPt[1] = $Focus * sinh($lambda) * sin($mu);
     } 
    elsif (@Pt == 3) {
      my($lambda,$mu,$theta) = ($Pt[0], deg2rad($Pt[1]), deg2rad($Pt[2]));
      $nPt[0] = $Focus * cosh($lambda) * cos($mu) * cos($theta);
      $nPt[1] = $Focus * sinh($lambda) * sin($mu); 
      $nPt[2] = $Focus * cosh($lambda) * cos($mu) * sin($theta);
    } 
    else {
      warn ('>>Error: Too many coordinates');
    }
  } 
  else {
    warn ('>>Error: Unknown coordinate system');
  }
  return @nPt;
}

=head2 B<Cylindrical_To_RC>

  @RCpt = Cylindrical_To_RC(@pt);

=cut
#GBS
sub Cylindrical_To_RC (@) {
  return Convert_To_RC(2,undef,@_);
}

=head2 B<Spherical_To_RC>

  @RCpt = Spherical_To_RC(@pt);

=cut
#GBS
sub Spherical_To_RC (@) {
  return Convert_To_RC(3,undef,@_);
}

=head2 B<Prolate_To_RC>

  @RCpt = Prolate_To_RC($focus,@pt);

=cut
#GBS
sub Prolate_To_RC ($@) {
  my ($focus, @pt) = @_;
  return Convert_To_RC(4,$focus,@pt);
}

=head2 B<Oblate_To_RC>

  @RCpt = Oblate_To_RC($focus,@pt);

=cut
#GBS
sub Oblate_To_RC ($@) {
  my ($focus, @pt) = @_;
  return Convert_To_RC(5,$focus,@pt);
}


1;

__END__
