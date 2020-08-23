package CmUtils::ViewNormal;

use strict;
use warnings;
use CmUtils::File qw/readExnode/;
use CmUtils::Vector;
use Carp;

=head1 CmUtils::ViewNormal

Sets the view of a named window normal to a plane defined by 
4 nodes in an F<.exnode> file.  

=cut

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK);

  $VERSION     = 1.02;
  @ISA         = qw( CmUtils::Exporter );
  @EXPORT      = qw( &viewNormal );
  @EXPORT_OK   = qw( &computeNormal &updateViewNormal );
}

=head1 VERSION

1.03 (17 December 2003)
=head1 CHANGES

=head2 1.03

Added B<updateViewNormal> which uses the Cmiss direct interface modules to 
interrogate the scene viewer and try and maintain a consistent interest point
when changing slices.  The routine compiles conditionally so that this 
Digitise module will still compile and fall back to the old behaviour if the
function does not initialise.

=head2 1.02

Added B<computeNormal> to compute normals from a exnode file without changing
windows - in essence prefilling the cache.  This node file may contain multiple
node groups, and each one is computed.

=head2 1.01

Added caching of nodal information in a hash C<%info>.

=head1 VARIABLES

=over 4

=item $CmUtils::ViewNormal::Scale

Scaling for viewing image.  By default, Scale is 100.

=item $CmUtils::ViewNormal::ViewRange

Range for near and far clipping planes.  By default, ViewRange is 10.

=back

=cut

our $Scale = 100;
our $ViewRange = 10;
my %info;

=head1 B<viewNormal(filename [,window] [,scale] [,viewrange])>

By default, the window to set is "Window 1".  For example

  use CmUtils::ViewNormal;
  viewNormal("axial10");

sets window 1 so that axial10.exnode (which may have axial10.rgb drawn in it)
is normal to the view.  Used so far to give a view normal to echo images.


The F<Scale> and F<ViewRange> can also be overridden for this call only. 
Neither can be zero;

=cut

sub viewNormal {
  my $filename = shift;
  croak "Must provide a filename to 'viewNormal'" unless $filename;
  my $window = shift || 1;
  my $scale = shift || $Scale;
  my $viewrange = shift || $ViewRange;
  $filename =~ s/\.exnode//;  # Remove ",exnode" if it exists.
  unless (exists $info{$filename}) {
    my $group                = readExnode($filename);
    my @nodes                = checkNodes($group);
    $info{$filename}{normal} = calculateNormal(@nodes);
    $info{$filename}{centre} = calculateCentre(@nodes);
    $info{$filename}{up}     = calculateUp(@nodes);
  }
  rotateView($window, $info{$filename}, $scale, $viewrange);
}

=head1 B<updateViewNormal(filename [,window] [,scale] [,viewrange])>

By default, the window to set is "Window 1".  For example

  use CmUtils::ViewNormal;
  updateViewNormal("axial10");

sets window 1 so that axial10.exnode (which may have axial10.rgb drawn in it)
is normal to the view.  Used so far to give a view normal to echo images.

The F<Scale> and F<ViewRange> can also be overridden for this call only. 
Neither can be zero;

This version is a replacement for viewNormal which maintains the same
viewing range when you change planes and ensures that the centre of the
new view is the intersection of the previous viewing direction and this
plane.

=cut

my $updateViewNormal__;
my $fallbackViewNormal__;

sub updateViewNormal {
  my $filename = shift;
  croak "Must provide a filename to 'updateViewNormal'" unless $filename;
  my $window = shift || 1;
  my $scale = shift || $Scale;
  my $viewrange = shift || $ViewRange;
  $filename =~ s/\.exnode//;  # Remove ".exnode" if it exists.
  unless (exists $info{$filename}) {
    my $group                = readExnode($filename);
    my @nodes                = checkNodes($group);
    $info{$filename}{normal} = calculateNormal(@nodes);
    $info{$filename}{centre} = calculateCentre(@nodes);
    $info{$filename}{up}     = calculateUp(@nodes);
  }

  # Conditional compilation to see if the Cmiss modules are available.
  if (! defined $fallbackViewNormal__) {
    if (defined $updateViewNormal__ or ($updateViewNormal__ = eval q[
      package CmUtils::ViewNormal;
      use strict;
      use Cmiss::Graphics_window;
      use Cmiss::Scene_viewer;
      use CmUtils::Vector;

      return sub {
        my $filename = shift;
        croak "Must provide a filename to 'updateViewNormal__'" unless $filename;
        my $window = shift;
        croak "Must provide a window name to 'updateViewNormal__'" unless $window;
        my $scale = shift;
        croak "Must provide a scale to 'updateViewNormal__'" unless $scale;
        my $viewrange = shift;
        croak "Must provide a viewrange to 'updateViewNormal__'" unless $viewrange;
        my $info = shift;

        my $scene_viewer = Cmiss::Graphics_window::get_scene_viewer_by_name($window, 0);
        my @view = Cmiss::Scene_viewer::get_lookat_parameters($scene_viewer);

        my $eyepoint = new CmUtils::Vector($view[0], $view[1], $view[2]);
        my $lookatpoint = new CmUtils::Vector($view[3], $view[4], $view[5]);

        # Find point on line and lookat which intersects the viewing plane.
        # t = (normal . centre) - (normal . line_origin) / (normal . line_direction)
        my $ratio = dotp($$info{$filename}{normal}, norm($lookatpoint - $eyepoint));
        my $t;
        if ($ratio > 0.0001 || $ratio < -0.0001) {
          $t = (dotp($$info{$filename}{normal}, $$info{$filename}{centre}) -
            dotp($$info{$filename}{normal}, $lookatpoint)) / $ratio;
        } else {
          $t = 0.0;
        }

        # print ("ratio $ratio t $t normal $$info{$filename}{normal} lookat $lookatpoint eye $eyepoint\n");
        my $newCentre = $lookatpoint + norm($lookatpoint - $eyepoint) * $t;

        my $eye    = $newCentre + $$info{$filename}{normal} * $scale;
        my @eye    = list($eye);
        my @centre = list($newCentre);
        my @up     = list($$info{$filename}{up});
        my ($near, $far) = ($scale - $viewrange, $scale + $viewrange);
        $near = 5.0 if $near < 5.0;

        Cmiss::Scene_viewer::set_lookat_parameters_non_skew($scene_viewer, @eye, @centre, @up);
        Cmiss::Scene_viewer::set_near_and_far_plane($scene_viewer, $near, $far);
      }
    ])) {
      return &$updateViewNormal__($filename, $window, $scale, $viewrange, \%info);
    } else {
      print ("Unable to initialise updateViewNormal, you probably have not installed the Cmiss perl modules correctly, falling back to standard viewNormal\n$@\n");
      $fallbackViewNormal__ = \&viewNormal;
      return &$fallbackViewNormal__($filename, $window, $scale, $viewrange);
    }
  } else {
    return &$fallbackViewNormal__($filename, $window, $scale, $viewrange);
  }
}

=head1 B<computeNormal(filename)>

Reads node file F<filename> and computes and stores normals for all node groups
in the file.

  use CmUtils::ViewNormal qw(computeNormal);
  computeNormal("file.exnode");

=cut

sub computeNormal {
  my $filename = shift;
  croak "Must provide a filename to 'computeNormal'" unless $filename;
  my @groups = readExnode($filename);
  foreach my $group (@groups) {
    unless (exists $info{$filename}) {
      my @nodes                = checkNodes($group);
      $info{$filename}{normal} = calculateNormal(@nodes);
      $info{$filename}{centre} = calculateCentre(@nodes);
      $info{$filename}{up}     = calculateUp(@nodes);
    }
  }
}

=head1 Other routines

The above routines call the following (unexported) routines.

=head2 B<checkNodes(group)>

Checks that nodal information in a group is valid, and puts the nodal
coordinates into the array as vectors.  Requires a single fieldset,
with 4 nodes of 3 coordinates each.

=cut

sub checkNodes {
  my ($group) = shift;
  my (@nodes);
  
  unless (1 == $group->getFieldSets()) { 
    croak "Too many FieldSets"
  };
  unless (3 == $group->getFieldSet(0)->numValues()) { 
    croak "Wrong number of values (need 3)"
  };
  unless (4 == $group->numberOfNodes()) { 
    croak "Wrong number of nodes (need 4)"
  };
  foreach my $n ($group->getNodes()) {
    push @nodes, new CmUtils::Vector($n->getValues());
  }
  
  return @nodes;
}

=head2 B<calculateNormal(nodelist)>

Calculates the normal to the plane by taking the cross-product of two vectors
in the plane.  Returns a normal vector as a CmUtils::Vector.

  $normal = calculateNormal(@nodes);

=cut

sub calculateNormal {
# normal vector is cross-product of two other vectors
# pass in three points
  my ($x1, $x2, $x3) = @_;
  my $n = norm(crossp($x3 - $x2, $x1 - $x2));
  return $n;
}

=head2 B<calculateCentre(nodelist)>

Calculates the centre of the nodes.

  $centre = calculateCentre(@nodes);

=cut

sub calculateCentre {
# calculate the centre of the 4 nodes
  my ($x1, $x2, $x3, $x4) = @_;
  my $c = ($x1 + $x2 + $x3 + $x4)/4;
  return $c;
}

=head2 B<calculateUp(nodelist)>

Calculates the up-vector of the image (node 3 - node 1).

  $up = calculateUp(@nodes);

=cut

sub calculateUp {
  my ($x1, $x2, $x3, $x4) = @_;
  my $u = $x3 - $x1;
  return $u;
}

=head2 B<rotateView(window, info, scale, viewrange)>

Rotates the view to align with a vector, centred on a point, with a given
up-vector.

  rotateView(1, $info);

rotates window 1 so that the interest point is $info->{centre}, and the view
direction is defined by $info->{normal} and an up-vector of $info->{up}.

=cut

sub rotateView {
  my ($window, $info, $scale, $viewrange) = @_;

  my $eye    = $info->{centre} + $info->{normal} * $scale;
  my @eye    = list($eye);
  my @centre = list($info->{centre});
  my @up     = list($info->{up});
  my ($near, $far) = ($scale - $viewrange, $scale + $viewrange);
  $near = 5.0 if $near < 5.0;
  cmiss::cmiss("gfx modify window $window view eye_point @eye view_angle 100 interest_point @centre up @up near $near far $far");
}

1;
