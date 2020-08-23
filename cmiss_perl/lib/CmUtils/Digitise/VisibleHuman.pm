package CmUtils::Digitise::VisibleHuman;

use strict;
use warnings;

=head1 CmUtils::Digitise::VisibleHuman

Perl subroutines for CMISS that are used in the digitizing process
and are specific to the Visible Human data set

=cut

use CmUtils::Exporter();
our( $VERSION, @ISA, @EXPORT, @EXPORT_OK );

$VERSION     = 1.00;
@ISA         = qw( CmUtils::Exporter );
@EXPORT      = qw( );
@EXPORT_OK   = qw( VH_CreateImageSet makemask editdata );

=head2 CONTAINS

 VH_CreateImageSet
 makemask 
 editdata

=cut

#-------------------------------------------------
sub VH_CreateImageSet {
#-------------------------------------------------

=head1 B<VH_CreateImageSet>

  $images = &VH_CreateImageSet( $pattern, $FirstImage, $LastImage, $ImageIncr,
                                $ImageDir, $FileExt );

This routine takes in the source image directory, the image file extension and
information about the image numbers, and returns a pointer to an array of the
images where each of the array elements is a hash of attributes that has been
set up specifically for the Visible Human image set.

=cut

  my( $pattern, $FirstImage, $LastImage, $ImageIncr, $ImageDir, $FileExt ) = @_;
  my $images;

  for( my $img = $FirstImage; $img <= $LastImage; $img = $img + $ImageIncr ) {
    my $name = sprintf( "$pattern", $img );
    my $image = {
      name      => $name,
      filename  => $ImageDir . $name,
      ext       => $FileExt,
      imgloaded => 0,
      visible2D => 0,
      glyph     => "sphere",
      colour    => "green",
    };
    
    push( @{$images}, $image );
  }

  return $images;
}

#-------------------------------------------------
sub makemask { 
#-------------------------------------------------

=head1 B<makemask>

  &makemask();

Make a mask that allows this editing of data points in the XY plane.
This subroutine is currently not called by any digitiser routines.

=cut

  #makes editing mask for editdata
  cmiss::cmiss("gfx define field edit.x edit_mask field coordinates edit_mask 1 0 0");
  cmiss::cmiss("gfx define field edit.y edit_mask field coordinates edit_mask 0 1 0");
  cmiss::cmiss("gfx define field edit.z edit_mask field coordinates edit_mask 0 0 1");
  cmiss::cmiss("gfx define field x_dirn coordinate_system rectangular_cartesian constant number_of_values 3 values 1 0 0");
  cmiss::cmiss("gfx define field y_dirn coordinate_system rectangular_cartesian constant number_of_values 3 values 0 1 0");
  cmiss::cmiss("gfx define field z_dirn coordinate_system rectangular_cartesian constant number_of_values 3 values 0 0 1");
}

#-------------------------------------------------
sub editdata { 
#-------------------------------------------------

=head1 B<editdata>

  &editdata( $dgroup );

Allows the data point group given by $dgroup to be edited in only the XY plane.
This subroutine is currently not called by any digitiser routines.

=cut

  my ( $dgroup ) = @_;
  
  cmiss::cmiss("gfx modify g_element $dgroup data_points coordinate coordinates glyph sphere size \"2*2*2\" centre 0,0,0 select_on material green selected_material purple");
  cmiss::cmiss("gfx modify g_element $dgroup data_points coordinate edit.x glyph arrow_solid size \"5*1*1\" centre 0,0,0 orientation x_dirn scale_factors \"0*0*0\" draw_selected material red selected_material red");
  cmiss::cmiss("gfx modify g_element $dgroup data_points coordinate edit.y glyph arrow_solid size \"5*1*1\" centre 0,0,0 orientation y_dirn scale_factors \"0*0*0\" draw_selected material bluey selected_material bluey");
  cmiss::cmiss("gfx data_tool edit");
  cmiss::cmiss("gfx data_tool no_create");
  cmiss::cmiss("gfx data_tool no_define");
}

1;
 
__END__
