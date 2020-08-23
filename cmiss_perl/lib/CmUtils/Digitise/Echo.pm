package CmUtils::Digitise::Echo;

use strict;
use warnings;

use File::Basename qw/basename/;

=head1 CmUtils::Digitise::Echo

Perl subroutines for CMISS that are used in the digitizing process
and are specific to sets of Echo images

=cut

use CmUtils::Exporter();
our( $VERSION, @ISA, @EXPORT, @EXPORT_OK );

$VERSION     = 1.00;
@ISA         = qw( CmUtils::Exporter );
@EXPORT      = qw( );
@EXPORT_OK   = qw( Echo_CreateImageSet );

=head2 CONTAINS

 Echo_CreateImageSet

=cut

#-------------------------------------------------
sub Echo_CreateImageSet {
#-------------------------------------------------

=head2 B<Echo_CreateImageSet>

  $images = &Echo_CreateImageSet( $ImageDir, $FileExt, [$pattern] );

This routine takes in the source image directory, the image file extension and
and an optimal pattern to match, and returns a pointer to an array of the
images where each of the array elements is a hash of attributes.

=cut

  my( $ImageDir, $FileExt, $pattern ) = @_;
  $pattern = qr/./ unless $pattern;
  $ImageDir .= "/" unless $ImageDir =~ /\/$/;

  my $images;

  foreach my $file ( glob "$ImageDir*$FileExt" ) {
    next unless $file =~ /$pattern/;
    my $name = basename($file, $FileExt);
    my $image = {
      name      => $name,
      filename  => $ImageDir . $name,
      ext       => $FileExt,
      imgloaded => 0,
      visible2D => 0,
      glyph     => "sphere",
    };
    
    push( @{$images}, $image );
  }

  return $images;
}

1;
 
__END__
