package CmUtils::Digitise;

# Original: March 2002
# Modified: May 2002 - added additional 3d window option
# Modified: GBS May 2002 - generalised for Visible Human / Echo / MRI / ...
# Modified: AJC Jan 2006 - quick hack to change from Gtk to Gtk2.


use strict;
use warnings;

use CmUtils;
use CmUtils::ViewNormal qw( viewNormal updateViewNormal );
use CmUtils::File qw( cmRead cmWrite );
use CmGtk2;
use File::Temp qw( :mktemp );

=head1 CmUtils::Digitise

CmUtils::Digitise contains a selection of routines for digitising a set of
images. 

=cut

use CmUtils::Exporter();
our( $VERSION, @ISA, @EXPORT, @EXPORT_OK );

$VERSION   = 2.00;
@ISA       = qw( CmUtils::Exporter );
@EXPORT    = qw( );
@EXPORT_OK = qw( cmguiInit loadExdata loadExdataGroup saveExdata splitExdata
                 loadImage loadImageSet visibleon visibleoff update_view
                 trace checkImageSet setDataGroups digitiseSlices setGlyphs);

=head2 CONTAINS

 cmguiInit *
 loadExdata 
 loadExdataGroup *
 saveExdata
 splitExdata *
 loadImage
 unloadImage
 loadImageSet *
 visibleon
 visibleoff
 update_view
 trace
 checkImageSet *
 setDataGroups *
 digitiseSlices *
 setGlyphs *

* = routines you may wish to call directly.  Short example of usage:

  use CmUtils::Digitise qw( cmguiInit checkImageSet digitiseSlices );
  cmguiInit( 600,600 );
  
  $images = [ glob( "*.rgb" ) ];
  checkImageSet( $images );
  digitiseSlices( $images );

=head2 Global variables

  %settings

List of settings and global information.

      <key>       <default>   <description>
User settings:
      {view3d}    0         whether 3D window is used
      {keep_view} 0         keeps the closest viewpoint when changing frames
      {dynamic}   1         images loaded dynamically
      {max_load}  10        maximum number of images in memory
      {saveall}   undef     whether to save all datagroups to named file
Set by "setDataGroups":
      {datanames} { data => "green" }
                            list of names to digitise, and associated colours
      {datalist}  ["data"]  sorted list of names
Used internally:
      {current}   0         index of current image
      {previous}  0         index of previous image
      {loaded}    undef     list of last n viewed images (= number loaded)
      {dataindex} 0         index of current data name in {datalist}
      {combo}     undef     combo widget with list of images
      {toggle}    undef     toggle button for digitising


=cut

our( %settings ) = (
  view3d    => 0,
  keep_view => 0,
  current   => 0,
  previous  => 0,
  dynamic   => 1,
  loaded    => undef,
  max_load  => 10,
  datanames => { data => "green" },
  datalist  => [ "data" ],
  dataindex => 0,
  saveall   => undef,
  combo     => undef,
  combosignalid => undef,
  toggle    => undef,
);

#-------------------------------------------------
sub cmguiInit {
#-------------------------------------------------

=head1 B<cmguiInit>

  &cmguiInit( $WinWidth2d, $WinHeight2d [, $View3d, $WinWidth3d, $WinHeight3d] );

Create a CMGUI 2d window of an appropriate size and also create five materials
of different colours for use.  Last 3 variables are optional, but if defined
create a second window use for 3d viewing of all datapoints.

=cut

  my( $WinWidth2d, $WinHeight2d, $View3d, $WinWidth3d, $WinHeight3d ) = @_;

  # Create a CMGUI window
  cmiss::cmiss( "gfx create window 1");
  cmiss::cmiss( "gfx modify window 1 layout width $WinWidth2d height $WinHeight2d" );

  # Create the 3d CMGUI window and scene if necessary
  if( $View3d || $settings{view3d} ) {
    cmiss::cmiss( "gfx create scene alldata" );
    cmiss::cmiss( "gfx modify scene alldata manual" );
    cmiss::cmiss( "gfx create window 2" );
    cmiss::cmiss( "gfx modify window 2 layout width $WinWidth3d height $WinHeight3d" );
    cmiss::cmiss( "gfx modify window 2 image scene alldata" );
    $settings{view3d} = 1;
  }

  # Create some coloured materials
  cmiss::cmiss( "gfx create material red ambient 1 0 0 diffuse 1 0 0 emission 0 0 0 specular 0.8 0.8 0.8 shininess 0.8" );
  cmiss::cmiss( "gfx create material green ambient 0 0.6 0 diffuse 0 0.6 0 emission 0 0 0 specular 0.8 0.8 0.8 shininess 0.8" );
  cmiss::cmiss( "gfx create material purple ambient 1 0 1 diffuse 1 0 1 emission 0 0 0 specular 1 1 1 shininess 0.74" );
  cmiss::cmiss( "gfx create material bluey ambient 0 0.2 0.4 diffuse 0 0.5 1 emission 0 0 0 specular 0.5 0.5 0.5 shininess 0.8" );
  cmiss::cmiss( "gfx create material gold ambient 1 0.7 0 diffuse 1 0.7 0 emission 0 0 0 specular 1 1 0.8 shininess 0.8" );

  # Set the initial viewing position
  cmiss::cmiss( "gfx modify window 1 layout 2d ortho_axes z -y eye_spacing 0.25" );
}

#-------------------------------------------------
sub loadExdata {
#-------------------------------------------------

=head1 B<loadExdata>

  &loadExdata( $image );

This routine loads the data points associated with a given image.

=cut

  my $image = shift;

  my $name = $image->{name};
  my $filename = "$name.exdata";
  return unless -r $filename;

  print( "Loading data: $name\n" );
  cmiss::cmiss( "gfx read data $name" );
  my @datagroups = cmRead( $filename );

  foreach my $group ( @datagroups ) {
    my $dgroup = $group->name();
    if ( $group->numberOfNodes() == 0 ) {
      # Destroy data group if it contains no data points
      cmiss::cmiss( "gfx destroy dgroup $dgroup" );
    } else {
      # draw points to 3d window
      unless ( $image->{datagroups}{$dgroup} ) {
        my( $dgroupname ) = ($dgroup =~ /.*_(.*)$/);
        if ( $settings{view3d} ) {
          cmiss::cmiss( "gfx draw group $dgroup as $dgroup scene alldata" );
          cmiss::cmiss( "gfx modify g_elem $dgroup general clear scene alldata" );
          cmiss::cmiss( "gfx modify g_elem $dgroup data_point glyph cross size $image->{glyphsize} material $settings{datanames}{$dgroupname} scene alldata" );
        }
        $image->{datagroups}{$dgroup} = 1;
      }
    }
  }
  
  return;
}

#-------------------------------------------------
sub loadExdataGroup {
#-------------------------------------------------

=head1 B<loadExdataGroup>

  &loadExdataGroup( $filename, $images);

This routine loads all data groups in a given file, then destroys all groups
not containing data points.

=cut

  my ($filename, $images) = @_;
  $filename =~ s/\.exdata$//;

  return unless -r "$filename.exdata";

  unless ( defined @{$images} ) {
    print "No images defined in loadExdataGroup -- must pass a pointer to an array!\n";
    return 0;
  }

  print "Loading data: $filename\n";
  cmiss::cmiss( "gfx read data $filename" );
  my @datagroups = cmRead( "$filename.exdata" );
  foreach my $group (@datagroups) {
    my $dgroup = $group->name();
    if ( $dgroup =~ /_/ ) { # data group
      if ( $group->numberOfNodes() == 0 ) {
        # Destroy data group if it contains no data points ...
        cmiss::cmiss( "gfx destroy dgroup $dgroup" );
      } else {
        # Find corresponding image in $images array
        my ($iname, $dgroupname) = ($dgroup =~ /(.*)_(.*)$/);

        next unless grep { $_ eq $dgroupname } @{$settings{datalist}};

        my $image = undef;
        foreach my $im (@{$images}) {
          if ($im->{name} eq $iname) {
            $image = $im;
            last;
          }

        }
        unless (defined $image) {
          print "Cannot find image for group $dgroup\n";
          next;
        }
        # ... and draw to 3d window
        unless ( $image->{datagroups}{$dgroup} ) {
          if ( $settings{view3d} ) {
            cmiss::cmiss( "gfx draw group $dgroup as $dgroup scene alldata" );
            cmiss::cmiss( "gfx modify g_elem $dgroup general clear scene alldata" );
            cmiss::cmiss( "gfx modify g_elem $dgroup data_point glyph cross size $image->{glyphsize} material $settings{datanames}{$dgroupname} scene alldata" );
#print ( "gfx modify g_elem $dgroup data_point glyph cross size $image->{glyphsize} material $settings{datanames}{$dgroupname} scene alldata" );
	    
          }
          $image->{datagroups}{$dgroup} = 1;
        }
      }
    } else { # image group
      cmiss::cmiss( "gfx set vis $dgroup off" );
    }
  }

  return;
}

#-------------------------------------------------
sub saveExdata {
#-------------------------------------------------

=head1 B<saveExdata>

  &saveExdata( $image );

Saves exdata from current image. 

=cut

  my $image = shift;

  if ( $settings{saveall} ) {
    cmiss::cmiss( "gfx write data $settings{saveall}" );
  } else {
    my @datagroups = keys %{$image->{datagroups}};
    if ( @datagroups ) {
      cmiss::cmiss( "gfx write data $image->{name} group @datagroups" );
    }
  }
}

#-------------------------------------------------
sub splitExdata {
#-------------------------------------------------

=head1 B<splitExdata>

  &splitExdata( $datafile, @groupnames );

Splits the exdata file F<$datafile> into separate files for each specified
I<groupname>.

=cut

  my ( $datafile, @groupnames ) = @_;
  
  $datafile .= ".exdata" unless $datafile =~ /\.exdata/;

  my @datagroups = cmRead( $datafile );

  foreach my $gname ( @groupnames ) {
    my $g = CmUtils::Objects::NodeGroup->new();
    $g->name( $gname );
    foreach my $gp ( @datagroups ) {
      $g->mergeGroup( $gp ) if $gp->name() =~ /$gname/;
    }
    cmWrite( "$gname.exdata", $g );
  }
}

#-------------------------------------------------
sub loadImage {
#-------------------------------------------------

=head1 B<loadImage>

  &loadImage( $image );

Load in a single image along with an exnode and an exelem file. This routine
then textures the element with the image.

=cut

  my $image = shift;

  my $name = $image->{name};
  print "Loading image: $name\n";

  my $filename = $image->{filename};
  my $IMGfilename;
  if ( $image->{imgtype} ) {
    $IMGfilename = "$image->{imgtype}$filename$image->{ext}";
  } else {
    $IMGfilename = "$filename$image->{ext}";
  }

  cmiss::cmiss( "gfx read node $filename" );
  cmiss::cmiss( "gfx read elem $filename" );
  cmiss::cmiss( "gfx modify g_element $name lines invisible" );

  cmiss::cmiss( "gfx create texture $name image $IMGfilename width 1.0 height 1.0 decal" );
  cmiss::cmiss( "gfx create material $name ambient 1 1 1 diffuse 1 1 1 emission 1 1 1 alpha 1.0 texture $name" );
  cmiss::cmiss( "gfx modify g_element $name surface material $name texture_coordinates xi" );
  $image->{imgloaded} = 1;
  $image->{visible2D} = 1;

  # ... and draw square to 3d window
  if ( $settings{view3d} ) 
  {
     cmiss::cmiss("gfx draw group $name scene alldata");
     cmiss::cmiss( "gfx modify g_element $name line material default scene alldata" );
  }

  return;
}

#-------------------------------------------------
sub unloadImage {
#-------------------------------------------------

=head1 B<unloadImage>

  &unloadImage( $image );

This routine unloads an image from memory, destroying material and texture.

=cut

  my $image = shift;

  my $name = $image->{name};
  print "Unloading image: $name\n";

  cmiss::cmiss( "gfx modify g_element $name lines delete" );
  cmiss::cmiss( "gfx modify g_element $name surface material $name texture_coordinates xi delete" );

  cmiss::cmiss( "gfx destroy material $name" );
  cmiss::cmiss( "gfx destroy texture $name" );
  cmiss::cmiss( "gfx destroy egroup $name" );
  $image->{imgloaded} = 0;
  $image->{visible2D} = 0;

  return;
}

#-------------------------------------------------
sub loadImageSet { 
#-------------------------------------------------

=head1 B<loadImageSet>

  &loadImageSet( $images );

This subroutine loops though an image list and calls loadImage for each image
unless dynamically loaded, and loads all exdata.

=cut

  my $images = shift;
  
  foreach my $image ( @{$images} ) {
    if ($settings{dynamic} == 0) {
      loadImage( $image );
    }
    loadExdata( $image );
  }

  return;
}

#-------------------------------------------------
sub visibleon {
#-------------------------------------------------

=head1 B<visibleon>

  &visibleon( $image(s), [$pattern] );

Turn the visibility on for all graphical elements that contain the string given
by pattern.

=cut

  my ($im, $pattern) = @_;
  my @images;
  if (ref($im) eq "ARRAY") {
    @images = grep { $_->{name} =~ /$pattern/ } @{$im};
  } else {
    @images = ( $im );
  }
  foreach my $image ( @images ) {
    my $gname = $image->{name};
    cmiss::cmiss( "gfx set vis $gname on" ) if $image->{imgloaded};

    # LKC do the same for the 3d window (scene alldata)
    cmiss::cmiss( "gfx set vis $gname on scene alldata" ) if $image->{imgloaded};

    foreach my $dname ( keys %{$image->{datagroups}} ) {
      cmiss::cmiss( "gfx set vis $dname on" );
    }
    $image->{visible2D} = 1;
  }
}

#-------------------------------------------------
sub visibleoff {
#-------------------------------------------------

=head1 B<visibleoff>

  &visibleoff( $image(s), [$pattern] );

Turn the visibility off on all graphical elements that match the supplied
pattern.

=cut

  my ($im, $pattern) = @_;
  $pattern = qr// unless defined $pattern;
  my @images;
  if (ref($im) eq "ARRAY") {
    @images = grep { $_->{name} =~ /$pattern/ } @{$im};
  } else {
    @images = ( $im );
  }
  foreach my $image ( @images ) {
    my $gname = $image->{name};
    cmiss::cmiss( "gfx set vis $gname off" ) if $image->{imgloaded};

    # LKC do the same for the 3d window (scene alldata)
    cmiss::cmiss( "gfx set vis $gname off scene alldata" ) if $image->{imgloaded};

    foreach my $dname ( keys %{$image->{datagroups}} ) {
      cmiss::cmiss( "gfx set vis $dname off" );
    }
    $image->{visible2D} = 0;
  }
}
  
#-------------------------------------------------
sub update_view  {
#-------------------------------------------------

=head1 B<update_view>

  &update_view( $images );

Update the 2D and 3D windows.

=cut

  my $images = shift;

  $settings{previous} = (${$settings{loaded}}[0] || 0);
  if ( $settings{previous} != $settings{current} ) {
    my $image = ${$images}[$settings{previous}];
    visibleoff( $image );
    foreach my $dgroup ( keys %{$image->{datagroups}} ) {
      # destroy data group if there were no points defined in it

		# generate a unique temporary filename, open it and close it
		# the file is opened at the same time it is generated to 
		# ensure no race conditions
		my $template = "/tmp/tempXXXXXX";
		my $suffix = ".exdata";
		my ($fh, $filename);
		($fh, $filename) = mkstemps($template,$suffix);
		close($fh);

		# Now we can use the filename to write out temporary data
		# without fear of another process (using the 
		# Digitise.pm module) writing to the same file
      cmiss::cmiss( "gfx write data $filename group $dgroup" );
		# read it back in
      my $group = cmRead( "$filename" );

      if ( $group->numberOfNodes() == 0 ) {
        cmiss::cmiss( "gfx destroy dgroup $dgroup" );
        delete $image->{datagroups}{$dgroup};
      } elsif ( $settings{view3d} ) {
        my ($dgroupname) = ($dgroup =~ /.*_(.*)$/);
        cmiss::cmiss( "gfx modify g_elem $dgroup general clear scene alldata" );
        cmiss::cmiss( "gfx modify g_elem $dgroup data_point glyph cross size $image->{glyphsize} material $settings{datanames}{$dgroupname} scene alldata" );
      }

		# remove temporary file
      unlink "$filename";
    }
  }

  my $image = ${$images}[$settings{current}];
  # View the new slice
  loadImage( $image ) unless $image->{imgloaded};
  # make sure all data groups are defined, and draw to 2d/3d windows
  foreach my $dgroupname ( @{$settings{datalist}} ) {
    my $dgroup = "$image->{name}_$dgroupname";
    unless ( $image->{datagroups}{$dgroup} ) {
      cmiss::cmiss( "gfx create dgroup $dgroup" );
      $image->{datagroups}{$dgroup} = 1;
      if ( $settings{view3d} ) {
        cmiss::cmiss( "gfx draw group $dgroup as $dgroup scene alldata" );
      }
    }
    cmiss::cmiss( "gfx modify g_element $dgroup data_points glyph $image->{glyph} size $image->{glyphsize} material $settings{datanames}{$dgroupname} select_on selected_material purple" );
    if ( $settings{view3d} ) {
      cmiss::cmiss( "gfx modify g_element $dgroup data_points glyph sphere size $image->{glyphsize} material gold select_on selected_material purple scene alldata" );
    }
  }
  $image->{visible2D} = 1;

  # unload images if over the limit
  @{$settings{loaded}} = grep { $_ != $settings{current} } @{$settings{loaded}};
  unshift @{$settings{loaded}}, $settings{current};
  while ( $#{$settings{loaded}} >= $settings{max_load} ) {
    my $old = pop @{$settings{loaded}};
    unloadImage( ${$images}[$old] );
  }

  visibleon( $image );

  if( $settings{toggle}->get_active() ) {
    trace( $images );
  } else {
    cmiss::cmiss( "gfx mod win 1 set transform_tool" );
    cmiss::cmiss( "gfx modify win 1 image view_all" );
  }


##  $settings{combo}->set_popdown_strings( 
##    map{ $_->{name} . (( scalar keys %{$_->{datagroups}} ) ? "  *" : "") } 
##    @{$images} );
##  $settings{combo}->list->select_item( $settings{current} );

# gtk2 translation of above code

  my $combobox = $settings{combo};

# Stop signals while we update the combobox so these changes don't call us back
 $combobox->signal_handler_block($settings{combosignalid});

# Remove all entries
  my $ii = scalar @{$images};
  while ( $ii > 0 ) {
       $combobox->remove_text(0);
       $ii--;
  }
# Add them back in with or without *
  foreach my $string (map{ $_->{name} . (( scalar keys %{$_->{datagroups}} ) ? "  *" : "")} @{$images}) {
      $combobox->append_text($string);
  }

  $combobox->set_active($settings{current});

# Resume signalling
  $combobox->signal_handler_unblock($settings{combosignalid});

}

#-------------------------------------------------
sub trace {
#-------------------------------------------------

=head1 B<trace>

  &trace( $images );

Set up the viewing window and clipping planes with the correct spacings
so the digitised points lie on the correct slice. Start up the data point
tool with the data creation options set.

=cut

  my( $images ) = @_;
  my $image = ${$images}[$settings{previous}];

  saveExdata( $image );

  $image = ${$images}[$settings{current}];
  if ($settings{keep_view})
  {
	 updateViewNormal( $image->{filename} );
	 updateViewNormal( $image->{filename}, 2, 150, 140 ) if $settings{view3d};
  }
  else
  {
	 viewNormal( $image->{filename} );
# LKC why would we ever want to reset the view in this window? Remove for now.
#	 viewNormal( $image->{filename}, 2, 150, 140 ) if $settings{view3d};
  }

  my $dgroupname = $settings{datalist}[$settings{dataindex}];
  my $dgroup = "$image->{name}_$dgroupname";

  cmiss::cmiss( "gfx data_tool open_dialog coordinate_field coordinates edit create group $dgroup" );
  cmiss::cmiss( "gfx mod win 1 set data_tool" );
}

#-------------------------------------------------
sub checkImageSet  {
#-------------------------------------------------

=head1 B<checkImageSet>

  &checkImageSet( $images );

Checks that the image set $images has all necessary information to load, and
fills in default values for undefined keys.  On input, $images must either be
an array of hashes e.g

  foreach $f ( "image1", "image2" ) {
    push @{$images}, {name=>$f, filename=>$f, ext=>".rgb", glyph=>"cross"};
  }
  checkImageSet( $images );

or an array of filenames (with or without extensions) e.g.

  $images = [ glob "*.rgb" ];
  checkImageSet( $images );

Other modules including CmUtils::Digitise::VisibleHuman and
CmUtils::Digitise::Echo have routines for easily creating the $images
structure.

Returns $images as a hash structure:

    <key>       <default> <description>
    {name}                Basename of image
    {filename}  {name}    Image file name including directory (no ext)
    {ext}       {.rgb}    Image file extension (including '.' if needed)
    {imgtype}   undef     Required if no extension, or if different type needed
    {imgloaded} 0         If image is loaded into cmgui
    {visible2D} 0         If image is visible in 2D window
    {glyph}     "sphere"  Glyph used for drawing data on this image
    {glyphsize} 3         Size of glyph used

=cut

  my $images = shift;

# LKC 12-MAY-2003
#  unless ( @{$images} ) {
#
  unless ( defined @{$images} ) {
    print "No images defined in checkImageSet -- must pass a pointer to an array!\n";
    return 0;
  }
  my $i = 0;
  foreach my $image ( @{$images} ) {
    if ( ref($image) ne "HASH" ) { #it is a filename
      if ( $image =~ /^(.*)(\..*)$/ ) {
        my ($name, $ext) = ($image =~ /^(.*)(\..*)$/);
        undef $image;
        $image = { name => $name, ext => $ext };
      } else {
        my $name = $image;
        undef $image;
        $image = { name => $name };
      }
    }
    unless ( defined $image->{name} ) {
      print "Image ($i) name not defined -- must at least do this!\n";
      return 0;
    }
    unless ( defined $image->{filename} ) {
      print "Image filename not defined -- setting to $image->{name}\n";
      $image->{filename} = $image->{name};
    }
    if ( defined $image->{ext} ) {
      unless (-r "$image->{filename}$image->{ext}") {
        print "Cannot locate $image->{filename}$image->{ext}\n";
        return 0;
      }
      if ( $image->{ext} eq "" && !$image->{imgtype} ) {
        print "File without extension has no type defined\n";
        return 0;
      }
      unless ( -r "$image->{filename}.exnode" ) {
        print "Cannot locate $image->{filename}.exnode\n";
        return 0;
      }
      unless ( -r "$image->{filename}.exelem" ) {
        print "Cannot locate $image->{filename}.exelem\n";
        return 0;
      }
    } else { # no extension
      if ( -r $image->{filename} ) { #test for file w/o extension
        $image->{ext} = "";
        unless ( $image->{imgtype} ) {
          print "File without extension has no type defined\n";
          return 0;
        }
      } else {
        print "Image extension not defined -- trying (";
        foreach my $ext ( qw(.rgb .jpg .tif .png .gif) ) {
          print " $ext";
          if ( -r "$image->{filename}$ext" ) {
            $image->{ext} = $ext;
            print " ... $image->{filename}$ext found";
            last;
          }
        }
        print " )\n";
        unless ( defined $image->{ext} ) {
          print "Cannot locate an image file for $image->{name}\n";
          return 0;
        }
      }
    } # extension
    if ( $image->{imgtype} ) {
      $image->{imgtype} =~ s/:$//;
      $image->{imgtype} .= ":";
    }
    unless ( defined $image->{imgloaded} ) {
      $image->{imgloaded} = 0;
    }
    unless ( defined $image->{visible2D} ) {
      $image->{visible2D} = 0;
    }
    unless ( defined $image->{glyph} ) {
      $image->{glyph} = "sphere";
    }
    unless ( defined $image->{glyphsize} ) {
      $image->{glyphsize} = 3;
    }

    $i++;
  }

  return 1;
}
  
#-------------------------------------------------
sub setDataGroups  {
#-------------------------------------------------

=head1 B<setDataGroups>

  &setDataGroups( @names_and_colours );

Sets new data groups.   e.g. setDataGroups(qw(Epi red Endo blue));

=cut

  my %names = @_;
  $settings{datanames} = \%names;
  @{$settings{datalist}} = sort keys %{$settings{datanames}};
}

#-------------------------------------------------
sub setGlyphs {
#-------------------------------------------------

=head1 B<setGlyphs>

  &setGlyphs( $glyph, $glyphSize, $images );

Sets the glyph type and size for the image set.
On input, $images must be an array of hashes.  This may have been generated using
the VH_CreateImageSet routine or another method.

The hash structure for each image object is as follows:

    <key>       <default> <description>
    {name}                Basename of image
    {filename}  {name}    Image file name including directory (no ext)
    {ext}       {.rgb}    Image file extension (including '.' if needed)
    {imgtype}   undef     Required if no extension, or if different type needed
    {imgloaded} 0         If image is loaded into cmgui
    {visible2D} 0         If image is visible in 2D window
    {glyph}     "sphere"  Glyph used for drawing data on this image
    {glyphsize} 3         Size of glyph used

=cut
  my $glyph = shift;
  my $glyphSize = shift;
  my $images = shift;

  if (! defined $glyph) {
	 die "setGlyphs called with an undefined \$glyph parameter.\nMaybe you used 'my' in a comfile without brackets";
  }

  if (! defined $glyphSize) {
	 die "setGlyphs called with an undefined \$glyphSize parameter.\nMaybe you used 'my' in a comfile without brackets";
  }

  foreach my $image ( @{$images} ) {
    $image->{glyph} = $glyph;
    $image->{glyphsize} = $glyphSize;
  }
  return 1;
}


#-------------------------------------------------
sub digitiseSlices  {
#-------------------------------------------------

=head1 B<digitiseSlices>

  &digitiseSlices( $images );

This routine takes care of the digitisation book keeping. If image viewing
is all that is needed then this is unnecessary although it is handy if you
want to go backwards and forwards through the image stack. 

=cut

  my $images = shift;
  
  unless (checkImageSet( $images )) {
    print "Invalid image set provided -- please fix\n";
    return;
  }
  loadImageSet( $images );

  CmGtk2->init();

  #Override the combobox style so we always get a list
  my $combobox_style = <<EOS
    style "combobox_style"
    {
      GtkComboBox::appears-as-list = 1
    }
    class "GtkComboBox" style "combobox_style"
EOS
	 ;
  Gtk2::Rc->parse_string($combobox_style);

  # Add a button to change between viewing and digitizing
  my $toggle = Gtk2::ToggleButton->new( " Digitise " );
  $toggle->signal_connect( "clicked" => \&toggle_signal, $images );
  $settings{toggle} = $toggle;

  # Add a list of all images 
  my $combo = Gtk2::ComboBox->new_text();

#  my $style = $combo->get_style();

  #$combo->set_property('appears-as-list' => 'true');
  foreach my $str (map{ $_->{name} . (( scalar keys %{$_->{datagroups}} ) ? "  *" : "")} @{$images}) { 
      $combo->append_text($str);
  }
  $combo->set_active(0);
  $settings{combosignalid} = $combo->signal_connect( "changed" => \&combo_signal, $images );
  $settings{combo} = $combo;

  
  # Add a list of all data groups
  # create combobox for datalist and show first in list as default.
  my $combo2 = Gtk2::ComboBox->new_text();
  foreach my $tempstr (@{$settings{datalist}} ) {
      $combo2->append_text($tempstr);
  }
  $combo2->set_active(0);
  $combo2->signal_connect( "changed" => \&data_signal, $images );
     
  # Add next/previous buttons
  my $next = new Gtk2::Button( "   Next   " );
  $next->signal_connect( "clicked" => \&next_signal, $images );
  my $prev = new Gtk2::Button( " Previous " );
  $prev->signal_connect( "clicked" => \&prev_signal, $images );  

  # Create a box to put all the widgets in
  my $box = new Gtk2::HBox();
  $box->set_spacing( 15 );
  $box->add( $prev );
  $box->add( $next );
  $box->add( $combo );
  $box->add( $combo2 );
  $box->add( $toggle );

  # Create a window to put the box in
  my $window = new Gtk2::Window( "toplevel" );
  $window->set_title( "Image Digitising" );
  $window->signal_connect( delete_event => sub { return 0;} );
  $window->set_border_width( 15 );
  $window->add( $box );
  $window->show_all();

  # View the images
  visibleoff( $images );
  update_view( $images );
  cmiss::cmiss("gfx modify win 1 image view_all");
  cmiss::cmiss("gfx modify win 2 image view_all") if $settings{view3d};


  # Start the Gtk2 main loop
  CmGtk2->init();
  CmGtk2->main();
}

#-------------------------------------------------
sub toggle_signal  {
#-------------------------------------------------

  my( $button, $images ) = @_;

  update_view( $images );
}

#-------------------------------------------------
sub combo_signal  {
#-------------------------------------------------

  my( $button, $images ) = @_;

  my $text = $button->get_active_text();
  my $i = 0;
  while ( $i < $#{$images} ) {
    last if $text =~ /${$images}[$i]->{name}/;
    $i ++;
  }
  $settings{current} = $i;

  update_view( $images );
}

#-------------------------------------------------
sub data_signal  {
#-------------------------------------------------

  my( $button, $images ) = @_;

  my $text = $button->get_active_text();
  my $i = 0;
  while ( $i < $#{$settings{datalist}} ) {
    last if $text eq ${$settings{datalist}}[$i];
    $i ++;
  }
  $settings{dataindex} = $i;

  update_view( $images );
}

#-------------------------------------------------
sub next_signal  {
#-------------------------------------------------

  my( $button, $images ) = @_;

  $settings{current}++;
  if( $settings{current} > $#{$images} ) {
    $settings{current} = 0;
  }

  update_view( $images );
}

#-------------------------------------------------
sub prev_signal  {
#-------------------------------------------------

  my( $button, $images ) = @_;

  $settings{current}--;
  if( $settings{current} < 0 ) {
    $settings{current} = $#{$images};
  }

  update_view( $images );
}

1;

__END__
