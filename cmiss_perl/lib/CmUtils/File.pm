package CmUtils::File;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::File::Exnode;
use CmUtils::File::Ipnode;
use CmUtils::File::Ipdata;
use CmUtils::File::Ipfiel;
use CmUtils::File::Ipfibr;
use CmUtils::File::Ascii;

=head1 CmUtils::File

A collection of routines to read and write node-based files for CMISS and 
CMGUI.  All file features should be supported, including multiple
node-groups, derivatives and versions.  

=head1 USAGE

  use CmUtils::File qw(readExnode writeIpnode);

will make available the two routines.  The read routines return one or more
node-groups.

  use CmUtils::File qw(:Exnode :Ipnode);

which will make the read and write Exnode and Ipnode routines available.

  $group = readExnode("file.exnode");
  writeIpnode("newfile.ipnode", $group);

Generic routines are also available, which determine the file type from the
file extensions, or from a I<format> option.  These routines (B<cmRead> and
B<cmWrite>) are exported by default.

  use CmUtils::File;
  $group = cmRead('file.ipnode');
  $options->{format}='exdata';
  cmWrite('newfile', $group, $options);

=head1 TAGS

Available tags to use are

I<:all> for all routines

I<:Exnode> for readExnode and writeExnode

I<:Exdata> for readExdata and writeExdata

I<:Ipnode> for readIpnode and writeIpnode

I<:Ipdata> for readIpdata and writeIpdata

I<:Ipfiel> for readIpfiel and writeIpfiel

I<:Ipfibr> for readIpfibr and writeIpfibr

I<:Ascii>  for readAscii  and writeAscii

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

=head1 VERSION

0.6 (6 November 2000)

=cut

  $VERSION     = 0.6;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&cmRead &cmWrite);

=head1 SUBROUTINES

The following generic subroutines are exported by default:

    cmRead
    cmWrite

The following subroutines are imported from the File modules, and are
exported on request:

    readExnode
    writeExnode
    readExdata
    writeExdata
    readIpnode
    writeIpnode
    readIpdata
    writeIpdata    
    readIpfiel
    writeIpfiel
    readIpfibr
    writeIpfibr
    readAscii
    writeAscii

=cut

  @EXPORT_OK   = qw(
    &readExnode &writeExnode
    &readExdata &writeExdata
    &readIpnode &writeIpnode
    &readIpdata &writeIpdata    
    &readIpfiel &writeIpfiel
    &readIpfibr &writeIpfibr
    &readAscii  &writeAscii
    &cmRead     &cmWrite
  );
  %EXPORT_TAGS = (
    all    => [qw(&readExnode &writeExnode &readExdata &writeExdata 
                  &readIpnode &writeIpnode &readIpdata &writeIpdata
                  &readIpfiel &writeIpfiel &readIpfibr &writeIpfibr
                  &readAscii  &writeAscii  &cmRead     &cmWrite)],
    Exnode => [qw(&readExnode &writeExnode)],
    Exdata => [qw(&readExdata &writeExdata)],
    Ipnode => [qw(&readIpnode &writeIpnode)],
    Ipdata => [qw(&readIpdata &writeIpdata)],
    Ipfiel => [qw(&readIpfiel &writeIpfiel)],
    Ipfibr => [qw(&readIpfibr &writeIpfibr)],
    Ascii  => [qw(&readAscii  &writeAscii )],
  );
}

=head2 B<cmRead(filename,[options])>

Generic file reading routine for CMISS Ip and Ex format files.  File format
is determined from the I<options> hash if possible, otherwise it should be 
specified in the I<filename>.

=over 4

=item USAGE

  $group = cmRead('file.ipnode');

  $o->{format}='exdata';
  $g = cmRead('file', $o);

=item OPTIONS

=over 8

=item format

File format to read.

=back

=back

=cut

use File::Basename qw/fileparse/;

sub cmRead {
  my ($filename, $options) = @_;
  my ($name, $path);
  
  # Pick up file format from options
  my $format = delete $options->{format};
  # otherwise from file name
  unless ($format) {
    ($name, $path, $format) = fileparse($filename, '\..*')
      or croak("Input file format could not be determined for $filename");
  }
  for ($format) {
    /exnode/i && do{ goto &readExnode };
    /exdata/i && do{ goto &readExdata };
    /ipnode/i && do{ goto &readIpnode };
    /ipdata/i && do{ goto &readIpdata };
    /ipfiel/i && do{ goto &readIpfiel };
    /ipfibr/i && do{ goto &readIpfibr };
    /asci?i?/i && do{ goto &readAscii };
    /txt/i    && do{ goto &readAscii };
    croak("Unknown file format <$format>");
  }
}

=head2 B<cmWrite(filename,group,[options])>

Generic file writing routine for CMISS Ip and Ex format files.  File format
is determined from the I<options> hash if possible, otherwise it should be 
specified in the I<filename>.

=over 4

=item USAGE

  cmWrite('file.ipnode', $group);

  $o->{format}='exdata';
  cmWrite('file', $g, $o);

=item OPTIONS

=over 8

=item format

File format to read.

=item ...

Other options as defined by the file writing routines.

=back

=back

=cut

sub cmWrite {
  my ($filename, $group, $options) = @_;
  my ($name, $path);
  
  # Pick up file format from options
  my $format = delete $options->{format};
  # otherwise from file name
  unless ($format) {
    ($name, $path, $format) = fileparse($filename, '\..*')
      or croak("Output file format could not be determined for $filename");
  }
  for ($format) {
    /exnode/i && do{ goto &writeExnode };
    /exdata/i && do{ goto &writeExdata };
    /ipnode/i && do{ goto &writeIpnode };
    /ipdata/i && do{ goto &writeIpdata };
    /ipfiel/i && do{ goto &writeIpfiel };
    /ipfibr/i && do{ goto &writeIpfibr };
    /asci?i?/i && do{ goto &writeAscii };
    /txt/i    && do{ goto &writeAscii };
    croak("Unknown file format <$format>");
  }
}

=head1 TODO

ipxi?

=cut

1;

__END__
