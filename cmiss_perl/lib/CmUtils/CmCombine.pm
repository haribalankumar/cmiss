eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
##!/usr/local/perl5.6/bin/perl
package CmUtils::CmCombine;
require 5.006;

use strict;
use warnings;
BEGIN{push @INC, "$ENV{CMISS_ROOT}/cmiss_perl/lib"}

use CmUtils::Pod::Usage qw(pod2usage);
use Getopt::Long;
use File::Basename qw(fileparse);
use Carp;

use CmUtils::File;
use CmUtils::Objects::NodeGroup;

=head1 CmUtils::CmCombine

cmCombine - combines several CMISS format files into one

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

=head1 VERSION

1.0 (28 February 2001)

=cut

  $VERSION     = 1.0;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&cmCombine);
  @EXPORT_OK   = qw(&combineExdata &combineIpdata &combineExnode &combineIpnode);
  %EXPORT_TAGS = (
    all => [qw(&combineExdata &combineIpdata &combineExnode &combineIpnode)],
  );
}

my $format;

my $podoptions = {
  -exitval => -1,
  -verbose => 1,
  -input   => "$ENV{CMISS_ROOT}/cmiss_perl/lib/CmUtils/CmCombine.pm",
  -output  => \*STDERR,
};

unless (caller) {
  # if run as another program name, set format as appropriate
  for ($0) {
    /exdata/i && do { $format = "exdata" };
    /ipdata/i && do { $format = "ipdata" };
    /exnode/i && do { $format = "exnode" };
    /ipnode/i && do { $format = "ipnode" };
  }
  $podoptions->{-exitval} = 2;
  exit( cmCombine(@ARGV) );
}

=head1 USAGE

The routines in this package can be called in several ways.

1/  From the shell prompt, call as cmCombine, or one of the aliases for
specific filetypes. Prefix options with '-'. e.g.

    cmCombine -output file.ipdata old*ipdata        # creates file.ipdata

    In order for this to work cmCombine and other names must be
    set up as links to CmCombine.pm.

2/  From a Perl script, or from inside CMISS, load this module (with optional
flags): 

    use CmUtils::CmCombine;                           # loads cmCombine()
    use CmUtils::CmCombine qw(:all);                  # loads all routines
    use CmUtils::CmCombine qw(:combineExnode);        # combineExnode()

Options and filenames can then be specified in several ways.

a)  as a string equivalent to a command line

    combineExnode("-renumber -output file.exnode file1.exnode file2.exnode"); 
    # create file.exnode

b)  as options and filenames

    cmCombine("-output","file.ipdata",@ipfiles);

c)  as a hash array of options, followed by a list of files

    $options->{-renumber} = 1;
    $options->{-offset} = 1000;
    $options->{-reduce} = 5;
    $options->{-output} = "newfile.exnode";
    cmCombine($options, @files);

    $options = {-output => $newfile};
    combineIpdata($options, @ipfiles);

=cut

sub cmCombine {
  my %hashopts = ();

  local @ARGV;
  if (@_ && ref($_[0]) =~ /HASH/) {
    # if first argument is a hashref, grab options (removing dashes), 
    # and then filenames
    %hashopts = map { s/^-(\D)/$1/;$_ } %{shift(@_)};
    @ARGV = @_;
  } elsif (@_ == 1) {
    # if only one argument, split on whitespace
    @ARGV = split ' ',$_[0];
  } else {
    # grab all options and files together
    @ARGV = @_;
  }

=head1 SYNOPSIS of cmCombine

=over 12

=item B<cmCombine>

[B<-help>]
[B<-man>]
[B<-verbose>]
[B<-name>S< >I<GroupName>]
[B<-format>S< >I<format>]
[B<-renumber>S< >[B<-offset>S< >I<Offset>]]
[B<-output>S< >I<FileName>]
[B<-reduce>S< >I<Number>]
I<filename(s)>

=back

=cut
  my %options = ();
  my @opt_specs = (
  	"help",
  	"man",
    "verbose",
  	"name=s",
    "format=s",
  	"renumber",
    "output=s",
  	"offset=i",
    "reduce=i",
  );
  
=head1 OPTIONS AND ARGUMENTS

=over 8

=item B<-help>

=item B<-man>

=item B<-verbose>

=item B<-name>S< >I<GroupName>

New group name to give file.

=item B<-format>S< >I<format>

File format.  Automatically determined from filename if not specified.

=item B<-renumber>

Renumber nodes sequentially.

=item B<-offset>S< >I<Offset>

Start numbering from I<Offset>+1.

=item B<-output>S< >I<FileName>

New filename - otherwise prints to STDOUT.

=item B<-reduce>S< >I<Number>

Reduce points by a factor. (2=> every other point etc)

=item I<filename(s)>

File(s) to combine.  Must specify at least one file.

=back

=cut

  GetOptions(\%options, @opt_specs)  ||  return pod2usage($podoptions);
  %options = (%options, %hashopts);
  return pod2usage(%{$podoptions}, -verbose => 1) if ($options{help});
  return pod2usage(%{$podoptions}, -verbose => 2) if ($options{man});
  
  @ARGV = grep { -f } @ARGV;

  # Need at least one filename
  return pod2usage($podoptions) if (@ARGV == 0);

  my $offset  = defined($options{offset}) ? $options{offset}+1 : 1;
  my $name    = defined($options{name})   ? $options{name}     : 
    join ",", map { m!([^/]*)\.! } @ARGV;
  my $reduce  = defined($options{reduce}) ? $options{reduce}   : 1;
  $format     = exists $options{format} ? delete $options{format}  : $format;
  unless ($format) {
    my ($name, $path);
    ($name, $path, $format) = fileparse($ARGV[0], '\..*')
      or croak("Input file format could not be determined for $ARGV[0]");
    $format =~ s/\.//;
  }

  my $outfile = '-';
  if (defined($options{output})) {
    $outfile = $options{output};
    $outfile .= ".$format" unless $outfile =~ /\.(ip|ex)/;
    print STDERR "Creating $outfile\n";
  }
  
  my $group = CmUtils::Objects::NodeGroup->new();
  $group->name($name);
  my $fieldset;
  my $num = 0;
  my $opt;
  while (my $file = shift @ARGV) {
    $opt->{format} = $format;
    my $g = cmRead($file,$opt);
    foreach my $node ($g->getNodes()) {
      $node->name($offset) if $options{renumber};
      unless ($group->isNode($node->name()) || ($num++ % $reduce)) {
        $group->addNode($node);
        $offset++;

        my $same = 0;
        my $i = 0;
        foreach my $f ($group->getFieldSets()) {
          last if ($same = $node->fieldSet()->isSame($f));
          $i++;
        }
        unless ($same) {
          $i = $group->addFieldSet($node->fieldSet());
        }
        $node->fieldSet($group->getFieldSet($i));
      }
    }
  }
  $opt->{format} = $format;
  cmWrite($outfile, $group, $opt);
}

sub combineExdata {
  $format = "exdata";
  cmCombine(@_);
}

sub combineIpdata {
  $format = "ipdata";
  cmCombine(@_);
}

sub combineExnode {
  $format = "exnode";
  cmCombine(@_);
}

sub combineIpnode {
  $format = "ipnode";
  cmCombine(@_);
}

1;
__END__

=head1 DESCRIPTION

B<cmCombine> combines a set of CMISS format files into a single file.  All
files must be of the same format.
If the option I<renumber> is used, then nodes are sequentially renumbered, 
from I<offset> or zero.
A new group name can be given with the I<name> option.
The I<reduce> option allows a specified reduction in the data.

=head1 SUBROUTINES

The following subroutines are available if cmCombine is used as a module.
Parameters may include an I<options> hash, followed by a list of filenames.

  cmCombine($options, @files);

=head2 B<cmCombine>

=head2 B<combineExdata>

=head2 B<combineIpdata>

=head2 B<combineExnode>

=head2 B<combineIpnode>

=cut

