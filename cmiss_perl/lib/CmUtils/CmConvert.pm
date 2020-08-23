eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
##!/usr/local/perl5.6/bin/perl
package CmUtils::CmConvert;
require 5.006;

use strict;
use warnings;
BEGIN{push @INC, "$ENV{CMISS_ROOT}/cmiss_perl/lib"};

use CmUtils::Pod::Usage qw(pod2usage);
use CmUtils::Utils      qw(list_to_string);
use File::Basename      qw(fileparse);
use Getopt::Long;
use Carp;

use CmUtils::File;

=head1 CmUtils::CmConvert

cmConvert - converts files between CMISS formats

=cut
BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

=head1 VERSION

1.01 (23 March 2001)

=cut

  $VERSION     = 1.01;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&cmConvert);
  @EXPORT_OK   = qw(
    &ipExdata &exIpdata &ipExnode &exIpnode
    &toExdata &toIpdata &toExnode &toIpnode
  );
  %EXPORT_TAGS = (
    all => [qw(&ipExdata &exIpdata &ipExnode &exIpnode
               &toExdata &toIpdata &toExnode &toIpnode)],
  );
}

my ($fromformat, $toformat) = ();

my $podoptions = {
  -exitval => -1,
  -verbose => 1,
  -input   => "$ENV{CMISS_ROOT}/cmiss_perl/lib/CmUtils/CmConvert.pm",
  -output  => \*STDERR,
};

unless (caller) {
  # if run as another program name, set input and output formats as appropriate
  for ($0) {
    /ipexdata/i && do {
      $fromformat = "ipdata";
      $toformat   = "exdata";
    };
    /ipexnode/i && do {
      $fromformat = "ipnode";
      $toformat   = "exnode";
    };
    /exipdata/i && do {
      $fromformat = "exdata";
      $toformat   = "ipdata";
    };
    /exipnode/i && do {
      $fromformat = "exnode";
      $toformat   = "ipnode";
    };
    /toexdata/i && do { $toformat = "exdata" };
    /toexnode/i && do { $toformat = "exnode" };
    /toipdata/i && do { $toformat = "ipdata" };
    /toipnode/i && do { $toformat = "ipnode" };
  }
  $podoptions->{-exitval} = 2;
  exit ( cmConvert(@ARGV) );
}

=head1 USAGE

The routines in this package can be called in several ways.

1/  From the shell prompt, call as cmConvert, or one of the aliases for
specific conversions. Prefix options with '-'. e.g.

    cmConvert -toformat exdata -force file.ipdata   # creates file.exdata
    toIpnode file.exnode                            # creates file.ipnode
    cmConvert file.ipnode                           # creates file.exnode

	In order for this to work cmConvert and other names must be
	set up as links to CmConvert.pm.

2/  If data is piped in (or out) then the file formats must be explicitly
specified. e.g.

    cat file.ipnode | cmConvert -fromformat ipnode -toformat exnode > file.exnode

3/  From a Perl script, or from inside CMISS, load this module (with optional
flags): 

    use CmUtils::CmConvert;                           # loads cmConvert()
    use CmUtils::CmConvert qw(:all);                  # loads all routines
    use CmUtils::CmConvert qw(ipExnode toIpdata);     # ipExnode() & toIpdata()

Options and filenames can then be specified in several ways.

a)  as a string equivalent to a command line

    ipExnode("-force -unique file.ipnode");   # create file.exnode

b)  as options and filenames

    cmConvert("-toformat","ipdata",$file);

c)  as a hash array of options, followed by a list of files

    $options->{-force} = 1;
    $options->{-toformat} = "ipdata";
    cmConvert($options, @files);

    $options = {-force => 1, -unique => 1, x0 => 5.7};
    toIpdata($options, $exdatafile);
    
    %options = (-force => 1, -unique => 1, x0 => 5.7);
    toIpdata(\%options, "file.exdata");

=cut

sub cmConvert {
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

=head1 SYNOPSIS of cmConvert

=over 12

=item B<cmConvert>

[B<-help>]
[B<-man>]
[B<-verbose>]
[B<-force>]
[B<-offset>S< >I<offset>]
[B<-unique>S< >[B<-thousand>|B<-log>]]
[B<-x0>S< >I<x0>]
[B<-y0>S< >I<y0>]
[B<-z0>S< >I<z0>]
[B<-sort>]
[B<-fields>S< >I<fieldlist>]
[B<-allfields>]
[B<-nodes>S< >I<nodelist>]
[B<-fromformat>S< >I<format>]
[B<-toformat>S< >I<format>]
[I<filename(s)>]

=back

=cut
  my (%options);
  my (@opt_specs) = (
    "help",
    "man",
    "verbose",
    "force",
    "offset=i",
    "unique",
    "thousand",
    "log",
    "x0=i",
    "y0=i",
    "z0=i",
    "sort",
    "fields=s",
    "allfields",
    "nodes=s",
    "fromformat=s",
    "toformat=s",
  );

=head1 OPTIONS AND ARGUMENTS

=over 8

=item B<-help>

=item B<-man>

=item B<-verbose>

=item B<-force>

Force output file to be overwritten if it exists.

=item B<-offset>S< >I<offset>

Node number offset to begin numbering from.

=item B<-unique>

Causes all nodes for all files to have unique numbers.

=item B<-thousand>

Starts each new numbering from next multiple of 1000.

=item B<-log>

Starts each new numbering according to power of ten of previous set.

=item B<-x0> I<x0> B<-y0> I<y0> B<-z0> I<z0>

Translation of position.

=item B<-sort>

Will sort multiple files by decreasing size.

=item B<-allfields>

Write all fields (usually "coordinates" and "weights") to file.  Default
is to only write "coordinate" field.

=item B<-fields>S< >I<fieldlist>

Write the given field(s) only to the file. Default is to only write
"coordinate" field.

=item B<-nodes>S< >I<nodelist>

Select a subset of nodes to convert.  Any nodes lying in the range given will
be output if they exist.  Nodelist must be defined without spaces, e.g.

  -nodes 4,7,10..14,32

=item B<-fromformat>S< >I<format>

Input file format.  If not specified, it will be determined from the input 
filename.  Must be specified if file comes from STDIN.

=item B<-toformat>S< >I<format>

Output file format.  By default, is the corresponding alternate format (i.e.
F<exdata> from F<ipdata>).

=item I<filename(s)>

Input file(s) to convert.  Also accepts piped input.

=back

=cut

  GetOptions(\%options, @opt_specs)  ||  return pod2usage($podoptions);
  %options = (%options, %hashopts);
  return pod2usage(%{$podoptions}, -verbose => 1) if ($options{help});
  return pod2usage(%{$podoptions}, -verbose => 2) if ($options{man});
  
  # Need at least one filename, or have data piped from STDIN
  return pod2usage($podoptions) if ((@ARGV == 0) && (-t STDIN));
  @ARGV = ("-")  unless (@ARGV > 0);

  my $offset = exists $options{offset} ? $options{offset} : 0;
  my $x0     = exists $options{x0}     ? $options{x0}     : 0;
  my $y0     = exists $options{y0}     ? $options{y0}     : 0;
  my $z0     = exists $options{z0}     ? $options{z0}     : 0;
  
  $fromformat = exists $options{fromformat} ? $options{fromformat} : $fromformat;
  $toformat   = exists $options{toformat}   ? $options{toformat}   : $toformat;

  # sort files by decreasing size
  @ARGV = sort {-s $b <=> -s $a} @ARGV if (defined $options{sort});

  foreach my $file (@ARGV){
    my $ext;

    # read input file
    my $readopt;
    $readopt->{format} = $fromformat;
    my $group = cmRead($file,$readopt);

    # determine new filename
    my $newfile = '-';
    if ($file eq '-') {
      $ext = $toformat;
    } else {
      my ($name, $path, $ext) = fileparse ($file, qr/\..*/);
      $ext =~ s/^\.//;
      if ($toformat) {
        $ext = $toformat;
      } else {
        $ext =~ s/^ex/ip/ or $ext =~ s/^ip/ex/;
      }
      $newfile = $path . $name . "." . $ext;
      if ((! defined $options{force}) && (-e $newfile)) {
        warn "Output file $newfile already exists.  Skipping...\n";
        next;
      } else {
        print STDERR "Converting $file to \u$ext format\n";
        if ($options{verbose}) {
          $group->listNodes(\*STDERR);
          $group->listFields(\*STDERR);
        }
      }
    }

    # modify nodal values as required
    my $n;
    my @nodes = defined $options{nodes} ? eval $options{nodes} : ();
    foreach my $node ($group->getNodes(@nodes)) {
      $n = $node->name();
      $node->name($n + $offset);
      foreach my $field ($node->fieldSet()->getFields()) {
        foreach my $component ($field->getComponents()) {
          my $index = $component->valIndex();
          my $value = $node->value($index);
          for ($component->name()) {
            /^x$/ && do { $node->value($index, $value - $x0) };
            /^y$/ && do { $node->value($index, $value - $y0) };
            /^z$/ && do { $node->value($index, $value - $z0) };
          }
        }
      }
    }

    # write output file
    my $writeopt;
    if (defined $options{nodes}) {
      $writeopt->{nodes} = [ $group->getNodes(@nodes) ];
    }
    if ($options{allfields}) {
      $writeopt->{fields} = [ $group->getAllFieldNames() ];
    } elsif (my $flist = $options{fields}) {
      $flist =~ s/,/|/g;
      $writeopt->{fields} = [ grep { /$flist/i } $group->getAllFieldNames() ];
    } else {
      # look for fields called coordinate, coordinates, Coordinates etc
      $writeopt->{fields} = [ grep { /^coordinate/i } $group->getAllFieldNames() ];
    }
    unless (@{$writeopt->{fields}}) {
      croak "No fields found to write.\n";
    }
    $writeopt->{format} = $ext;
    if ($options{verbose}) {
      print STDERR "Writing fields: ",
        (join ", "=> @{$writeopt->{fields}}),
        "\n";
      print STDERR "Writing nodes : ",
        list_to_string($group->getNodeNames(@nodes)),
        "\n";
    }
    cmWrite($newfile, $group, $writeopt);
        
    # update incrementing variables
    if (defined $options{unique}) {
      $offset += $n;
      if (defined $options{log}) {
        my $power = 10.0**int((log($n)/log(10))+1);
        $offset = (int($offset/$power)+1)*$power;
      }
      if (defined $options{thousand}) {
        $offset = (int($offset/1000)+1)*1000;
      } 
    }
    
  }

}

=head1 DESCRIPTION

B<cmConvert> converts file between different CMISS formats, including F<ipnode>,
F<ipdata>, F<exnode>, F<exdata>, F<ipfiel> and F<ipfibr>.
Will also read from STDIN and write to STDOUT.  If multiple files 
are specified and I<sort> is set, then files will be sorted by decreasing size.
If an output file already exists, it will not be written, unless the I<force>
option is used.

File input and output formats are specified using the I<fromformat> and
I<toformat> options, or are determined from the input file(s).

If the variable I<offset> is set, then the numbering will 
begin from I<offset>, else it is the same as the input file.
If multiple files are parsed, the I<unique> option allows numbering to
increment for each file in sequence so that all written nodes for
all files have unique numbers.  The I<thousand> option causes this numbering
to jump to the next multiple of 1000.  The I<log> option instead jumps to 
next appropriate power of 10.

Can also set a translation using I<x0>, I<y0> and I<z0>.  Translates the nodes
so this specified point is at the origin.

If symbolic links are made to this program, it can also be run as one of the
following names, with appropriate input and output formats preset.

  ipExdata
  exIpdata
  ipExnode
  exIpnode
  toExdata
  toIpdata
  toExnode
  toIpnode

=cut

=head1 SUBROUTINES

=head2 B<cmConvert>

Used as a module, the B<cmConvert> routine is available by default, with
additional routines exportable if needed, individually, or together with
qualifier I<:all>.

=head2 B<ipExdata>

Converts from F<ipdata> to F<exdata> formats.

=cut
sub ipExdata {
  $fromformat = "ipdata";
  $toformat   = "exdata";
  cmConvert(@_);
}

=head2 B<ipExnode>

Converts from F<ipnode> to F<exnode> formats.

=cut
sub ipExnode {
  $fromformat = "ipnode";
  $toformat   = "exnode";
  cmConvert(@_);
}

=head2 B<exIpdata>

Converts from F<exdata> to F<ipdata> formats.

=cut
sub exIpdata {
  $fromformat = "exdata";
  $toformat   = "ipdata";
  cmConvert(@_);
}

=head2 B<exIpnode>

Converts from F<exnode> to F<ipnode> formats.

=cut
sub exIpnode {
  $fromformat = "exnode";
  $toformat   = "ipnode";
  cmConvert(@_);
}

=head2 B<toExdata>

Converts to F<exdata> formats.

=cut
sub toExdata {
  $toformat   = "exdata";
  cmConvert(@_);
}

=head2 B<toIpdata>

Converts to F<ipdata> formats.

=cut
sub toIpdata {
  $toformat   = "ipdata";
  cmConvert(@_);
}

=head2 B<toExnode>

Converts to F<exnode> formats.

=cut
sub toExnode {
  $toformat   = "exnode";
  cmConvert(@_);
}

=head2 B<toIpnode>

Converts to F<ipnode> formats.

=cut
sub toIpnode {
  $toformat   = "ipnode";
  cmConvert(@_);
}

=head1 CHANGES

1.01 - Added options for specifying fields to use.  Also verbose option.

1.0 - Release

0.7 - Added new Pod documentation system to give help in CMISS or perldl etc

0.6 - Added documentation and examples describing usage

=cut
1;
__END__
