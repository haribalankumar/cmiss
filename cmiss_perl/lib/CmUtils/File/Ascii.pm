package CmUtils::File::Ascii;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::File::Utils;

use CmUtils::Objects::NodeGroup;
use CmUtils::Objects::FieldSet;
use CmUtils::Objects::Field;
use CmUtils::Objects::Component;
use CmUtils::Objects::Node;

=head1 CmUtils::File::Ascii

Routines for reading and writing F<ascii> files.  An F<ascii> file contains
raw coordinate values, one node per line.  Nodes are not numbered. e.g.

  1.0 1.0 1.0
  2.0 1.0 0.0
  ...

F<Ascii> files may have an extension of F<.asc> or F<.txt>.

=head1 VERSION

1.1 (16 May 2002)

=head1 CHANGES

1.1 (16 May 2002) :
  Added F<.txt> extension.

1.0 (15 May 2002) : Created module for converting to/from raw ascii file. 
  Prompted by output from laser scanners.

=head1 SUBROUTINES

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 1.1;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readAscii &writeAscii);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readAscii>

=over 4

=item PARAMETERS

(filename)

=item USAGE

  $group = readAscii("newfile");

=item FUNCTION

Returns a NodeGroup read from F<filename> (or F<filename.asc>, or
F<filename.txt>).  Creates a single FieldSet relevant to an F<ascii> file.  It
has a I<coordinates> Field with I<x>, I<y> and I<z> Components.

=back

=cut

sub readAscii {
  my ($filename) = shift;

  if ($filename ne '-' && $filename !~ /\.(asc|txt)/) {
    $filename .= ".asc";
  }
  open ASCFILE, "<$filename" or croak "Cannot open file $filename";

  my $group = CmUtils::Objects::NodeGroup->new();
  $group->name("data");
  my $fieldset = CmUtils::Objects::FieldSet->new();
  $group->addFieldSet($fieldset);

  my $numvalues;
  my $name = 1;
  while (<ASCFILE>) {
    my (@values) = split;
    $numvalues = scalar @values;

    my $node = CmUtils::Objects::Node->new();
    $node->fieldSet($fieldset);
    $node->name($name++);
    $node->addValues(@values);
    carp("Created new node ". $node->name()) if $node->debug();
    $group->addNode($node);
  }
  close ASCFILE;

  my $index = 0;
  my $field = CmUtils::Objects::Field->new();
  $field->name("coordinates");
  $field->type("coordinate");
  $field->coordsys("rectangular cartesian");
  foreach my $cname ((qw/x y z/,"c1".."c9")[0..($numvalues)-1]) {
    my $component = CmUtils::Objects::Component->new();
    $component->name($cname);
    $component->valIndex($index++);
    $field->addComponent($component);
  }
  $fieldset->addField($field);
  $fieldset->numValues($index);

  return $group;
}

=head2 B<writeAscii>

=over 4

=item PARAMETERS

(filename, group, options)

=item USAGE

  writeAscii("newfile", $group, \%options);

=item FUNCTION

Writes a NodeGroup to an F<ascii> file.  Options are passed through
the I<options> hash to control the subset of the NodeGroup which will be
printed.

=item OPTIONS

=over 8

=item nodes

A array of node structures which will be output.  A number of routines are
available to get nodes from a group (e.g. B<getNodes>.  See CmUtils::Objects::NodeGroup).
If not set, all nodes will be printed.

=item fields

An field name which will be printed.  If not set, the "coordinates" field
will be printed if it exists.

=back

=cut

sub writeAscii {
  my ($filename, $group, $options) = @_;
  
  if ($filename ne '-' && $filename !~ /\.(asc|txt)/) {
    $filename .= ".asc";
  }

  my @nodes = $options->{nodes} 
    ? @{$options->{nodes}} 
    : $group->getNodes();

  open ASCFILE, ">$filename" or croak "Cannot open file $filename";
  
  foreach my $node (@nodes) {
    my $fieldset = $node->fieldSet();
    my @fields = $options->{fields} 
      ? $fieldset->getFields(@{$options->{fields}})
      : $fieldset->getFields("coordinates");
      
    # Determine whether the field is valid for ascii files  
    for my $i (0 .. $#fields) {
      if ($fields[$i]->coordsys() eq 'string') {
        if ($node->debug()) {
          carp "Ascii file format cannot handle string values. 
            Deflating fields";
        }
        # Remove field from the @fields array
        splice (@fields, $i, 1); 
        last;
      }
    }
    my $field = $fields[0];
    foreach my $component ($field->getComponents()) {
      printf ASCFILE "%16.8f ", $node->value($component->valIndex());
    }
    print ASCFILE "\n";
  }

  close ASCFILE;
}

1;

__END__
