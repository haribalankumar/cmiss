package CmUtils::File::Ipdata;
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

=head1 CmUtils::File::Ipdata

Routines for reading and writing F<ipdata> files.

=head1 VERSION

0.5 (12 October 2000)
0.6 (08 August 2001)

=head1 CHANGES

version 0.6
- added necessary field validity check for writing Ipdata files JMB

=head1 SUBROUTINES

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 0.5;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readIpdata &writeIpdata);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readIpdata>

=over 4

=item PARAMETERS

(filename)

=item USAGE

  $group = readIpdata("newfile");

=item FUNCTION

Returns a NodeGroup read from F<filename> (or F<filename.ipdata>).  Creates a
single FieldSet relevant to an F<ipdata> file.  It has a I<coordinates> Field
with I<x>, I<y> and I<z> Components, and a I<weights> Field  with the same
Components.

=back

=cut

sub readIpdata {
  my ($filename) = shift;

  if ($filename ne '-' && $filename !~ /\.ip/) {
    $filename .= ".ipdata";
  }

  open IPFILE, "<$filename" or croak "Cannot open file $filename";

  my $group = CmUtils::Objects::NodeGroup->new();
  my $fieldset = CmUtils::Objects::FieldSet->new();
  $group->addFieldSet($fieldset);

  my $first=1;
  my $numvalues;
  while (<IPFILE>) {
    s/^\s+//;
    s/\s+$//;
    if ($first) {
      $first=0;
      $group->name($_);
      next;
    }
    my($name, @values) = split;
    $numvalues = scalar @values;
    
    my $node = CmUtils::Objects::Node->new();
    $node->fieldSet($fieldset);
    $node->name($name);
    $node->addValues(@values);
    if ($node->debug()) {
      carp "Created new node ". $node->name();
    }
    
    $group->addNode($node);
  }
  close IPFILE;

  my $index = 0;
  foreach my $fieldname ("coordinates", "weights") {
    my $field = CmUtils::Objects::Field->new();
    $field->name($fieldname);
    $field->type("coordinate");
    $field->coordsys("rectangular cartesian");
    foreach my $cname ((qw/x y z/,"c1".."c9")[0..($numvalues/2)-1]) {
      my $component = CmUtils::Objects::Component->new();
      $component->name($cname);
      $component->valIndex($index++);
      $field->addComponent($component);
    }
    $fieldset->addField($field);
  }
  
  $fieldset->numValues($index);

  return $group;
}

=head2 B<writeIpdata>

=over 4

=item PARAMETERS

(filename, group, options)

=item USAGE

  writeIpdata("newfile", $group, \%options);

=item FUNCTION

Writes a NodeGroup to an F<ipdata> file.  Options are passed through
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
will be printed if it exists, with weights of 1.0.  If two fields are
specified, the second one will be written as the weights.

=back

=cut

sub writeIpdata {
  my ($filename, $group, $options) = @_;
  
  if ($filename ne '-' && $filename !~ /\.ip/) {
    $filename .= ".ipdata";
  }

  my @nodes = $options->{nodes} 
    ? @{$options->{nodes}} 
    : $group->getNodes();

  open IPFILE, ">$filename" or croak "Cannot open file $filename";
  
  printf IPFILE "%s\n", $group->name();
  
  foreach my $node (@nodes) {
    printf IPFILE "%10s ", $node->name();
    my $fieldset = $node->fieldSet();
    my @fields = $options->{fields} 
      ? $fieldset->getFields(@{$options->{fields}})
      : $fieldset->getFields("coordinates");
      
    # Determine whether the field is valid for ipdata files  
    LINE : for my $i (0 .. $#fields) {
      if ($fields[$i]->coordsys() eq 'string') {
        if ($node->debug()) {
          carp "Ipdata file format cannot handle string values. 
            Deflating fields";
        }
        
        # Deflate the @fields array
        splice (@fields, $i, 1); 
        last LINE;
      }
    }
    if (@fields > 2) {
      # restrict to two fields at most
      @fields = @fields[0,1];
    }
    my $numvalues;
    foreach my $field (@fields) {
      $numvalues = scalar $field->getComponents();
      foreach my $component ($field->getComponents()) {
        printf IPFILE "%16.8f ", $node->value($component->valIndex());
      }
    }
    if (@fields == 1) {
      print IPFILE " 1.0" x $numvalues;
    }
    print IPFILE "\n";
  }

  close IPFILE;
}

1;

__END__
