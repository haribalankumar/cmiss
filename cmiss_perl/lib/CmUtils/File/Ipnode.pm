package CmUtils::File::Ipnode;
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

=head1 CmUtils::File::Ipnode

Routines for reading and writing F<ipnode> files.

=head1 VERSION

0.5 (12 October 2000)
0.6 (08 August 2001)
0.7 (11 December 2003)

=head1 CHANGES

version 0.6
- added necessary field validity check for writing Ipnode files JMB 

version 0.7
- GBS moved CMISS::File:: routines to CmUtils::File:: and CmUtils::Objects::

=head1 SUBROUTINES

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 0.5;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readIpnode &writeIpnode);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readIpnode>

=over 4

=item PARAMETERS

(filename)

=item USAGE

  $group = readIpnode("newfile");

=item FUNCTION

Returns a NodeGroup read from F<filename> (or F<filename.ipnode>).  Creates
a new FieldSet for each type of versioning in the file, but only one FieldSet
if there are no versions.

Assigns a focus to the NodeGroup if it is present (for prolate coordinates).

Creates a Field named I<coordinates> and reads in the number of coordinates
(components), and derivatives.  Does not cope well with malformed files :-(

=back

=cut

sub readIpnode {
  my ($filename) = shift;

  if ($filename ne '-' && $filename !~ /\.ip/) {
    $filename .= ".ipnode";
  }

  open IPFILE, "<$filename" or croak "Cannot open file $filename";
  my $group = CmUtils::Objects::NodeGroup->new();
  my $genericfs = CmUtils::Objects::FieldSet->new();
  my $multiversion;
  my @multiversion;
  my $line;

  while (<IPFILE>) {
    s/^\s+//;
    s/\s+$//;
    SWITCH: { for ($_) {
      /^Heading:\s*(.*)$/ && do {
        $group->name($1);
        last SWITCH;
      };

      /^Specify the focus position [1.0]:\s*(.*)$/ && do {
        $group->property("Focus",$1);
        last SWITCH;
      };

      /^The number of nodes .*:\s*(\d+)$/ && do {
        my @coords = qw(x y z);
        my @derivs = qw(d/ds1 d/ds2 d2/ds1ds2 d/ds3 d2/ds1ds3 d2/ds2ds3 d3/ds1ds2ds3);
        
        my $field = CmUtils::Objects::Field->new();
        $field->name("coordinates");
        $field->type("coordinate");
        $field->coordsys("rectangular cartesian");
      
        $line = <IPFILE> =~ /Number of coordinates.*\]:(.*)/;
        my $numcomponents = strtonum($1);
        
        foreach (0..$numcomponents-1) {
          $line = <IPFILE> =~ /prompting for different versions.*\]\? (.*)/;
          $multiversion[$_] = ($1 =~ /[Yy]/);
        }
        $multiversion = grep { $_ > 0 } @multiversion;
        
        my $index = 0;
        foreach (0..$numcomponents-1) {
          $line = <IPFILE> =~ /The number of derivatives.*\]:(.*)/;
          my $numderivs = strtonum($1);
          my $component = CmUtils::Objects::Component->new();
          $component->name($coords[$_]);
          $component->valIndex($index++);
          $component->derivatives($numderivs);
          if ($numderivs) {
            $component->setDerivNames(@derivs[0..$numderivs-1]);
          }
          $field->addComponent($component);
          $index += $numderivs;
        }
        $genericfs->addField($field);
        $genericfs->numValues($index);

        $group->addFieldSet($genericfs) unless $multiversion;

        last SWITCH;
      };

      /^Node number .*\]:\s*(.*)$/ && do {
        my $node = CmUtils::Objects::Node->new();
        $node->name(strtonum($1));
        my $index = 0;
        my @mv = @multiversion;
      
        my $fs = $genericfs->copy();
        foreach my $field ($fs->getFields()) {
          foreach my $component ($field->getComponents()) {
            $component->valIndex($index);
            if (shift @mv) {
              $line = <IPFILE> =~ /number of versions.*\]:(.*)/;
              $component->versions(strtonum($1));
            }
            foreach (1..$component->versions()) {
              if ($component->versions() > 1) {
                <IPFILE>; # For version number ##:
              }
              $line = <IPFILE> =~ /coordinate is.*\]:(.*)/;
              $node->addValues(strtonum($1));
              $index++;
              foreach (1..$component->derivatives()) {
                $line = <IPFILE> =~ /derivative wrt.*\]:(.*)/;
                $node->addValues(strtonum($1));
                $index++;
              }
            }
          }
        }
        $fs->numValues($index);
        if ($node->debug()) {
          carp "Created new node ". $node->name();
        }
      

        if ($multiversion) {
          my $same = 0;
          my $i = 0;
          foreach my $fset ($group->getFieldSets()) {
            last if ($same = $fs->isSame($fset));
            $i++;
          }
          if ($same) {
            $node->fieldSet($group->getFieldSet($i));
          } else {
            $group->addFieldSet($fs);
            $node->fieldSet($fs);
          }
          
        } else {
          $node->fieldSet($genericfs);
        }

        $group->addNode($node);
        last SWITCH;
      };

    }}
  }

  close IPFILE;
  return $group;
}

=head2 B<writeIpnode>

=over 4

=item PARAMETERS

(filename, group, options)

=item USAGE

  $options{nodes} = $group->getNodes(10..20);
  writeIpnode($filename, $group, \%options);

=item FUNCTION

Writes a NodeGroup to an F<ipnode> file.  Options are passed through
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
will be printed if it exists.  Only one field may be printed to an F<ipnode>
file.

=item noDerivatives

Set this option to suppress printing of nodal derivatives.

=item noVersions

Set this option to suppress printing of multiple versions - only the first
will be printed.  (not yet implemented)

=back

=back

=cut

sub writeIpnode {
  my ($filename, $group, $options) = @_;

  if ($filename ne '-' && $filename !~ /\.ip/) {
    $filename .= ".ipnode";
  }

  if (ref($group) =~ /ARRAY/) {
    $group = $group->[0]; # can only print one group to an ipnode file
  }

  my @nodes = $options->{nodes} ? @{$options->{nodes}} : $group->getNodes();
  my $fieldname = $options->{fields} 
    ? (ref($options->{fields}) =~ /ARRAY/ 
      ? $options->{fields}->[0] # can only print one field to an ipnode file
      : $options->{fields})
    : "coordinates";
  
  # Determine whether the field name is valid for ipnode files
  LINE : foreach my $fs ($group->getFieldSets()) {
    if ($fs->getField($fieldname)->coordsys() eq 'string') {
      carp "Ipnode file format cannot handle string values. Reverting to default"; 
      $fieldname = 'coordinates';
      last LINE;
    }
  } 
    
  my $noDeriv = defined $options->{noDerivatives} ? $options->{noDerivatives} : 0;
  my $noVersions = defined $options->{noVersions} ? $options->{noVersions} : 0;
  my @indexlist = ();
  
  my @multiversion=();
  my @deriv = ('1','2','1 & 2','3','1 & 3','2 & 3','1, 2 & 3');
  
  foreach my $fs ($group->getFieldSets()) {
    if (my $field = $fs->getField($fieldname)) {
      my $nj = 1;
      foreach my $component ($field->getComponents()) {
        if ($component->versions() > 1) {
            $multiversion[$nj] = 1;
        }
        $nj++;
      }
    }
  }
  
  open IPFILE, ">$filename" or croak "Cannot open file $filename";

  printf IPFILE " CMISS Version 1.21 ipnode File Version 2\n Heading: %s\n\n", 
    $group->name();
  if ($group->property("Focus")) {
    printf IPFILE " Specify the focus position [1.0]: %f\n", $group->property("Focus");
  }
  printf IPFILE " The number of nodes is [    1]: %5d\n", scalar @nodes;
  
  my $first = 1;
  foreach my $node (@nodes) {
    my $field = $node->fieldSet()->getField($fieldname);
    if ($first) {
      $first = 0;
      my @components = $field->getComponents();
      printf IPFILE " Number of coordinates [3]: %d\n", scalar @components;
      foreach my $nj (1 .. scalar @components) {
        print IPFILE " Do you want prompting for different versions of nj=$nj [N]? ";
        print IPFILE $multiversion[$nj] ? "Y\n" : "N\n";
      }
      my $i = 1;
      foreach my $component (@components) {
        my $derivatives = $noDeriv ? 0 : $component->derivatives();
        print IPFILE " The number of derivatives for coordinate $i is [0]: $derivatives\n";
        $i++;
      }
    }

    printf IPFILE "\n Node number [%5d]: %5d\n", $node->name(), $node->name();
    my $nj = 1;
    
    foreach my $component ($field->getComponents()) {
      if ($multiversion[$nj]) {
        printf IPFILE " The number of versions for nj=%1d is [1]: %2d\n",
          $nj, $component->versions();
      }
      # $component->versions() is 1 if not multi-version => all OK
      foreach my $v (0..$component->versions()-1) {
        my $index = $component->valIndex() + ($component->derivatives()+1)*$v;
        if ($component->versions() > 1) {
          printf IPFILE " For version number%2d:\n", $v+1;
        }
        printf IPFILE " The Xj(%1d) coordinate is [ 0.00000E+00]: %25.14G\n",
          $nj, $node->value($index);
        unless ($noDeriv) {
          foreach my $d (0..$component->derivatives()-1) {
            $index++;
            my $dirn = 'direction';
            if ($deriv[$d] =~ /&/) {
              $dirn = 'directions';
            }
            printf IPFILE " The derivative wrt %s %s is [ 0.00000E+00]: %25.14G\n",
              $dirn, $deriv[$d], $node->value($index);
          }
        }
      }
      $nj++;
    }
  }

  close IPFILE;
}

1;

__END__
