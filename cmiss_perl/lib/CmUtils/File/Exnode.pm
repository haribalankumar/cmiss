package CmUtils::File::Exnode;
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

=head1 CmUtils::File::Exnode

Routines for reading and writing F<exnode>/F<exdata> files.

=head1 VERSION

0.5 (12 October 2000)
0.6 (08 August 2001)

=head1 CHANGES

version 0.6
- added the ability to read and write string values for Exnode files JMB


=head1 SUBROUTINES

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 0.5;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readExnode &writeExnode &readExdata &writeExdata);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readExnode>

=over 4

=item PARAMETERS

(filename)

=item USAGE

  @groups = readExnode($filename);
  $group = readExnode("newfile");

=item FUNCTION

Returns a list of NodeGroups read from F<filename> (or F<filename.exnode>).
If a scalar is assigned the return value, it will contain the first NodeGroup
read.

=back

=cut

sub readExnode ($) {
  my ($filename) = shift;
  
  # Read from an .exnode file ('-' means STDIN)
  if ($filename ne '-' && $filename !~ /\.ex/) {
    $filename .= ".exnode";
  }

  # Group will contain the list of ExNode groups
  my @Group;
  # Current group, fieldset and number of values to read
  my ($group, $curFieldSet, $numvalues);
  
  # Component field type array
  my @typecomponent;
  
  open EXFILE, "<$filename" or croak "Cannot open file $filename";
  
  while (<EXFILE>) {
    s/^\s+//;
    s/\s+$//;
    SWITCH: { for ($_) {
      /^Group name\s*:\s*(.*)$/ && do {
        # we have a new group -> create a new NodeGroup
        $group = CmUtils::Objects::NodeGroup->new();
        $group->name($1);
        push @Group, $group;
        if ($group->debug()) {
          carp "Created new group ". $group->name();
        }
        last SWITCH;
      };

      /^#Fields\s*=\s*(.*)$/ && do {
        my $numfields = strtonum($1);
        # Read in FieldSet header
        my $fieldset = CmUtils::Objects::FieldSet->new();
        my $totalvalues = 0; 
      
        foreach (1..$numfields) {
      	  my $field = CmUtils::Objects::Field->new();
          my $line;
          defined ($line = <EXFILE>) or croak "Error reading file";
          # Remove all whitespace around commas before processing
          $line =~ s/\s*,\s*/,/g;
          if($line =~ /\)\s*(\w+),(\w+),(.+),#Components\s*=\s*(\d+)/) {
            $field->name($1);
            $field->type($2);
            $field->coordsys($3);
            my $numcomponents = $4;
            foreach (1..$numcomponents) {
              
              my $component = CmUtils::Objects::Component->new();
              defined ($line = <EXFILE>) or croak "Error reading file";
              # Remove all whitespace before processing
              $line =~ s/\s//g;
              $line =~ /^(\w+)\./;
              $component->name($1);
              if ($line =~ /Valueindex=(\d+)/) { $component->valIndex($1-1) };
              if ($line =~ /#Versions=(\d+)/) { $component->versions($1) };
              if ($line =~ /#Derivatives=(\d+)/) { $component->derivatives($1) };
              if ($component->derivatives()) {
                my ($dnames) = ($line =~ /\((.*)\)/);
                my @dnames = split ',', $dnames;
                $component->setDerivNames(@dnames);
              }
              
              if ($component->debug()) {
                carp "Created new component ". $component->name();
              }
              my $numcompvalues = ($component->derivatives()+1) * $component->versions();
              $totalvalues += $numcompvalues;
              $field->addComponent($component);

              # Associate the component with the field type
              push (@typecomponent, $field->coordsys()) foreach 1 .. $numcompvalues; 
            }
          }
          if ($field->debug()) {
            carp "Created new field ". $field->name();
          }
          $fieldset->addField($field);
        }
        $fieldset->numValues($totalvalues);

        # Only add FieldSet to Group if it doesn't already exist
        my $same = 0;
        my $i = 0;
        foreach my $f ($group->getFieldSets()) {
          last if ($same = $fieldset->isSame($f));
          $i++;
        }
        if ($same && $fieldset->debug()) {
          carp "Fieldset already exists in this group - not adding";
        }
        unless ($same) {
          $i = $group->addFieldSet($fieldset);
        }

        # Set current FieldSet and number of values to read
        $curFieldSet = $group->getFieldSet($i);
        $numvalues = $curFieldSet->numValues();
        last SWITCH;
      };

      /^Node\s*:\s*(.*)$/ && do {
        # Read a new node into the current Group with the current FieldSet
        my $node = CmUtils::Objects::Node->new();        
        $node->name($1);
        $node->fieldSet($curFieldSet);
        while ((scalar $node->getValues())<$numvalues) {
          defined (my $line = <EXFILE>) or croak "Error reading file";
          
          my @line = split ' ', $line;
          for my $i (0 .. $#line) {
          
            # Determine the field associated with the current value index
            if ($typecomponent[scalar $node->getValues()] eq 'string') {
          
              # The node value is defined as a string
              $node->addStrings($line[$i]);
            } else {
          
              # The node value is defined as a number
              $node->addValues($line[$i]);
            }
          }  
        }
        if ($node->debug()) {
          carp "Created new node ". $node->name();
        }

        $group->addNode($node);
        last SWITCH;
      };

    }}
  }
  close EXFILE;
  
  return wantarray ? @Group : $Group[0];
}

=head2 B<readExdata>

=over 4

=item PARAMETERS

(filename)

=item USAGE

  @groups = readExdata($filename);

=item FUNCTION

Returns a list of NodeGroups read from F<filename> (or F<filename.exdata>).
If a scalar is assigned the return value, it will contain the first NodeGroup
read.  Calls B<readExnode>.

=back

=cut

sub readExdata ($) {
  my ($filename) = shift;
  
  if ($filename ne '-' && $filename !~ /\.ex/) {
    $filename .= ".exdata";
  }

  return readExnode ($filename);
}

=head2 B<writeExnode>

=over 4

=item PARAMETERS

(filename, group(s), options)

=item USAGE

  $options{nodes} = $group->getNodes(10..20);
  writeExnode($filename, $group, \%options);

=item FUNCTION

Writes one or more NodeGroups to an F<exnode> file.  Options are passed through
the I<options> hash to control the subset of the NodeGroups which will be
printed.

=item OPTIONS

=over 8

=item nodes

A array of node structures which will be output.  A number of routines are
available to get nodes from a group (e.g. B<getNodes>.  See CmUtils::Objects::NodeGroup).
If not set, all nodes will be printed.

=item fields

An array of field names which will be printed.  If not set, all fields will be
printed.

=item noDerivatives

Set this option to True to suppress printing of nodal derivatives.

=item noVersions

Set this option to suppress printing of multiple versions - only the first
will be printed.  (not yet implemented)

=back

=back

=cut

sub writeExnode {
  my ($filename, $groups, $options) = @_;

  # make sure filename is 'exnode' if already specified
  if ($filename ne '-' && $filename !~ /\.ex/) {
    $filename .= ".exnode";
  }

  # list of groups to write out - may be one, or an array
  my @group;
  if (ref($groups) =~ /ARRAY/) {
    @group = @{$groups};
  } else {
    @group = ($groups);
  }
  my $noDeriv = $options->{noDerivatives} ? $options->{noDerivatives} : 0;

  open EXFILE, ">$filename" or croak "Cannot open file $filename";

  foreach my $group (@group) {  
    my @nodes = $options->{nodes} ? @{$options->{nodes}} : $group->getNodes();
          
    printf EXFILE "Group name: %s\n", $group->name();
    
    my $fieldset = '';
    my @indexlist = ();

    foreach my $node (@nodes) {
      if ($node->fieldSet() ne $fieldset) {
        # write a new FieldSet if different to that for the previous node
        $fieldset = $node->fieldSet();
        my @fields = $options->{fields} 
          ? $fieldset->getFields(@{$options->{fields}})
          : $fieldset->getFields();
              
        printf EXFILE "#Fields=%d\n", scalar @fields;
        my $num=1;
        my $lastindex = 1;
        @indexlist = ();
      
        foreach my $field (@fields) {
          my @components = $field->getComponents();
          printf EXFILE "%d) %s, %s, %s, #Components=%d\n", 
            $num, $field->name(), $field->type(), $field->coordsys(), 
            scalar @components;
          foreach my $component (@components) {
            my $derivatives = $noDeriv ? 0 : $component->derivatives();
            printf EXFILE " %s.  Value index=%d, #Derivatives=%d", 
              $component->name(), $lastindex, $derivatives;
            
            if ($derivatives > 0) {
              print EXFILE " (";
              print EXFILE join (",", $component->getDerivNames());
              print EXFILE ")";
            }
            if ($component->versions() > 1) {
              printf EXFILE ", #Versions=%d", $component->versions();
            }
            print EXFILE "\n";
          
            my $i = $component->valIndex();
            for (1..$component->versions()) {
              my @line = ();
              push @line, $i;
              $lastindex += 1;
              if ($derivatives > 0) {
                push @line, ($i+1..$i+$derivatives);
                $lastindex += $derivatives;
              }
              $i+=$component->derivatives()+1;
              push @indexlist, \@line;
            }
          }

          $num++;
        }
      }
      # print out the node
      printf EXFILE "Node: %6s\n", $node->name();
      my $numprint = 0;
      foreach my $line (@indexlist) {
        foreach (@{$line}) {
        
          # Determine whether we are dealing with a string or a number
          if ($node->value($_) =~ /\D/) {
            printf EXFILE "%s ", $node->value($_);
          } else {
            printf EXFILE "%25.16G ", $node->value($_);
          }  
        }
        $numprint += scalar @{$line};
        print EXFILE "\n" unless scalar @{$line} == 1;
      }
      print EXFILE "\n" if $numprint = scalar @indexlist;
    }
  }
  
  close EXFILE;
}

=head2 B<writeExdata>

=over 4

=item PARAMETERS

(filename, group(s), options)

=item USAGE

  $options{nodes} = $group->getNodes(10..20);
  writeExdata($filename, $group, \%options);

=item FUNCTION

Calls B<writeExnode>, but writes to an F<exdata> file.

=back

=cut

sub writeExdata {
  my ($filename) = shift;
  
  if ($filename ne '-' && $filename !~ /\.ex/) {
    $filename .= ".exdata";
  }

  return writeExnode ($filename, @_);
}

1;

__END__
