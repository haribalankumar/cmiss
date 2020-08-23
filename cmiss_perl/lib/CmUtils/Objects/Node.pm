package CmUtils::Objects::Node;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::debug;
use CmUtils::Objects::Component;

=head1  CmUtils::Objects::Node

Defines a Node, which is a point with Field values defined by a FieldSet.

=head1 Node structure

=over 4

=item NAME

Name of the Node

=item FIELDSET

The FieldSet, used by this Node, which describes the nature of the values stored

=item VALUES

A list of values as defined by the FieldSet

=back

=head1 VERSION

0.6 (20 February 2001)
0.7 (08 August 2001)
0.8 (12 May 2004)

=head1 CHANGES

version 0.6:
- added values structure

version 0.7:
- added addStrings() method to handle string values JMB

version 0.8:
- added setValuesHash() method to allow values, versions and
derivatives to be added using a hash. Glenn Ramsey

=head1 METHODS

=cut

BEGIN {  
  our ($VERSION, @ISA);

  $VERSION     = 0.6;
  @ISA         = qw(CmUtils::debug);
}

=head2 B<new()>

Returns a new Node.

=cut

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self = {
    NAME     => undef,
    FIELDSET => undef,
    VALUES   => []
  };

  bless $self, $class;
  return $self;
}

=head2 B<copy()>

Returns an identical copy of an existing Node.

=cut

sub copy {
  my $self = shift;
  my $class = ref($self) || $self;

  my $new = $self->new();
  $new->name($self->name());
  $new->fieldSet($self->fieldSet());
  $new->setValues($self->getValues());

  return $new;
}

=head2 B<name([string/number])>

Optional argument is name for Node.  Returns the name of the Node.

=cut

sub name {
  my $self = shift;
  croak "usage: thing->name(name)"    unless @_ <= 1;
  if (@_) { $self->{NAME} = shift }
  return $self->{NAME};
}

=head2 B<fieldSet([FieldSet])>

Optional argument defines the FieldSet for Node.  Returns the FieldSet
of the Node.

=cut

sub fieldSet {
  my $self = shift;
  croak "usage: thing->fieldSet(fieldset)"    unless @_ <= 1;
  if (@_) { $self->{FIELDSET} = shift }
  return $self->{FIELDSET};
}

=head2 B<addValues(number(s))>

Adds a list of values to the Nodal value array.  Returns the last defined index
of the values.

=cut

sub addValues {
  use CmUtils::File::Utils qw/strtonum/;

  my $self = shift;
  # "usage: thing->addValues(values)"

  #push @{$self->{VALUES}}, strtonum($_) foreach @_;
  foreach (@_) {
    if( defined($_) ){
      push @{$self->{VALUES}}, strtonum($_);
    } else {
      carp "Attempt to add undefined values to node ".$self->name().", using 0 instead\n";
      push @{$self->{VALUES}}, 0;
    }
  }
  return $#{$self->{VALUES}};
}

=head2 B<addStrings(number(s))>

Adds a list of strings to the Nodal value array. Returns the last defined index
of the values.

=cut

sub addStrings {
  my $self = shift;
  push @{$self->{ VALUES }}, $_ foreach @_;
  return $#{$self->{ VALUES }};
}

=head2 B<clearValues()>

Clears the list of values of the Node.

=cut

sub clearValues {
  my $self = shift;
  # "usage: thing->getValues()"
  $self->{VALUES} = [];
}

=head2 B<setValues()>

Sets the list of values of the Node (clears, then adds the new values).

=cut

sub setValues {
  my $self = shift;
  $self->clearValues();
  return $self->addValues(@_);
}

=head2 B<getValues()>

Returns an array of the values of the Node.

=cut

sub getValues {
  my $self = shift;
  # "usage: thing->getValues()"
  return @{$self->{VALUES}};
}

=head2 B<getValuesHash()>

Returns the values of the Node in a hash format.
Hash of field names
  Hash of component names
    Array of versions
      Hash of
        Value
        Array of derivatives (if derivs are defined)

$node->{fieldname}{componentname}[versionnumber]{value}
$node->{fieldname}{componentname}[versionnumber]{derivatives}[derivnumber]
e.g.

  $node = $group->getNode(5)->getValuesHash();
  print $node->{coordinates}{x}[0]{value};
  @ver3derivs = @{$node->{coordinates}{z}[3]{derivatives}};


=cut

sub getValuesHash {
  my $self = shift;

  my @values = $self->getValues();
  my $valstruct;
  foreach my $f ($self->fieldSet()->getFields()) {
    my $field;
    foreach my $c ($f->getComponents()) {
      my $index = $c->valIndex();
      my $component;
      foreach my $version (1..$c->versions()) {
        my $vals;
        $vals->{value} = $values[$index];
        @{$vals->{derivatives}} = @values[$index+1..$index+$c->derivatives()]
          if ($c->derivatives() > 0);
        push @{$component}, $vals;
        $index += $c->derivatives()+1;
      }
      $field->{$c->name()} = $component;
    }
    $valstruct->{$f->name()} = $field;
  }
  return $valstruct;
}

=head2 B<setValuesHash(values, options, group)>

Sets all of the values of a Node from hash and updates its
fieldset if necessary.

This function can be used to add versions and derivatives to the Node.

It will fail if the fieldset in the values hash does not match
an existing fieldset. If you want to add versions or derivatives then
this should be explicitly specified in the options hash using the
newfieldset option.

An alternative implementation would be to have it automatically add
a new fieldset but that may cause a bug that is hard to track down
if that wasn't the intended behviour.

options:
to add a new fieldset 
newfieldset => 1 or 0

Warning: has not been tested extensively

=cut

sub setValuesHash {

  use Data::Dumper;

  my $self = shift;
  croak "usage: thing->setValuesHash(values, options, group)"    unless @_ <= 3;
  my $data = shift;
  my $options = shift;
  my $group = shift; # TODO: the nodegroup should be a member variable
  
  #
  # Determine if this FieldSet is the same as the Node's fieldset.
  #
  #my @coords = qw(x y z);
  my @derivs = qw(d/ds1 d/ds2 d2/ds1ds2 d/ds3 d2/ds1ds3 d2/ds2ds3 d3/ds1ds2ds3);
  
  # First construct a fieldSet from the hash and put the values
  # into an array.

  my $fs = CmUtils::Objects::FieldSet->new();
  my @values;
  my $index = 0;
  # $values contains a Hash of field names
  foreach my $fname (keys %{$data}) {
    my $f = CmUtils::Objects::Field->new();
    #TODO: these should be in the options hash
    $f->name("coordinates");
    $f->type("coordinate");
    $f->coordsys("rectangular cartesian");
    # $fname contains a Hash of component names
    #foreach my $cname (keys %{$data->{coordinates}}) {
    foreach my $cname (sort keys %{$data->{coordinates}}) {
      my $c = CmUtils::Objects::Component->new();
      my $numderivs = 0;
      $c->name($cname);
      $c->versions(scalar @{$data->{coordinates}{$cname}});
      if(exists($data->{coordinates}{$cname}[0]{derivatives})){
        # TODO: this assumes that all components have the same number of derivs
        $numderivs = scalar @{$data->{coordinates}{$cname}[0]{derivatives}};
        $c->derivatives($numderivs);
        if ($numderivs > 0) {
          $c->setDerivNames(@derivs[0..$numderivs-1]);
        }
       }
      $c->valIndex($index);
      # $cname contains an Array of versions      
      foreach my $version (0..@{$data->{coordinates}{$cname}}-1) {
        # $version contains a Hash of
        #         Value
        #         Array of derivatives (if derivs are defined)
        my $val = $data->{coordinates}{$cname}[$version]{value};
        push @values, $val;
        $index++;
        if(exists($data->{coordinates}{$cname}[$version]{derivatives})){
          my @derivs = @{$data->{coordinates}{$cname}[$version]{derivatives}};
          push @values, @derivs;
          $index += scalar @derivs;
        }
      }    
      $f->addComponent($c);
    }
    $fs->addField($f);
  }

  #now compare this to the node's field set
  my $isSame = $self->fieldSet()->isSame($fs);
 
  my $newfieldset = defined $options->{newfieldset} ? $options->{newfieldset} : 0;
 
  if($isSame){
    croak "Internal error in setValuesHash: wrong number of values"
      if (scalar $self->getValues() != scalar @values);
    $self->setValues(@values);
  } else {
    if($options->{newfieldset}>=1){
      $self->fieldSet($fs);
      $group->addFieldSet($fs);
      $self->setValues(@values);
    } else {
      carp "FieldSet does not match the Node's fieldSet,
             use the \"newfieldset\" option to override";
      print "--------Existing fieldset---------\n";
      print Dumper($self->fieldSet());
      print "----------New fieldset---------\n";
      print Dumper($fs);
      
    }
  }
}


=head2 B<value(index,[value])>

Setting and/or returning a single value of a Node by index.

=cut

sub value {
  my $self = shift;
  croak "usage: thing->value(number,[value])"    unless @_ <= 2;
  my $number = shift;
  if (@_) { $self->{VALUES}[$number] = shift }
  return $self->{VALUES}[$number];
}


=head2 B<merge(Node[,Versions])>

Merges the values of the given node with this node. If an optional
versions list is given then it only merges the versions in the list. 

=cut

sub merge {
  my $node0 = shift;
  croak "usage: thing->merge(Node[,Versions])"  unless @_ >= 1;
  my $node1 = shift;
  
  #my @versions = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_; 
  my $vref = shift;
  my @versions = @{$vref}; 
  
  #check that the fieldsets are the same
  my $field0 = $node0->fieldSet();
  my $field1 = $node1->fieldSet();
  
  if(!$field0->isSame($field1)){
    croak "Nodes to be merged must have the same field.\n"; 
  }
  
  #print "Node::merge merging Node ".$node0->name()."\n";
  my @values0 = $node0->getValues();
  my @values1 = $node1->getValues();

  foreach my $f ($node0->fieldSet()->getFields()) {
    foreach my $c ($f->getComponents()) {
      my $index = $c->valIndex();
      my $component;
            
      # if the versions array is defined then use it, else use all versions
      if(!@versions){
         carp "node ".$node0->name()." versions not defined using them all.";
         foreach (1..$c->versions()){
            push @versions, $_;
         }
      }
      
      #print Dumper( \@versions );
      
      # convert version array to hash so we can check if elements exist 
      my %vhash = ();
      for (@versions) { $vhash{$_} = 1 }
      
      foreach my $version (1..$c->versions()) {
      #foreach my $version (@versions) {
        if(exists $vhash{$version}){
          #print @versions." $c->versions() ".$c->versions()."\n";
          if($version <= $c->versions()){
            print "Node::merge() updating version $version\n";
            my $position0 = $values0[$index];
            my $position1 = $values1[$index];
            if($position0 != $position1){
              print "Warning: merging node ".$node0->name()." positions of component ".$c->name()." are not identical\n";
              print "Node0 position=$position0 Node1 position=$position1 \n";
            }
            my $deriv=1;
            for my $i ($index+1..$index+$c->derivatives())
            {
              if($values0[$i]!=0){
                print "Warning: merging node ".$node0->name()." derivative $deriv of component ".$c->name()." was not zero.\n";
                print "Node0 deriv=$values0[$i] Node1 derive=$values1[$i] \n";
              }
              $values0[$i]=$values1[$i];
              $deriv++; 
            }
          } else {
            croak "version doesn't exist (version > numversions)";
          }
        }
        $index += $c->derivatives()+1;
      } # foreach 
    }
  }
  $node0->setValues(@values0);
}



=head2 B<list([FH])>

List the contents of this Node to the filehandle I<FH> (default to STDOUT).
Set the "verbose" flag ($s->verbose(1)) to see extended information.

=cut

sub list {
  my $self = shift;
  my $FH = (@_) ? shift : \*STDOUT;

  my @values = $self->getValues();
  printf $FH <<EOF, $self->name(), scalar @values, $self;
   Node     %s     %X
     Number of values       %d
EOF
  foreach my $f ($self->fieldSet()->getFields()) {
    foreach my $c ($f->getComponents()) {
      my $index = $c->valIndex();
      if ($self->verbose()) {
        printf $FH "     %s\n", join ", " => 
          @values[$index..$index+($c->derivatives()+1)*$c->versions()-1];
      } else {
        printf $FH " %s\n", $values[$index];
      }
    }
  }
}



1;

__END__
