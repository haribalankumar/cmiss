package CmUtils::Objects::NodeGroup;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::debug;

=head1 CmUtils::Objects::NodeGroup

The NodeGroup class defines structures and methods for creating and
manipulating groups of nodes.  The NodeGroup structure has 4 fields:

=over 4

=item NAME

Name of the NodeGroup

=item FOCUS

Value of the focus for NodeGroups based on prolate-spheroidal coordinates

=item FIELDSETS

A set of FieldSets which describe the Nodes in this NodeGroup

=item NODES

The set of nodes belonging to this NodeGroup

=back

=cut

BEGIN {  
  our ($VERSION, @ISA);

  $VERSION     = 0.6;
  @ISA         = qw(CmUtils::debug);
}

=head1 VERSION

0.6 (7 September 2001)

=head1 CHANGES

0.6 GBS: Added B<distance> method to compute distance between two nodes.

=head1 METHODS

=head2 B<new()>

Returns a new NodeGroup.

=cut

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self = {
    NAME      => '',
    PROPERTY  => {},
    FIELDSETS => [],
    NODES     => {},
  };
  bless $self, $class;
  return $self;
}

=head2 B<name([string/number])>

Optional argument is name for NodeGroup.  Returns the name of the NodeGroup.

=cut

sub name {
  my $self = shift;
  croak "usage: thing->name(name)" unless @_ <= 1;
  if (@_) { $self->{NAME} = shift }
  return $self->{NAME};
}

=head2 B<property(name,[value])>

Sets a problem-specific/dependent property (named NAME) of the NodeGroup.
Optional argument is the value of the property.  Returns the value.

=cut

sub property {
  my $self = shift;
  my $name = shift;
  
  if (@_) { $self->{PROPERTY}{$name} = shift }
  return $self->{PROPERTY}{$name};
}

=head2 B<getPropertyNames([Name(s)])>

Returns the sorted names of all (or NAMED) Properties defined for a NodeGroup.

=cut

sub getPropertyNames {
  my $self = shift;

  return (@_) 
    ? grep { exists $self->{PROPERTY}{$_} } @_ 
    : sort keys %{$self->{PROPERTY}};
}

=head2 B<addFieldSet(FieldSet)>

Adds a FieldSet to the list of FieldSets of this NodeGroup.  Returns the 
index of the last FieldSet.

=cut

sub addFieldSet {
  my $self = shift;
  croak "usage: thing->addFieldSet(fs)" unless @_ == 1;
  my $field = shift;
  push @{$self->{FIELDSETS}}, $field;
  if ($self->debug()) {
    carp "Added fieldset to group ". $self->name();
  }
  return $#{$self->{FIELDSETS}};
}

=head2 B<getFieldSets()>

Returns the list of defined FieldSets.

=cut

sub getFieldSets {
  my $self = shift;
  croak "usage: thing->getFieldSets()" unless @_ == 0;
  return @{$self->{FIELDSETS}};
}

=head2 B<getFieldSet(number)>

Returns a specific FieldSet by number.  FieldSets do not have names :-(

=cut

sub getFieldSet {
  my $self = shift;
  croak "usage: thing->getFieldSet(num)" unless @_ == 1;
  return $self->{FIELDSETS}[shift];
}

=head2 B<getAllFieldNames([FieldSet(s)])>

Returns the sorted names of all Fields in DEFINED or, by default, 
ALL FieldSets of a NodeGroup.

=cut

sub getAllFieldNames {
  my $self = shift;

  my @fieldsets = @_ ? @_ : @{$self->{FIELDSETS}};
  my %fields;
  foreach my $f (@fieldsets) {
    @fields{$f->getFieldNames()} = ();
  }
  return sort keys %fields;
}

=head2 B<addNode(Node)>

Adds a Node to the list of Nodes in the NodeGroup.  Node will NOT be added
(replaced) if it already exists.  Returns the total number of defined Nodes
in the NodeGroup.

=cut

sub addNode {
  my $self = shift;
  croak "usage: thing->addNode(node)" unless @_ == 1;
  my $node = shift;
  unless (exists $self->{NODES}{$node->name()}) {
    $self->{NODES}{$node->name()} = $node;
    my $fs = $node->fieldSet();
    my $same = 0;
    foreach my $fset ($self->getFieldSets()) {
      last if ($same = $fs->isSame($fset));
    }
    if (!$same) {
      $self->addFieldSet($fs);
    }    
  } else {    
    carp "Warning: Node ".$node->name() ." already exists in NodeGroup '". $self->name()."'";
  }
  return $self->numberOfNodes();
}

=head2 B<addNodes(Node(s))>

Adds one or more Nodes to the list of Nodes in the NodeGroup.  Nodes will 
NOT be added (replaced) if they already exist.  Returns the total number of 
defined Nodes in the NodeGroup.

=cut

sub addNodes {
  my $self = shift;
  croak "usage: thing->addNodes(nodes)" unless @_ >= 1;
  $self->addNode($_) foreach @_;
  return $self->numberOfNodes();
}

=head2 B<updateNode(Node)>

Updates the Node in the NodeGroup.  Node will NOT be changed unless it already
exists.  Returns the total number of defined Nodes in the NodeGroup.

=cut

sub updateNode {
  my $self = shift;
  croak "usage: thing->updateNode(node)" unless @_ == 1;
  my $node = shift;
  if (exists $self->{NODES}{$node->name()}) {
    delete $self->{NODES}{$node->name()};
    $self->{NODES}{$node->name()} = $node;
  } else {
    carp $node->name() ." does not exist in ". $self->name();
  }
  return $self->numberOfNodes();
}

=head2 B<mergeNode(Node)>

Merges the nodes by copying the non-zero derivatives of the given Node
into the existing node. 

=cut

sub mergeNode {
  my $self = shift;
  croak "usage: thing->mergeNode(Node[,versions])" unless @_ >= 1;
  my $newnode = shift;
  my $vref = shift;
  my @versions = @{$vref};
  #print "NodeGroup::mergeNode()\n";
  #print Dumper(\@versions);
  if (exists $self->{NODES}{$newnode->name()}) {
    my $thisnode = $self->{NODES}{$newnode->name()};
    $thisnode->merge($newnode, $vref);    
  } else {
    carp $newnode->name() ." does not exist in ". $self->name();
  }
  return $self->numberOfNodes();
}

=head2 B<deleteNode(NodeName)>

Deletes Node named NodeName in the NodeGroup.  Returns the total number of 
defined Nodes in the NodeGroup.

=cut

sub deleteNode {
  my $self = shift;
  croak "usage: thing->deleteNode(name)" unless @_ == 1;
  my $name = shift;
  if (exists $self->{NODES}{$name}) {
    delete $self->{NODES}{$name};
  }
  return $self->numberOfNodes();
}

=head2 B<deleteNodes(NodeName(s))>

Deletes one or more Nodes named NodeName(s) in the NodeGroup.
Returns the total number of defined Nodes in the NodeGroup.

=cut

sub deleteNodes {
  my $self = shift;
  croak "usage: thing->deleteNodes(name(s))" unless @_ >= 1;
  my @names = $self->getNodeNames(@_);
  $self->deleteNode($_) foreach @names;
  return $self->numberOfNodes();
}

=head2 B<isNode(NodeName)>

Returns True if the Node named NodeName is in the NodeGroup, otherwise returns
False.

=cut

sub isNode {
  my $self = shift;
  croak "usage: thing->isNode(name)" unless @_ == 1;
  my $name = shift;
  return exists $self->{NODES}{$name};
}

=head2 B<getNodeNames([NodeName(s)])>

If NodeName(s) are specified, returns a list of names of all Nodes with these 
names that exist in the NodeGroup, in the order they are specified.  If no 
argument is given, returns a list of names of ALL Nodes, sorted in 
numerical/alphabetical order.

=cut

sub getNodeNames {
  my $self = shift;
  return (@_) 
    ? grep { exists $self->{NODES}{$_} } @_ 
    : sort {$a <=> $b || $a cmp $b} keys %{$self->{NODES}};
}

=head2 B<getNodes([NodeName(s)])>

If NodeName(s) are specified, returns a list of all Nodes with these names that
exist in the NodeGroup, in the order they are specified.  If no argument is 
given, returns a list of ALL Nodes, sorted in numerical/alphabetical order.

=cut

sub getNodes {
  my $self = shift;
  my @keys = $self->getNodeNames(@_);

  my @nodes = ();
  push @nodes, $self->{NODES}{$_} foreach @keys;

  return @nodes;
}

=head2 B<getNodesByFieldSet([NodeName(s)])>

If NodeName(s) are specified, returns a list of all Nodes with these names that
exist in the NodeGroup.  If no argument is given, returns a list of ALL Nodes.
Returned Nodes are sorted by FieldSet, then by name
(numerically/alphabetically).

=cut

sub getNodesByFieldSet {
  my $self = shift;
  my @nodes = $self->getNodes(@_);
  return 
    map {$_->[0]} 
    sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] || $a->[2] cmp $b->[2]} 
    map {[$_, $_->fieldSet(), $_->name()]} 
    @nodes;
}

=head2 B<getNode(NodeName)>

Returns the Node named NodeName from the NodeGroup.

=cut

sub getNode {
  my $self = shift;
  croak "usage: thing->getNode(name)" unless @_ == 1;
  my $name = shift;
  return $self->{NODES}{$name};
}

=head2 B<numberOfNodes()>

Returns the number of Nodes defined in the NodeGroup.

=cut

sub numberOfNodes {
  my $self = shift;
  return scalar keys %{$self->{NODES}};
}

=head2 B<maximumNodeNumber()>

Returns the number of Nodes defined in the NodeGroup.

=cut

sub maximumNodeNumber {
  my $self = shift;
  my $max = 0;
  do { $max = $_ if $_ > $max } foreach keys %{$self->{NODES}};
  return $max;
}

=head2 B<listFields([FH])>

List the Field information of this NodeGroup to the filehandle I<FH>
(default to STDOUT).  Set the "verbose" flag ($s->verbose(1)) to see
extended information.

=cut

sub listFields {
  my $self = shift;
  my $FH = (@_) ? shift : \*STDOUT;

  printf $FH " NodeGroup                %s\n", $self->name();
  if ($self->verbose()) {
    if (my @properties = $self->getPropertyNames()) {
      print $FH "   User Defined Properties:\n";
      foreach my $p (@properties) {
        printf $FH "     %-20s = %s\n", $p, $self->property($p);
      }
    }
  }
  foreach my $fieldset ($self->getFieldSets()) {
    $fieldset->list($FH);
  }
}

=head2 B<listNodes([FH])>

List the Node information of this NodeGroup to the filehandle I<FH>
(default to STDOUT).  Set the "verbose" flag ($s->verbose(1)) to see
extended information.

=cut

sub listNodes {
  use CmUtils::File::Utils qw(strtonum);
  use CmUtils::Utils qw(list_to_string);
  
  my $self = shift;
  my $FH = (@_) ? shift : \*STDOUT;

  printf $FH " NodeGroup                %s\n", $self->name();
  printf $FH "   Number of Nodes          %d\n", $self->numberOfNodes();
  if ($self->numberOfNodes()) {
    my @nodes = $self->getNodeNames();
    if (strtonum($nodes[0])) { #numeric node names
      printf $FH "   Nodes                    %s\n", list_to_string(@nodes);
    }
    if ($self->verbose()) {
      foreach my $node ($self->getNodes()) {
        $node->list($FH);
      }
    }
  }
}

=head2 B<mergeGroup(Group[,Options])>

Merges another Group into this NodeGroup.  By default, nodes are added as long
as they do not already exist.
The Options hashref can alter this by either
forcing an overwrite of existing nodes (key "overwrite"), or assigning a new
number to the new nodes (key "newnumber").

  $options{newnumber} = 1;
  $group->mergeGroup($another_group,$options);
  
The "merge" key specifies that only the derivatives will be overwritten and if a hash
of versions is given in the "versions" key then only the specfied versions will
be overwritten.

  $options{merge} = 1;
  $options{versions}{10} = [1,3]; # for node 10 overwrite versions 1 and 3
  $group->mergeGroup($another_group,$options);  

=cut

sub mergeGroup {

  use Data::Dumper;
  my $self = shift;
  croak "usage: thing->mergeGroup(group[,options])" unless @_ >= 1;
  my ($newgroup, $options) = @_;
  croak "invalid new Group specified" unless ref($newgroup) =~ /CmUtils::Objects::NodeGroup/;
  my %options = (newnumber => 0, overwrite => 0);
  if ($options && ref($options) =~ /HASH/) {
    %options = (%options, %{$options});
  }
  #print Dumper( $options{versions} );
  my $number = $options{newnumber} ? $self->maxNodeNumber() : 0;
  foreach my $node ($newgroup->getNodes()) {
    my $fs = $node->fieldSet();
    my $same = 0;
    my $i = 0;
    foreach my $fset ($self->getFieldSets()) {
      last if ($same = $fs->isSame($fset));
      $i++;
    }
    if ($same) {
      $node->fieldSet($self->getFieldSet($i));
    } else {
      $self->addFieldSet($fs);
      $node->fieldSet($fs);
    }
    if ($options{newnumber}) {
      $number++;
      $node->name($number);
      $self->addNode($node);
    } else {
      if($self->isNode($node->name())) {
        if($options{overwrite}) {
          $self->updateNode($node);
        } elsif($options{merge}) {
          #print "mergeGroup() node: ".$node->name()."\n";
          #print Dumper( $options{versions}{$node->name()} );
          my @versions = @{$options{versions}{$node->name()}};
          #print Dumper( \@versions );
          $self->mergeNode($node, \@versions);
        }
        #else ignore it
        #TODO: probably should warn the user that it was ignored
      } else {
        $self->addNode($node);
      }
    }
  }
}

=head1 Transformation METHODS

The following routines define useful transformations of existing NodeGroups. 
By default each routine processes the whole NodeGroup, but a list of NodeNames
may be specified to restrict processing to just a portion of it.  NodeNames can
be specified as either an array or a list, and are always the last argument(s)
to the method call.

Where operations apply to coordinate values, they are restricted to fields
matching C</coordinate/i> with components matching C</[xyz]/i>.

=head2 B<reduce(Factor[,NodeName(s)])>

Deletes nodes so that only 1/Factor of the original number remain.  For
example:

  $group->reduce(2);

keeps only every other node.  Useful for subsampling a dataset.

=cut

sub reduce {
  croak "usage: group->reduce(factor[,nodes])" unless @_ >= 2;
  my $group  = shift;
  my $factor = shift;
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  my $index  = 0;
  foreach my $node ($group->getNodeNames(@nodes)) {
    if ($index % $factor) {
      $group->deleteNode($node);
    }
    $index++;
  }
}

=head2 B<renumber([Offset[,Increment[,NodeNames]]])>

Renumbers the NodeGroup starting from Offset (default 1) by an Increment
(default 1).

B<Note:> Nodes may be lost if a portion of the group is renumbered to the same
as other existing nodes.

=cut

sub renumber {
  croak "usage: group->renumber([offset[,increment[,nodes]]])" unless @_ >= 1;
  my $group  = shift;
  my $offset = shift || 1;
  my $increment = shift || 1;
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  @nodes = $group->getNodeNames() unless @nodes;
  my @nodelist = $group->getNodes(@nodes);
  $group->deleteNodes(@nodes);
  foreach my $node (@nodelist) {
    $node->name($offset);
    $offset += $increment;
  }
  $group->addNodes(@nodelist);
}

=head2 B<offset_number([Offset[,NodeNames]])>

Renumbers the NodeGroup by adding Offset to each node number.

B<Note:> Nodes may be lost if a portion of the group is renumbered to the same
as other existing nodes.

=cut

sub offset_number {
  croak "usage: group->offset_number(offset[,nodes])" unless @_ >= 2;
  my $group  = shift;
  my $offset = shift;
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  @nodes = $group->getNodeNames() unless @nodes;
  my @nodelist = $group->getNodes(@nodes);
  $group->deleteNodes(@nodes);
  foreach my $node (@nodelist) {
    my $n = $node->name();
    $node->name($n + $offset);
  }
  $group->addNodes(@nodelist);
}

=head2 B<translate(Translation[,NodeName(s)])>

Translates nodal x/y/z coordinates by Translation.  The translation may be
specified by either:
  * a hashref,   with keys of one or more of x,y,z
  * an arrayref, with values in order of x,y,z
  * a scalar,    where the value applies equally to each direction

  $group->translate({x => 3, z => -7});
  $group->translate([4,5,6]);

=cut

sub translate {
  croak "usage: group->translate(translation[,nodes])" unless @_ >= 2;
  my $group  = shift;
  my $trans = shift;
  my %trans = (x => 0, y => 0, z => 0);
  if (ref($trans) =~ /HASH/) {
    %trans = (%trans, %{$trans});
  } elsif (ref($trans) =~ /ARRAY/) {
    my $coord = "x";
    %trans = (%trans, map { $coord++ => $_ } @{$trans});
  } else {
    %trans = (x => $trans, y => $trans, z => $trans);
  }
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  foreach my $node ($group->getNodes(@nodes)) {
    foreach my $field ($node->fieldSet()->getFieldsMatching('coordinate')) {
      foreach my $component ($field->getComponentsMatching('[xyz]')) {
        my $index = $component->valIndex();
        my $value = $node->value($index);
        $node->value($index, $value + $trans{lc $component->name()});
      }
    }
  }
}

=head2 B<setorigin(Origin[,NodeName(s)])>

Translates nodal x/y/z coordinates so that Origin is the origin.  The origin
may be specified by either:
  * a hashref,   with keys of one or more of x,y,z
  * an arrayref, with values in order of x,y,z

  $group->setorigin([-10,10,5]);

=cut

sub setorigin {
  croak "usage: group->setorigin(origin[,nodes])" unless @_ >= 2;
  my $group  = shift;
  my $origin = shift;
  my %trans = (x => 0, y => 0, z => 0);
  if (ref($origin) =~ /HASH/) {
    while (my ($coord, $value) = each %{$origin}) {
      $trans{$coord} = -$value;
    }
  } elsif (ref($origin) =~ /ARRAY/) {
    my $coord = "x";
    %trans = (%trans, map { $coord++ => -$_ } @{$origin});
  }
  $group->translate(\%trans,@_);
}

=head2 B<scale(Scale[,NodeName(s)])>

Scales nodal x/y/z coordinates by the given Scale.  The values to scale by
may be specified by either:
  * a hashref,   with keys of one or more of x,y,z
  * an arrayref, with values in order of x,y,z
  * a scalar,    where the value applies equally to each direction

  $group->scale({x => 5});
  $group->scale(2);

=cut

sub scale {
  croak "usage: group->scale(scale[,nodes])" unless @_ >= 2;
  my $group  = shift;
  my $scale = shift;
  my %scale = (x => 1, y => 1, z => 1);
  if (ref($scale) =~ /HASH/) {
    %scale = (%scale, %{$scale});
  } elsif (ref($scale) =~ /ARRAY/) {
    my $coord = "x";
    %scale = (%scale, map { $coord++ => $_ } @{$scale});
  } else {
    %scale = (x => $scale, y => $scale, z => $scale);
  }
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  foreach my $node ($group->getNodes(@nodes)) {
    foreach my $field ($node->fieldSet()->getFieldsMatching('coordinate')) {
      foreach my $component ($field->getComponentsMatching('[xyz]')) {
        my $index = $component->valIndex();
        my $value = $node->value($index);
        $node->value($index, $value * $scale{lc $component->name()});
      }
    }
  }
}

=head2 B<multmatrix4(matrix[,NodeName(s)])>

Perform a general transformation on nodal coordinates by multiplying by a
standard 16 element transformation Matrix.  The matrix must be specified as an
arrayref with 16 elements (or a 4x4 array of arrays).  The matrix is given
row-first.

  $mat = [ 
    [ 1, 0, 0, 5],
    [ 0, 0,-1,10],
    [ 0, 1, 0, 0],
    [ 0, 0, 0, 1] ];
  $group->multmatrix4($mat);

=cut

sub multmatrix4 {
  croak "usage: group->multmatrix4(matrix[,nodes])" unless @_ >= 2;
  eval "use PDL";
  my $group  = shift;
  my $matrix = pdl( shift )->reshape(4,4);
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  foreach my $node ($group->getNodes(@nodes)) {
    my @xyz;
    foreach my $field ($node->fieldSet()->getFieldsMatching('coordinate')) {
      foreach my $component ($field->getComponentsMatching('[xyz]')) {
        push @xyz, $node->value($component->valIndex());
      }
    }
    push @xyz, 1;
    my $xyz = pdl(@xyz);
    my $new = $xyz x $matrix ;
    my @new = $new->list();
    foreach my $field ($node->fieldSet()->getFieldsMatching('coordinate')) {
      foreach my $component ($field->getComponentsMatching('[xyz]')) {
        $node->value($component->valIndex(), shift @new);
      }
    }
  }
}

=head2 B<rotateradians(Angles[,NodeName(s)]) rotatedegrees(Angles[,NodeName(s)])>

Rotates nodal x/y/z coordinates by the given Angles.  The values to rotate by
may be specified by either:
  * a hashref,   with keys of one or more of x,y,z
  * an arrayref, with values in order of x,y,z
  * a scalar,    where the rotation is assumed to be in the x-direction
Angles are given in radians or degrees depending on the function called

  $group->rotateradians({z => $pi/2});
  $group->rotatedegrees([45,-30,10]);

Both these routines called the internal routine C<_rotate> which converts
the angles to a 4x4 transformation matrix and calls the C<multmatrix4> routine.
CmUtils::RotationTransform is used to do the conversion.

=cut

sub rotateradians {
  croak "usage: group->rotateradians(angles[,nodes])" unless @_ >= 2;
  eval "use CmUtils::RotationTransform qw(angle2trans)";
  my $group  = shift;
  $CmUtils::RotationTransform::Transpose = 0;
  $group->_rotate(\&angle2trans,@_);
}

sub rotatedegrees {
  croak "usage: group->rotatedegrees(angles[,nodes])" unless @_ >= 2;
  eval "use CmUtils::RotationTransform qw(degangle2trans)";
  my $group  = shift;
  $CmUtils::RotationTransform::Transpose = 0;
  $group->_rotate(\&degangle2trans,@_);
}

sub _rotate {
  my $group  = shift;
  my $sub    = shift;
  my @rotate = (0,0,0);
  my $rotate = shift;
  if (ref($rotate) =~ /HASH/) {
    push @rotate, $rotate->{x} || 0;
    push @rotate, $rotate->{y} || 0;
    push @rotate, $rotate->{z} || 0;
  } elsif (ref($rotate) =~ /ARRAY/) {
    push @rotate, @{$rotate};
  } else {
    push @rotate, $rotate, 0, 0;
  }
  my @matrix = $sub->(@rotate);
  $group->multmatrix4(\@matrix, @_);
}

=head2 B<centroid([NodeName(s)])>

Computes the centroid of the xyz coordinate components of the NodeGroup.

  @centroid = $group->centroid();

=cut

sub centroid {
  my $group = shift;
  my @nodes = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  my (%sum, $num);
  foreach my $node ($group->getNodes(@nodes)) {
    foreach my $field ($node->fieldSet()->getFieldsMatching('coordinate')) {
      foreach my $component ($field->getComponentsMatching('[xyz]')) {
        $sum{lc $component->name()} += $node->value($component->valIndex());
        $num++;
      }
    }
  }
  my @centroid = map { $sum{$_}/$num } sort keys %sum;
  return wantarray ? @centroid : "@centroid";
}

=head2 B<distance(node1, node2)>

Computes the distance between two nodes ina nodegroup (assuming xyz
coordinates).  e.g. nodes 1 and 5:

  $distance = $group->distance(1, 5);

=cut

sub distance {
  return unless @_ == 3;
  my ($group, $node1, $node2) = @_;
  my $values1 = $group->getNode($node1)->getValuesHash();
  my $values2 = $group->getNode($node2)->getValuesHash();
  my $dist = 0;
  foreach my $coord (qw(x y z)) {
    $dist += ($values1->{coordinates}{$coord}[0]{value} - 
      $values2->{coordinates}{$coord}[0]{value}) ** 2;
  }
  return sqrt $dist;
}

1;

__END__
