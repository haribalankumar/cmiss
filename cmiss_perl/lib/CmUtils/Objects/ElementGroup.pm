package CmUtils::Objects::ElementGroup;
require 5.006;
use strict;
use warnings;
use Carp;
use CmUtils::Objects::Element;
use CmUtils::Objects::Node;
use CmUtils::debug;
use Tie::RefHash;
use Data::Dumper;

=head1 CMISS::ElementGroup

The ElementGroup class defines structures and methods for creating and
manipulating groups of elements.  The ElementGroup structure has fields:

The name ElementGroup is no longer a good name since group labels have been
added. It would be more appropriate to name it TriCubicMesh or similar.

=over 4

=item NAME

Name of the ElementGroup

=item ELEMENTS

The set of Elements belonging to this ElementGroup

=item NODEGROUP

The NodeGroup containing the nodes in this ElementGroup

=item CONNECTIONS

A hash containing:

  {nmap}[0..versions]
	An array with one entry for each version of the node containing a
	3*3 array with the nodes that this node is connected to and a flag
	indicating if the direction is collapsed.
	
	{node}
	A ref to the Node itself
	
	{elements}[0..versions][0..7]
	An hash indicating which elements each version of this node occurs in.
  The array index is the position of the element in relation to the node.
  It is labelled the same way as the node positions within an element. 

  TODO: might be better arranged like this:
  [0..versions]
    {nmap}[][]
    {elements}[][]
  
=item GROUPS

A hash containing lists of element names. The key is the group name.

=item VERSIONMAP

A hash containing mapping for nodal mesh dofs. The structure is:

  {nodename}[direction]{version}{derivative number}


=back

=head2 Acknowledgements

A lot of the methods here were copied from CmUtils::Objects::Node
and modified to suit.

=head2 Warning

A lot of the functions that were adapted from other code have not
been tested. Use with caution. 

=head1 VERSION

$Id: ElementGroup.pm,v 1.12 2012-07-27 06:34:22 gram003 Exp $

=head1 CHANGES

=head1 METHODS

=cut

sub log2 {
  my $n = shift;
  return log($n) / log(2);
}

BEGIN {
  our ( $VERSION, @ISA );
  $VERSION = 0.1;
  @ISA     = qw(CmUtils::debug);
}

# These tables define the way that the nodes are connected in a hexahedral
# element. They are used by the addVersions() function and its helper
# functions.
# Use the position index of the node (0-7) within the element to select the
# row in the table. The column (0-2) gives the relevant information for
# the xi directions (xi dir = (col index) + 1).
# E.g. for node 3 connections[3][n] shows it is connected to nodes 2, 1, and
# 7. The Xi 1 direction node is node 2 and occurs at Xi1=xi[3][0].
our @xi = (
  [ 0, 0, 0 ],
  [ 1, 0, 0 ],
  [ 0, 1, 0 ],
  [ 1, 1, 0 ],
  [ 0, 0, 1 ],
  [ 1, 0, 1 ],
  [ 0, 1, 1 ],
  [ 1, 1, 1 ]
);
our @nmap = (
  [ 1, 2, 4 ],
  [ 0, 3, 5 ],
  [ 3, 0, 6 ],
  [ 2, 1, 7 ],
  [ 5, 6, 0 ],
  [ 4, 7, 1 ],
  [ 7, 4, 2 ],
  [ 6, 5, 3 ]
);
our @emap = ( 7, 3, 5, 1, 6, 2, 4, 0 );

=head1 VERSION

0.1 (12 May 2004)

=head1 METHODS

=head2 B<new(CmUtils::Objects::NodeGroup)>

Returns a new ElementGroup. The argument is a NodeGroup array.

=cut

sub new {
  croak "usage: thing->new(NodeGroup)" unless @_ >= 2;
  my $proto     = shift;
  my $class     = ref($proto) || $proto;
  my $nodegroup = shift;

# check that the argument has the correct type
# not sure if this is really required because perl will find the error eventually
#	if(defined($nodegroups[0])) {
#		carp "CmUtils::Objects::ElementGroup::new(): the array argument must contain CmUtils::Objects::NodeGroup objects"
#			unless((ref($nodegroups[0]) =~ /CmUtils::Objects::NodeGroup/));
#	}
  my $self = {
    NAME         => '',
    ELEMENTS     => {},
    NODEGROUP    => $nodegroup,
    CONNECTIONS => undef,
    GROUPS       => {},
    VERSIONMAP   => {},
    RETRIES      => {}
  };
  
  bless $self, $class;
  
  # allow use of references as hash keys 
  tie %{$self->{VERSIONMAP}}, 'Tie::RefHash::Nestable';
  
  return $self;
}

=head2 B<name([string/number])>

Optional argument is name for ElementGroup.  Returns the name of the ElementGroup.

=cut

sub name {
  my $self = shift;
  croak "usage: thing->name(name)" unless @_ <= 1;
  if (@_) { $self->{NAME} = shift }
  return $self->{NAME};
}

sub getGroup {
  my $self  = shift;
  my $group = shift;
  return $self->{GROUPS}{$group};
}

sub getGroups {
  my $self   = shift;
  my @keys   = $self->getGroupNames(@_);
  my @groups = ();
  push @groups, $self->{GROUPS}{$_} foreach @keys;
  return @groups;
}

sub getGroupNames {
  my $self = shift;
#  return (@_)
#    ? grep { exists $self->{GROUPS}{$_} }
#    @_
#
#    #: sort {$a <=> $b || $a cmp $b} keys %{$self->{GROUPS}};
#    : keys %{ $self->{GROUPS} };

  return keys %{ $self->{GROUPS} };
}

sub getElementGroup {
  my $self = shift;
  my $groupname = shift;
  
  #print "getElementGroup() groupname=$groupname\n";
  
  my $egr = $self->new($self->{NODEGROUP});
  my @keys  = @{$self->getGroup($groupname)};
  my @elems = ();
  push @elems, $self->{ELEMENTS}{$_} foreach @keys;
  $egr->addElements(\@elems,undef);
  return $egr;    
}

=head2 B<addELement(Element [,groupname])>

Adds one Element to the list of Elements in the ElementGroup.  Elements will
NOT be added (replaced) if they already exist.  Returns the total number of
defined Elements in the ElementGroup.

=cut

sub addElement {
  my $self = shift;
  croak "usage: thing->addElement(element [,groupname])" unless @_ >= 1;
  my $elem  = shift;
  my $groupname = shift;

  #tie %{$self->{VERSIONMAP}}, 'Tie::RefHash::Nestable';

  unless ( exists $self->{ELEMENTS}{ $elem->name() } ) {
    #print "addElement() adding element ".$elem->name()."\n";
    $self->{ELEMENTS}{ $elem->name() } = $elem;
  } else {
    carp "element ".$elem->name() . " already exists in group" . $self->name();
  }
  if ( defined $groupname ) {
  
    #print "addElement() setting to group name to $groupname \n";
    
    $elem->group($groupname);
    push @{ $self->{GROUPS}{$groupname} }, $elem->name();
  }

  # add mapping for this element to the group's map
    
  print "get_mapping: element: ".$elem->name()."\n" if $self->debug();

  my $nodemap;
  tie %{$nodemap}, 'Tie::RefHash::Nestable';
  $nodemap = $elem->getMapping();
      
  foreach my $node (keys %{$nodemap}){
    print "get_mapping: node: ".$node." name=".$node->name()."\n" if $self->debug();
    foreach my $component (0..@{$nodemap->{$node}}-1){
      print "get_mapping: component: ".$component."\n" if $self->debug();
      foreach my $version (sort keys %{$nodemap->{$node}[$component]}){
        print "get_mapping: version ".$version."\n" if $self->debug();
        foreach my $deriv (sort keys %{$nodemap->{$node}[$component]{$version} }){
          print "get_mapping: ".$node->name()." $component $version $deriv \n" if $self->debug();
#            if(exists $self->{VERSIONMAP}{$node}[$component]{$version}{$deriv}){
#              print "warning overwriting mapping for ".$node->name().":$component:$version:$deriv\n";
#            } 
          
          $self->{VERSIONMAP}->{$node}[$component]{$version}{$deriv} = $nodemap->{$node}[$component]{$version}{$deriv};
        }
      }
    }
  }
    
  $self->{HASVERSIONS} = undef;
  
  print "element ".$elem->name()."\n";
  #update the CONNECTIONS data structure
  for my $i (0..7)
  {
    my $node = $elem->node($i);
    $self->{CONNECTIONS}{$node->name()}{elements}[$elem->nodeVersion($i)-1][$emap[$i]] = $elem->name();
    #$self->{CONNECTIONS}{$node->name()}{elements}[$emap[$i]] = $elem->name();
    #print "node ".$node->name()." version ".$elem->nodeVersion($i)."\n";
  }
  
  return $self->numberOfElements();
}

=head2 B<addELements(Element(s)[,groupname])>

Adds one or more Elements to the list of Elements in the ElementGroup.  Elements will
NOT be added (replaced) if they already exist.  Returns the total number of
defined Elements in the ElementGroup.

=cut

sub addElements {
  my $self  = shift;
  croak "usage: thing->addElements(elements [,group])" unless @_ >= 1;
  my $elements = shift;
  #print Dumper( $elements );
  my $groupname = shift;
  #print "addElements() groupname=$groupname\n" if defined $groupname;
  $self->addElement( $_, $groupname ) foreach @{$elements};
  return $self->numberOfElements();
}

=head2 B<updateElement(Element)>

Updates the Element in the ElementGroup.  Element will NOT be changed unless it already
exists.  Returns the total number of defined Elements in the ElementGroup.

=cut

sub updateElement {
  my $self = shift;
  croak "usage: thing->updateElement(elem)" unless @_ == 1;
  my $elem = shift;
  if ( exists $self->{ELEMENTS}{ $elem->name() } ) {
    delete $self->{ELEMENTS}{ $elem->name() };
    $self->{ELEMENTS}{ $elem->name() } = $elem;
  } else {
    carp $elem->name() . " does not exist in " . $self->name();
  }
  return $self->numberOfElements();
}

=head2 B<updateElementLabel(ElementName ,GroupName])>

Change the group of the named element.

=cut

sub updateElementLabel {
  my $self = shift;
  croak "usage: thing->updateElementLabel(element,groupname)" unless @_ >= 2;
  my $elemname  = shift;
  my $groupname = shift;

  if ( exists $self->{ELEMENTS}{ $elemname } ) {
    
    my $elem = $self->{ELEMENTS}{ $elemname };

    my $currentgroup = $elem->group();
    if(!$currentgroup){
      $currentgroup = "undef";
    }

    print " found ".$elemname." cur grp=".$currentgroup." new group=".$groupname."\n";
    #print Dumper $elem;

    # find the current group that the element is in and remove it from that group
      
    if(defined $currentgroup && exists $self->{GROUPS}{$currentgroup}){
      # my @a = qw(a b c);
      # @a = grep { $_ ne 'b' } @a; # remove 'b' from the array
      @{$self->{GROUPS}{$currentgroup}} = grep { $_ != $elem } @{$self->{GROUPS}{$currentgroup}};
    }
    
    # add the element to the desired group
    $elem->group($groupname);
    push @{ $self->{GROUPS}{$groupname} }, $elem->name();
     
  } else {
    carp "element ".$elemname . " not found";
  }
}

=head2 B<deleteElement(ElementName)>

Deletes Element named ElementName in the ElementGroup.  Returns the total number of
defined Elements in the ElementGroup.

=cut

sub deleteElement {
  my $self = shift;
  croak "usage: thing->deleteELement(name)" unless @_ == 1;
  my $name = shift;
  if ( exists $self->{ELEMENTS}{$name} ) {
    delete $self->{ELEMENTS}{$name};
  } else
  {
    carp "can't delete element $name : it doesn't exist\n"; 
  }
  # map needs to be recreated if an element is deleted
  $self->{VERSIONMAP} = undef; 
  $self->{HASVERSIONS} = undef;
  
  return $self->numberOfElements();
}

=head2 B<deleteElements(ElementName(s))>

Deletes one or more Elements named ElementName(s) in the ElementGroup.
Returns the total number of defined Elements in the ElementGroup.

=cut

sub deleteElements {
  my $self = shift;
  croak "usage: thing->deleteElements(name(s))" unless @_ >= 1;
  my @names = $self->getElementNames(@_);
  $self->deleteElement($_) foreach @names;
  return $self->numberOfElements();
}

=head2 B<deleteElementGroup(GroupName)>

Deletes a named (sub)group of the elements in this ElementGroup.
Returns the total number of defined Elements in the ElementGroup.

=cut

sub deleteElementGroup {
  my $self = shift;
  croak "usage: thing->deleteElementGroup(GroupName)" unless @_ == 1;
  my $groupname = shift;
  if ( exists $self->{GROUPS}{$groupname} ) {
    if($self->debug())
    {
      carp "deleting group $groupname\n";
    }
    my $egr = $self->getElementGroup($groupname);
    $self->deleteElement($_) foreach $egr->getElementNames();
    delete $self->{GROUPS}{$groupname};
  } else {
    carp "group: $groupname doesn't exist\n";
  }
  
  return $self->numberOfElements();
}

=head2 B<isElement(ElementName)>

Returns True if the Element named ElementName is in the ElementGroup, otherwise returns
False.

=cut

sub isElement {
  my $self = shift;
  croak "usage: thing->isElement(name)" unless @_ == 1;
  my $name = shift;
  return exists $self->{ELEMENTS}{$name};
}

=head2 B<getElements([ElementName(s)])>

If ElementName(s) are specified, returns a list of names of all Elements with these
names that exist in the ElementGroup, in the order they are specified.  If no
argument is given, returns a list of names of ALL Elements, sorted in
numerical/alphabetical order.

=cut

sub getElementNames {
  my $self = shift;
  return (@_)
    ? grep { exists $self->{ELEMENTS}{$_} } @_
    : sort { $a <=> $b || $a cmp $b } keys %{ $self->{ELEMENTS} };
}

=head2 B<getElements([ElementName(s)])>

If ElementName(s) are specified, returns a list of all Elements with these names that
exist in the ElementGroup, in the order they are specified.  If no argument is
given, returns a list of ALL Elements, sorted in numerical/alphabetical order.

=cut

sub getElements {
  my $self  = shift;
  my @keys  = $self->getElementNames(@_);
  my @elems = ();
  push @elems, $self->{ELEMENTS}{$_} foreach @keys;
  return @elems;
}

=head2 B<getElement(ElementName)>

Returns the Element named ElementName from the ElementGroup.

=cut

sub getElement {
  my $self = shift;
  croak "usage: thing->getElement(name)" unless @_ == 1;
  my $name = shift;
  return $self->{ELEMENTS}{$name};
}

=head2 B<numberOfElements()>

Returns the number of Elements defined in the ElementGroup.

=cut

sub numberOfElements {
  my $self = shift;
  return scalar keys %{ $self->{ELEMENTS} };
}

=head2 B<maximumElementNumber()>

Returns the number of Elements defined in the ElementGroup.

=cut

sub maximumElementNumber {
  my $self = shift;
  my $max  = 0;
  do { $max = $_ if $_ > $max }
    foreach keys %{ $self->{ELEMENTS} };
  return $max;
}

=head2 B<listElements([FH])>

Lists the Element information of this ElementGroup to the filehandle I<FH>
(default to STDOUT).  Set the "verbose" flag ($s->verbose(1)) to see
extended information.

=cut

sub listElements {
  use CmUtils::File::Utils qw(strtonum);
  use CmUtils::Utils qw(list_to_string);
  my $self = shift;
  my $FH = (@_) ? shift: \*STDOUT;
  printf $FH " ElementGroup                %s\n",   $self->name();
  printf $FH "   Number of Elements          %d\n", $self->numberOfElements();
  if ( $self->numberOfElements() ) {
    my @elems = $self->getElementNames();
    if ( strtonum( $elems[0] ) ) {    #numeric element names
      printf $FH "   Elements                    %s\n", list_to_string(@elems);
    }
    if ( $self->verbose() ) {
      foreach my $elem ( $self->getElements() ) {
        $elem->list($FH);
      }
    }
  }
}

sub getNodeGroup {
  my $self  = shift;
  return $self->{NODEGROUP};
}

# TODO:
# =head2 B<mergeGroup(Group[,Options])>
#
# Merges another Group into this ElementGroup.  By default, elements are added as long
# as they do not already exist.  The Options hashref can alter this by either
# forcing an overwrite of existing nodes (key "overwrite"), or assigning a new
# number to the new elements (key "newnumber").
#
#   $options{newnumber} = 1;
#   $group->mergeGroup($another_group,$options);
#
# =cut
#
# sub mergeGroup {
#   my $group = shift;
#   croak "usage: thing->mergeGroup(group[,options])" unless @_ >= 1;
#   my ($newgroup, $options) = @_;
#   croak "invalid new Group specified" unless ref($newgroup) =~ /CmUtils::Objects::ElementGroup/;
#   my %options = (newnumber => 0, overwrite => 0);
#   if ($options && ref($options) =~ /HASH/) {
#     %options = (%options, %{$options});
#   }numberOfElements
#   my $number = $options{newnumber} ? $group->maxElementNumber() : 0;
#   foreach my $elem ($newgroup->getElements()) {
#     my $fs = $elem->fieldSet();
#     my $same = 0;
#     my $i = 0;
#     foreach my $fset ($group->getFieldSets()) {
#       last if ($same = $fs->isSame($fset));
#       $i++;
#     }
#     if ($same) {
#       $elem->fieldSet($group->getFieldSet($i));
#     } else {
#       $group->addFieldSet($fs);
#       $elem->fieldSet($fs);
#     }
#     if ($options{newnumber}) {
#       $number++;
#       $elem->name($number);
#       $group->addElement($elem);
#     } else {
#       if($group->isElement($elem->name())) {
#         $group->updateElement($elem) if $options{overwrite};
#       } else {
#         $group->addElement($elem);
#       }
#     }
#   }
# }

=head2 B<renumber([Offset[,Increment[,ElementNames]]])>

Renumbers the ElementGroup starting from Offset (default 1) by an Increment
(default 1).

B<Note:> Elements may be lost if a portion of the group is renumbered to the same
as other existing elements.

=cut

sub renumber {
  croak "usage: group->renumber([offset[,increment[,elements]]])"
    unless @_ >= 1;
  my $group     = shift;
  my $offset    = shift || 1;
  my $increment = shift || 1;
  my @elems     = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  @elems = $group->getElementNames() unless @elems;
  my @elemlist = $group->getElements(@elems);
  $group->deleteElements(@elems);

  foreach my $elem (@elemlist) {
    $elem->name($offset);
    $offset += $increment;
  }
  $group->addElements(\@elemlist);
}

=head2 B<offset_number([Offset[,ElementNames]])>

Renumbers the ElementGroup by adding Offset to each element number.

B<Note:> Elements may be lost if a portion of the group is renumbered to the same
as other existing elements.

=cut

sub offset_number {
  croak "usage: group->offset_number(offset[,elements])" unless @_ >= 2;
  my $group  = shift;
  my $offset = shift;
  my @elems  = map { ref($_) =~ /ARRAY/ ? @{$_} : $_ } @_;
  @elems = $group->getElementNames() unless @elems;
  my @elemlist = $group->getElements(@elems);
  $group->deleteElements(@elems);
  foreach my $elem (@elemlist) {
    my $n = $elem->name();
    $elem->name( $n + $offset );
  }
  $group->addElements(@elemlist);
}

=head2 B<get_connected_nodes(Element, NodeIndex)>

Used internally by addVersions.

Creates a 3x3 array which defines the nodes that are connected to
this node. Cols are the Xi direction, rows 0 & 1 are the inverse of the Xi value.
Row 2 is a flag indicating if the direction is collapsed.

E.g. [a,b,c]
     [d,e,f]
     [0,0,1]
     
Means that node a is in the +Xi1 direction, node d in the -Xi1 direction, etc, 
and that the node is collapsed in xi3. 

=cut

sub get_connected_nodes {
  my $self = shift;
  my $elem = shift;
  my $i    = shift;             # node index of this occurrence of the node
  #my $k    = shift;             # node index of second occurrence of the node
  my $node = $elem->node($i);

  # Determine which nodes the node connects to - use thelookup table
  # TODO: direction is not considered here - it might be needed
  # and the @xi array can be used to determine that
  # print "connected nodes for node ".$node->name()." at index $i:\n";
  #my $kk;

  # If this is the first occurrence of the collapsed node then
  # find the next occurrance so that later we can work out which
  # xi dir is collapsed.

  #if ( $k == $i ) {
  #  for ( $kk = $i - 1 ; $kk >= 0 ; $kk-- ) {
  #   last if ( $node == $elem->node($kk) );
  #  }
  #} else {
  #  $kk = $k;
  #}

  # format of @connectednodes is 3*3, rows 0 & 1 are Xi, col is direction
  # the third row contains flags that indicate if the direction is collapsed
  my $connectednodes = [ [ 0, 0, 0 ], [ 0, 0, 0 ], [ 0, 0, 0 ] ];
  #my $collapsedxi;
  for ( my $j = 0 ; $j < 3 ; $j++ ) {

    # Retrieve the index numbers for the connected nodes
    # from the lookup table and get the node refs for the
    # nodes at those indices. Store these in 2x3 array.
    my $pos   = $nmap[$i][$j];
    my $cnode = $elem->node($pos);
    my $dir   = $xi[$i][$j];

    #print "row=$row i=$i k=$k j=$j\n";
    $connectednodes->[$dir][$j] = $cnode->name();

    # test if this is the collapsed xi dir.
    #if(($k!=$kk) && ( $xi[$kk][$j] != $xi[$i][$j] )) {
    if($node == $cnode){
      #$collapsedxi = $j;
      $connectednodes->[2][$j] = 1;
      
      if($self->debug()){print "get_connected_nodes: collapsed xi=".($j+1)." at index $pos \n";}

#$connectednodes->[$row][$j] = 0;
#$connectednodes->[1-$row][$j] = $cnode->name(); # 1-$row gives the binary inverse
    }
  }

  #$connectednodes->[2] = $collapsedxi;
  return $connectednodes;
}

sub compare_connections {
  my $self            = shift;
  my $connectednodes  = shift;
  my $existingversion = shift;
  #my $collapsedxi = $connectednodes->[2];
  my @e               = ( [ 0, 0, 0 ], [ 0, 0, 0 ], $existingversion->[2] );
  my @c               = ( [ 0, 0, 0 ], [ 0, 0, 0 ], $connectednodes->[2] );
  my @m               = ( [ " ", " ", " " ], [ " ", " ", " " ], [ " ", " ", " " ]);
  my $match           = 0;
  my $nomatch         = 0;
  my $nomatchxi       = 0;
  my $collapsedxi     =-1;
  my @matchxi;
  
  for ( my $xj = 0 ; $xj < 3 ; $xj++ ) {
    for ( my $d = 0 ; $d < 2 ; $d++ ) {

      #if($self->debug()) {print "[$d][$xj] ";}
      # only consider the defined nodes in the connectednodes we are comparing
      #if(defined($connectednodes->[$d][$xj])) {
      $c[$d][$xj] = $connectednodes->[$d][$xj];
      if ((1==$existingversion->[2][$xj]) || (1==$connectednodes->[2][$xj])){
        $collapsedxi = $xj;
      }
      if ( ( 0 != $existingversion->[$d][$xj] )
        && ( 0 != $connectednodes->[$d][$xj] ) )
      {
        if ( $existingversion->[$d][$xj] == $connectednodes->[$d][$xj] ) {

          #if($d!=$collapsedxi){
            $match++;
            $m[$d][$xj] = "Y";
            push @matchxi, $xj;  
          #}
        } else {
          unless((1==$existingversion->[2][$xj]) || (1==$connectednodes->[2][$xj])){
            $nomatch++;
            $m[$d][$xj] = "N";
            $nomatchxi = $xj;
          }
        }
      }
      $e[$d][$xj] = $existingversion->[$d][$xj];
    }
#    my $d = 2;
#    if($existingversion->[$d][$xj] == $connectednodes->[$d][$xj]){
#          $nomatch++;
#          $m[$d][$xj] = "N";      
#    }

    #				if($self->debug()) {print "\n";}
  }
  if ( $self->debug() ) {
    print "  match        connected         existing\n";
    for ( my $xj = 0 ; $xj < 3 ; $xj++ ) {
      for ( my $d = 0 ; $d < 3 ; $d++ ) {
        print "[" . $m[$d][$xj] . "]";
      }
      print "   ";
      for ( my $d = 0 ; $d < 3 ; $d++ ) {
        printf "[%3d]", $c[$d][$xj];
      }
      print "   ";
      for ( my $d = 0 ; $d < 3 ; $d++ ) {
        printf "[%3d]", $e[$d][$xj];
      }
      print "\n";
    }
    print "match=$match nomatch=$nomatch \n";
  }
  return ( $match, $nomatch, $nomatchxi, $collapsedxi, @matchxi);
}

=head2 B<get_version(Node, ConnectedNodes)>

Used internally by addVersions

=cut

sub get_version {
  my $self           = shift;
  my $node           = shift;
  my $connectednodes = shift;


  # If there are existing versions then try to find a matching
  # version for this node in the version list.
  #This doesn't consider Xi orientation. If the mesh has
  # consistent Xi directions then this will probably do the
  # correct thing.
  # Actually if the Xi directions are consistent but the directions
  # are reversed then this is OK and this algorithm needs to
  # work in that case - I think it will because it doesn't consider
  # the orientation.
  
  my $versionnum = 1;
  my @possible;
  my $numcollapsed=0;
  my $collapsedxi;
  my $createversion = 1;
  my $match;
  my @matchxi;
  my $nomatch;
  my $nomatchxi;
  
  my $nodename = $node->name();

  #print $node." ".$node->name." dump follows...........\n";
  #print Dumper($self->{CONNECTIONS});
  if ( exists $self->{CONNECTIONS}{ $nodename }{nmap} ){    #AoA of Nodes
     #print "Node already has ".scalar @{$self->{CONNECTIONS}{$node->name()}{nmap}}." versions \n";

    # A matching node id is defined as one where at least 2 out of 3 of the
    # adjacent nodes are the same.
    my @versarrays = @{ $self->{CONNECTIONS}{ $nodename }{nmap} };

    #my @versarrays = @{$self->{CONNECTIONS}{$node->name()}{nmap}};
    foreach my $existingversion (@versarrays) {
      if ( $self->debug() ) {
        print "Trying to match version ". $versionnum . " of ". @versarrays ."\n";
          #. " collapsed in Xi". (1+$connectednodes->[2]) . ":\n";

        #			  print Dumper($existingversion);
      }
      my @matches = $self->compare_connections( $connectednodes, $existingversion);
      $match   = shift @matches; #$matches[0];
      $nomatch = shift @matches; #$matches[1];
      $nomatchxi = shift @matches; #$matches[2];
      $collapsedxi = shift @matches; #$matches[3];
      @matchxi = @matches;
      
      if ( $nomatch == 0 ) {    # it might match
                                #if($nomatch<=1) { # it might match
        if ( $match >= 2 ) {

          # this version matches - use it and don't create a new one
          if ( $self->debug() ) { print "version $versionnum matches\n" }
          $createversion = -1;

          #update the version if the match had an empty space
          if ( 2 == $match ) {
            for ( my $xj = 0 ; $xj < 3 ; $xj++ ) {
              # increment counter if this dir is collapsed
              #$numcollapsed++ if (1==$connectednodes->[2][$xj]);
              for ( my $d = 0 ; $d < 2 ; $d++ ) {
                # only consider the defined nodes in the connectednodes we are comparing
                if ( 0 != ( $connectednodes->[$d][$xj] ) ) {
                  if ( 0 == ( $existingversion->[$d][$xj] ) ) {
                    $existingversion->[$d][$xj] = $connectednodes->[$d][$xj];
                  } else {
                    if ( $existingversion->[$d][$xj] != $connectednodes->[$d][$xj] ){
                      #check that it matches, it should unless this is the collapsed node
                      if($self->debug()){print "non matching versions detected in node ".$nodename." version $versionnum \n";}
                    }
                  }
                }
              }
            }
          }
          last;
        } elsif ( $match <= 1 ) {

          #print "version $versionnum might match - but we can't tell yet\n";
          $createversion = 0;
          push @possible, $versionnum;
        }
      #} elsif((1==$nomatch) && (1==$match) && (1==$numcollapsed) && ($nomatchxi != $collapsedxi)){
      #} elsif((1==$nomatch) && (1==$match) && ($nomatchxi != $collapsedxi)){
      #} elsif((1==$nomatch) && ($match>=1) && ($nomatchxi != $collapsedxi)){

      }
      if((1==$createversion) && ($match>=1) && ($nomatchxi != $collapsedxi)){
      
        # update the mapping data structure
        
        # if node is collapsed and there is one match not in the collasped dir and
        # there is one non-match then map the non-match dir to version one.

        my $maptover = 1;                 # always map to version 1
        #my $node     = $elem->node($i);
        my $vers = (1 == $createversion) ? 1+$versionnum : $versionnum;
        foreach my $m (@matchxi){
          #my $dindex = 2**$matchxi;
          my $dindex = 2**$m;
          if ( $self->debug() ) {
              print "mapping out node:vers:deriv $nodename:$vers:$dindex \n";
          }
          my $xj = 1;
          foreach my $comp ($node->fieldSet()->getField("coordinates")->getComponents()) {
            # index in cm is 1 based so need to add 1 to this value
            # always map out the position (index=0)
            #print "comp:".$comp->name()." ver ";
            
            #my $dindex = 2**$dir;
  
            #foreach my $ver (@mvers) {
  
              $self->{VERSIONMAP}->{ $nodename }[ $xj - 1 ]{$vers}{0} = [ $nodename, $maptover, $xj, 1, 1 ];
              $self->{VERSIONMAP}->{ $nodename }[ $xj - 1 ]{$vers}{$dindex} =[ $nodename, $maptover, $xj, $dindex+1, 1 ];
            #}
            $xj++;
          }
        }
          
      }
      $versionnum++;
    } # foreach @versarrays
  }

  # createversion:
  #  1 - new version created
  # -1 - existing version matched
  #  0 - possible match

  if ( 1 == $createversion ) {
    
    # the current connection array becomes the new version 
    push @{$self->{CONNECTIONS}{$nodename}{nmap}}, $connectednodes;
    $self->{CONNECTIONS}{$nodename}{node} = $node;
    if ( $self->debug() ) {
      print "created new version "
        . ( @{ $self->{CONNECTIONS}{$nodename}{nmap} } )
        . " for node "
        . $node->name() . " :\n";
    }
        
  } elsif(0 == $createversion) {
    $versionnum = 0;
  }
  
    
  # first value of array:
  # yes = matched an existing version (in 2 or more positions)
  # maybe = possible match versions
  # no = created a new version
  
  my @results;
  push @results, $versionnum;
  if(0==$versionnum){
    push @results, @possible;
  }
 
  return wantarray ? @results : $versionnum;
}

=head2 B<addVersions(NodeGroup)>

Finds elements with collapsed nodes and adds versions to the relevant nodes.

=over

=head4 B<Implementation notes>
	
Builds a list of the number of versions each node requires.

The following 2 tables are used to find the indices of nodes that
are connected to the other nodes in the element. 

	my @xi=
	(	[0,0,0],
		[1,0,0],
		[0,1,0],
		[1,1,0],
		[0,0,1],
		[1,0,1],
		[0,1,1],
		[1,1,1]);
	
	my @connections=
	(	[1,2,4],
		[0,3,5],
		[3,0,6],
		[2,1,7],
		[5,6,0],
		[4,7,1],
		[7,4,2],
		[6,5,3]);

If a node occurs at index $pos then the nodes which it is connected
to in the $xj direction are given by:

$connections[pos,$xj-1] - gives the index of the other node.
$xi[$pos,$xj-1] - gives the sign of direction of the connection,
0=+ve, 1=-ve.

E.g. $connections[3,1] = 1, means that the node at index 3 is connected
to the node at index 1 in the Xi2 direction. $xi[3,1] = 1 means that
the other node is connected in the +ve Xi2 direction.
  		
=back

=cut

sub addVersions {
  my $self = shift;
  my $ngroup;

  # if a nodegroup parameter was given then use it, otherwise use the
  # instance version
  if ( @_ > 0 ) {
    carp "ElementGroup::addVersions(): The NODEGROUP parameter is no longer required\n";
    $ngroup = shift;
  } else {
    $ngroup = $self->{NODEGROUP};
  }

  # TODO: make this work for 2D (4 noded) elements
  # first need to build the list of versions and connected nodes
  # this will tell us how many versions each node has
  # then we need to determine which of the versions each element node has
  # If we can't find a match for a node, then we put the element back
  # in the list.
  #my @retry = values %{$self->{ELEMENTS}};
  my @retry = sort { $a <=> $b } keys %{ $self->{ELEMENTS} };

  #my @retry;
  #foreach my $elem ( values %{$self->{ELEMENTS}}) {
  while (@retry) {

    #my $elem = shift @retry;
    my $key  = shift @retry;
    my $elem = $self->{ELEMENTS}{$key};
    if ( $self->debug() ) {
      print "\n\n==================================================================\n";
      print "\nProcessing Element: " . $elem->name() . "\n";
      for my $name ( $elem->getNodeNames() ) {
        print "$name ";
      }
      print "\n";
    }

    #	don't try to add versions to elements with fewer than 8 nodes
    next unless ( 8 == $elem->getNodes() );
    #for ( my $i = 7 ; $i >= 0 ; $i-- ) {
    for my $i (0..7) {
      my $node = $elem->node($i);

      if($self->debug()){
        print "\nIndex $i Processing Node: ".$node->name()." Occurrence: ".$elem->occurrence($i)."\n";
      }
      my $versionnum=1;

      # if the node occurs more than once
      #if ( $elem->occurrence($i) > 1 ) {

        # find the next occurrence
        #my $debug = $self->debug();
        #$self->debug(0);


#        for ( my $k = $i ; $k >= 0 ; $k-- ) {
#          if ( $node == $elem->node($k) ) {
#            if ( $self->debug() ) { print "node ". $node->name(). " has occurrence "
#                . $elem->occurrence($k). " at index $k\n";
#            }

            #if($self->debug()) {$elem->list();}
            my $connectednodes = $self->get_connected_nodes( $elem, $i);
 
            #if($self->debug()){print "connected nodes: ".Dumper($connectednodes)."\n";}
            #$versionnum = $self->get_version( $node, $connectednodes );
            my @results = $self->get_version( $node, $connectednodes );
            
            $versionnum = shift @results;
            
            if ( 0 == $versionnum ) {

        # we couldn't tell, put the element back in the list try again after all
        # the others have been processed
              if ( $self->debug() ) {
                print "Couldn't match any existing version for node "
                  . $node->name() . "\n";
              }

              if($self->{RETRIES}{$key}++<3){              
                #push @retry, $elem;
                push @retry, $key;
              } else {
              
                $versionnum = shift @results; 
                carp "Element ".$elem->name()." Node ".$node->name()." Local Node ".$i
                  ." could not find a double match... using version $versionnum.\n";
              }
              
            }
            
            if($versionnum > 0)
            {
              if($self->debug()){print "Setting versionnum for node at index $i to $versionnum \n";} 
              $elem->nodeVersion( $i, $versionnum );
#              $versionnum = 1;
            }
          }
#        }# for $k

#      } else {
#
## This detects nodes that require extra versions because of inconsistent xi dirs.
## Not %100 sure how/why/if it works.
#        if ( $self->debug() ) {
#          print "Occurrence <=1 Node "
#            . $node->name()
#            . " has occurrence "
#            . $elem->occurrence($i)
#            . " at index $i \n";
#        }
#        my $connectednodes = $self->get_connected_nodes( $elem, $i);
#
##
##				#if($self->debug()){print "connected nodes: ".Dumper($connectednodes)."\n";}
##
#        $versionnum = $self->get_version( $node, $connectednodes);
#        if ( $versionnum < 1 ) {
#            carp "Simple Element ".$elem->name()." Node ".$node->name()." Local Node ".$i
#              ." could not find a match, using version 1.\n";
#             $versionnum = 1; 
#        }
#        $elem->nodeVersion( $i, $versionnum );
#
#        #
#        if ( $self->debug() ) { print "version: $versionnum \n"; }
#      }
#    }

    #$elem->list();
  }    # foreach my $elem

  #	# Reprocess the elements that we couldn't match on the first pass.
  #	foreach my $elem ( @retry) {
  #		#print "Reprocessing element: ".$elem->name()."\n";
  #
  #		for(my $i=7; $i>=0; $i--) {
  #			my $node = $elem->node($i);
  #			my $versionnum;
  #
  #			# if the node occurs more than once
  #			if($elem->occurrence($i)>1) {
  #				# find the next occurrence
  #				for(my $k=$i; $k>=0; $k--) {
  #				 	if($node == $elem->node($k)) {
  #						#print "node has occurrence ".$elem->occurrence($k)." at index $k\n";
  #						my $connectednodes = $self->get_connected_nodes($elem, $i, $k);
  #						#print "connected nodes: ".Dumper($connectednodes)."\n";
  #						$versionnum = $self->get_version($node, $connectednodes);
  #
  #						if(-1 == $versionnum){
  #								croak "couldn't match on the 2nd try";
  #						}
  #
  #					}
  #					# Fill in the VERSIONS field in the element. The version number is
  #					# the index of the version in the CONNECTIONS hash.
  #					$elem->nodeVersion($k, $versionnum);
  #					$versionnum = 1;
  #				}
  #			}
  #		}
  #		#$elem->list();
  #	} # foreach my $elem
  #$self->listNodeVersions();
  #now add versions to each node that needs them
  my $debug = $self->debug();
  $self->debug(0);
  foreach my $nodename ( keys %{ $self->{CONNECTIONS} } ) {

    #get the node
    my $node = $self->{CONNECTIONS}{$nodename}{node};
    croak "Node $nodename not found" unless ( ref $node );
    my $values = $node->getValuesHash();

   #$values->{fieldname}{componentname}[versionnumber]{value}
   #$values->{fieldname}{componentname}[versionnumber]{derivatives}[derivnumber]
    my $nversions = scalar @{ $self->{CONNECTIONS}{$nodename}{nmap} };

    #copy the version 1 values into the other versions
    foreach my $i ( 1 .. $nversions - 1 ) {
      foreach my $component ( keys %{ $values->{coordinates} } ) {
        $values->{coordinates}{$component}[$i]{value} = $values->{coordinates}{$component}[0]{value};
        if ( exists( $values->{coordinates}{$component}[0]{derivatives} ) ) {
          $values->{coordinates}{$component}[$i]{derivatives} = $values->{coordinates}{$component}[0]{derivatives};
        }
      }
    }
    my $options = { newfieldset => 1 };
    $self->debug($debug);
    if ( $self->debug() ) {
      print "Adding values to node " . $node->name() . "\n";
    }
    $self->debug(0);
    #$self->debug(0);
    $node->setValuesHash( $values, $options, $ngroup );

    #$self->debug(1);
  }
  $self->debug($debug);

#  # Update the version list for each element.
#  foreach my $key ( sort { $a <=> $b } keys %{ $self->{ELEMENTS} } ) {
#
#    #$elem->list();
#    my $elem = $self->{ELEMENTS}{$key};
#    if ( $self->debug() ) {
#      print "\nElement " . $elem->name() . "\n";
#      for my $name ( $elem->getNodeNames() ) {
#        print "$name ";
#      }
#      print "\n";
#     }
#
#    #check each node in the element
#    my @names = $elem->getNodeNames();
#    for my $i ( 0 .. @names - 1 ) {
#      my $nodename = $names[$i];
#      if ( $self->debug() ) { print "node: $nodename \n"; }
#      my $node = $self->{CONNECTIONS}{$nodename}{node};
#
#      #if the node has versions then choose the correct one
#      next unless ( defined @{ $self->{CONNECTIONS}{$nodename}{nmap} } );
#      my @connected = @{ $self->{CONNECTIONS}{$nodename}{nmap} };
#      for my $n ( 0 .. @connected - 1 ) {
#        my $connectednodes = $self->get_connected_nodes( $elem, $i, $i );
#        my @matches = $self->compare_connections( $connected[$n], $connectednodes );
#        my $match   = $matches[0];
#        my $nomatch = $matches[0];
#        if ( $self->debug() ) { print "matches: $match nomatches $nomatch \n"; }
#        if ((0==$nomatch) && ($match >= 2))    {
#          $elem->nodeVersion( $i, $n + 1 );
#          last;
#        }
#      }
#    }
#  }

  $self->{HASVERSIONS}=1;

}

=head2 B<get_mapping> used internally by writeIpmap

=over

=head4 B<Implementation notes>

Find collapsed nodes using the occurrences array.

The number of index positions between the 2 occurrences indicates which
direction is collapsed, eg:
 1 1 2 2 5 6 7 8 indices 1 & 2 are 1 index apart and are collapsed in Xi1
 1 2 3 4 5 6 5 6 indices 5 & 6 are 2 indices apart and collapsed in Xi2
 1 2 3 4 1 2 7 8 indices 5 & 6 are 4 indices apart and collapsed in Xi3

distance apart is 2^(xj-1)

Any matching nodes that are more than 4 indices apart do not need to be
mapped.

algorithm:

for each collapsed element having more than 1 collapsed node
  for each collapsed node (determine using occurrence>1)
    for each other collasped node (determine using occurrence>1)
      find the direction of the connection to this node
      add an entry to the mapping list
    end for
  end for
end for

=back	

=cut

sub get_mapping {
  my $self = shift;
  
#  if(!defined $self->{HASVERSIONS}){
#    $self->addVersions();
#  }

# tie %{$self->{VERSIONMAP}}, 'Tie::RefHash::Nestable';
  
#  foreach my $elem ($self->getElements()){
#    print "get_mapping: element: ".$elem->name()."\n" if $self->debug();
#  
#    my $nodemap;
#    tie %{$nodemap}, 'Tie::RefHash::Nestable';
#    $nodemap = $elem->getMapping();
#        
#    foreach my $node (keys %{$nodemap}){
#
#
#      print "get_mapping: node: ".$node." name=".$node->name()."\n" if $self->debug();
#      foreach my $component (0..@{$nodemap->{$node}}-1){
#        print "get_mapping: component: ".$component."\n" if $self->debug();
#        foreach my $version (sort keys %{$nodemap->{$node}[$component]}){
#          print "get_mapping: version ".$version."\n" if $self->debug();
#          foreach my $deriv (sort keys %{$nodemap->{$node}[$component]{$version} }){
#            print "get_mapping: ".$node->name()." $component $version $deriv \n" if $self->debug();
#            if(exists $self->{VERSIONMAP}{$node}[$component]{$version}{$deriv}){
#              print "warning overwriting mapping for ".$node->name().":$component:$version:$deriv\n";
#            } 
#            
#            $self->{VERSIONMAP}->{$node}[$component]{$version}{$deriv} = $nodemap->{$node}[$component]{$version}{$deriv};
#          }
#        }
#      }
#    }
#    
#  }
  #print "get_mapping VERSIONMAP\n";  
  #print Dumper($self->{VERSIONMAP});
  return $self->{VERSIONMAP};
}

sub addMapping {
  my $self = shift;
  my $vmap = shift;
  
  # re tieing clobbers the existing hash      
  #tie %{$self->{VERSIONMAP}}, 'Tie::RefHash::Nestable';
  
#  my $map;
#  tie %{$map}, 'Tie::RefHash::Nestable';
#  $map = $self->{VERSIONMAP};
    
#  my $versionmap;
#  tie %{$versionmap}, 'Tie::RefHash::Nestable';
#  $versionmap = $self->get_mapping();
  
  print "ElementGroup->addMapping VERSIONMAP:\n";
  #print Dumper($self->{VERSIONMAP});
  
  for my $nodename (keys %{$vmap})
  {
    print "addMapping for node $nodename\n";
    my $fromnode = $self->{NODEGROUP}->getNode($nodename);
    #for my $component (0..(@{$map->{$nodename}}-1))
    my $xj = 1;
    foreach my $component ( $fromnode->fieldSet()->getField("coordinates")->getComponents() )
    {
      print "component ".$component->name()." Xj=$xj\n";
      for my $version (keys %{$vmap->{$nodename}[$xj-1]}){
        print "version $version\n";
        for my $deriv (keys %{$vmap->{$nodename}[$xj-1]{$version}}){
          print "deriv $deriv\n";
          
          my $tonodename = $vmap->{$nodename}[$xj-1]{$version}{$deriv}[0];
          my $tonode = $self->{NODEGROUP}->getNode($tonodename);
          my $toversion = $vmap->{$nodename}[$xj-1]{$version}{$deriv}[1];
          
          #my $tocomponent = $xj;
          my $tocomponent = $vmap->{$nodename}[$xj-1]{$version}{$deriv}[2];
          
          my $toderiv = $vmap->{$nodename}[$xj-1]{$version}{$deriv}[3];
 
          my $factor = $vmap->{$nodename}[$xj-1]{$version}{$deriv}[4];
          
          print "to: $tonodename $tocomponent $toversion $toderiv $factor\n";

          if(exists $self->{VERSIONMAP}->{$fromnode}[$xj-1]{$version}{$deriv}){
            print "warning overwriting mapping for ".$fromnode->name().":".($xj).":$version:$deriv\n";
          } 
                  
         $self->{VERSIONMAP}->{$fromnode}[$xj-1]{$version}{$deriv} = [ $tonode,  $toversion, $tocomponent, $toderiv, $factor ];
          #print "fromnode ".$self->{VERSIONMAP}->{$fromnode}[$xj-1]{$version}{$deriv}[0]->name()."\n";

          #print Dumper($self->{VERSIONMAP}->{$fromnode}[$xj-1]{$version}{$deriv});

#          for my $value (1..@{$map->{$nodename}{$version}{$deriv}}-1){
#            print "$nodename $version $xj $deriv = ".@{$map->{$nodename}{$version}{$deriv}}[$value]."\n";
#            push $self->{VERSIONMAP}{$fromnode}[$xj-1]{$version}{$deriv},
#                $map->{$nodename}{$version}{$deriv}[$value];
#          }
        }
      }
      $xj++;
    }
  } 
  
#print "ElementGroup->addMapping VERSIONMAP 2:\n";
#print Dumper($self->{VERSIONMAP});
  
  #print Dumper($map);
#  if(exists $self->{VERSIONMAP}{$fromnode}[$fromcomponent]{$fromversion}{$fromderiv}){
#    carp "warning overwriting mapping for $fromnode:$fromcomponent:$fromversion:$fromderiv";
#  } 
#  push @{$self->{VERSIONMAP}{$fromnode}[$fromcomponent]{$fromversion}{$fromderiv}}, @values;
}

# checks for consistent Xi directions
# potentially could also automatically correct inconsistencies
sub checkXi {

  my $self = shift;
  my $ngroup = $self->{NODEGROUP};

  # TODO: make this work for 2D (4 noded) elements
  # first need to build the list of versions and connected nodes
  # this will tell us how many versions each node has
  # then we need to determine which of the versions each element node has
  # If we can't find a match for a node, then we put the element back
  # in the list.
  #my @retry = values %{$self->{ELEMENTS}};
  my @retry = sort { $a <=> $b } keys %{ $self->{ELEMENTS} };

  #my @retry;
  #foreach my $elem ( values %{$self->{ELEMENTS}}) {
  while (@retry) {

    #my $elem = shift @retry;
    my $key  = shift @retry;
    my $elem = $self->{ELEMENTS}{$key};
    if ( $self->debug() ) {
      print "\n\n==================================================================\n";
      print "\nProcessing Element: " . $elem->name() . "\n";
      for my $name ( $elem->getNodeNames() ) {
        print "$name ";
      }
      print "\n";
    }

    #	don't try to add versions to elements with fewer than 8 nodes
    next unless ( 8 == $elem->getNodes() );
    for my $i (0..7) {
      my $node = $elem->node($i);

      if($self->debug()){
        print "\nIndex $i Processing Node: ".$node->name()." Occurrence: ".$elem->occurrence($i)."\n";
      }
      my $versionnum=1;

      # if the node occurs more than once
      #if ( $elem->occurrence($i) > 1 ) {

        # find the next occurrence
        #my $debug = $self->debug();
        #$self->debug(0);


#        for ( my $k = $i ; $k >= 0 ; $k-- ) {
#          if ( $node == $elem->node($k) ) {
#            if ( $self->debug() ) { print "node ". $node->name(). " has occurrence "
#                . $elem->occurrence($k). " at index $k\n";
#            }

            #if($self->debug()) {$elem->list();}
            my $connectednodes = $self->get_connected_nodes( $elem, $i);
 
            #if($self->debug()){print "connected nodes: ".Dumper($connectednodes)."\n";}
            #$versionnum = $self->get_version( $node, $connectednodes );
            my @results = $self->get_version( $node, $connectednodes );
            
            $versionnum = shift @results;
            
            if ( 0 == $versionnum ) {

        # we couldn't tell, put the element back in the list try again after all
        # the others have been processed
              if ( $self->debug() ) {
                print "Couldn't match any existing version for node "
                  . $node->name() . "\n";
              }

              if($self->{RETRIES}{$key}++<3){              
                #push @retry, $elem;
                push @retry, $key;
              } else {
              
                $versionnum = shift @results; 
                carp "Element ".$elem->name()." Node ".$node->name()." Local Node ".$i
                  ." could not find a double match... using version $versionnum.\n";
              }
              
            }
            
            if($versionnum > 0)
            {
              if($self->debug()){print "Setting versionnum for node at index $i to $versionnum \n";} 
              $elem->nodeVersion( $i, $versionnum );
#              $versionnum = 1;
            }
          }
#        }# for $k

#      } else {
#
## This detects nodes that require extra versions because of inconsistent xi dirs.
## Not %100 sure how/why/if it works.
#        if ( $self->debug() ) {
#          print "Occurrence <=1 Node "
#            . $node->name()
#            . " has occurrence "
#            . $elem->occurrence($i)
#            . " at index $i \n";
#        }
#        my $connectednodes = $self->get_connected_nodes( $elem, $i);
#
##
##				#if($self->debug()){print "connected nodes: ".Dumper($connectednodes)."\n";}
##
#        $versionnum = $self->get_version( $node, $connectednodes);
#        if ( $versionnum < 1 ) {
#            carp "Simple Element ".$elem->name()." Node ".$node->name()." Local Node ".$i
#              ." could not find a match, using version 1.\n";
#             $versionnum = 1; 
#        }
#        $elem->nodeVersion( $i, $versionnum );
#
#        #
#        if ( $self->debug() ) { print "version: $versionnum \n"; }
#      }
#    }

    #$elem->list();
  }    # foreach my $elem

  #	# Reprocess the elements that we couldn't match on the first pass.
  #	foreach my $elem ( @retry) {
  #		#print "Reprocessing element: ".$elem->name()."\n";
  #
  #		for(my $i=7; $i>=0; $i--) {
  #			my $node = $elem->node($i);
  #			my $versionnum;
  #
  #			# if the node occurs more than once
  #			if($elem->occurrence($i)>1) {
  #				# find the next occurrence
  #				for(my $k=$i; $k>=0; $k--) {
  #				 	if($node == $elem->node($k)) {
  #						#print "node has occurrence ".$elem->occurrence($k)." at index $k\n";
  #						my $connectednodes = $self->get_connected_nodes($elem, $i, $k);
  #						#print "connected nodes: ".Dumper($connectednodes)."\n";
  #						$versionnum = $self->get_version($node, $connectednodes);
  #
  #						if(-1 == $versionnum){
  #								croak "couldn't match on the 2nd try";
  #						}
  #
  #					}
  #					# Fill in the VERSIONS field in the element. The version number is
  #					# the index of the version in the CONNECTIONS hash.
  #					$elem->nodeVersion($k, $versionnum);
  #					$versionnum = 1;
  #				}
  #			}
  #		}
  #		#$elem->list();
  #	} # foreach my $elem
  #$self->listNodeVersions();
  #now add versions to each node that needs them
  my $debug = $self->debug();
  $self->debug(0);
  foreach my $nodename ( keys %{ $self->{CONNECTIONS} } ) {

    #get the node
    my $node = $self->{CONNECTIONS}{$nodename}{node};
    croak "Node $nodename not found" unless ( ref $node );
    my $values = $node->getValuesHash();

   #$values->{fieldname}{componentname}[versionnumber]{value}
   #$values->{fieldname}{componentname}[versionnumber]{derivatives}[derivnumber]
    my $nversions = scalar @{ $self->{CONNECTIONS}{$nodename}{nmap} };

    #copy the version 1 values into the other versions
    foreach my $i ( 1 .. $nversions - 1 ) {
      foreach my $component ( keys %{ $values->{coordinates} } ) {
        $values->{coordinates}{$component}[$i]{value} = $values->{coordinates}{$component}[0]{value};
        if ( exists( $values->{coordinates}{$component}[0]{derivatives} ) ) {
          $values->{coordinates}{$component}[$i]{derivatives} = $values->{coordinates}{$component}[0]{derivatives};
        }
      }
    }
    my $options = { newfieldset => 1 };
    $self->debug($debug);
    if ( $self->debug() ) {
      print "Adding values to node " . $node->name() . "\n";
    }
    $self->debug(0);
    #$self->debug(0);
    $node->setValuesHash( $values, $options, $ngroup );

    #$self->debug(1);
  }
  $self->debug($debug);

#  # Update the version list for each element.
#  foreach my $key ( sort { $a <=> $b } keys %{ $self->{ELEMENTS} } ) {
#
#    #$elem->list();
#    my $elem = $self->{ELEMENTS}{$key};
#    if ( $self->debug() ) {
#      print "\nElement " . $elem->name() . "\n";
#      for my $name ( $elem->getNodeNames() ) {
#        print "$name ";
#      }
#      print "\n";
#     }
#
#    #check each node in the element
#    my @names = $elem->getNodeNames();
#    for my $i ( 0 .. @names - 1 ) {
#      my $nodename = $names[$i];
#      if ( $self->debug() ) { print "node: $nodename \n"; }
#      my $node = $self->{CONNECTIONS}{$nodename}{node};
#
#      #if the node has versions then choose the correct one
#      next unless ( defined @{ $self->{CONNECTIONS}{$nodename}{nmap} } );
#      my @connected = @{ $self->{CONNECTIONS}{$nodename}{nmap} };
#      for my $n ( 0 .. @connected - 1 ) {
#        my $connectednodes = $self->get_connected_nodes( $elem, $i, $i );
#        my @matches = $self->compare_connections( $connected[$n], $connectednodes );
#        my $match   = $matches[0];
#        my $nomatch = $matches[0];
#        if ( $self->debug() ) { print "matches: $match nomatches $nomatch \n"; }
#        if ((0==$nomatch) && ($match >= 2))    {
#          $elem->nodeVersion( $i, $n + 1 );
#          last;
#        }
#      }
#    }
#  }

  $self->{HASVERSIONS}=1;

}


=head2 B<swapXi(dir1, dir2)>

Swaps Xi directions for all elements in the ElementGroup.

=cut

sub swapXi {
  my $self = shift;
  croak "usage: thing->swapXi( direction1, direction2 )" unless @_ == 2;
  my $dir1 = shift;
  my $dir2 = shift;
  foreach my $elem ( $self->getElements() ) {
    $elem->swapXi( $dir1, $dir2 );
  }
}

sub listNodeVersions {
  my $self     = shift;
  my $FH       = (@_) ? shift: \*STDOUT;
#  my %elems = %{ $self->{ELEMENTS} };    # HoA
#    printf $FH <<EOF, $node, scalar @{ $versions{$node} };
#   Node     %s Number of versions       %d
#   Connected nodes are:
#EOF
#    foreach my $version ( @{ $versions{$node} } ) {
#      foreach my $n ( @{$version} ) {
#        print $FH " $n ";
#      }
#      print $FH "\n";
#    }
#  }
  foreach my $elemname  ( sort {$a <=> $b} keys %{ $self->{ELEMENTS} }) {
    print " element ".$elemname."\n";
    my $elem = $self->{ELEMENTS}{$elemname};
    for my $i (0..7)
    {
      my $node = $elem->node($i);
      #$self->{CONNECTIONS}{$node->name()}{elements}[$elem->nodeVersion($i)-1][$emap[$i]] = $elem->name();
      #$self->{CONNECTIONS}{$node->name()}{elements}[$emap[$i]] = $elem->name();
      print "node ".$node->name()." version ".$elem->nodeVersion($i)."\n";
    }
  }
}

# make node names (which are numbers) contiguous and remove nodes
# that are not in any elements from the nodegroup
sub compact {
  my $self = shift;
  
#  if($self->debug()){
#    $self->{NODEGROUP}->verbose(1);
#    $self->{NODEGROUP}->listNodes();
#  }
  
  # make a list of all nodes in the elements
  my $nodegroup = CmUtils::Objects::NodeGroup->new();
  $nodegroup->name($self->{NODEGROUP}->name());
  
  foreach my $elem ($self->getElements()){
    foreach my $node ($elem->getNodes()){
       $nodegroup->addNode($node) if(!defined $nodegroup->getNode($node->name()));
    }
  }
  
  my @a = $nodegroup->getNodeNames();
  if(@a){
    $nodegroup->renumber();
  }
  
  $self->{NODEGROUP} = $nodegroup;
#  if($self->debug()){
#    $self->{NODEGROUP}->verbose(1);
#    $self->{NODEGROUP}->listNodes();
#    $self->{NODEGROUP}->verbose(0);
#  }

  return $nodegroup; 
}

1;
__END__
