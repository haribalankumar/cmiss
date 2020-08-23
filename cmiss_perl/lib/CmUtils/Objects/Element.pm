package CmUtils::Objects::Element;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::Objects::Node;
use CmUtils::debug;
use Tie::RefHash;
use Data::Dumper;


=head1  CmUtils::Objects::Element

TODO: This is a work in progress until this message is removed.

Defines an Element, which is an ordered collection of nodes. The node
references are stored an array in the same order as cmiss uses.

        	Xi1	Xi2	Xi3
node # 0	0	0	0
node # 1	1	0	0
node # 2	0	1	0
node # 3	1	1	0
node # 4	0	0	1
node # 5	1	0	1
node # 6	0	1	1
node # 7	1	1	1

=head1 Element structure

=over 4

=item NAME

Name of the Element

=item NODE

An ordered list of the CmUtils::Objects::Nodes in the element

=item OCCURRENCE

A list of the number of each occurrence of the node. Corresponds
with NODE:
e.g. If NODE contains           1 2 3 2 5 2 7 8
then OCCURRENCE will contain    1 1 1 2 1 3 1 1

=item VERSION

A list of the number of each occurrence of the node. Corresponds
with NODE:
e.g. If NODE contains       1 2 3 2 5 2 7 8
then OCCURRENCE contains    1 1 1 2 1 3 1 1
and VERSION might be        1 5 1 4 1 9 1 1 .

The actual version numbers will be 1 if there is only one version
of the node or the real version number of the node if there is
more than one version. Also the version number may be greater than
1 even if none of the nodes in the element are shared.

=item COLLAPSED

Indicates the number of collapsed nodes in the element.

=item MAP

A list of the mesh DOF mappings for the element

=back

=head1 VERSION

0.1 (14 May 2004) Glenn Ramsey

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

# @nmap describes the connectivity of nodes by index.
# eg $nmap[0][0] gives the index of the node connected
# to the node at index 0 in direction 0 (xi1).

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
  
=head2 B<new()>

Returns a new Element.

=cut

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self  = {
    NAME        => undef,
    NODE        => [],
    OCCURRENCE  => [],
    VERSION     => [],
    GROUP       => undef,
    HASVERSIONS => 0,
    VERSIONMAP  => {},
    COLLAPSED   => 0
  };

  # options are
  #my $options = shift;

  # 	$options = {
  # 		VERSION => undef, # add versions to the node
  # 	}

  bless $self, $class;
  # allow use of references as hash keys 
  tie %{$self->{VERSIONMAP}}, 'Tie::RefHash::Nestable';
  return $self;
}

=head2 B<copy()>

Returns an identical copy of an existing Element.

=cut

sub copy {
  my $self = shift;

  my $class = ref($self) || $self;

  my $new = $self->new();
  $new->name( $self->name() );
  $new->setNodes( $self->getNodes() );
  $new->setOccurrences( $self->getOccurrences() );
  $new->setVersions( $self->getVersions() );
  $new->{COLLAPSED} = $self->collapsed();
  return $new;
}

=head2 B<name([string/number])>

Optional argument is name for Element.  Returns the name of the Element.

=cut

sub name {
  my $self = shift;
  croak "usage: thing->name(name)" unless @_ <= 1;
  if (@_) { $self->{NAME} = shift }
  return $self->{NAME};
}

=head2 B<group([string/number])>

Optional argument is name for Element group.  Returns the name of the Element group.

=cut

sub group {
  my $self = shift;
  croak "usage: thing->group(name)" unless @_ <= 1;
  if (@_) { $self->{GROUP} = shift }
  return $self->{GROUP};
}


=head2 B<collapsed()>

Returns the number of collapsed nodes in the element.

=cut

sub collapsed {
  my $self = shift;
  croak "usage: thing->collapsed()" unless @_ <= 1;
  return $self->{COLLAPSED};
}

=head2 B<addNodes(Nodes)>

Adds nodes to the Element. The parameter is a list of Nodes and must
contain exactly 8 nodes.
 
=cut

sub addNodes {
  use CmUtils::File::Utils qw/strtonum/;

  my $self = shift;

  # "usage: thing->addNodes(values)"

  if ( $self->debug() ) {
    carp "Adding nodes to Element " . $self->name() . "\n";
  }
  # There is no reason to have this restriction here
  #croak "Only elements with 8 nodes are supported." if ( @_ != 8 );

  # count the number of occurrences of this node in the element
  # @_ should contain refs to CMISS::Node objects
  foreach my $node (@_) {

    if ( $self->debug() ) {
      carp "Node ".$node->name()."\n";
    }
    push @{ $self->{NODE} }, $node;

    my $occ = 1;
    # iterate backwards through the node array
    for ( my $j = scalar @{ $self->{NODE} } - 1 ;$j > 0 ; $j-- )
    {
      if ( $self->{NODE}[ $j - 1 ] == $node ) {

        #if(0==$occ){$occ++};
        $occ++;

        if ( $self->debug() ) {
          carp "Node "
            . $node->name()
            #. " :ref "
            #. $node
            . " occurrence "
            . $occ . "\n";
        }
      }
    }

    push @{ $self->{OCCURRENCE} }, $occ;
    push @{ $self->{VERSION} },    1;      #use version 1 by default

    # increment the count of collapsed nodes if necessary
    $self->{COLLAPSED}++ if ( $occ > 1 );

#    if ( $self->debug() ) {
#      carp "Element "
#        . $self->name()
#        . " added node "
#        . $node->name()
#        . " :ref "
#        . $node
#        . " occurrence "
#        . $#{ $self->{OCCURRENCE} } . "\n";
#    }
  }

  return $#{ $self->{NODE} };
}

=head2 B<clearValues()>

Clears the list of nodes in the Element.

=cut

sub clearNodes {
  my $self = shift;
  $self->{NODE}       = [];
  $self->{OCCURRENCE} = [];
  $self->{COLLAPSED}  = [];
}

=head2 B<setNodes()>

Sets the list of node names of the Element (clears, then adds the new values).

=cut

sub setNodes {
  my $self = shift;
  $self->clearNodes();
  return $self->addNodes(@_);
}

=head2 B<getNodes()>

Returns an array of the Nodes in the Element.

=cut

sub getNodes {
  my $self = shift;

  # "usage: thing->getNodes()"
  return @{ $self->{NODE} };
}

=head2 B<getNodeNames()>

Returns an array of the Node names in the Element.

=cut

sub getNodeNames {
  my $self = shift;

  # "usage: thing->getNodeNames()"
  my @nodenames;
  foreach my $node ( @{ $self->{NODE} } ) {
    push @nodenames, $node->name();
  }
  return @nodenames;
}

=head2 B<node(index,[value])>

Sets and/or return a single Node in a Element by index.

=cut

sub node {
  my $self = shift;
  croak "usage: thing->node(index)" unless @_ <= 1;
  my $number = shift;
  return $self->{NODE}[$number];
}

=head2 B<getOccurrences()>

Returns the occurrences array.

=cut

sub getOccurrences {
  my $self = shift;

  # "usage: thing->getNodes()"
  return @{ $self->{OCCURRENCE} };
}

=head2 B<occurrence(index,[value])>

Sets and/or returns a single item in the Element's occurrence array by index.

=cut

sub occurrence {
  my $self = shift;
  croak "usage: thing->occurrence(number,[value])" unless @_ <= 2;
  my $number = shift;
  if (@_) { $self->{OCCURRENCE}[$number] = shift }
  return $self->{OCCURRENCE}[$number];
}

=head2 B<nodeVersion(index,[value])>

Sets and/or returns a single item in the Element's version array by index.

=cut

sub nodeVersion {
  my $self = shift;
  croak "usage: thing->nodeVersion(number,[value])" unless @_ <= 2;
  my $number = shift;
  if (@_) { $self->{VERSION}[$number] = shift }
  return $self->{VERSION}[$number];
}

=head2 B<getNodeVersions(index,[value])>

Returns the Element's version array.

=cut

sub getNodeVersions {
  my $self = shift;

  # "usage: thing->getNodeVersions()"
  return @{ $self->{VERSION} };
}

=head2 B<setNodeVersions(values)>

Sets the Element's version array.

=cut

sub addNodeVersions {
  my $self = shift;
  my $ngroup = shift;
  croak "usage: thing->addNodeVersions(nodegroup, values)" unless @_ >= 1;
  #TODO check that the given array has eight values
  croak "Parameter array \"values\" must contain exactly 8 elements" unless @_ == 8; 
  
  if ( $self->debug() ) {print "\n"; $self->list();}

  my $nodeindex=0;
  foreach my $version (@_) {
    $self->nodeVersion($nodeindex, $version);
    
    if($version>1){
      $self->{HASVERSIONS}=1; 
      my $node = $self->node($nodeindex);
      #croak "Node $nodename not found" unless ( ref $node );
      my $values = $node->getValuesHash();
  
      #copy the version 1 values into this version, relies on version 1 being declared first 
      #foreach my $i ( 1 .. $version - 1 ) {
        foreach my $component ( keys %{ $values->{coordinates} } ) {
          $values->{coordinates}{$component}[$version-1]{value} =
              $values->{coordinates}{$component}[0]{value};
          if ( exists( $values->{coordinates}{$component}[0]{derivatives} ) ) {
            $values->{coordinates}{$component}[$version-1]{derivatives} =
                $values->{coordinates}{$component}[0]{derivatives};
            #check that all verisions have values
            foreach my $v (1..$version-1) {
              if(!exists $values->{coordinates}{$component}[$v]{derivatives}){
                 $values->{coordinates}{$component}[$v]{derivatives}= $values->{coordinates}{$component}[0]{derivatives};
              }
            }  
          } 
#          # set all the derivatives 
#          for my $v ($version-1..0) {
#            $values->{coordinates}{$component}[$v]{derivatives} = [0,0,0,0,0,0,0];
#          }
        }
      #}
      my $options = { newfieldset => 1 };

      if ( $self->debug() ) {print "addNodeVersions: Adding values to node ".$node->name()." for version $version \n";}
      # switch off debugging of this call to the Node class
      my $debug = $self->debug();
      $self->debug(0);
      $node->setValuesHash( $values, $options, $ngroup );
      $self->debug($debug);
      
      #always map out the position
      my $maptover = 1;                 # always map to version 1
      my $xj       = 1;
      foreach my $comp ($node->fieldSet()->getField("coordinates")->getComponents()) {
        # index in cm is 1 based so need to add 1 to this value
        # always map out the position (index=0)
 
        if ( $self->debug() ) {print "node ".$node->name()." mapping out component ".$comp->name()." version $version position\n";}
        $self->{VERSIONMAP}->{$node}[$xj-1]{$version}{0} = [ $node, $maptover, $xj, 1 ];
        # for debugging map out all dof's
#        foreach my $deriv ( 1 .. $comp->derivatives() ) {
#          print "node ".$node->name()." mapping out component ".$comp->name()." version $version deriv $deriv\n";
#          $self->{VERSIONMAP}{$node}[$xj-1]{$version}{$deriv} = [ $node, $maptover, $xj, $deriv+1 ];
#        }
               
        $xj++;
        
      }
    }
    $nodeindex++;
    
  }

  if ( $self->debug() ) {print "\n"; $self->list();}
    
  # Find the mapping for this element

  if($self->{HASVERSIONS}) {   
    my @occs      = $self->getOccurrences();
    my @vers      = $self->getNodeVersions();
    my @nodenames = $self->getNodeNames();
    
    if ( $self->debug() ) { print "addNodeVersions: Element " . $self->name() . "\n"; }
    
    for ( my $index1 = scalar @occs - 1 ; $index1 >= 0 ; $index1-- ) {

        #my $name1 = $elem->node($index1)->name();
      my $name1 = $nodenames[$index1];
      if ( $occs[$index1] > 1 ) {    #node is collapsed

        if ( $self->debug() ) { print "addNodeVersions: found a collapsed node $name1 at index $index1 \n";}
      
        if ( $self->debug() ) { print "name:i:o[i]:v[i] $name1:$index1:$occs[$index1]:$vers[$index1] "."\n";}

        # find the other collapsed node
        for ( my $index2 = $index1 - 1 ; $index2 >= 0 ; $index2-- ) {
        
          if ( $self->debug() ) { print "addNodeVersions:  index1: $index1 index2: $index2 \n"; }
          
          # ignore if indices are more than 4 apart - they aren't connected 
          last if ( ( $index1 - $index2 ) > 4 );
          my $name2 = $nodenames[$index2];
          
            #print "name1: $name1 name2: $name2\n";
#          if ( ( $index1 != $index2 )
#            && ( $occs[$index2] > 1 )
#            && ( $name2 != $name1 ) )
          if ( ( $occs[$index2] > 1 )
            && ( $name2 != $name1 ) )
          {
          
            # found a collapsed node
            if ( $self->debug() ) { print "addNodeVersions: found another collapsed node $name2 at index $index2\n";}
            # if this node is connected to the first then map the derivative
            # find the direction in which it is connected
            # only use exact powers of 2 ie 2**0=1, 2**1=2, 2**2=4
            my $connecteddir = log2( $index1 - $index2 );

            if ( $self->debug() ) {print "addNodeVersions: determining if it is connected to node $name1 ...\n";}
            
            if ( $connecteddir != POSIX::floor($connecteddir) ){
            
              if ( $self->debug() ) {print "addNodeVersions: nope.\n";}
            
              last;
            }
            
            if ( $self->debug() ) {print "addNodeVersions: yes it is! connected in xi".($connecteddir+1)." (dir=$connecteddir)\n";}
            
            if ( ( $nodenames[ $nmap[$index1][$connecteddir] ] == $name2 )
              && ( $nodenames[ $nmap[$index2][$connecteddir] ] == $name1 ) )
            {
            
              # Find the dir in which the node is collapsed so
              # that can be mapped out also.
              # This can be done by looking through the @nmap array and comparing names 
              my $collapseddir;
              for ($collapseddir=0;$collapseddir<2;$collapseddir++) {
                last if($nodenames[ $nmap[$index1][$collapseddir]]==$name1);
              }
              if ( $self->debug() ) {print "collapsed in xi".($collapseddir+1)."(dir=$collapseddir) \n";}
              
              
              if ( $self->debug() ) {
                print " $nodenames[$nmap[$index2][$connecteddir]]:$vers[$index1] "
                     ."<->  $nodenames[$nmap[$index1][$connecteddir]]:$vers[$index2] \n";
              }
              my @indices =
                ( $index1, $index2 );    # list of indices of nodes to map out
              my @names = ( $name1, $name2 );

                # only map the derivatives of nodes within themselves
              foreach my $n ( 0 .. @names - 1 ) {
                my $i = $indices[$n];    #index in element of this node
                my @mvers;
                push @mvers, $vers[$i] if($vers[$i]>1);    # the versions we are mapping out

                # Map the version of this occurrence to the version of the first occurrence.
                # Find the index of the first occurrence of this node.
                for my $k ( 0 .. @nodenames - 1 ) {
                
                  if ( $self->debug() ) { print " $k $nodenames[$k]-$names[$n]";}
                
                  if ( $nodenames[$k] eq $names[$n] ) {
                    push @mvers, $vers[$k] if($vers[$k]>1);

                    if ( $self->debug() ) { print ":$vers[$k] "; }
                  }
                }
                if ( $self->debug() ) { print "\n"; }

                #my $maptover = $vers[$index1]; # the version to map to
                my $maptover = 1;                 # always map to version 1
                my $node     = $self->node($i);
                
                
                for my $dir ($connecteddir,$collapseddir){
                #for my $dir ($connecteddir){
             
                  my $xj = 1;
                  
                  foreach my $comp ($node->fieldSet()->getField("coordinates")->getComponents()) {
                    # index in cm is 1 based so need to add 1 to this value
                    # always map out the position (index=0)
                    #print "comp:".$comp->name()." ver ";
                    #my $dindex = 2**$dir;
                    #my @dmap2=(
                    my @dmap; #1->1, 2->2, 4->3
                    my $dindex;
                    foreach my $version (@mvers) {
  
                      # position is already mapped out
                      #$self->{VERSIONMAP}->{ $node->name() }[ $xj - 1 ]{$version}{0} =
                      #  [ $node->name(), $maptover, $xj, 1 ];

                      if ( $self->debug() ){print "node ".$node->name()."\n";}
                      @dmap=(1,2,4); #1->1, 2->2, 4->3
                      $dindex=$dmap[$dir];
                      if ( $self->debug() ){print "mapping out component ".$comp->name()." version $version deriv $dindex\n";}
                      $self->{VERSIONMAP}->{ $node }[$xj - 1]{$version}{$dindex} =
                        [ $node, $maptover, $xj, 1 + $dindex ];
                      # also map out the cross derivatives
  
                      @dmap=(3,3,5);#3->1&2, 3->2&1, 5->2&3
                      $dindex=$dmap[$dir];
                      if ( $self->debug() ){print "mapping out component ".$comp->name()." version $version deriv $dindex\n";}
                      $self->{VERSIONMAP}->{ $node }[ $xj - 1 ]{$version}{$dindex} =
                        [ $node, $maptover, $xj, 1 + $dindex ];
                        
                      @dmap=(5,6,6);
                      $dindex=$dmap[$dir];
                      if ( $self->debug() ){print "mapping out component ".$comp->name()." version $version deriv $dindex\n";}
                      $self->{VERSIONMAP}->{ $node }[ $xj - 1 ]{$version}{$dindex} =
                        [ $node, $maptover, $xj, 1 + $dindex ];
  
                      # The triple cross derivative never appears in a face fitting problem
                      $dindex=7;                    
                      $self->{VERSIONMAP}->{ $node }[ $xj - 1 ]{$version}{$dindex} =
                        [ $node, $maptover, $xj, 1 + $dindex ];
                    }
                    $xj++;
                  }
                }
              }    # foreach $nodename
            } # if ( ( $nodenames[ $nmap[$index1][$dir] ] == $name2 ) ...
          }
        }
      }

        #$index1++;
    }
  #}    # if($elem->collapsed() > 0)
    if ( $self->debug() ) { print "\n"; }
  } #if($self->{HASVERSIONS})
  
  #print Dumper($self->{VERSIONMAP});
    
}

sub getMapping {
  my $self = shift;
  return $self->{VERSIONMAP};
}


=head2 B<setOptions(values)>

Set the options hash.

=cut

sub setOptions {
  my $self = shift;
  $self->clearOptions();
  return $self->addOption(@_);
}

=head2 B<clearOptions(index,[value])>

Clear the options hash.

=cut

sub clearOptions {
  my $self = shift;
  $self->{OPTIONS} = [];
}

=head2 B<getOptions(index,[value])>

Get the options hash.

=cut

sub getOptions {
  my $self = shift;
  return @{ $self->{OPTIONS} };
}

=head2 B<swapXi(Xidir1,Xidir2)>

Swaps the Xi directions in the element.
Takes 2 integer arguments that indicate the Xi directions to operate on.
e.g.
   1 2 swap Xi1 and Xi2
   1 3 swap Xi1 and Xi3
   2 1 same as 1 2
   1 1 reverse Xi1 direction
   2 2 reverse Xi1 direction

Warning: not tested very much. 

=cut

sub swapXi {
  # TODO: make this work for 2D and 1D elements
  my $self = shift;
  croak "usage: thing->swapXi( direction1, direction2)" unless @_ == 2;
  my $d1 = shift;
  my $d2 = shift;
  croak "direction indices must be 1,2 or 3"
    unless ( ( 1 <= $d1 ) && ( $d1 <= 3 ) && ( 1 <= $d2 ) && ( $d2 <= 3 ) );
  $d1--;
  $d2--;

  if ( $d1 != $d2 ) {
    my @first  = ( 2, 3, 5 );
    my @second = ( 7, 6, 4 );

    my $index1 = $first[$d1];
    my $index2 = $first[$d2];
    swap_nodes( $index1, $index2 );

    $index1 = $second[$d1];
    $index2 = $second[$d2];
    swap_nodes( $index1, $index2 );

  }
  else {

    # both indices the same - Xi dir reversal
    # Xi1 0<>1, 2<>3, 4<>5, 6<>7
    # Xi2 0<>2, 1<>3, 4<>6, 5<>7
    # Xi3 0<>4, 1<>5, 2<>5, 3<>7
    my $offset = 2**$d1;
    my $step   = 2**( $d1 + 1 );

    #print "$offset $step\n";
    #return;
    for ( my $i = 0 ; $i < 8 ; $i += $step ) {
      for ( my $j = 0 ; $j < $offset ; $j++, $i++ ) {

        #print "index1=".$i." index2=". ($i+$offset) ."\n";
        $self->swap_nodes( $i, $i + $offset );
      }
      $i = $i - $offset;
    }

  }

}

#used internally by swapXi
sub swap_nodes {
  my $self = shift;
  croak "usage: thing->swap_nodes(index1, index2)" unless @_ == 2;
  my $index1 = shift;
  my $index2 = shift;

  my $tmp;
  $tmp = @{ $self->{NODE} }[$index1];
  @{ $self->{NODE} }[$index1] = @{ $self->{NODE} }[$index2];
  @{ $self->{NODE} }[$index2] = $tmp;

}

=head2 B<list([FH])>

List the contents of this Element to the filehandle I<FH> (default to STDOUT).
Set the "verbose" flag ($s->verbose(1)) to see extended information.

=cut

sub list {
  my $self = shift;
  my $FH   = (@_) ? shift: \*STDOUT;

  my @nodes = $self->getNodes();
  printf $FH <<EOF, $self->name(), scalar @nodes;
   Element     %s
     Number of values       %d
     Index Name Ocrnce Vrsn
EOF
  my @occ  = $self->getOccurrences();
  my @vers = $self->getNodeVersions();
  for ( my $i = 0 ; $i < @nodes ; $i++ ) {
    printf $FH "     %4d     %s    %4d    %4d\n", $i, $nodes[$i]->name(),
      $occ[$i], $vers[$i];
  }
}

1;

__END__
