package CmUtils::File::Txtelem;
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
use CmUtils::Objects::ElementGroup;
use CmUtils::Objects::Element;

=head1 CmUtils::File::Ipelem

Routines for writing F<ipelem> files.

=head1 VERSION

0.1 (14 May 2004) Glenn Ramsey

=head1 CHANGES

=head1 SUBROUTINES

=cut

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 0.5;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readTxtelem &writeTxtelem);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readTxtelem>

=over 4

=item PARAMETERS$#

(filename)

=item USAGE

  @groups = readTxtelem(Filename, NodeGroup [, Mapfile]);
  $group = readTxtelem("newfile", NodeGroup [, Mapfile]);

=item FUNCTION

Returns a list of ElementGroups read from F<filename> (or F<filename.txtelem>).
If a scalar is assigned the return value, it will contain the first ElementGroup
read.
The second parameter is a ref to a NodeGroup that contains all the nodes
named in the element file. The function will produce warnings if nodes are
not defined.
The optional third parameter is the name of a map file that contains a mapping
relating the node names used in the Txtelem file to the node names used in the
Nodegroup.

=back

=cut


sub readTxtelem {

	my ($elemfile, $nodegroup, $mapfile) = @_;
	# TODO: check that the params are valid
	
	# Read from a .txtelem file ('-' means STDIN)
	if ($elemfile ne '-' && $elemfile !~ /\.txt/) {
		$elemfile .= ".txtelem";
	}

	my $debug = 0;

  #read the map file into a hash
  my %nmap;
  if(defined $mapfile){
  	if ($mapfile ne '-' && $mapfile !~ /\.txt/) {
  		$mapfile .= ".txtmap";
  	}
    open MAPFILE, "<$mapfile" or croak "Cannot open file $mapfile";
    
    my $key;
    my $val;
    while (<MAPFILE>) {
      my (@values) = split;
      $key = shift @values;
      $val = shift @values;
      $nmap{$key} = $val;   
    }
    close MAPFILE;
  }
  #read the element file into an element group
  
  #my %elemgroups;
  my $egroup = CmUtils::Objects::ElementGroup->new($nodegroup);
  my $element_name = 1;
  
  #TODO: create a separate element group for each label so
  # that the ElementGroup object can make com file fragments for the
  # element groups.
  
  open ELEMFILE, "<$elemfile" or croak "Cannot open file $elemfile";
  
  while (<ELEMFILE>) {
    next if /^#/; #ignore comments  
    my (@values) = split; # should contain 8 nodes followed by a label
    next if 9 != @values; #ignore malformed lines
    
    my $label = pop @values; # last value is the label
    my @noderefs; # nodes in this element
    my @versions; #versions of the nodes in this element
    # translate the node names and retrieve the node objects
    foreach my $name (@values){
      #extract the version number
      my @tokens = split(/:/,$name);
      $name=shift @tokens;
      my $v;
      if(scalar @tokens){
        $v=shift @tokens;
      } else {
        $v=1;    
      }
      if($debug){
        print "node $name version=$v\n";
      }
      push @versions, $v;
      
      my $nodename = $name;
      $nodename = $nmap{$name} if(exists $nmap{$name});
      my $node = $nodegroup->getNode($nodename);
      push @noderefs, $node;
    }
    
  	my $elem=CmUtils::Objects::Element->new();
  	$elem->name($element_name);
  	$elem->addNodes(@noderefs);  	
  	$elem->addNodeVersions($nodegroup, @versions);
  	$egroup->addElement($elem, $label);
  	    
    $element_name++;
  }
  close ELEMFILE;
   
	return $egroup;
}


sub writeTxtelem {
  croak ( "writeTxtelem is not implemented yet, perhaps you could implement it?");
  
# NOT IMPLEMENTED

}

sub readTxtmap {

	my ($mapfile, $nodegroup) = @_;

	# Read from a .txtmap file ('-' means STDIN)
	if ($mapfile ne '-' && $mapfile !~ /\.txt/) {
		$mapfile .= ".txtelem";
	}

}

1;

__END__
