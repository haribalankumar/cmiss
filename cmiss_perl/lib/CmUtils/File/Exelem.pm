package CmUtils::File::Exelem;
require 5.006;

#use strict 'vars';
#use strict 'subs';
use strict;
use warnings;
use Carp;

use CmUtils::File::Utils;

use CmUtils::Objects::Element;
use CmUtils::Objects::ElementGroup;

=head1 CmUtils::File::Exelem

Routines for reading (and maybe writing in the future) simple F<exelem> files.

It is intended to read exelem files that cmgui produces when you interactively
define a mesh by using element_creator. Currently it will only read linear
elements. 

=head1 VERSION

0.1 (12 May 2004) Glenn Ramsey

=head1 SUBROUTINES

=cut

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = "0.1";
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readExelem);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readExnode>

=over 4

=item PARAMETERS$#

(filename)

=item USAGE

  @groups = readExelem($filename, @nodegroups);
  $group = readExelem("newfile", @nodegroups);

=item FUNCTION

Returns a list of ElementGroups read from F<filename> (or F<filename.exelem>).
If a scalar is assigned the return value, it will contain the first ElementGroup
read.
The second parameter is a ref to a NodeGroup that contains all the nodes
named in the element file. The function will produce warnings if nodes are
not defined.

=back

=cut

sub readExelem {

	my ($filename, @nodegroups) = @_;
	# TODO: check that the params are valid
	
	# Read from an .exnode file ('-' means STDIN)
	if ($filename ne '-' && $filename !~ /\.ex/) {
		$filename .= ".exelem";
	}

	# Group will contain the list of Exelem groups
	my @Group;
	# Current group
	my $group;
	
	my $ngr;

	my $debug = 0;

	open EXFILE, "<$filename" or croak "Cannot open file $filename";

	my $read_node_names = 0;
	# the reader loop has 2 states
	#	1. looking for a top level element
	# 2. looking for the node names

	# the current element
	#my $elem;

	my $number_of_nodes;
	my $element_name;

	while (<EXFILE>) {
  	#ignore whitespace
  	s/^\s+//;
  	s/\s+$//;
  	if(!$read_node_names) {
  		SWITCH: { for ($_) {
  			/c.Hermite/ && do {
  				croak "Only linear elements are supported.";
  				last SWITCH;
  			};
  			/^Group name\s*:\s*(.*)$/ && do {
  				# we have a new group -> create a new ElementGroup
  				my $grname = $1;
  				for(my $i=0; $i<@nodegroups;$i++){
  					$ngr = $nodegroups[$i] if($grname eq $nodegroups[$i]->name()); 
  				}
  				croak "Node group ".$grname." not found." if(!defined $ngr);

  				$group = CmUtils::Objects::ElementGroup->new($ngr);
  				#$group->debug(1);  
  				$group->name($grname);
  				push @Group, $group;
  				if ($group->debug()) {
  					carp "Created new element group ". $group->name();
  				}
  				last SWITCH;
  			};
  			/#Nodes=\s*(\d+)/ && do {
  				$number_of_nodes = $1;
  				last SWITCH;
  			};
  			/Element:\s*([0-9]+).*/ && ($1 > 0) && do {
  			    #print "Element found in line:$_\n";
  				#if(8 == $number_of_nodes) {
  					$read_node_names = 1;
  					$element_name = $1;
  				#}
  				last SWITCH;
  			};
  		}}
  	}  else {
  		# eat lines until we find "Nodes:"
  		#print "line: $_\n";
  		#while ((defined ($line = <EXFILE>)) && !($line =~ /Nodes:/)){print "rejected line:".$line;}
  		next unless($_ =~ /Nodes:/);
  		my $line;
  		my @pos; # array of node names
  		if ((defined ($line = <EXFILE>)) && (my @nodenames = ($line =~ /(\d+)/g) ))
  		#if ((defined ($line = $_)) && (my @nodenames = ($line =~ /(\d+)/g) ))
  		{
  			if($debug){
  				carp "found ". scalar @nodenames ." nodes in $line\n";
  			}
  			# find the node ref for each node name
  			my @noderefs;
  			foreach my $name (@nodenames) {
  				# try to find the node in any of the NodeGroups
  				for(my $i=0; $i<@nodegroups;$i++){
  					$ngr = $nodegroups[$i];
  					last if $ngr->isNode($name);
  				}
  
  				# if the node is not found then fail
  				unless($ngr->isNode($name)) {
  					croak "Node $name not found."
  				} else {
  					my $thisnode = $ngr->getNode($name);
  					push @noderefs, $thisnode;
  					#if($debug){
  					#	carp "Element $element_name Node name $name has ref $thisnode";
  					#}
  				}
  			}
  
  			#finally create the element
  			my $elem=CmUtils::Objects::Element->new();
  			#$elem->debug(1);
  			$elem->name($element_name);
  			$elem->addNodes(@noderefs);
  			$group->addElement($elem);
  
  			# go back into 'looking for an element' state
  			$read_node_names = 0;
  		}
  	}
	}
	close EXFILE;
	return wantarray ? @Group : $Group[0];
}

sub writeExelem {
  carp "writeExelem() is not implemented because both cm and cmgui can write exelem files.";
}
