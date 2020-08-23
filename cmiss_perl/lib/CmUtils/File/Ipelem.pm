package CmUtils::File::Ipelem;
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

=head1 CmUtils::File::Ipnode

Routines for writing (and maybe reading in the future) F<ipelem> files.

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
  @EXPORT      = qw(&readIpelem &writeIpelem);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readIpelem>

=over 4

=item PARAMETERS$#

(filename)

=item USAGE

  @groups = readIpelem($filename, NodeGroup);
  $group = readIpelem("newfile", NodeGroup);

=item FUNCTION

Returns a list of ElementGroups read from F<filename> (or F<filename.ipelem>).
If a scalar is assigned the return value, it will contain the first ElementGroup
read.
The second parameter is a ref to a NodeGroup that contains all the nodes
named in the element file. The function will produce warnings if nodes are
not defined.

=back

=cut

sub readIpelem {

	my ($filename, $nodegroup) = @_;
	# TODO: check that the params are valid
	
	# Read from an .ipelem file ('-' means STDIN)
	if ($filename ne '-' && $filename !~ /\.ip/) {
		$filename .= ".ipelem";
	}

	my $group;

	my $debug = 0;

	open IPFILE, "<$filename" or croak "Cannot open file $filename";

	my $element_name;

	while (<IPFILE>) {
	#strip whitespace from beginning and end of line
	s/^\s+//;
	s/\s+$//;
		SWITCH: { for ($_) {
		
			/^Heading\s*:\s*(.*)$/ && do {
				# we have a new group -> create a new ElementGroup
				$group = CmUtils::Objects::ElementGroup->new($nodegroup);
				
				$group->debug(1) if ($debug!=0);

				$group->name($1);
				if ($group->debug()) {
					carp "Created new element group ". $group->name();
				}
				last SWITCH;
			};

      /^Element number .*\]:\s*(\d+)/ && do {
					$element_name = $1;
					my @nodenames;
    			#create the element
    			my $elem=CmUtils::Objects::Element->new();
    			$elem->debug() if ($debug!=0);
    			if($debug){
    				carp "Created element ".$elem." for $element_name\n";
    			}
			    $elem->name($element_name);
      		my $line;
      		while ((defined ($line = <IPFILE>)) && !($line =~ /Enter the/)){}
      		my @pos; # array of node names
      		my ($tmp, $nodeline) = split(/:/, $line, 2);
      		if(@nodenames = ($nodeline =~ /(\d+)/g) ){
      			if($debug){
      				carp "found ". scalar @nodenames ." nodes in $line\n";
      			}
      			# find the node ref for each node name
      			my @noderefs;
      			foreach my $name (@nodenames) {
      				# try to find the node in any of the NodeGroups
#      				for(my $i=0; $i<@nodegroups;$i++){
#      					$nodegroup = $nodegroups[$i];
#      					last if $nodegroup->isNode($name);
#      				}
      				# if the node is not found then fail
      				unless($nodegroup->isNode($name)) {
      					croak "Node $name not found."
      				} else {
      					my $thisnode = $nodegroup->getNode($name);
      					push @noderefs, $thisnode;
      					#if($debug){
      					#	carp "Element $element_name Node name $name has ref $thisnode";
      					#}
      				}
      			}
        		$elem->addNodes(@noderefs);
      		}
      		
			    if($debug){$elem->list();}     		
			
#    			# now read in the version numbers if they exist
      		while ((defined ($line = <IPFILE>)) && !($line =~ /^$/)) {
      			if($debug){print "got line: $line \n";}
        		my $pos = 0;
        		if(my @versdata = ($line =~ /(\d+)/g) ){
        		  my $occ = shift @versdata;
        		  my $nname = shift @versdata;
        		  my $njj = shift @versdata;
        		  my $version = pop @versdata;
        	
      				if($debug){print "got numbers: $occ $nname $njj $version \n";}
      				
      				# find the position index of this occurrence of the node in the element
        		  my $index = 0;
        		  my $nocc = 0;
        		  #$elem->occurrence($index);
        		  foreach my $nn (@nodenames){
        		    if($nname == $nn){
        		      $nocc++;
        		      if($nocc==$occ){ 
          		      $elem->nodeVersion($index, $version);
        		      }
        		    }
        		    $index++;
        		  } 
        		}
      		}
      		if($debug){$elem->list();}
        		
    		  $group->addElement($elem);					
				  last SWITCH;
			};
		}} # SWITCH			
	} # while <IPFILE>
	close IPFILE;
	return $group;
}

=head2 B<writeIpelem>

=over 4

=item PARAMETERS

(filename, elementgroups, nodegroups, options)

=item USAGE

  #$options{nodes} = $group->getElements(10..20);
  writeIpelem($filename, $elemgroups, $nodegroups, \%options);

=item FUNCTION

Writes an ElementGroup to an F<ipelem> file.  Options are passed through
the I<options> hash to control the subset of the ElementGroup which will be
printed.

Multiple groups are not supported by the backend and if multiple groups
are passed in then they are concatenated.

=item OPTIONS

=over 8

=item elements

An array of element structures which will be output.  A number of routines are
available to get nodes from a group (e.g. B<getELements>. 
See CmUtils::Obects::ElementGroup). If not set, all elements will be output.

=item comfile

If the comfile option is set then the function will also produce a cmiss
command file that contains commands for loading and grouping the elements.

=item noVersions

Set this option to disable outputting of multiple versions.

=back

=back

=cut

sub writeIpelem {
  use Data::Dumper;
  
	croak "not enough parameters" unless @_ >=3;
  #my ($filename, $elemgroups, $options) = @_;
  my ($filename, $elemgroup, $options) = @_;
  my $elemgroups = [ $elemgroup ];
  
	# TODO: validate the parameters

	if ($filename ne '-' && $filename !~ /\.ip/) {
		$filename .= ".ipelem";
	}

	# strip off the extension TODO: there must be a one liner for this...
	my $comfilename = $filename;
	$comfilename =~ s/(.*)\..*/$1\.com/;

	# TODO: validate the element and node groups arrays to make sure that they match

	# Use the name of the first group in the file header. not sure if this is actually
	# used for anything
	my $egroupname = $elemgroups->[0]->name();
	
	#my $egroupname = $elemgroup->name();

	# count the number of elements
	# TODO: This could be done afterwards and then inserted into the file.
	# That would avoid an extra loop. But CPU time is cheaper than debugging
	# and this is easier to understand...
	my $num_elements = 0;
	foreach my $egroup (@{$elemgroups}){
		$num_elements += $egroup->numberOfElements();
	}
	open IPFILE, ">$filename" or croak "Cannot open file $filename";
	if ($options->{comfile}) {
		open COMFILE, ">$comfilename" or croak "Cannot open file $comfilename";
	}

	print IPFILE <<END_HEADER;
 CMISS Version 2.0  ipelem File Version 2
 Heading: $egroupname

END_HEADER

	print IPFILE " The number of elements is [1]: ".$num_elements."\n\n";

	foreach my $egroup (@{$elemgroups}) {
		#carp $egroup->getElementNames($egroup->getElementNames());
		print $egroup->name()."\n";
		my @tmp = $egroup->getElementNames();
		if ((@tmp > 0) && ($options->{comfile})) {
			#print COMFILE "fem group elements ".join(',', $egroup->getElementNames())." as ".$egroup->name().";\n";
		  foreach my $group ($egroup->getGroupNames()){
			   print COMFILE "fem group elements ".join(',', @{$egroup->getGroup($group)})." as ".$group.";\n";
		  }
		}
		foreach my $elem ($egroup->getElements()) {
			print IPFILE " Element number [     1]:      ".$elem->name()."\n";
			my @pos = $elem->getNodeNames();
			my $numnodes = scalar @pos;
#			croak "Only elements with 8 nodes are supported" unless(@pos == 8);
#			print IPFILE <<END_ELEMENT;
# The number of geometric Xj-coordinates is [3]: 3
# The basis function type for geometric variable 1 is [1]:  1
# The basis function type for geometric variable 2 is [1]:  1
# The basis function type for geometric variable 3 is [1]:  1
# Enter the 8 global numbers for basis 1:  $pos[0] $pos[1] $pos[2] $pos[3] $pos[4] $pos[5] $pos[6] $pos[7]
#END_ELEMENT
      my $numxj = scalar keys %{$elem->node(0)->getValuesHash()->{coordinates}};
      print IPFILE " The number of geometric Xj-coordinates is [3]: $numxj\n";
      for(my $i=0; $i<$numxj; $i++){
        print IPFILE " The basis function type for geometric variable ".($i+1)." is [1]:  1\n";
      }
      print IPFILE " Enter the ".scalar @pos." global numbers for basis 1:";
      for(my $i=0; $i<scalar @pos; $i++){
        print IPFILE " ".$pos[$i];
      }
      print IPFILE "\n";
      
			if(!$options->{noVersions}){
			 	my @occurrences = $elem->getOccurrences();
			 	my @versions = $elem->getNodeVersions();
			 	my $index = 0;
			 	foreach my $num (@occurrences){
			 		my $nodename = $pos[$index];
			 		my $node = $elem->node($index);
			 		#
			 		# TODO: here we assume that all components have the same number of versions as 'x'
			 		#
			 		my $numversions=$node->fieldSet()->getField('coordinates')->getComponent('x')->versions();
			 		if($numversions>1){
				 			for(my $k=1;$k<4;$k++) {
			 				#print OUTPUT_FILE " The version number for occurrence $num of node $nodename, njj=$k is [ 1]:  $num\n";
			 				# The backend wants the numerical values formatted exactly as below.
			 				printf(IPFILE " The version number for occurrence %2d of node %5d, njj=%d is [ 1]: %2d\n",
			 						 $num, $nodename, $k, $versions[$index]);
			 			}
			 		}
			 		$index++;
			 	}
			 	#print IPFILE "\n";
			}
			print IPFILE "\n";
		}
	}

	close IPFILE;
	close COMFILE if ($options->{comfile});
}

1;

__END__
