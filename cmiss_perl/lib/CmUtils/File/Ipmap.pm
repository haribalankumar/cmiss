package CmUtils::File::Ipmap;
require 5.006;

use strict;
#use diagnostics;
use warnings;
use Carp;

use CmUtils::File::Utils;

use CmUtils::Objects::NodeGroup;
use CmUtils::Objects::FieldSet;
use CmUtils::Objects::Field;
use CmUtils::Objects::Component;
use CmUtils::Objects::Node;
use Tie::RefHash;

=head1 CmUtils::File::Ipmap

Routines for reading and writing F<ipmap> files.

=head1 VERSION

=head1 CHANGES

=head1 SUBROUTINES

=cut

BEGIN {  
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 0.1;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw( &writeIpmap );
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}


=head2 B<writeIpmap(filename)>

Determines the version mapping and creates an ipmap file. The elements must
have the version info already initialised either by readIpelem or addVersions().

=cut

sub writeIpmap {
  croak "usage: writeIpmap(filename, elementgroup [, options])" unless @_ >= 2;
  
  my ($filename, $egroup, $options) = @_;

  if ($filename ne '-' && $filename !~ /\.ip/) {
    $filename .= ".ipmap";
  }

  if (ref($egroup) =~ /ARRAY/) {
    $egroup = $egroup->[0]; # can only print one group to an ipmap file
  }
  
  use Data::Dumper;
  
  #if($self->debug()){print "writeIpmap: $filename \n"};
    
  my @prompt = (
    "Is the nodal position mapped out [N]? ",
    "Is the derivative wrt direction 1 is mapped out [N]? ",
    "Is the derivative wrt direction 2 is mapped out [N]? ",
    "Is the derivative wrt directions 1 & 2 is mapped out [N]? ",
    "Is the derivative wrt direction 3 is mapped out [N]? ",
    "Is the derivative wrt directions 1 & 3 is mapped out [N]? ",
    "Is the derivative wrt directions 2 & 3 is mapped out [N]? ",
    "Is the derivative wrt directions 1, 2 & 3 is mapped out [N]? "
  );

  #get a list containing the mapping info
  my @nodes;
  
  my $debug = 0;
  
  my $nodemap;
  # allow use of references as hash keys 
  tie %{$nodemap}, 'Tie::RefHash::Nestable';
  $nodemap = $egroup->get_mapping();
	  
  
  print "Ipmap::writeIpmap nodemap:\n";
  #print Dumper($nodemap);

  @nodes = sort {$a->name() <=> $b->name()} keys %{$nodemap}; 
  
  open IPFILE, ">$filename" or croak "Cannot open file $filename";

  print IPFILE <<END_HEADER;
 CMISS Version 2.0  ipmap File Version 1
 Heading:

 Define node position mapping [N]? N
END_HEADER
  print IPFILE " The number of nodes with special mappings is [    1]: "
    . ( scalar @nodes ) . "\n\n";
    
  foreach my $node ( @nodes ) {

    if ( !defined $node ) {
      croak "Node ".$node->name()." not found in NodeGroup '".$egroup->getNodeGroup()->name()."' ! (undefined ref)\n";
    }
    next if( !defined $nodemap->{$node});
    
    print IPFILE " Node number [    1]: ".$node->name()."\n";
    my $xj = 1;
    foreach my $c ( $node->fieldSet()->getField("coordinates")->getComponents() )
    {
      print IPFILE " For the Xj($xj) coordinate:\n";
      foreach my $version ( 1 .. $c->versions() ) {
        print IPFILE " For version number $version:\n";
        foreach my $deriv ( 1 .. $c->derivatives() ) {
          print IPFILE " $prompt[$deriv] ";
#          if(34==$node->name()){
#            #print Dumper($nodemap->{$node});
#            print $node->name()." $xj $version $deriv ".$nodemap->{$node}[$xj - 1]{$version}{$deriv}."\n";
#            #$nodemap->{$node}[$xj - 1]{$version}{$deriv}
#          }
          #if ( ( 1 != $version ) && ( defined $nodemap->{$node}[$xj - 1]{$version}{$deriv} ) )
          if (defined $nodemap->{$node}[$xj - 1]{$version}{$deriv} )
          {
          
            print IPFILE " Y \n";
            print IPFILE " Enter node, version, direction, derivative numbers to map to [1,1,1,1]: ";
            print IPFILE $nodemap->{$node}[$xj - 1]{$version}{$deriv}[0]->name() . " ";
            for my $n ( 1 .. 3 ) {
              print IPFILE $nodemap->{$node}[$xj - 1]{$version}{$deriv}[$n] . " ";
            }
            print IPFILE "\n";
            print IPFILE " Enter the mapping coefficient [1]: 1.00000E+00 \n";
          } else {
            print IPFILE " N\n";
          }
        }
      }
      $xj++;
    }
    print IPFILE "\n"; 
  }
}

1;

__END__

