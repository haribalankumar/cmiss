package CmUtils::File::Txtmap;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::File::Utils;

=head1 CmUtils::File::Ipelem

Routines for writing F<txtmap> files.

=head1 VERSION

0.1 () Glenn Ramsey

=head1 CHANGES

=head1 SUBROUTINES

=cut

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION     = 0.5;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(&readTxtmap &writeTxtmap);
  @EXPORT_OK   = qw();
  %EXPORT_TAGS = ();
}

=head2 B<readTxtmap>

=over 4

=item PARAMETERS$#

(filename)

=item USAGE

  @groups = readTxtmapFilename, Mapfile);
  $group = readTxtmap("newfile",Mapfile);

=item FUNCTION

The optional third parameter is the name of a map file that contains a mapping
relating the node names used in the Txtelem file to the node names used in the
Nodegroup.

=back

=cut


sub readTxtmap {
	use Data::Dumper;
	my ($file, $nodemapfile) = @_;
	# TODO: check that the params are valid
	
	# Read from a .txtmap file ('-' means STDIN)
	if ($file ne '-' && $file !~ /\.txt/) {
		$file .= ".txtmap";
	}

  my $map;
  
	my $debug = 0;

  #read the map file into a hash
  my %nmap;
  if(defined $nodemapfile){
  	if ($nodemapfile ne '-' && $nodemapfile !~ /\.txt/) {
  		$nodemapfile .= ".txtmap";
  	}
    open NODEMAPFILE, "<$nodemapfile" or croak "Cannot open file $nodemapfile";
    
    my $key;
    my $val;
    while (<NODEMAPFILE>) {
      my (@values) = split;
      $key = shift @values;
      $val = shift @values;
      $nmap{$key} = $val;   
    }
    close NODEMAPFILE;
  }
    
  open FILE, "<$file" or croak "Cannot open file $file";
  
  while (<FILE>) {
    next if /^#/; #ignore comments  
    my (@values) = split; # should contain 2 mapping groups
    next if 2 != @values; #ignore malformed lines
    my @tokens = split(/:/, $values[0]);
    
    my $fromnode = shift @tokens;
    $fromnode = 1.0 * $nmap{$fromnode} if (exists $nmap{$fromnode});
    
    my $fromversion = 1.0 * shift @tokens;
    
    my $fromderiv = 1.0 * shift @tokens;
    
    my @maptovalues = (split(/:/, $values[1]));
    my $tonode = shift @maptovalues;
    $tonode = 1.0 * $nmap{$tonode} if (exists $nmap{$tonode});

    my $toversion = 1.0 * shift @maptovalues;
    my $toderiv = 1.0 * shift @maptovalues;

    push @{ $map->{$fromnode}{$fromversion}{$fromderiv} },
        ($tonode, $toversion, $toderiv+1);

  }
  close FILE;
   
	return $map;
}


sub writeTxtmap {
  croak ( "writeTxtmap is not implemented yet, perhaps you could implement it?");
  
# NOT IMPLEMENTED

}


1;

__END__
