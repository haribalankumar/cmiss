eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
##!/product/cmiss/bin/mips-irix/perl
package CmUtils::CmAddVersions;

use strict;
BEGIN {
  use Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK);

  $VERSION   = 1.00;
  @ISA       = qw(Exporter);
  @EXPORT    = qw(&CmAddVersions);
  @EXPORT_OK = qw();
}

use CmUtils::File;

exit( &CmAddVersions(@ARGV) ) unless caller;

=head1 CmUtils::CmAddVersions

cmAddVersions - add versions to an ipnode file for specified nodes.

=head1 VERSION

1.0 (17 December 2003) GBS

=head1 USAGE

  > cmAddVersions nodefile versionsfile outputfile

  CmUtils::CmAddVersions("initial.ipnode", "versionchanges", "new.ipnode");

All three parameters are required.  The Versions file contains a list of pairs
where the first number is the node number to add versions to, and the second
number is how many versions this node should have.

If a node listed already contains versions, the original information is ignored,
and new versions are created from the first version of the node.

=cut

sub CmAddVersions {
  my $usage = "Usage:\n$0 nodefile versionfile newfile\n";
  die "$usage" unless (@_ == 3);
  
  my ($nodefile, $versionfile, $newfile) = @_;
  -f $nodefile or die "$usage  >>Node file doesn't exist\n";
  -f $versionfile or die "$usage  >>Version file doesn't exist\n";
  
  my $n = cmRead($nodefile);
  my %versions = map {split} `cat $versionfile`;
  
  my ($end, @fs, @index);
  
  $fs[1] = $n->getFieldSet(0);
  foreach my $field ( $fs[1]->getFields() ) {
    foreach my $comp ( $field->getComponents() ) {
      push @index, $comp->valIndex();
      $end = $comp->valIndex()+($comp->derivatives()+1)*$comp->versions();
    }
  }
  push @index, $end;
  foreach my $nodenum ( keys %versions ) {
    my $v = $versions{$nodenum};
    unless ($fs[$v]) {
      $fs[$v] = $fs[1]->copy();
      foreach my $field ( $fs[$v]->getFields() ) {
        foreach my $comp ( $field->getComponents() ) {
          $comp->valIndex($comp->valIndex()*$v);
          $comp->versions($v);
        }
      }
      $n->addFieldSet($fs[$v]);
    }
    my $node = $n->getNode($nodenum);
    $node->fieldSet($fs[$v]);
    my @vals = $node->getValues();
    my @newvals;
    foreach my $i (0 .. $#index-1) {
      push @newvals, @vals[$index[$i]..$index[$i+1]-1] foreach (1 .. $v);
    }
    $node->setValues(@newvals);
    print "  Changed node $nodenum\n";
  }
  
  cmWrite($newfile, $n, {format => 'ipnode'});
}

1;
__END__
