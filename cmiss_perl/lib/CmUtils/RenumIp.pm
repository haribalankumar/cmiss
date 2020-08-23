eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
##!/usr/local/perl5.6/bin/perl
package CmUtils::RenumIp;

use strict;
use warnings;
BEGIN{push @INC, "$ENV{CMISS_ROOT}/cmiss_perl/lib"};

use CmUtils::InplaceEdit;

=head1 CmUtils::RenumIp

Renumbers a pair of F<.ipnode> and F<.ipelem> files with specified node and
element offsets.

=cut

exit( renumIp(@ARGV) ) unless caller;

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT);
  $VERSION     = 1.00;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(renumIp);
}

my ($nodeoffset, $nodefile, $elemoffset, $elemfile, $verbose);

sub nodeChange {
  s/^(\s*Node number .*:\s*)(\d+)/$1 . ($2+$nodeoffset)/e;
}

sub elemChange {
  s/^(\s*Element number .*:\s*)(\d+)/$1 . ($2+$elemoffset)/e;
  /Enter.*:/g && s/\G(\s*)(\d+)/$1 . ($2+$nodeoffset)/eg;
}

=head1 renumIp

Takes the following arguments:

=over 4

=item nodefile

=item nodeoffset

=item elemfile

=item elemoffset

=item verbose

=back

Only the last argument is optional.  Arguments may either be passed as a list,
or as a hash array with the above names as the hash keys.

=cut

sub renumIp {
  local @ARGV;
  if (@_ && $_[0] =~ / /) {
    @ARGV = split ' ',$_[0];
  } elsif (@_ && ref($_[0]) =~ /HASH/) {
    @ARGV = @{$_[0]}{qw/nodefile nodeoffset elemfile elemoffset verbose/};
  } else {
    @ARGV = @_;
  }
  	  
  if (@ARGV < 2) {
    print STDERR <<EOF;
Usage:
renumIp NODEFILE NODEOFFSET ELEMFILE ELEMOFFSET [VERBOSE]
    NODEFILE    node file name
    NODEOFFSET  node offset
    ELEMFILE    element file name
    ELEMOFFSET  element offset
    VERBOSE
EOF
    exit(-1);
  }

  ($nodefile, $nodeoffset, $elemfile, $elemoffset, $verbose) = @ARGV;
  print STDERR <<EOF if defined $verbose;
Renumbering files $nodefile and $elemfile:
  Node offset    : $nodeoffset
  Element offset : $elemoffset
EOF
  
  inplaceEdit($nodefile, \&nodeChange);
  if(defined $elemfile)
	{
    inplaceEdit($elemfile, \&elemChange);
  }
}

1;
__END__
