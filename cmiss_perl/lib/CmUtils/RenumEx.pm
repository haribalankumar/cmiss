eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
##!/usr/local/perl5.6/bin/perl
package CmUtils::RenumEx;

use strict;
use warnings;
BEGIN{push @INC, "$ENV{CMISS_ROOT}/cmiss_perl/lib"};


use CmUtils::InplaceEdit;

=head1 CmUtils::RenumEx

Renumbers a pair of F<.exnode> and F<.exelem> files with specified node and
element offsets.

=cut

exit( renumEx(@ARGV) ) unless caller;

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT);
  $VERSION     = 1.00;
  @ISA         = qw(CmUtils::Exporter);
  @EXPORT      = qw(renumEx);
}

my ($nodeoffset, $nodefile, $elemoffset, $faceoffset, $lineoffset, $elemfile, $verbose);
my ($face, $index, $node, @offset, $shape);
my $onlynodes=0;

sub nodeChange {
  /^\s*Node:\s+/g && s/\G(\d+)/($1+$nodeoffset)/e;
}

sub elemChange {
  # Node offset
	# LKC added \s* before the colon
	#  if (/^\s*Nodes:/) {
  if (/^\s*Nodes\s*:/) {
    $node = 1;
  } elsif (defined ($node)) {
    s/(\s*)(\d+)/$1 . ($2+$nodeoffset)/eg;
    $node = undef;  
  }
  # Determine first shape occurance
  if (/^\s*Shape.\s+Dimension/) {
    if (defined ($index)) {
      $index--;
    } else {
      $shape = 1;
    }  
  } elsif (defined ($shape)) {
    # Determine offset index
		# LKC added \s* before the colon
    #if (/^\s*Element:(.*)$/) {
		if (/^\s*Element\s*:(.*)$/) {
      my @values = (split ' ', $1);
      $index = 0;
      while ($values[$index] == 0) {
        $index++;
      }
      $shape = undef;
    }
  }
  # Element offset
	# LKC added \s* before the colon
  # /^\s*Element:\s+/g && s/\G(\s*)(\d+)/$1 . ($2 ? $2+$offset[$index] : 0)/eg;
  /^\s*Element\s*:\s+/g && s/\G(\s*)(\d+)/$1 . ($2 ? $2+$offset[$index] : 0)/eg;
	# LKC added \s* before the colon
  # if (/^\s*Faces:/) {
  if (/^\s*Faces\s*:/) {	
    $face = 1;
  } elsif (defined ($face)) {
    if (/^\s*\d+/) {
      s/(\s*)(\d+)/$1 . ($2 ? $2+$offset[$index + 1] : 0)/eg;
    } else {   
      $face = undef;
    }  
  }
}

=head1 renumEx

Takes the following arguments:

=over 4

=item nodefile

=item nodeoffset

=item elemfile

=item elemoffset

=item faceoffset

=item lineoffset

=item verbose

=back

Only the last argument is optional.  Arguments may either be passed as a list,
or as a hash array with the above names as the hash keys.

=cut

sub renumEx {
  local @ARGV;
  if (@_ && $_[0] =~ / /) {
    @ARGV = split ' ',$_[0];
  } elsif (@_ && ref($_[0]) =~ /HASH/) {
    @ARGV = @{$_[0]}{qw/nodefile nodeoffset elemfile elemoffset faceoffset lineoffset verbose/};
  } else {
    @ARGV = @_;
  }

# Ability to offset only nodes or both nodes and elements
  if (@ARGV != 2 && @ARGV < 6 ) {
    print STDERR <<EOF;
Usage:
renumEx NODEFILE NODEOFFSET ELEMFILE ELEMOFFSET FACEOFFSET LINEOFFSET [VERBOSE]
    NODEFILE    node file name
    NODEOFFSET  node offset
    ELEMFILE    element file name
    ELEMOFFSET  element offset
    FACEOFFSET  face offset
    LINEOFFSET  line offset
    VERBOSE
EOF
    exit(-1);
  }

  if(@ARGV == 2)
  {
    $onlynodes=1;
  }


  ($nodefile, $nodeoffset, $elemfile, $elemoffset, $faceoffset, $lineoffset, $verbose) = @ARGV;
  @offset = ($elemoffset, $faceoffset, $lineoffset);
  
  # Initialization
  $index = $node = $shape = undef;
  
  print STDERR <<EOF if defined $verbose;
Renumbering files $nodefile and $elemfile:
  Node offset    : $nodeoffset
  Element offset : $elemoffset
  Face offset    : $faceoffset
  Line offset    : $lineoffset
EOF
  
  inplaceEdit($nodefile, \&nodeChange);
  if(!$onlynodes)
  {
    inplaceEdit($elemfile, \&elemChange);
  }
}

1;
__END__
