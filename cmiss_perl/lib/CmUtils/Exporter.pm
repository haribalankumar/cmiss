package CmUtils::Exporter;
require 5.006;

=head1 CmUtils::Exporter

Base routines to automatically add and provide help for CMISS modules.

=cut

use strict qw(vars);
eval 'use warnings';

BEGIN {
  use Exporter ();

=head1 VERSION

1.04 (13 September 2001)

=head1 CHANGES

1.04 - Added support of :DEFAULT flag

1.03 - Moved storage of defined subroutines back into this module

=cut
  our $VERSION = 1.04;
  our @ISA     = qw(Exporter);
}

# %SUBS is a private array containing a list of defined subroutines and is 
#     defined by the import routine
our %SUBS = ();

# Redefine the "import" routine to add defined subroutines to the SUBS hash
sub import {
  my @imports;
  my ($pkg, @list) = @_;
  if (@list) {
    foreach (@list) {
      if (/^:(\w+)/) {
        if ($1 eq 'DEFAULT') {
          push @imports, map {s/^&//;$_} @{"${pkg}::EXPORT"};
        } else {
          push @imports, map {s/^&//;$_} @{${"${pkg}::EXPORT_TAGS"}{$1}};
        }
      } else {
        push @imports, map {s/^&//;$_} $_;
      }
    }
  } else {
    @imports = map {s/^&//;$_} @{"${pkg}::EXPORT"};
  }
  foreach (@imports) {
    $SUBS{$_} = $pkg unless exists $SUBS{$_};
    print qq/Defining subroutine "$_" from $pkg\n/ if caller =~ /cmiss/;
  }
  goto &Exporter::import;
}

1;
