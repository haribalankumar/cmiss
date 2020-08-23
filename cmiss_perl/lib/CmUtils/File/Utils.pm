package CmUtils::File::Utils;
require 5.006;

use strict;
use warnings;

=head1 CmUtils::File::Utils

Useful routines for parsing files.

=cut

BEGIN {  
  use Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK);

=head1 VERSION

1.01 (28 March 2001)

Moved 'list_to_string' to CmUtils::Utils as it is more globally useful than
just in parsing Cmiss files.

=cut

  $VERSION     = 1.01;
  @ISA         = qw(Exporter);
  @EXPORT      = qw(strtonum compare_arrays);
  @EXPORT_OK   = qw();
  
}

=head2 strtonum(string)

Returns a numeric representation of a string, or 'undef' if the string is 
not a number.  Handles Fortran strings with 'D'.

=cut

sub strtonum {
  use POSIX qw(strtod);
  my $str = shift;
  $str =~ tr/dD/eE/; # Cope with Fortran numbers!!!
  $str =~ s/^\s+//;
  $str =~ s/\s+$//;
  $! = 0;
  my($num, $unparsed) = strtod($str);
  if (($str eq '') || ($unparsed != 0) || $!) {
    return undef;
  } else {
    return $num;
  }
}

=head2 compare_arrays($array1, $array2)

Returns True (1) if two arrays (passed by reference) are identical 
(same elements and same number of elements), otherwise returns False (0).

=cut

sub compare_arrays {
  # Check for two arguments, which must both be arrays
  return 0 unless @_ == 2;
  my ($first, $second) = @_;
  return 0 unless ref($first) eq 'ARRAY' && ref($second) eq 'ARRAY';
  
  # no warnings;  # silence spurious -w undef complaints
  return 0 unless @$first == @$second; #size
  for (my $i = 0; $i < @$first; $i++) { #elements
    return 0 if $first->[$i] ne $second->[$i];
  }
  return 1;
}

1;
__END__
