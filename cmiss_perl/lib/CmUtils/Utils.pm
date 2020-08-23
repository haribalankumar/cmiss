package CmUtils::Utils;
require 5.006;

use strict;
use warnings;

=head1 CmUtils::Utils

Useful routines for use in cmiss perl scripts and com files.

=cut

BEGIN {  
  use CmUtils::Exporter ();

=head1 VERSION

=over 4

=item 1.03 (10 April 2001)

runningCmgui() and runningCm() functions now give correct information.

=item 1.02 (2 April 2001)

Added runningCmgui() and runningCm() functions.

=item 1.01 (28 March 2001)

Moved list_to_string() from CMISS::File::Utils as it is more globally useful
than just in parsing Cmiss files.

Added doCmiss().

=back

=cut

  our $VERSION     = 1.03;
  our @ISA         = qw(CmUtils::Exporter);
  our @EXPORT      = qw(); # nothing by default
  our @EXPORT_OK   = qw(
    doCmiss 
    list_to_string
    runningCmgui
    runningCm
  );
}

=head1 SUBROUTINES

=head2 B<doCmiss(command_string(s))>

Executes a set of CMISS commands sequentially.  Obvious method for calling is:

  doCmiss(<<EOF);
    fem define node;p
    fem define elem;p
EOF

which will execute both commands.  Alternate method is:

  doCmiss("fem define node;p","fem define elem;p");
  doCmiss(@command_list);

which may be useful if a series of commands has been defined in an array

=cut

sub doCmiss {
  foreach (@_) {
    # break each multiline string up into separate command lines
    foreach (/^.*?$/mg) {
      $_ = substr($_,0,index($_,'#')) if /#/;  # strip out comments
      chomp;                            # cmiss doesn't like newlines
      next if /^\s*$/;                  # don't execute blank lines
      cmiss::cmiss($_);
    }
  }
}

=head2 B<list_to_string(@list)>

Converts a LIST of numbers to the shortest string describing that list.  e.g.

  list_to_string(4,5,15,6,7,8,9,12,14,8..10)

is converted to

  4..10,12,14,15

List is sorted prior to condensing.  Duplicate entries are ignored.  May also
be called with an array reference.

=cut

sub list_to_string {
  my @list = sort {$a <=> $b} (ref $_[0] eq 'ARRAY') ? @{shift()} : @_;
  return "" unless @list;
  
  my $old = shift @list;
  my $string = "$old";
  my $start = $old;
  while (my $i = shift @list) {
    next if $i == $old;
    unless ($i == $old+1) {
      $string .= ($old == $start) ? ",$i" : 
        ($old == $start+1) ? ",$old,$i" : "..$old,$i";
      $start = $i;
    }
    $old=$i;
  }
  $string .= ($old == $start) ? "" : ($old == $start+1) ? ",$old" : "..$old";

  return $string;
}

=head2 B<runningCmgui()>

Returns the name of the cmgui process if Cmgui is running, undefined otherwise.

=cut

sub runningCmgui {
  return grep { /^cmgui/ } getexecname();
}

=head2 B<runningCm()>

Returns the name of the cm process if Cm is running, undefined otherwise.

=cut

sub runningCm {
  return grep { /^(cm\b|cm[^g])/ } getexecname();
}

sub getexecname {
  my @plist = `ps -o args=`;  # all current processes with full name of program
  return map {
    /cm/m && do {             # if this process matches cm, cmiss, cmgui etc
      $_ = (split)[0];        #   extract only the program name (first in list)
      $_ = readlink while -l; #   follow if a symbolic link
      s{\S*/}{};              #   remove path
      $_                      # and return...
    } || ();
  } @plist;
}

1;

__END__
