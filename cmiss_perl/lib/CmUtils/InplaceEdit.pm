eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
##!/usr/local/perl5.6/bin/perl
package CmUtils::InplaceEdit;

use strict;
use warnings;

=head1 inplaceEdit

Replicates the "perl -pi -e" command, without a separate process.  
Edits a file "inplace"

=cut

use Carp;
use POSIX qw(tmpnam);
use File::Basename qw/basename/;

use CmUtils::Exporter ();
our ($VERSION, @ISA, @EXPORT, @EXPORT_OK);

$VERSION     = 1.00;
@ISA         = qw(CmUtils::Exporter);
@EXPORT      = qw(inplaceEdit);
@EXPORT_OK   = qw();

=head2 USAGE

  inplaceEdit ($filename, $code);

where $filename is the file to edit, and $code is the action to perform on
each line before it is written e.g.

  $code = sub { s/\s+/ /g }

or

  sub dostuff { s/\s+/ /g }
  inplaceEdit ($filename, \&dostuff);

or

  inplaceEdit ($filename, sub { s/\s+/ /g });

=cut

sub inplaceEdit {
  my $filename = shift;
  my $code = shift;
  
  my $tmpfile = basename(tmpnam());
  
  rename $filename, $tmpfile or croak "Cannot rename file";
  
  open NEW, ">$filename" or croak "Cannot open file $filename";
  open OLD, "<$tmpfile" or croak "Cannot open file $tmpfile";
  
  while (<OLD>) {
    $code->($_);
    print NEW;
  }
  close NEW;
  close OLD;
  
  unlink $tmpfile;
}

1;

__END__

