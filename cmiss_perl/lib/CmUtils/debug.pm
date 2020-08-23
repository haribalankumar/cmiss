package CmUtils::debug;
require 5.006;

use strict;
use warnings;
use Carp;

=head1 CmUtils::debug

Base class for the CMISS modules related to nodes, files etc.  Defines shared
methods and variables for debugging and help.  Any class or structure derived
from CMISS have these methods available.  This includes:

  CmUtils::Objects::Component
  CmUtils::Objects::Field
  CmUtils::Objects::FieldSet
  CmUtils::Objects::Node
  CmUtils::Objects::NodeGroup

=head1 METHODS

=head2 B<debug([number])>

Activates or deactivates a shared debug flag for the modules derived from
this class.  Returns the state of this flag.

=cut

our $Debugging = 0;

sub debug {
  my $self = shift;
  if (@_) {$Debugging = shift};
  return $Debugging;
}

=head2 B<verbose([number])>

Activates or deactivates a shared verbose flag for the modules derived from
this class.  Set/test this for printing out information.

=cut

our $Verbose = 0;

sub verbose {
  my $self = shift;
  if (@_) {$Verbose = shift};
  return $Verbose;
}

=head2 B<listMethods()>

Print out the list of available methods for this object to STDERR.
Gets list from POD docs, not from Package variables.

=cut

sub listMethods {
  no strict 'refs';
  my $self = shift;
  use CmUtils::Pod::Usage qw(pod2usage);
  my $name = ref $self;
  my $podoptions = {
    -exitval => -1,
    -verbose => -1,
    -topic   => 'METHODS',
    -input   => $name,
    -output  => \*STDERR,
    -headings_only => 2,
  };
  pod2usage($podoptions);
  my @parents = @{"${name}::ISA"};
  my $i=0;
  while (my $parent = $parents[$i]) {
    push @parents, @{"${parent}::ISA"};
    $i++;
  }
  foreach my $parent (@parents) {
    $podoptions = {
      -exitval => -1,
      -verbose => -1,
      -topic   => 'METHODS',
      -input   => $parent,
      -output  => \*STDERR,
      -headings_only => 2,
    };
    pod2usage($podoptions);
  }
}

=head2 B<help([method])>

Print out the POD documentation for this object to STDERR.  If METHOD is
given, only the help for this method is printed.

=cut

sub help {
  no strict 'refs';
  my $self = shift;
  my $name = ref $self;
  my $topic = @_ ? ".*@_.*" : $name;
  
  use CmUtils::Pod::Usage qw(pod2usage);
  my $podoptions = {
    -exitval => -1,
    -verbose => -1,
    -topic   => $topic,
    -input   => $name,
    -output  => \*STDERR,
  };
  pod2usage($podoptions);
  my @parents = @{"${name}::ISA"};
  my $i=0;
  while (my $parent = $parents[$i]) {
    push @parents, @{"${parent}::ISA"};
    $i++;
  }
  foreach my $parent (@parents) {
    $podoptions = {
      -exitval => -1,
      -verbose => -1,
      -topic   => $topic,
      -input   => $parent,
      -output  => \*STDERR,
    };
    pod2usage($podoptions);
  }
  $self->listMethods() unless @_;
}


=head1 SEE ALSO

L<CmUtils::File> for routines related to handling CMISS nodal files (exnode,
ipnode, ipdata, ipfiel).

=cut

1;

__END__
