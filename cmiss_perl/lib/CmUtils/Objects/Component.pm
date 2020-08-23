package CmUtils::Objects::Component;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::debug;

=head1  CmUtils::Objects::Component

Components are the basic elements of Fields.

=head1 Component structure

=over 4

=item NAME

Name of the Component

=item INDEX

Index number for the value of this Component in the Node array.  Derivative
values follow this index, and then each version is repeated as required.

=item VERSIONS

Number of versions of this component.  Default is 1.

=item DERIVATIVES

Number of derivatives of this component.

=item DERIVNAMES

Names of the derivatives of this component.

=back

=head1 VERSION

0.51 (15 November 2000)
0.7 (11 December 2003)

=head1 METHODS

=cut

BEGIN {  
  our ($VERSION, @ISA);

  $VERSION     = 0.7;
  @ISA         = qw(CmUtils::debug);
}

=head2 B<new()>

Returns a new Component.

=cut

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self = {
    NAME        => '',
    INDEX       => 0,
    VERSIONS    => 1,
    DERIVATIVES => 0,
    DERIVNAMES  => [],
  };
  bless $self, $class;
  return $self;
}

=head2 B<copy()>

Returns an identical copy of an existing Component.

=cut

sub copy {
  my $self = shift;
  my $class = ref($self) || $self;
  
  my $new = $self->new();
  $new->name($self->name());
  $new->valIndex($self->valIndex());
  $new->versions($self->versions());
  $new->derivatives($self->derivatives());
  $new->setDerivNames($self->getDerivNames());

  bless $new, $class;
  return $new;
}

=head2 B<isSame(Component)>

Returns True if the Component is the same as the owner Component, otherwise
returns False.

=cut

sub isSame ($$) {
  use CmUtils::File::Utils qw(compare_arrays);

  my ($c1,$c2) = @_;
  return (
    $c1->name()        eq $c2->name()        &&
    $c1->valIndex()    == $c2->valIndex()    &&
    $c1->versions()    == $c2->versions()    &&
    $c1->derivatives() == $c2->derivatives() &&
    compare_arrays(\@{$c1->{DERIVNAMES}},\@{$c2->{DERIVNAMES}})
  );
}

=head2 B<name([string])>

Optional argument is name for Component.  Returns the name of the Component.

=cut

sub name {
  my $self = shift;
  croak "usage: thing->name(name)"    unless @_ <= 1;
  if (@_) { $self->{NAME} = shift }
  return $self->{NAME};
}

=head2 B<valIndex([number])>

Optional argument sets the index number for the value of this Component.
Returns the index value of the Component.

=cut

sub valIndex {
  my $self = shift;
  croak "usage: thing->valIndex(index)" unless @_ <= 1;
  if (@_) { $self->{INDEX} = shift }
  return $self->{INDEX};
}

=head2 B<versions([number])>

Optional argument sets the number of versions of this Component.
Returns the number of versions of the Component.

=cut

sub versions {
  my $self = shift;
  croak "usage: thing->versions(numversions)" unless @_ <= 1;
  if (@_) { $self->{VERSIONS} = shift }
  return $self->{VERSIONS};
}

=head2 B<derivatives([number])>

Optional argument sets the number of derivatives of this Component.
Returns the number of derivatives of the Component.

=cut

sub derivatives {
  my $self = shift;
  croak "usage: thing->derivatives(numderivatives)" unless @_ <= 1;
  if (@_) { $self->{DERIVATIVES} = shift }
  return $self->{DERIVATIVES};
}

=head2 B<addDerivName(string)>

Add a derivative name to the end of the list of derivative names for
this Component.  Returns the list of derivative names of the Component.

=cut

sub addDerivName {
  my $self = shift;
  croak "usage: thing->addDerivName(dname)" unless @_ == 1;
  if (@_) { push @{$self->{DERIVNAMES}}, shift }
  return @{$self->{DERIVNAMES}};
}

=head2 B<setDerivNames(list)>

Sets the list of derivative names for this Component to be the passed list
of strings.  Clears list if no arguments are passed. 
Returns the list of derivative names of the Component.

=cut

sub setDerivNames {
  my $self = shift;
  $self->{DERIVNAMES} = [];
  $self->addDerivName($_) foreach @_;
  return $self->getDerivNames();
}

=head2 B<getDerivNames()>

Returns the list of derivative names of the Component.

=cut

sub getDerivNames {
  my $self = shift;
  return @{$self->{DERIVNAMES}};
}

=head2 B<list(FH)>

List the contents of this Component to the filehandle I<FH>.  Set the 
"verbose" flag ($s->verbose(1)) to see extended information.

=cut

sub list {
  my $self = shift;
  my $FH = (@_) ? shift : \*STDOUT;

  printf $FH <<EOF, $self->name();
       Component          %s
EOF
  if ($self->verbose()) {
    printf $FH <<EOF, $self->valIndex(), $self->versions(), $self->derivatives();
         Value at index     %d
         # Versions         %d
         # Derivatives      %d
EOF
    if ($self->derivatives()) {
      printf $FH <<EOF, join ', ' => $self->getDerivNames();
         Derivative names   %s
EOF
    }
  }
  
}

1;

__END__
