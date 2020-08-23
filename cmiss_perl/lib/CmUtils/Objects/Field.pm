package CmUtils::Objects::Field;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::debug;

=head1  CmUtils::Objects::Field

Fields are a collection of one or more Components.

=head1 Field structure

=over 4

=item NAME

Name of the Field

=item TYPE

Field type

=item COORDSYS

Name of the coordinate system associated with the Field

=item COMPONENTS

A set of Components which are part of this Field

=back

=head1 VERSION

0.51 (20 February 2001)

=head1 METHODS

=cut

BEGIN {  
  our ($VERSION, @ISA);

  $VERSION     = 0.51;
  @ISA         = qw(CmUtils::debug);
}

=head2 B<new()>

Returns a new Field.

=cut

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self = {
    NAME          => '',
    TYPE          => '',
    COORDSYS      => '',
    COMPONENTS    => {},
  };
  bless $self, $class;
  return $self;
}

=head2 B<copy()>

Returns an identical copy of an existing Field.

=cut

sub copy {
  my $self = shift;
  my $class = ref($self) || $self;
  
  my $new = $self->new();
  $new->name($self->name());
  $new->type($self->type());
  $new->coordsys($self->coordsys());
  foreach my $component ($self->getComponents()) {
    $new->addComponent($component->copy());
  }

  bless $new, $class;
  return $new;
}

=head2 B<name([string])>

Optional argument is name for Field.  Returns the name of the Field.

=cut

sub name {
  my $self = shift;
  croak "usage: thing->name(name)"    unless @_ <= 1;
  if (@_) { $self->{NAME} = shift }
  return $self->{NAME};
}

=head2 B<type([string])>

Optional argument is type of Field.  Returns the type of the Field.

=cut

sub type {
  my $self = shift;
  croak "usage: thing->type(type)"    unless @_ <= 1;
  if (@_) { $self->{TYPE} = shift }
  return $self->{TYPE};
}

=head2 B<coordsys([string])>

Optional argument is coordinate system of Field.  Returns the coordinate
system of the Field.

=cut

sub coordsys {
  my $self = shift;
  croak "usage: thing->coordsys(coordsys)"    unless @_ <= 1;
  if (@_) { $self->{COORDSYS} = shift }
  return $self->{COORDSYS};
}

=head2 B<addComponent(Component)>

Adds a Component to the list of Components of this Field.  Returns the 
number of Components.

=cut

sub addComponent {
  my $self = shift;
  croak "usage: thing->addComponent(component)"    unless @_ == 1;
  my $component = shift;
  unless (exists $self->{COMPONENTS}{$component->name()}) {
    $self->{COMPONENTS}{$component->name()} = $component;
    if ($self->debug()) {
      carp "Added " . $component->name() . " to " . $self->name();
    }
  } else {
    carp $component->name() ." already exists in ". $self->name();
  }
  return $self->numberOfComponents();
}

=head2 B<updateComponent(Component)>

Updates a Component in the list of Components of this ComponentSet.  Returns the 
index of the last Component.

=cut

sub updateComponent {
  my $self = shift;
  croak "usage: thing->updateComponent(component)" unless @_ == 1;
  my $component = shift;
  if (exists $self->{COMPONENTS}{$component->name()}) {
    delete $self->{COMPONENTS}{$component->name()};
    $self->{COMPONENTS}{$component->name()} = $component;
    if ($self->debug()) {
      carp "Redefined " . $component->name() . " in " . $self->name();
    }
  } else {
    carp $component->name() ." does not exist in " . $self->name();
  }
  return $self->numberOfComponents();
}

=head2 B<getComponentNames([ComponentName(s)])>

If ComponentName(s) are specified, returns a list of names of all 
Components with these names that exist in the Field, in the order they 
are specified.  If no argument is given, returns a list of names of ALL
Components, sorted in alphabetical order.

=cut

sub getComponentNames {
  my $self = shift;
  return (@_) 
    ? grep { exists $self->{COMPONENTS}{$_} } @_ 
    : sort keys %{$self->{COMPONENTS}};
}

=head2 B<getComponents([ComponentName(s)])>

If ComponentName(s) are specified, returns a list of all Components with these
names that exist in the Field, in the order they are specified.  If no
argument is given, returns a list of ALL Components, sorted in
alphabetical order.

=cut

sub getComponents {
  my $self = shift;
  my @keys = $self->getComponentNames(@_);

  my @components = ();
  push @components, $self->{COMPONENTS}{$_} foreach @keys;

  return @components;
}

=head2 B<getComponentsMatching([Pattern])>

If Pattern are specified, returns a list of all Components matching /Pattern/i
that exist in the Field, in the order they are specified.  If no argument is
given, returns a list of ALL Components, sorted in alphabetical order.

=cut

sub getComponentsMatching {
  my $self = shift;
  my $match = shift || '';
  my @keys = grep { /$match/i } $self->getComponentNames();

  my @components = ();
  push @components, $self->{COMPONENTS}{$_} foreach @keys;

  return @components;
}

=head2 B<getComponent(ComponentName)>

Returns the Component named ComponentName from the Field.

=cut

sub getComponent {
  my $self = shift;
  croak "usage: thing->getComponent(name)" unless @_ == 1;
  my $name = shift;
  return $self->{COMPONENTS}{$name};
}

=head2 B<deleteComponent(ComponentName)>

Deletes the Component named ComponentName from the Field.  Returns the new
number of Components.

=cut

sub deleteComponent {
  my $self = shift;
  croak "usage: thing->deleteComponent(name)" unless @_ == 1;
  my $name = shift;
  if (exists $self->{COMPONENTS}{$name}) {
    delete $self->{COMPONENTS}{$name};
  }
  
  return $self->numberOfComponents();
}

=head2 B<numberOfComponents()>

Returns the number of Components defined in the Field.

=cut

sub numberOfComponents {
  my $self = shift;
  return scalar keys %{$self->{COMPONENTS}};
}

=head2 B<isSame(Field)>

Returns True if the Field is the same as the owner Field, otherwise
returns False.

=cut

sub isSame ($$) {
  use CmUtils::File::Utils qw/compare_arrays/;
  
  my ($f1,$f2) = @_;
  my $true = 1;
  ( $f1->name()     eq $f2->name()     &&
    $f1->type()     eq $f2->type()     &&
    $f1->coordsys() eq $f2->coordsys()) or return 0;
  my @f1 = sort keys %{$f1->{COMPONENTS}};
  my @f2 = sort keys %{$f2->{COMPONENTS}};
  compare_arrays(\@f1,\@f2) or return 0;
  foreach my $comp (@f1) {
    $f1->{COMPONENTS}{$comp}->isSame($f2->{COMPONENTS}{$comp}) or return 0;
  }
  return 1;
}

=head2 B<list(FH)>

List the contents of this Field to the filehandle I<FH> (default to STDOUT).
Set the "verbose" flag ($s->verbose(1)) to see extended information.

=cut

sub list {
  my $self = shift;
  my $FH = (@_) ? shift : \*STDOUT;

  printf $FH <<EOF, $self->name();
     Field                %s
EOF
  if ($self->verbose()) {
    printf $FH <<EOF, $self->type(), $self->coordsys(), ($self->numberOfComponents());
       Field type           %s
       Coordinate System    %s
       # Components         %d
EOF
  }
  foreach my $component ($self->getComponents()) {
    $component->list($FH);
  }
}

1;

__END__
