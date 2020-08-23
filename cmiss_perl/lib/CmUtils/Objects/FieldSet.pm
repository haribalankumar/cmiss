package CmUtils::Objects::FieldSet;
require 5.006;

use strict;
use warnings;
use Carp;

use CmUtils::debug;

=head1  CmUtils::Objects::FieldSet

FieldSets are a collection of one or more Fields.  Each Node has a FieldSet
associated with it which defines the values stored at the Node.

=head1 FieldSet structure

=over 4

=item FIELDS

A set of Fields which are part of this FieldSet

=item NUMVALS

For speed in reading in F<exnode> files, this is the total number of values
defined by all Fields in this FieldSet

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

Returns a new FieldSet.

=cut

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self = {
    FIELDS  => {},
    NUMVALS => 0,
  };
  bless $self, $class;
  return $self;
}

=head2 B<copy()>

Returns an identical copy of an existing FieldSet.

=cut

sub copy {
  my $self = shift;
  my $class = ref($self) || $self;
  
  my $new = $self->new();
  $new->numValues($self->numValues());
  foreach my $field ($self->getFields()) {
    $new->addField($field->copy());
  }

  bless $new, $class;
  return $new;
}

=head2 B<addField(Field)>

Adds a Field to the list of Fields of this FieldSet.  Returns the 
index of the last Field.

=cut

sub addField {
  my $self = shift;
  croak "usage: thing->addField(field)" unless @_ == 1;
  my $field = shift;
  unless (exists $self->{FIELDS}{$field->name()}) {
    $self->{FIELDS}{$field->name()} = $field;
    if ($self->debug()) {
      print "Added " . $field->name() . " to FieldSet";
    }
  } else {
    print $field->name() ." already exists in ". $self->name();
  }
  return $self->numberOfFields();
}

=head2 B<updateField(Field)>

Updates a Field in the list of Fields of this FieldSet.  Returns the 
index of the last Field.

=cut

sub updateField {
  my $self = shift;
  croak "usage: thing->updateField(field)" unless @_ == 1;
  my $field = shift;
  if (exists $self->{FIELDS}{$field->name()}) {
    delete $self->{FIELDS}{$field->name()};
    $self->{FIELDS}{$field->name()} = $field;
    if ($self->debug()) {
      carp "Redefined " . $field->name() . " in FieldSet";
    }
  } else {
    carp $field->name() ." does not exist in FieldSet";
  }
  return $self->numberOfFields();
}

=head2 B<deleteField(FieldName)>

Deletes Field named FieldName in the FieldSet.  Returns the total number of 
defined Fields in the FieldSet.

=cut

sub deleteField {
  my $self = shift;
  croak "usage: thing->deleteField(name)" unless @_ == 1;
  my $name = shift;
  delete $self->{FIELDS}{$name};
  
  return $self->numberOfFields();
}

=head2 B<getFieldNames([FieldName(s)])>

If FieldName(s) are specified, returns a list of names of all Fields with 
these names that exist in the FieldSet, in the order they are specified.  
If no argument is given, returns a list of names of ALL Fields, sorted in
alphabetical order.

=cut

sub getFieldNames {
  my $self = shift;
  return (@_) 
    ? grep { exists $self->{FIELDS}{$_} } @_ 
    : sort keys %{$self->{FIELDS}};
}

=head2 B<getFields([FieldName(s)])>

If FieldName(s) are specified, returns a list of all Fields with these names
that exist in the FieldSet, in the order they are specified.  If no argument is 
given, returns a list of ALL Fields, sorted in alphabetical order.

=cut

sub getFields {
  my $self = shift;
  my @keys = $self->getFieldNames(@_);

  my @fields = ();
  push @fields, $self->{FIELDS}{$_} foreach @keys;

  return @fields;
}

=head2 B<getFieldsMatching([Pattern])>

If Pattern is specified, returns a list of all Fields matching /Pattern/i that
exist in the FieldSet.  If no argument is given, returns a list of  ALL Fields,
sorted in alphabetical order.

=cut

sub getFieldsMatching {
  my $self = shift;
  my $match = shift || '';
  my @keys = grep { /$match/i } $self->getFieldNames();

  my @fields = ();
  push @fields, $self->{FIELDS}{$_} foreach @keys;

  return @fields;
}

=head2 B<getField(FieldName)>

Returns the Field named FieldName from the FieldSet.

=cut

sub getField {
  my $self = shift;
  croak "usage: thing->getField(name)" unless @_ == 1;
  my $name = shift;
  return $self->{FIELDS}{$name};
}

=head2 B<numberOfFields()>

Returns the number of Fields defined in the FieldSet.

=cut

sub numberOfFields {
  my $self = shift;
  return scalar keys %{$self->{FIELDS}};
}

=head2 B<isSame(FieldSet)>

Returns True if the FieldSet is the same as the owner FieldSet, otherwise
returns False.

=cut

sub isSame ($$) {
  use CmUtils::File::Utils qw/compare_arrays/;
  
  my ($fs1,$fs2) = @_;
  my @fs1 = sort keys %{$fs1->{FIELDS}};
  my @fs2 = sort keys %{$fs2->{FIELDS}};
  compare_arrays(\@fs1,\@fs2) or return 0;
  foreach my $field (@fs1) {
    $fs1->{FIELDS}{$field}->isSame($fs2->{FIELDS}{$field}) or return 0;
  }
  return 1;
}

=head2 B<numValues([number])>

Optional argument sets the number of values in this FieldSet.  Returns the
number of values in the FieldSet.

=cut

sub numValues {
  my $self = shift;
  croak "usage: thing->numValues(num)" unless @_ <= 1;
  if (@_) { $self->{NUMVALS} = shift }
  return $self->{NUMVALS};
}

=head2 B<list(FH)>

List the contents of this FieldSet to the filehandle I<FH> (default to STDOUT).
Set the "verbose" flag ($s->verbose(1)) to see extended information.

=cut

sub list {
  my $self = shift;
  my $FH = (@_) ? shift : \*STDOUT;

  printf $FH <<EOF, $self->numValues();
   FieldSet:
     Number of values     %d
EOF
  if ($self->verbose()) {
    printf $FH <<EOF, ($self->numberOfFields());
     Number of Fields     %d
EOF
  }
  foreach my $field ($self->getFields()) {
    $field->list($FH);
  }
}

1;

__END__
