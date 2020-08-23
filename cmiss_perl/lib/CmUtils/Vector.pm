package CmUtils::Vector;

use strict;
use CmUtils::Exporter ();
our ($VERSION, @ISA, @EXPORT, @EXPORT_OK);

=head1 CmUtils::Vector

Provides simple 3d vector functionality including, addition, subtraction,
scalar multiplication, cross product, norm, plot(print)  
and converting to a perl list.

=cut

$VERSION     = 1.00;
@ISA         = qw( CmUtils::Exporter );
@EXPORT      = qw( crossp dotp norm list );
@EXPORT_OK   = qw();

=head1 VERSION

1.00 (23 October 2002)

=cut

use overload
    '+' => \&myadd,
    '*' => \&mymult,
    '/' => \&mydiv,
    '-' => \&mysub,
    '""' => \&myprint;

=head1 B<new(array)>

Constructor to create a new Vector object.  Takes an array of
3 numbers to populate the vector coordinates.

    $v1 = new CmUtils::Vector( @array );
    $v2 = new CmUtils::Vector( (1,2,3) );

=cut

sub new
{
    my($class, @v) = @_;
    my $self = \@v;        
    bless $self, $class;
    return $self;
}


=head1 B<list()>

Returns the vector components as an array

    @array = $v1->list();

=cut

sub list
{
    my $self = shift;
    return @$self;
}

=head1 B<myprint()>

Evaluates the vector in a string context;

    @v1->myprint();

=cut

sub myprint
{
    my $self = shift;
    return join ", ", @$self;
}    

=head1 B<myadd>

Implementation of vector addition for overloaded operator +

=cut

sub myadd
{
    my $self = shift;
    my $x = shift;
    my($i, @z);

    for( $i = 0; $i < 3; $i++) {
	$z[$i] = $$self[$i] + $$x[$i];
    }

    return new CmUtils::Vector($z[0],$z[1],$z[2]);
}

=head1 B<mymult>

Implementation of vector scalar multiplication for overloaded operator *

=cut

sub mymult
{
    my $self = shift;
    my $x = shift;
    my($i, @z);

    for( $i = 0; $i < 3; $i++) {
	$z[$i] = $$self[$i] * $x;
    }

    return new CmUtils::Vector($z[0],$z[1],$z[2]);
}

=head1 B<mydiv>

Implementation of vector scalar division for overloaded operator /

=cut

sub mydiv
{
    my $self = shift;
    my $x = shift;
    my($i, @z);

    for( $i = 0; $i < 3; $i++) {
	$z[$i] = $$self[$i] / $x;
    }

    return new CmUtils::Vector($z[0],$z[1],$z[2]);
}

=head1 B<mysub>

Implementation of vector subtraction for overloaded operator -

=cut

sub mysub
{
    my $self = shift;
    my $x = shift;
    my($i, @z);

    for( $i = 0; $i < 3; $i++) {
	$z[$i] = $$self[$i] - $$x[$i];
    }

    return new CmUtils::Vector($z[0],$z[1],$z[2]);
}

=head1 B<crossp($v1,$v2)>

Returns the vector cross product of two vectors

    $v3 = crossp($v1,$v2);

=cut

sub crossp
{
    my $self = shift;
    my $x = shift;
    my($i, @z);

    $z[0] = $$self[1] * $$x[2] - $$self[2] * $$x[1];
    $z[1] = -1 * ($$self[0] * $$x[2] - $$self[2] * $$x[0]);
    $z[2] = $$self[0] * $$x[1] - $$self[1] * $$x[0];

    return new CmUtils::Vector($z[0],$z[1],$z[2]);
}

=head1 B<dotp($v1,$v2)>

Returns the vector dot product of two vectors

    $v3 = dotp($v1,$v2);

=cut

sub dotp
{
    my $self = shift;
    my $x = shift;
    my($product);

    $product = $$self[0] * $$x[0] + $$self[1] * $$x[1] + $$self[2] * $$x[2];

    return ($product);
}

=head1 B<norm($v1)>

Returns the norm of a vector

    $norm = norm($v1);

=cut

sub norm
{
   my $self = shift;
   my (@z, $i, $l);
   $l = $self->length;

   for($i = 0; $i < 3; $i++) {
       $z[$i] = $$self[$i] / $l;
   }
   return new CmUtils::Vector($z[0], $z[1], $z[2]);
}

=head1 B<length()>

Returns the length of a vector

    $length = $v1->length;

=cut

sub length
{
    my $self =shift;
    my $length;
    $length = ($$self[0]**2 + $$self[1]**2 + $$self[2]**2)**0.5;
    return $length;
}

=head1 B<plot()>

Prints out the vector components, useful in debugging

    @v1->plot();

=cut

sub plot
{
    my $self = shift;
    print "($self)\n";
}    

1;


