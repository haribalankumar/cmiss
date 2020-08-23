#!/usr/bin/perl


################################################
# Usage is fe2fz.pl fe04.f > fz04.f            #
#                                              #
# Created : 28/07/99 By : Matthew Stevenson    #
#                                              # 
# Modified :  6/3/01 By : Karl Tomlinson       #
# Modified :  27/01/2003 By : David Nickerson  #
# Modified :  08/09/2003 By : David Nickerson  #
#                                              #
# By :                                         #
################################################

          
use strict;

if( $#ARGV != 0 ){
  print STDERR "Usage : fe2fz.pl fe_file.f > ...\n";
  exit(1);
}

if( $#ARGV == 0){
  open(STDIN,"<$ARGV[0]") or die "can't open $ARGV[0]: $!\n ";
}

my @routines;
my @functions;

while(<STDIN>){
  if( m/^      SUBROUTINE (\w*?)(\n|\()/i ){
    push(@routines,$1);
  }
  if(m/^      (\S*) ?FUNCTION (.*?)(\n|\().*/i ){
    push(@functions,[$1,$2]);
  }
}

# foreach my $rout (@routines){
#   print $rout,"\n";
# }


my $module;
my $module_nofe;

$_= $ARGV[0];
m/(fe\w+\/\w+\.f)/;
$module = $1;
m/fe(\w+\/\w+\.f)/;
$module_nofe=$1;

print "C**** CMISS Module fz$module_nofe\n\n";

foreach my $rout (@routines){
  print <<FORTRAN;
      SUBROUTINE $rout(A,*)
      DIMENSION A(*)
      CALL FLAG_ERROR(0,'Link with module '//
     &  '$module: '//
     &  'need sub $rout')
      RETURN 1
      END

FORTRAN
#    print "      INCLUDE 'cmiss\$reference:cbdi02.cmn'\n";
}

foreach my $func (@functions){
    my ($type,$name) = @$func;
    my $result =
      $type =~ /^real/i ? # produce a NaN.
	'0.0/0.0'
      : $type =~ /^integer/i ?
	-1
      : $type =~ /^logical/i ?
	'.FALSE.'
      : # character?
	"'*'";
  print <<FORTRAN;
      $type FUNCTION $name()
      CALL FLAG_ERROR(0,'Link with module '//
     &  '$module: '//
     &  'need func $name()')
      $name=$result
      RETURN
      END

FORTRAN
}
