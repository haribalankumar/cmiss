eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

use strict;

BEGIN {
  die "Environment variable CMISS_ROOT not defined\n"
    unless exists $ENV{CMISS_ROOT};
}
use lib "$ENV{CMISS_ROOT}/cmiss_perl/lib";
use CmUtils::IpFileFix qw(fixfiles);

fixfiles(
  filter => sub {
    my $file = shift;

    {
      # KFA: Extended the exponential form with an orthotropic law as per Costa.
      $file->add_after(
        target => "   (3) Tong & Fung skin function\n",
        insert => "   (4) W=C1*exp(Q) where Q=C2*Ef^2+C3*Es^2+C4*En^2+2*C5*Efs*Esf+2*C6*Efn*Enf+2*C7*Esn*Ens\n",
        unless => "   (4) W=C1*exp(Q) where Q=C2*Ef^2+C3*Es^2+C4*En^2+2*C5*Efs*Esf+2*C6*Efn*Enf+2*C7*Esn*Ens"
      );
    }

    {
      # D Malcolm: Adding GRAD(v) parameter to ipmate for time-dependent
      #            advection-diffusion problems
      my $questionx = " Grad Mass flow in x direction is [1]:";
      my $questiony = " Grad Mass flow in y direction is [1]:";
      my $questionz = " Grad Mass flow in z direction is [1]:";
      my $questiontail = "   (1) Constant spatially\n   (2) Piecewise constant (defined by elements)\n   (3) Piecewise linear   (defined by nodes)\n   (4) Defined by grid points (CQ)\n   (5) Defined elsewhere by Gauss pt array (YG)\n 1\n The value is [0.0000E+00]: 0.00000E+00\n";
      
      my $one   = qr/( Mass flow in x direction is \[1\]:\n(\s*\(\d+\).*\n){5}\s*1\s*\n The value is \[.*\]:\s*.*\n)/;
      my $two   = qr/( Mass flow in y direction is \[1\]:\n(\s*\(\d+\).*\n){5}\s*1\s*\n The value is \[.*\]:\s*.*\n)/;
      my $three = qr/( Mass flow in z direction is \[1\]:\n(\s*\(\d+\).*\n){5}\s*1\s*\n The value is \[.*\]:\s*.*\n)/;
      my $four  = qr/ Grad Mass flow in x direction is \[1\]:\n(\s*\(\d+\).*\n){5}\s*1\s*\n The value is \[.*\]:\s*.*\n/;
      my $five  = qr/ Grad Mass flow in y direction is \[1\]:\n(\s*\(\d+\).*\n){5}\s*1\s*\n The value is \[.*\]:\s*.*\n/;
      my $six   = qr/ Grad Mass flow in z direction is \[1\]:\n(\s*\(\d+\).*\n){5}\s*1\s*\n The value is \[.*\]:\s*.*\n/;
      
      if (my ($oldstring) = $file->match(pattern => $one) and not($file->match(pattern => $four)))  {
        my $newstring = "$oldstring\n$questionx\n$questiontail";
        print "\033[32mX direction\033[0m ... ";
        $file->replace(target => $oldstring, replace => $newstring, unless => $questionx);
      }
      if (my ($oldstring) = $file->match(pattern => $two) and not($file->match(pattern => $five)))  {
        my $newstring = "$oldstring\n$questiony\n$questiontail";
        print "\033[32mY direction\033[0m ... ";
        $file->replace(target => $oldstring, replace => $newstring, unless => $questiony);
      }
      if (my ($oldstring) = $file->match(pattern => $three) and not($file->match(pattern => $six)))  {
        my $newstring = "$oldstring\n$questionz\n$questiontail";
        print "\033[32mZ direction\033[0m ... ";
        $file->replace(target => $oldstring, replace => $newstring, unless => $questionz);
      }
      1;
    }
    {
      # D Malcolm: Changing material parameters to be defined wrt fibre
      # ipmate for time-dependent advection-diffusion problems
      my $one   = qr/( Diffusion coeffient \(x\))/;
      my $two   = qr/( Diffusion coeffient \(y\))/;
      my $three = qr/( Diffusion coefficient \(z\))/;
      my $four  = qr/( Mass flow in x direction)/;
      my $five  = qr/( Mass flow in y direction)/;
      my $six   = qr/( Mass flow in z direction)/;
      my $seven = qr/( Grad Mass flow in x direction)/;
      my $eight = qr/( Grad Mass flow in y direction)/;
      my $nine  = qr/( Grad Mass flow in z direction)/;

      my $one_sub   = " Diffusion coefficient (fibre_1)";
      my $two_sub   = " Diffusion coefficient (fibre_2)";
      my $three_sub = " Diffusion coefficient (fibre_3)";
      my $four_sub  = " Mass flow (fibre_1)"; 
      my $five_sub  = " Mass flow (fibre_2)";
      my $six_sub   = " Mass flow (fibre_3)";
      my $seven_sub = " Grad Mass flow (fibre_1)";
      my $eight_sub = " Grad Mass flow (fibre_2)";
      my $nine_sub  = " Grad Mass flow (fibre_3)";
    
      my $change=0;
    
      if (my ($oldstring) = $file->match(pattern => $one))  {
        print "\033[31m 1\033[0m";
        $file->replace(target => $oldstring, replace => $one_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $two))  {
        print "\033[31m 2\033[0m";
        $file->replace(target => $oldstring, replace => $two_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $three))  {
        print "\033[31m 3\033[0m";
        $file->replace(target => $oldstring, replace => $three_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $four))  {
        print "\033[31m 4\033[0m";
        $file->replace(target => $oldstring, replace => $four_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $five))  {
        print "\033[31m 5\033[0m";
        $file->replace(target => $oldstring, replace => $five_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $six))  {
        print "\033[31m 6\033[0m";
        $file->replace(target => $oldstring, replace => $six_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $seven))  {
        print "\033[31m 7\033[0m";
        $file->replace(target => $oldstring, replace => $seven_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $eight))  {
        print "\033[31m 8\033[0m";
        $file->replace(target => $oldstring, replace => $eight_sub);
        $change=$change+1;
      }
      if (my ($oldstring) = $file->match(pattern => $nine))  {
        print "\033[31m 9\033[0m";
        $file->replace(target => $oldstring, replace => $nine_sub);
        $change=$change+1;
      }
      if ($change==0){
        print "\033[32m No Change\033[0m";
      }
      1;
    }
  }
)

# $file->delete(string => "");
# $file->replace(target => "", replace => "");
# $file->add_before(target => "", insert => "", unless => "");
# $file->add_after(target => "", insert => "", unless => "");
# $file->match(pattern => qr//);
