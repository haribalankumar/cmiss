
#
# Leo Cheng, 30 November 1999
#
# Add a new question to ipinve files to split 
#  potential and activation problems
#
# Leo Cheng 4 April 2000
# Adding transmebrane jump and resting potential
#
# John Bodley 8 May 2000
# Modifying potential inverse approach calls
#
# Greg Sands 21 July 2000
# Added question regarding residual calculation
#
# John Bodley 13 October 2000
# Modifying potential inverse approach calls


BEGIN{ADDED_POTENTIAL_ACTIVATION=0;

    getline;     #CMISS Version...
    print $0;    
    getline;     #Heading ...
    print $0;
    getline;     #Blank
    print $0;
    getline;     #Blank
    print $0;
 
  }
{
  if($2=="inverse" && $3=="approach:")
  {
    ADDED_POTENTIAL_ACTIVATION=1;
  }

  if($1=="Specify" && $2=="whether:" && ADDED_POTENTIAL_ACTIVATION=="0")
  {
    printf(" Specify inverse approach: [1]\n");
    printf("   (1) Potential imaging\n");
    printf("   (2) Activation imaging\n");
    printf("   (3) Unused\n");
    printf("    1\n\n");

  }

  if($1=="(1)" && $2=="Zero-crossing")
  {
    print $0; # (1)
    getline;
    print $0; # (2)
    getline;
    print $0; # (3)
    getline;
    print $0; # answer
    getline;
    print $0; # space
    getline;

    if($3!="resting")
    {
      printf(" Enter the resting potential [0.0] : 0.00000D+00\n")
      printf(" Enter the transmembrane jump (b) [1.0] : 0.10000D+01\n\n")

      print $0;
      getline;      
    }
  }
    
  if($1=="(4)" && $2=="Greensite")
  { 
    getline;
    if($1==4)
    { 
      ICOUPLING = 2
      printf("    2\n")
      getline
    }
    else
    {
      ICOUPLING = 1
    }  
  }  
  
  if($1=="Enter" && $3=="regularisation" && ICOUPLING > 0)
  {
    printf(" Specify any additional coupling [1]\n");
    printf("   (1) None\n");
    printf("   (2) Greensite\n");
    printf("    %d\n\n", ICOUPLING);
  }
  
  # GBS 21 July 2000
  if($1 == "(2)" && $2=="Surface" && $3=="Laplacian")
  {
    print $0;
    getline;
    print $0; # (3)
    getline;
    print $0; # answer
    if($1=="2")
    {
      getline;
      print $0; # parameter
    }
    getline;
    print $0;
    getline;
    if ($3!="residual")
    {
      printf(" Specify the residual function calculation  [1]\n");
      printf("   (1) Normal\n");
      printf("   (2) Optimised\n");
      printf("   (3) Unused\n");
      printf("    1\n");
    }
  }
  
  # default
  print $0
}
END{}
