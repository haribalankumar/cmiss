
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



BEGIN{ADDED_POTENTIAL_ACTIVATION=0;

    getline;     #CMISS Version...
    print $0;    
    getline;     #Heading ...
    print $0;
    getline;     #Blank
    print $0;
    getline;     #Blank
    print $0;
    
    # Conversion of stabilisation methods
    convert[1] = 1;
    convert[2] = 4;
    convert[3] = 5;
    convert[4] = 2;
    convert[5] = 6;
    convert[6] = 7;
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
    
  if($1=="(3)" && $2=="Greensite")
  {
    getline;
    if($1==3) 
    {     
      printf("    2\n\n");
      printf(" Enter the type of inversion stabilisation scheme [3]\n");
      printf("   (1) None\n");
      printf("   (2) Truncated SVD\n");
      printf("   (3) Tikhonov\n");
      printf("   (4) Greensite\n");
      printf("    4\n\n");
      
      printf(" Specify any additional contraints [1]\n");
      printf("   (1) None\n");
      printf("  *(2) Surface Gradient Operator\n");
      printf("  *(3) Surface Laplacian Operator\n");
      printf("    1\n\n");
      
      printf(" Enter the regularisation scheme [1]\n");
      printf("   (1) Generalised Cross Validation (GCV) Criterion\n");
      printf("   (2) L-curve Criterion\n");
      printf("   (3) Picard Criterion\n");
      printf("   (4) Quasi-optimality Criterion\n");
      printf("    3\n\n");
      getline;
    }   
  } 

  if($1=="(3)" && $2=="Tikhonov")
  {
    print $0
    getline;
    if($1!="(4)")
    {
      printf("   (4) Greensite\n");
      ISTABILISE=$0
      print ISTABILISE 
    }
    else
    {
      print $0
    }
    getline; 
  }  
  
  if(ISTABILISE && $3=="regularisation") 
  { 
    if(ISTABILISE==2) # TSVD
    {
      NUM = 5
    }
    else if(ISTABILISE==3) # Tikhonov
    {
      NUM = 7
    }    
    for(i=1;i<=NUM;i++) getline;  
    IREGULARISE=$0
    getline;
      
    printf(" Enter the regularisation scheme [1]\n");
    printf("   (1) Generalised Cross Validation (GCV) Criterion\n");
    printf("   (2) L-curve Criterion\n");
    printf("   (3) Picard Criterion\n");
    printf("   (4) Quasi-optimality Criterion\n");
    printf("   (5) Optimal Criterion\n");
      
    if(ISTABILISE==3)
    { 
      printf("   (6) CRESO Criterion\n");
      printf("   (7) Zero-crossing Criterion\n");
    }
    printf("    %i\n", convert[int(IREGULARISE)]);
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
