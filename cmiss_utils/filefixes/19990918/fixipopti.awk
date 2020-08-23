#
# Leo Cheng, 7 November 1999
#
# Add a new question to limit the maximum number of minor
# iterations.
#
# Leo Cheng, 20 December 1999
#
# Adding Solution timing output
#
# Chris Bradley, 9 August 2000
#
# Fixing formats for activation time optimisation
#
# Chris Bradley, 10 August 2000
#
# Adding sum-of-squares and correlation coefficient weightings
#

BEGIN{FOUND_MAX_MAJOR=0;FOUND_ACTIVATION_OPTI=0}
{

if(($1=="(12)" && $2=="Activation") || ($1=="(9)" && $2=="Unused" ))
{
    if($1=="(12)" && $2=="Activation")
	{
	    FOUND_ACTIVATION_OPTI=1;
	}
  print$0
  getline;
  print$0  #the answer
  getline;
  if($1=="Specify" && $2=="output" && $4=="optimisation")
  {
  }
  else
  {
    printf(" Specify output from optimisation [1]: \n");
    printf("   (1) None\n");
    printf("   (2) Total Solution Time\n");
    printf("  *(3) Unused\n");
    printf("    2 \n");
  }
  print $0
  getline;
}


if($1=="Do" && $7=="maximum" && $10=="major")
{


  if($NF=="y" || $NF=="Y")
  {
      if($11=="iterations")
	  {
	      print $0;
	  }
      else
	  {
	      printf(" Do you wish to specify a maximum number of major iterations [N]? Y\n");
	  }
     getline;
     if($6=="[50]:")
	 {
	     print $0;
	 }
     else
	 {
	     printf(" Enter the maximum major iterations [50]: %i\n",$NF);
	 }
  }
  else
  {
      if($11=="iterations")
	  {
	      print $0;
	  }
      else
	  {
	      printf(" Do you wish to specify a maximum number of major iterations [N]? N\n");
	  }
  }

  getline;
  if($1=="Do" && $7=="maximum" && $10=="minor")
  {

      if($NF=="y" || $NF=="Y")
	  {
	      if($11=="iterations")
		  {
		      print $0;
		  }
	      else
		  {
		      printf(" Do you wish to specify a maximum number of minor iterations [N]? Y\n");
		  }
	      getline;
	      if($6=="[50]:")
		  {
		      print $0;
		  }
	      else
		  {
		      printf(" Enter the maximum minor iterations [50]: %i\n",$NF);
		  }
	  }
      else
	  {
	      if($11=="iterations")
		  {
		      print $0;
		  }
	      else
		  {
		      printf(" Do you wish to specify a maximum number of minor iterations [N]? N\n");
		  }
	  }
  }
  else
  {
    printf(" Do you wish to specify a maximum number of minor iterations [N]? N\n");
    print $0
  }
  getline;
}



if($1=="Specify" && $2=="whether" && $3=="optimising")
{


    print $0;
    printf("   (1) Material parameters\n");
    printf("   (2) Geometric parameters\n");
    printf("   (3) Micro-structure parameters\n");
    getline;	
    getline;	
    getline;	
}

if(FOUND_ACTIVATION_OPTI==1 && $1=="Enter" && $3=="width" && $6=="wavefront")
{
    print $0;
    getline;
    if($3=="weighting" && $5=="sum-of-squares")
	{
	    #Already done.
	}
    else
	{
	    printf(" Enter the weighting for sum-of-squares objective [1.0]: 1.0\n");
	    printf(" Enter the weighting for correlation coefficient objective [0.0]: 0.0\n");
	}
    
}

print $0;
  
}
END{}

