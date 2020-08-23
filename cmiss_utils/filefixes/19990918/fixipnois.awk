
#
# John Bodley, 29 May 2000
#
# Modifying noise mesurement from uV (in reality nV) to mV

BEGIN{
       getline;    #CMISS Version...
       print $0;
       getline;    #Heading ...
       print $0;
       getline;    #Blank
       print $0;
       getline;    #Blank
       print $0;
     } 
     {  
       if($0==" Specify absolute noise level (microV) [10]:")
       {
         getline;
         NOIS_LEV=$1
         getline;
         
         printf(" Specify absolute noise level (mV) [1]:\n")
         printf("    %f\n", NOIS_LEV/1e+6)
       }
       print $0     
     }
END  {  
     }   
