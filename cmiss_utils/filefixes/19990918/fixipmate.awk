BEGIN { fe50problem=0; found=0}

{
    	if( $1== "Stresses" && $2=="in" && $3=="constitutive" && $4=="law" && $5=="are" && $6=="referred" && $7=="to" )
	{
		fe50problem=1;
	}
    	if( $1== "Specify" && $2=="whether" && $3=="the" && $4=="spring" && $5=="stiffness" && $6=="for" && $8=="cavity" )
	{
		fe50problem=1;
	}
	
	if($1== "(5)" && $2=="Defined" && $3=="by" && $4=="Grid")
	{
		found=1;
	}
	else
	{
		found=0;
	}

	if(fe50problem)
	{
		if(  $1=="(4)" && $2=="Defined" && $3=="by" && $4=="Gauss" && $5=="points")
		{
    			print $0;
			if(found==0)
			{
				print "  (5) Defined by Grid points";
			}
		}
		else
		{
		print $0;	
		}
	}
	else
	{
		print $0;	
	}

}
