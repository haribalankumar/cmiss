BEGIN{INCOMPANDINEXT=0;}
{

# 25/05/2000 Tim Kirk made this file to update finite elasticity .ipequa
# files. It adds a 6th option for finite elasticity problem type



        if(INCOMPANDINEXT == 1)
	{
		if($1=="(6)" && $2=="Compressible")
		{
		#file already up to date
		INCOMPANDINEXT=0;
		}
		else
		{
		printf("   (6) Compressible + fluid for lung\n");
		INCOMPANDINEXT=0;
		}
	}
        if($1=="(5)" && $2=="Incomp" && $4=="inextensible")
	{
		INCOMPANDINEXT=1;
	}	

	print $0;

}
END{}
