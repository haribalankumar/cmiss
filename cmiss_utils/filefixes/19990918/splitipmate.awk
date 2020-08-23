BEGIN{found = 0;linecount = 1}
{
	if($1 == "Int." || $1 == "Ext.")
	{
		if($2 == "conductivity-fibre" || $2 == "conductivity-sheet" || $2 == "conductivity-cross" )
		{
			linecount = 0;
			print $0;
			count = 0;
			found = 1;
		}
	}
	else
	{
		if(found)
		{
			print $0;
			count++;
			if(count == 8)
			{
				found = 0;
			}
		}
	}

	if($1 == "Membrane" && $2 == "capacitance" && $3 == "Cm" && linecount)
	{
		linecount=0;
	}
	else
	{
		if(linecount)
		{
			print $0;
		}
	}
}
END{}