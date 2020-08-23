BEGIN{found = 0}
{
	if($1 == "CMISS" && $2 == "Version" && $3 == "1.21" && $4 == "ipmate")
	{
		printf(" CMISS Version 1.21 ipcell File Version 1\n");
	}
        else
	{
		if($3 == "material" && $4 == "parameters" && $5 == "time-varying")
		{
		}
		else
		{
			if($1 == "Int." || $1 == "Ext.")
			{
				if($2 == "conductivity-fibre" || $2 == "conductivity-sheet" || $2 == "conductivity-cross" )
				{
					count = 0;
					found = 1;
				}
			}
			else
			{
				if(found)
				{
					count++;
					if(count == 8)
					{
						found = 0;
					}
				}
				else
				{
					print $0;
				}

			}
		}
	}
}
END{}
