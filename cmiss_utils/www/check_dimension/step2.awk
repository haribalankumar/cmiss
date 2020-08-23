BEGIN	{ FS="(\\()|(,)|(\\))|(\\|)"
}
FILENAME ~ "variables.check$"	{	present[$1] = 1
					declaration[$1] = $0
					for (count=2;count<NF;count++)	{
						indices[$1,count-1] = $count
					}
}
FILENAME ~ "all_variables$"	{	if (present[$3]==1)	{
						correct = 1
						for(count=4;count<NF;count++)	{
							if ($count!=indices[$3,count-3])	{
								correct = 0
							}
						}
						if (correct==0)	{
							if(($2!="FEM")&&($2!="PARALLEL_ELEM_STIFF"))	{
								printf "%s|%s\n",$0,declaration[$3]
							}
						}
					}
} 
