BEGIN 				{	FS="(,[ \\t]*)|([ \\t]+)"}
function get_name(long_name)	{	where = match(long_name,/\(/)
					if (where>0)	{
						module = substr(long_name,0,where-1)
					}
					else	{
						module = long_name
					}
}
function print_vars(num_vars)	{	bracket = 0
					for (count=3;count<=num_vars;count++) 	{
						if (bracket==0)	{
							printf "%s|%s|",FILENAME,module
						}
						printf "%s",$count
						if (match($count,/\(/)) 	{
							bracket = 1
						}
						if (match($count,/\)/)) 	{
							bracket = 0
						}
						if (bracket==1)	{
							printf ","
						}
						else	{
							printf "\n"
						}
					}
}
$2!="'" && in_dec==1		{	in_dec = 0
}
(($2=="INTEGER" || $2=="REAL*8") && $3!="FUNCTION") && (0==match($1,/C/))	{ 	print_vars(NF)
					in_dec = 1
}
($2=="'" && in_dec==1)	&& (0==match($1,/C/))	{	print_vars(NF)
}
$2=="SUBROUTINE"		{	get_name($3)					
}
$2=="PROGRAM"			{	get_name($3)					
}
$3=="FUNCTION"			{	get_name($4)					
} 
