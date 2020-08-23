BEGIN	{ 	FS="|" 
		first = 1
}
{	if ($2!=prev_routine)	{
		prev_routine = $2
		if (first==0)	{
			printf "</UL>\n"
		}
# KAT 5Jan00: want calls/called-from info
		printf "<H3>Routine <A HREF=\"%sscripts/cmiss-lookup.cgi?name=%s\">%s</A> (%s)</H3>",ENVIRON["CMISS_WWW_URL"],$2,$2,$1
#		printf "<H3>Routine <A HREF=\"%sscripts/routine-browser.cgi?routine=%s\">%s</A> (%s)</H3>",ENVIRON["CMISS_WWW_URL"],$2,$2,$1
		printf "<UL>\n"
		first = 0
	}
	printf "<LI> %s does not match %s.\n",$3,$4
}
END	{	if (first==0)	{
			printf "</UL>\n"
		}
} 
