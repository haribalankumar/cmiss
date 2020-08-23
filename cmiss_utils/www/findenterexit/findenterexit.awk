BEGIN{
    current_subroutine=0;
}

$1=="SUBROUTINE" {
    if( current_subroutine )
    {
	printf "<LI>End statement not detected for subroutine %s\n",current_subroutine;
    }
    bracketpos=index($2,"(")
    if(bracketpos != 0)
    {
	current_subroutine=substr($2,1,bracketpos-1);
    }
    else
    {
	current_subroutine=$2
    }
    has_enters=0;
    called_exits=0;
    returned=0;
    routinename=0;
}

{
if( current_subroutine ) {
    if( $1=="END" || $1 == "end" ) {
	if( has_enters && !called_exits ) 
	{
	    printf "<LI>Subroutine %s does not have any exit calls\n",current_subroutine;
	}
	else if( !has_enters && called_exits )
	{
	    printf "<LI>Subroutine %s does not have an enters calls\n",current_subroutine;
	}
	current_subroutine = 0;
    } else if( $0 ~ /^ *[0-9]* +RETURN/ ) {
	if( called_exits && returned ) {
	    printf "<LI>Error in subroutine: %s -- RETURN without CALL EXITS\n",current_subroutine;
	} else if( has_enters && $1 ~ /^[0-9]+$/ ) {
	    printf "<LI>Error in subroutine: %s -- Label on RETURN provides path around CALL EXITS\n",current_subroutine;
	}
	returned = 1;
    } else if( $0 ~ / *PARAMETER\(ROUTINENAME/ ) {
	match($0,/'[^']*/);
	routinename = substr($0,RSTART+1,RLENGTH-1);
	if(routinename != current_subroutine)
	    {
		printf "<LI>Error in subroutine: %s -- ROUTINENAME=%s\n",current_subroutine,routinename;
	    }
    } else {
	callindex = 0;
	if($1=="CALL")
	    {
		callindex = 2;
	    }
	else if($2=="CALL" && substr($0,1,1) == " ") # not comment
	    {
		callindex = 3;
	    }
	if(callindex) {		
	    subroutine=substr($callindex,1,index($callindex,"(")-1);
	    if(subroutine=="ENTERS")
		{
		    has_enters=1;
		}
	    else if(subroutine=="EXITS") {
		if( called_exits && !returned ) {
		    printf "<LI>Error in subroutine: %s -- extra CALL EXITS without RETURN\n",current_subroutine;
		}
		called_exits = 1;
		returned = 0;
	    }
	    if(routinename)
		{
		    if(!index($0,"(ROUTINENAME"))
			{
			    printf "<LI>Error in subroutine: %s -- use ROUTINENAME\n",current_subroutine;
			}
		}
	    else
		{
		    split($callindex,name,"'");
		    if(name[2] != current_subroutine)
			{
			    printf "<LI>Error in subroutine: %s -- enter/exit name: %s\n",current_subroutine,name[2];
			}
		}
	}
    }
}
}
