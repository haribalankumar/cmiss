BEGIN { basis = 0; }

{
    if( basis == -2 )
	{
	    if( match($0,"^ Do you want to set cross derivatives to zero") )
		{
		    basis = 0;
		}
	    else if( match($0,"^ Enter the node position indices") )
		{
		    saveline = $0;
		    getline;
		    if( match($0,"^ Enter the derivative order indices \\[11211222\\]") )
			print " Do you want to set cross derivatives to zero [N]? N";
		    print saveline;
		    basis = 0;
		}
	}
    else if( basis == -1 )
	{
	    if( match($0,"^ low and high order schemes") )
		basis = -2;
	}
    else if ( basis == 9 )
	{
	    if( NF == 1 && $1 == "5" )
		basis = -1;
	    else
		basis = 0;
	}
    else if( basis )
	basis = basis + 1;
    else if( $0 == " For basis function type 1 the type of nodal interpolation is [1]:" )
	basis = 1;

    print $0;
}
