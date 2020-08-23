BEGIN 	{	in_item = 0;
		print "2 module line";
		print "Incorrect Subscripts";
	}
	{
		if($1=="subscript")
		{
			print identifier[0],identifier[1],$0;
		}
		parse_file_module($0,identifier);
	}
