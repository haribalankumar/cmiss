function checkunits() # KAT 15/3/00: mA -> uA
{
  if(match($0,/\([^ )]*\) \[/))
    {
      unitPosition = RSTART + 1;
      unitLength = RLENGTH - 4;
      unit = substr($0,unitPosition,unitLength);
      if(unit in unitReplacement)
	{
	  defaultPosition = RSTART + RLENGTH;
	  if(defaultLength = index(substr($0,defaultPosition),"]:"))
	    {
	      defaultLength--;
	      default = substr($0,defaultPosition,defaultLength);
	      value = substr($0,defaultPosition+defaultLength+2);
	      {
		if(multiplier[unit] != 1)
		  {
# change value
		    sub(/[Dd]/,"e",value);
		    sub(/[Dd]/,"e",default);
		    default *= multiplier[unit];
		    value *= multiplier[unit];
		    value = " " value;
		  }
# change units in prompt
		$0 = substr($0,1,unitPosition-1) unitReplacement[unit] ") [" default "]:" value;
	      }
	    }
	}
    }
}

BEGIN{
  FOUND_BIDOMAIN=0;
  FOUND_TRA_POT=0;
  FOUND_TRA_FLX=0;
  FOUND_EXIT=0;
  pulseQuestion = " Enter current for pulse ";
  unitReplacement["mA/mm^3"] = "uA/mm^3";
  multiplier["mA/mm^3"] = 1e3;
}
{
	if($1 =="Do" && $2 =="you" && $3 =="want" && $7 =="boundary" && $8 =="extracellular" && $9 =="potentials?")
	{
		FOUND_BIDOMAIN=1;
	}
	if($1 =="Do" && $2 =="you" && $3 =="want" && $5 =="time" && $6 =="dependent" && $7 =="extracellular")
	{
		FOUND_BIDOMAIN=0;
	}
	if($7 == "boundary" && $8 == "transmembrane" && $9 == "potentials?")
	{
		FOUND_TRA_POT=1;
	}
	if($7 == "boundary" && $8 == "transmembrane" && $9 == "fluxes?")
	{
		FOUND_TRA_FLX=1;
	}
	if(FOUND_TRA_POT)
	{
		print(" Specify the type of intracellular b.c. [1]:");
		print("   (1) Zero transmembrane flux");
		print("   (2) Zero intracellular flux");
		print("    1");
		if($11 == "N" || $11 == "n")
		{
			FOUND_TRA_POT=0;
		}
		else
		{
			FOUND_TRA_POT=0;
			FOUND_EXIT=1;
		}
	}
	else if(FOUND_TRA_FLX)
	{
		if($11 == "N" || $11 == "n")
		{
			FOUND_TRA_FLX=0;
		}
		else
		{
			FOUND_TRA_FLX=0;
			FOUND_EXIT=1;
		}
	}
	else
	{
		if(FOUND_EXIT == 0)
		{
		  if(substr($0,1,length(pulseQuestion))==pulseQuestion)
		    {
		      checkunits();
		    }
		  print $0;
		}
	}
	if($1 == "Enter" && $2 == "collocation" && $3 == "point" && $4 == "#s/name" && $5 == "[EXIT]:" && $6 == "0")
	{
		FOUND_EXIT=0;
	}
}
END{
	if(FOUND_BIDOMAIN)
	{
		print(" Do you want any time dependent extracellular bcs [N]? N\n");
	}

}




