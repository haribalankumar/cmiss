# KAT 14/3/00: Unit changes for capacitance and conductance
BEGIN{
  replacement["uMho/mm^2"] = "mS/mm^2";
  multiplier["uMho/mm^2"] = 1e-3;
  replacement["nf/mm^2"] = "uF/mm^2";
  multiplier["nf/mm^2"] = 1e-3;
  replacement["mMho/mm^2"] = "mS/mm^2";
  multiplier["mMho/mm^2"] = 1;
  firstLineEnd = " is [1]:";
  valuePrompt = " The value is [";
  valuePromptLength = length(valuePrompt);
}

{
  numLines = 0;
# collect lines in line
  line[++numLines] = $0;
  if(match(line[1],/^ [A-Z]/))
    { # first line of multiline question
      if(match(line[1],/\([^ )]*\) is \[1\]:/))
	{
	  unitPosition = RSTART + 1;
	  unitLength = RLENGTH - length(firstLineEnd) - 2;
	  unit = substr(line[1],unitPosition,unitLength);
	  if(unit in replacement)
	    {
	      for ( i = 0; i < 6; i++ ) # get options and answer
		{
		  getline line[++numLines];
		}
	      if(match(line[numLines],/ *1 */))
		{
		  getline line[++numLines]; # should contain value
		  if(substr(line[numLines],1,valuePromptLength) == valuePrompt && defaultLength = index(substr(line[numLines],valuePromptLength+1),"]:"))
		    {
		      if(multiplier[unit] != 1)
			{
# change value line
			  value = substr(line[numLines],valuePromptLength+defaultLength+2);
			  defaultLength--;
			  default = substr(line[numLines],valuePromptLength+1,defaultLength);
			  sub(/[Dd]/,"e",value);
			  sub(/[Dd]/,"e",default);
			  value *= multiplier[unit];
			  default *= multiplier[unit];
			  line[numLines] = valuePrompt default "]: " value;
			}
# change units in prompt
		      line[1] = substr(line[1],1,unitPosition-1) replacement[unit] ")" firstLineEnd;
		    }
		}
	    }
	}
    }
# output collected lines
  for ( i = 0; i++ < numLines; )
    {
      print line[i];
    }
}
