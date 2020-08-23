# KAT 14/3/00: Unit changes for conductance
BEGIN{
  replacement["Mho/mm"] = "mS/mm";
  multiplier["Mho/mm"] = 1e3;
  firstLineEnd = " is [1]:";
  valuePrompt = " The value is [";
  valuePromptLength = length(valuePrompt);
}

{
  numLines = 0;
# collect lines in line
  line[++numLines] = $0;
  if(match(line[1],/^ (In|Ex)t\. conductivity-(fibre|sheet|cross)/))
    { # first line of multiline question
      if(match(line[1],/\([^ )]*\)/) && substr(line[1],RSTART+RLENGTH) == firstLineEnd)
	{
	  unitPosition = RSTART + 1;
	  unitLength = RLENGTH - 2;
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
		      value = substr(line[numLines],valuePromptLength+defaultLength+2);
		      defaultLength--;
		      default = substr(line[numLines],valuePromptLength+1,defaultLength);
# change value line
		      sub(/[Dd]/,"e",value);
		      sub(/[Dd]/,"e",default);
		      value *= multiplier[unit];
		      default *= multiplier[unit];
		      line[numLines] = valuePrompt default "]: " value;
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
