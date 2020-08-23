# fixdescription.awk

# Fixes changes to cm commands from last global release

function create_string()
{
# KAT 14/3/00: Including spaces at beginning
  match($0,/^ */);
  OUTLINE = substr($0,1,RLENGTH) $1;
#    OUTLINE = $1;
  for(i = 2 ; i <= NF ; i++)
    {
      if( $i != "" )
	{
	  OUTLINE = OUTLINE " " $i;
	}
    }
  newlines[NUMBER_OF_NEW_LINES] = OUTLINE;
  NUMBER_OF_NEW_LINES++;
}


BEGIN{
  ADJUST_FOR_INTERPRETER = 1;
  split("",varlist);
  variableSeparators = "[ ;,.:=+-*/()[\\]^]"
}

{
  ORIGINAL_LINE = $0;
  NUMBER_OF_NEW_LINES = 0;
  CHANGED_MAIN_LINE = 0;
  ORIENTATION_FLAG = 0;

  gsub("\\.\\./images", "web_data/images", $0);
  gsub("generated/test_page.html", "", $0);
  gsub("examples/example_files", "examples", $0);
  if (! match ($0, "cmiss_input/draw.com"))
    {
      gsub("draw.com", "cmiss_input/draw.com", $0);
    }
  if (! match ($0, "cmiss_input/movie.com"))
    {
      gsub("movie.com", "cmiss_input/movie.com", $0);
    }
  if (! match ($0, "cmiss_input/draw_elec.com"))
    {
      gsub("draw_elec.com", "cmiss_input/draw_elec.com", $0);
    }
  gsub("cmiss_input/cmiss_input/movie.com", "cmiss_input/draw.com", $0);
  print $0;
}
