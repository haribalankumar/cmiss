# fixcom.awk

# Fixes changes to cm commands from last global release
# Originally written SAB 17 June 1999
# GRC 20 August 1999. Added window background texture distortion fix
# KAT 13 March 2000. ! comments changed to # for interpreter
# KAT 22 March 2000. ! user variables and if constructs for interpreter
# SAB 30 May 2000. Quote "5*5*5" so that the interpreter doesn't evaluate it.
# SAB 2 June 2000. Changes for the Perl interpreter.
# JMB 12 June 2000. Removed time flags from "fem evaluate zeroxing" 

function window_background_texture_distortion_fix()
{
# gfx mod win 1 back ... tex_distortion CX CY K1 MAX_PIXELS_PER_POLYGON ...
# Distortion centre and k1 are now stored with the texture in physical coords.
# Changes above to:
# gfx mod win 1 back ... max_pixels_per_polygon MAX_PIXELS_PER_POLYGON ...
  for (j = i ; j <= (NF-4) ; j++)
  {
    if (tolower(substr($j,0,5)) == "tex_d")
    {
      $j = "max_pixels_per_polygon";
      j++;
      $j="";
      j++;
      $j="";
      j++;
      $j="";
      j++;
      CHANGED_MAIN_LINE = 1;
    }
  }
}

function glyph_edit_mode_to_select_mode_fix()
{
# From gfx mod g_e NAME:
# Glyph_edit_mode no longer exists; use select_mode instead
# Convert edit_off, edit_position, edit_select, edit_vector to select_on.
# Note that select mode is now available to all settings, but glyph edit mode
# was only available to node_points, data_points and element_points.
  if ((tolower(substr($i,0,4)) == "elem")||
      (tolower(substr($i,0,4)) == "node")||
      (tolower(substr($i,0,4)) == "data"))
    {
      for (j = i+1 ; j <= NF ; j++)
	{
	  if(tolower(substr($j,0,6)) == "edit_o")
	    {
	      CHANGED_MAIN_LINE = 1;
	      $j = "no_select";
	    }
	  else if((tolower(substr($j,0,6)) == "edit_p")||
		  (tolower(substr($j,0,6)) == "edit_s")||
		  (tolower(substr($j,0,6)) == "edit_v"))
	    {
	      CHANGED_MAIN_LINE = 1;
	      $j = "select_on";
	    }
	}
    }
}

function discretization_evaluation_escape()
{
# From gfx mod g_e NAME:
# Quotes discrtisations using * so that the interpreter doesn't
# evaluate them 5*5*5
    for (j = i+1 ; j <= NF ; j++)
    {
      if ((tolower(substr($j,0,3)) == "siz") || 
	  (tolower(substr($j,0,3)) == "sca") || 
	  (tolower(substr($j,0,3)) == "dis") || 
	  (tolower(substr($j,0,3)) == "ele"))
      {
         j++;
	 if ((!match ($j, /^["'].*["']$/)) && (match ($j, /\*/)))
	 {
            CHANGED_MAIN_LINE = 1;
	    $j = "\"" $j "\"";
	 }
      }
    }
}

function detect_interpreter()
# KAT 13/3/00:  Checks to see if it looks like there are interpreter specific
# commands.  If there are then it can be assume that syntax does not need to
# be changed for the interpreter.
{
  if((abbrev($1,"assign",2) && (abbrev($2,"variable",2) || abbrev($2,"function",2)) && NF >= 4 && $4 !~ /^[^"]*#/) || tolower($0) ~ /^ *((if|select|while) *\(|repeat)/)
    ADJUST_FOR_INTERPRETER = 0;
}

function comment_char_change(following_char)
# KAT 13/3/00:  ! comments changed to # for interpreter.
{
  if(match($0,/^[^"!#]*("[^"]*"[^"!#]*)*!/))
    {
# != never used to work after `if'.
#  # ! had special meaning after `if'.  Check for another !
#        if(tolower($1) == "if" && $3 == "!=") match($0,/^[^!]*![^!]*!/);
#        if(RLENGTH > 0)
#  	{
      following_char = substr($0,RLENGTH+1,1);
      if(following_char == "!") # change ## to !!
	$0 = substr($0,1,RLENGTH-1) "##" substr($0,RLENGTH+2);
      else if(following_char == "#") # change !# to #!
	$0 = substr($0,1,RLENGTH-1) "#!" substr($0,RLENGTH+2);
      else # change ! to #
	$0 = substr($0,1,RLENGTH-1) "#" substr($0,RLENGTH+1);
#	  CHANGED_MAIN_LINE = 1;
#  	}
    }
}

function variable_substitution(string, # 1 parameter
			       comment,position)
# does variable substitution in string and returns the result
{
# truncate at comment
  if(match(string,/^[^"#]*("[^"#]*"[^"#]*)*#/))
    {
      comment = substr(string,RLENGTH);
      string = substr(string,1,RLENGTH-1);
    }

  for (variable in varlist)
    {
# repeat in case of multiple matches
      while(position = match(string,"(^|" variableSeparators ")" variable "($|" variableSeparators ")"))
	{
	  position = position + index(substr(string,position,RLENGTH),variable) - 1;
	  string = substr(string,1,position-1) "$" substr(string,position);
	  CHANGED_MAIN_LINE = 1;
	}
    }

  string = string comment;

  return string;
}

function change_variable_syntax(keyword1,varname,remainder) # no parameters
# KAT 13/3/00:  Changes `(de)assign varname' to `(de)assign variable varname'
# and replaces usage of varname with $varname.
{
  if(match($0,/^ *[Aa][Ss][Ss]?[Ii]?[Gg]?[Nn]? +/)) #abbrev($1,"assign",2))
    {
      keyword1 = substr($0,1,RLENGTH-1);
      remainder = substr($0,RLENGTH+1);
      if(match(remainder,/[ 	=]+/))
	{
	  varname = substr(remainder,1,RSTART-1);
	  remainder = variable_substitution(substr(remainder,RSTART+RLENGTH));
# quoting things that are not numbers but look like numbers
	  if(match(remainder,/^0[0-9]+[^ 	]*/))
	    remainder = "\"" substr(remainder,1,RLENGTH) "\"" substr(remainder,RLENGTH+1);
# checking for assignment already made
#!!! doesn't consider if constructs
	  if(varname in varlist && varlist[varname] == "")
	    {
	      $0 = substr(keyword1,1,match(keyword1,/[^ ]/)-1) varname " = " remainder;
	      CHANGED_MAIN_LINE = 1;
	    }
	  else if(varname ~ varname) #varname matches itself as a regexp
	    {
#store variable name
	      varlist[varname] = "";
#add in variable keyword
	      $0 = keyword1 " variable " varname " " remainder;
	      CHANGED_MAIN_LINE = 1;
	    }
	}
    }
  else if(match($0,/^ *[Dd][Ee][Aa][Ss]?[Ss]?[Ii]?[Gg]?[Nn]? /)) #abbrev($1,"deassign",3))
    {
      if($2 in varlist)
	{
	  if(varlist[$2] == "do")
	    { #it's a do variable.  `deassign' command no longer needed
	      $0 = "";
	    }
	  else
	    { #add in variable keyword
	      $0 = substr($0,1,RLENGTH) "variable" substr($0,RLENGTH);
	    }
	  delete varlist[$2];
	  CHANGED_MAIN_LINE = 1;
	}
    }
  else if(match($0,/^[Dd][Oo] +/)) #abbrev($1,"do",2))
    {
      remainder = substr($0,RLENGTH+1);
      if(match(remainder,/[ =]+/))
	{
	  varname = substr(remainder,1,RSTART-1);
	  if(varname in varlist && varlist[varname] == "")
	    { # the variable has already been assigned - deassign it
	      newlines[NUMBER_OF_NEW_LINES] = "deassign variable " varname;
	      NUMBER_OF_NEW_LINES++;
	    }
	  if(varname ~ varname) #varname matches itself as a regexp
	    {
#store variable name
	      varlist[varname] = "do";
	    }
	}
    }
  else
    {
      $0 = variable_substitution($0);
    }
}

function change_if_constructs(position) # no parameters
# KAT 13/3/00:  Adds () to if and elseif and changes <> to != and = to ==
{
  if(match($0,/^ *([Ee][Ll][Ss][Ee])?[Ii][Ff] /))
    {
# opening (
      $0 = substr($0,1,RLENGTH-1) "(" substr($0,RLENGTH);
      if(match($0,/\( +[$A-Za-z][A-Za-z0-9_]* *=/))
	{ # insert =
	  position = RSTART + RLENGTH - 1;
	  $0 = substr($0,1,position) substr($0,position);
	}
#        if(position = index($2,"="))
#  	{ # insert =
#  	  $2 = substr($2,1,position) substr($2,position);
#  	}
#        else if($3 == "=")
#  	{ # insert =
#  	  $3 = "==";
#  	}
      else if(match($0,/\( +[A-Za-z_$]+ +<>/))
	{ # insert =
	  position = RSTART + RLENGTH;
	  $0 = substr($0,1,position-3) "!=" substr($0,position);
	}
#        else if($3 == "<>")
#  	{
#  	  $3 = "!=";
#  	}
# closing )
      if(match($0,/ +#/))
	$0 = substr($0,1,RSTART) ")" substr($0,RSTART);
      else
	$0 = $0 " )";

      CHANGED_MAIN_LINE = 1;
    }
}

function detect_perl_interpreter()
# SAB 2/6/00:  Looks for Perl specific interpreter commands.  If it finds one then there is no need to process the rest of the file.
{
  if(abbrev($1,"use",3))
    ADJUST_FOR_PERL_INTERPRETER = 0;
}

function change_assign_command()
{
  if((abbrev($1,"assign",2)) && (abbrev($2,"variable",2)))
  {
	 if (NF > 3)
	 {
		 # Leave the change variable to fix the assign command.
		 $1 = $3;
		 $2 = "=";
		 for (i = 4 ; i <= NF ; i++)
		 {
		    if (abbrev($i,"unmodifiable",3) || abbrev($i,"modifiable",3))
			 {
				$i = "";
				NF = NF - 1;
			 }
			 $(i-1) = $i;
		 }
		 $NF = "";
		 NF = NF - 1;
	 }
	 else
	 {
		 $1 = "";
		 $2 = "";
		 $3 = "";
	 }
	 CHANGED_MAIN_LINE = 1;
  }
}

function change_deassign_command()
{
  if(abbrev($1,"deassign",2))
  {
	  for (i = 1 ; i <= NF ; i++)
	  {
		  $i = "";
	  }
	  CHANGED_MAIN_LINE = 1;
  }
}

function change_variable_assignment()
{
   match_string = "^((\\$[A-Za-z0-9_]+)|[0-9/<>\\+\\-*\\\\%!\\~\\|\\(\\)\\{\\}\\$\\^\\&\\.=%]|([A-Za-z_]+\\())+$";
   second_match_string = "[0-9]\\.\\.[0-9]";
	if (! match($1, "^\\$") && (match($1, "=") || match($2, "=")))
	{
		if (2 == split($1, fields, "="))
		{
		   quote = 0;
			if (match (fields[2], "^\""))
			{
				$1 = "$" fields[1] "=" fields[2];
			}
			else
			{
			  if (!match(fields[2], match_string))
			  {
			     quote = 1;
			  }
			  if (match(fields[2], second_match_string))
			  {
				 quote = 1;
			  }
			  i = 2;
			  while (!quote && i <= NF)
			  {
				 if (!match($i, match_string))
				 {
					quote = 1;
				 }
				 if (match($i, second_match_string))
				 {
					quote = 1;
				 }
				 i++;
			  }
			  if (quote)
			  {
				  $1 = "$" fields[1] "= \"" fields[2];
				  $NF = $NF "\"";
			  }
			  else
			  {
				  $1 = "$" fields[1] "=" fields[2];
			  }
			}
			CHANGED_MAIN_LINE = 1;			 
		}
		else if ($2 == "=")
		{
		   quote = 0;
			if (match ($3, "^\""))
			{
				$1 = "$" $1;
			}
			else
			{
			   i = 3;
			   while (!quote && i <= NF)
			   {
				   if (!match($i, match_string))
					{
				 	   quote = 1;
				   }
					if (match($i, second_match_string))
					{
					   quote = 1;
					}
					i++;
			   }
				$1 = "$" $1;
				if (quote)
				{
				   $3 = "\"" $3;
				   $NF = $NF "\"";
				}
			}
			CHANGED_MAIN_LINE = 1;			 
		}
		else if (substr($2,1,1) == "=")
		{
		   quote = 0;
			if (match (substr($2, 2), "^\""))
			{
				$1 = "$" $1;
			}
			else
			{
			   if (!match(substr($2, 2), match_string))
				{
			      quote = 1;
				}
				if (match(substr($2, 2), second_match_string))
				{
					quote = 1;
				}
				i = 3;
			   while (!quote && i <= NF)
			   {
				   if (!match($i, match_string))
				   {
				 	   quote = 1;
				   }
					if (match($i, second_match_string))
					{
					   quote = 1;
					}
					i++;
			   }
				$1 = "$" $1;
			   if (quote)
			   {
				   $2 = "=\"" substr($2, 2);
					$NF = $NF "\"";
				}
			}
			CHANGED_MAIN_LINE = 1;			 
		}
	}
}

function remove_doubleslash_separator()
{
	for(i = 1 ; i <= NF ; i++)
	{  
		if (match($i, "//"))
		{
			if (match($i, ";"))
			{
				QUOTE = 0;
			}
			else
			{
			   QUOTE = 1;
			}

				number_of_terms = split ($i, terms, "//");
				if (number_of_terms > 1)
				{
				   for (j = 1 ; j <= (number_of_terms - 1) ; j++)
					{
# Don't worry about the last one as it is the end of the token
			         if ((pos = match(terms[j], "\\$[^ ;/]*$")) && 
							(!match (terms[j + 1], "^[ ]")))
			         {
                     terms[j] = substr(terms[j], 1, pos-1) "${" substr(terms[j], pos + 1) "}";
				      }
					}
					   if (QUOTE)
					   {
                     if (match (terms[1], "\"$"))
					      {
                        $i = substr(terms[1], 1, -1); 
                     }
					      else
					      {
                        $i = "\"" terms[1];
					      }
					   }
					   else
					   {
                     $i = terms[1];
					   } 
		            for(j = 2 ; j <= number_of_terms - 1 ; j++)
	    	         {
                     if( terms[j] != "" )
                     {
	                     $i = $i terms[j];
			            }
		            }
					   if (QUOTE)
					   {
                     if (match (terms[number_of_terms], "^\""))
					      {
                        $i = $i substr(terms[number_of_terms], 2); 
                     }
					      else
					      {
			               $i = $i terms[number_of_terms] "\"";
					      }
					   }
					   else
					   {
                     $i = $i terms[number_of_terms];
					   } 
					CHANGED_MAIN_LINE = 1;			 
				} 
		}
   } 
}

function change_do_to_for()
{
	if (abbrev($1,"do",2))
	{
# Try and get f90 parser do only
		if (2 == split($2, fields, "="))
		{
			newlines[NUMBER_OF_NEW_LINES] = "for $" fields[1] " ( " fields[2] " ) ";
			NUMBER_OF_NEW_LINES++;
			$1 = "{";
			$2 = "";
			CHANGED_MAIN_LINE = 1;			 
		}
		else if ($3 == "=")
		{
			newlines[NUMBER_OF_NEW_LINES] = "for $" $2 " ( " $4 " ) ";
			NUMBER_OF_NEW_LINES++;
			$1 = "{";
			$2 = "";
			$3 = "";
			$4 = "";
			CHANGED_MAIN_LINE = 1;			 
		}
		else if (abbrev($3, "=", 1))
		{
			newlines[NUMBER_OF_NEW_LINES] = "for $" $2 " ( " substr($3, 2) " ) ";
			NUMBER_OF_NEW_LINES++;
			$1 = "{";
			$2 = "";
			$3 = "";
			CHANGED_MAIN_LINE = 1;			 
		}
	}
	if (abbrev($1, "enddo", 5))
	{
		$1 = "}";
		CHANGED_MAIN_LINE = 1;			 
	}
}

function change_if()
{
	if (if_bracket_required)
	{
		if (!match ($1, "^ *{"))
		{
			newlines[NUMBER_OF_NEW_LINES] = "{";
			NUMBER_OF_NEW_LINES++;
		}
		if_bracket_required = 0;
	}
	if (match ($1, "^[Ii][Ff]"))
	{
		$1 = "if" substr($1, 3);
		if_bracket_required = 1;
		#Only if case changed. CHANGED_MAIN_LINE = 1;			 
	}
	if (abbrev($1,"else",2))
	{
		$1 = "} else";
		if_bracket_required = 1;
		CHANGED_MAIN_LINE = 1;			 
	}
	if (match ($1, "^[Ee][Ll][Ss][Ee][Ii][Ff]"))
	{
		$1 = "} elsif" substr($1, 7);
		if_bracket_required = 1;
		CHANGED_MAIN_LINE = 1;			 
	}
	if (abbrev($1, "endif", 5))
	{
		$1 = "}";
		CHANGED_MAIN_LINE = 1;			 
	}
}

function change_while()
{
	if (while_bracket_required)
	{
		if (!match ($1, "^ *{"))
		{
			newlines[NUMBER_OF_NEW_LINES] = "{";
			NUMBER_OF_NEW_LINES++;
		}
		while_bracket_required = 0;
	}
	if (match ($1, "^[Ww][Hh][Ii][Ll][Ee]"))
	{
	   $1 = "while" substr($1, 6);
		while_bracket_required = 1;
		#Only if case changed. CHANGED_MAIN_LINE = 1;		
	}
	if (abbrev($1, "endwhile", 5))
	{
		$1 = "}";
		CHANGED_MAIN_LINE = 1;			 
	}
}

function change_output()
{
  if (sub("^ *output", "print", $1))
  {
	 if (match($NF,"\"$"))
		{
		  $NF = substr($NF, 1, RSTART-1) "\\n\"";
		}
	 else
		{
		  $NF = $NF ".\"\\n\"";
		}
     CHANGED_MAIN_LINE = 1;			 
  }
}

function change_optimise()
{
	if (abbrev($1, "optimis", 3))
	{
		$1 = "optimise";
		CHANGED_MAIN_LINE = 1;			 
	}
}

function change_open()
{
	if (abbrev($1, "open", 2))
	{
		$1 = "open";
		CHANGED_MAIN_LINE = 1;			 
	}
}

function change_format()
{
  for (i = 1 ; i <= NF ; i++)
	{
	  while (match($i, /format\([^\)]*\)/))
	  {
		  format_start = RSTART;
		  format_length = RLENGTH;
		  comma_position = match(substr($i, format_start + 7, format_length - 8), ",");
		  variable_name = substr ($i, format_start + 7, comma_position - 1);
		  format = substr($i, format_start + 7 + comma_position, format_length - 8 - comma_position);
		  if (substr(format, 1, 1) == "F")
		  {
			  format = substr(format, 2) "f";
		  }
		  else if (substr(format, 1, 1) == "I")
		  {
			  format = substr(format, 2) "d";
		  }
		  else
		  {
			  #This isn't exhaustive but lets only implement things as necessary.
			  format = substr(format, 2) substr(format, 1, 1);
		  }
		  sub(/format\([^\)]*\)/,
			  "sprintf(\"%" format "\"," variable_name ")", $i);
		  CHANGED_MAIN_LINE = 1;			 
	  }
	}
}

function create_string()
{
  OUTLINE = $1;
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

function abbrev(short, long, minchar)
# KAT 13/3/00:  Determines whether short is an abbreviation of length at least
# minchar for long.  Long must be lower case.
{
  if(length(short) < minchar)
    return(0);
  else
    return(tolower(short) == substr(long,1,length(short)));
}

BEGIN{
  ADJUST_FOR_INTERPRETER = 1;
  ADJUST_FOR_PERL_INTERPRETER = 1;
  split("",varlist);
  variableSeparators = "[ ;,.:=+-*/()[\\]^]";
  if_bracket_required = 0;
  while_bracket_required = 0;
}

{
  ORIGINAL_LINE = $0;
  NUMBER_OF_NEW_LINES = 0;
  CHANGED_MAIN_LINE = 0;
  ORIENTATION_FLAG = 0;

  if(match($0,/^ */))
	 {
		SPACES = substr($0,1,RLENGTH);
	 }
  else
	 {
		SPACES = "";
	 }


  if(ADJUST_FOR_INTERPRETER)
    {
      detect_interpreter();
      if(ADJUST_FOR_INTERPRETER)
	{
	  comment_char_change();
	  change_variable_syntax();
	  change_if_constructs();
	}
    }

  if(ADJUST_FOR_PERL_INTERPRETER)
    {
      detect_perl_interpreter();
      if(ADJUST_FOR_PERL_INTERPRETER && (!match($0, "^#")) && (!match($0, "^ *$")))
		  {
			 #Preserve comments without parsing them
			 keep_NF = NF;
			 for (i = 1 ; i <= NF ; i++)
			 {
				 if (match ($i, "^#"))
				 {
					NF = i - 1;
				 }
			 }
			 remove_doubleslash_separator();
			 change_assign_command();
			 change_deassign_command();
			 change_do_to_for();
			 change_if();
			 change_while();
			 change_variable_assignment();
			 change_output();
			 change_optimise();
			 change_open();
			 change_format();

			 NF = keep_NF;
		  }
    }  

  i = 1;
  if(tolower($i) == "gfx")
    {
      i++;
# modify
      if(tolower(substr($i,0,3)) == "mod")
	{
	  i++;
# g_element
	  if(tolower(substr($i,0,3)) == "g_e")
	    {
	      i+=2;
	      glyph_edit_mode_to_select_mode_fix();
              discretization_evaluation_escape();
	    }
# window background
	  else if(tolower(substr($i,0,3)) == "win")
	    {
	      i+=2;
	      if (tolower(substr($i,0,4)) == "back")
		{
		  i++;
		  window_background_texture_distortion_fix();
		}
	    }
	}
    }

  # fem evaluate zeroxing     
  i = 1;
  if(tolower($i) == "fem")
  {
    i++;
    if(tolower(substr($i,0,2)) == "ev")
    {
      i++;
      if(tolower(substr($i,0,3)) == "zer")
      {
        j = i;
        i++;
        if(tolower(substr($i,0,7)) == "tstep_s")
        {
          i+=2;          
          CHANGED_MAIN_LINE = 1;
        }
        if(tolower(substr($i,0,7)) == "tstep_e")
        {
          i+=2;
          CHANGED_MAIN_LINE = 1;
        }
        
        if(CHANGED_MAIN_LINE)
        {
          if(tolower(substr($i,0,5)) == "class")
          {
            j++;
            $j = "class";
            j++;
            i++;
            $j = $i;
          }
          NF = j;
        }  
      }
    }
  }          
    
  if(CHANGED_MAIN_LINE || NUMBER_OF_NEW_LINES)
    {
      create_string();
      printf("\n# Command updated by %s on %s\n# Old command: %s\n",
	     PROGRAM_NAME, DATE, ORIGINAL_LINE);
      for(m = 0 ; m < NUMBER_OF_NEW_LINES ; m++)
	{
	  print(SPACES newlines[m]);
	}
      printf("\n");
    }
  else
    {
      print($0); 
    }
}
