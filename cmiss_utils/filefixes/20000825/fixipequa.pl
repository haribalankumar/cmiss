#!/usr/bin/perl -w
#
# This is a perl script to update ipequa files with changing commands
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/
#
# Created:
#   Shane Blackett 2 October 2000
# Modified from fixcom.pl:
#   Scott Marsden 11 October 2000
# 
# The -c option is really useful for the example tree.
# All you do is specify a path like examples/1/13 and then
# this program will automatically find and checkout any files
# of the required extensions, then it will run the script on them.
# If the file isn't changed by this script it is automatically released.
# This means that you are left with a tree containing just those files
# which need to be updated, which you can then inspect and commit.

use strict;
use File::Compare;
use File::Find;
use Getopt::Std;

my $FILESTOCHANGE = q/(\.ipequa|_ipequa\.cmiss|\.irequa|_irequa\.cmiss)$/;

#
# Location of the log file recording script usage
#
my $logfile = "$ENV{CMISS_ROOT}/filefixes/20000825/fixipequa.log";

sub update_file
  {
	  # Given the input filename and output filename this routine parses
	  # the given file and makes any required changes.
	  my $input = shift;
	  my $output = shift;
	  my $line;
	  my $date = localtime(time());

    # @options_to_confirm is an array holding menu names and options needed
    # for $lines_to_insert to be inserted above $insert_above_this_line
    my @options_to_confirm = ( ['Specify whether','2'],
                              ['Specify equation','9'] );
    my $lines_to_insert = <<EOF;
 Specify whether solution is by [4]:
   (1)
   (2)
   (3)
   (4) Collocation
   (5)
   (6) Grid-based Finite Element
    4
EOF
    my $insert_above_this_line = 
    " Is the basis function type for dependent variable 1";

    # @add_option_to_menu adds another option to an existing menu
    # (checks that there is already enough options already ie 5
    my @add_option_to_menu = 
    ('Specify whether solution is by','6','Grid-based Finite Element');


# ----------------------------------------------------------------------
    my $want_to_add_menu = 1;
    my $want_to_add_menu_option = 1;
    my $in_menu_for_confirm = 0;
    my $in_menu_for_add = 0;
    my $options_confirmed = 0; # counts confirmed options
    my $correct_no_options = 0;
    my $menu_name;
    my $total_no_regions = 1; # default to only one region
    my $current_region = 1; 

    my  $prev_line1 = ''; # stores previous lines
    my  $prev_line2 = '';
    my  $prev_line3 = '';
   
	  open(INPUT, "<$input") || die "Could not open file $input: $!\n";
	  open(OUTPUT, ">$output") || die "Could not open file $output: $!\n";

	  while (defined($line = <INPUT>))
	  {
      $line =~ m| Total number of regions :  ([0-9]+)| # regions defined
        && ($total_no_regions = $1);
      
      if ($line =~ m| Region number :  ([0-9]+)|) # new region
      {
        $current_region = $1; # re-initialise all variables
        $want_to_add_menu = 1;
        $want_to_add_menu_option = 1;
        $in_menu_for_confirm = 0;
        $in_menu_for_add = 0;
        $options_confirmed = 0; 
        $correct_no_options = 0;        
      } 
      
      if ($current_region <= $total_no_regions) # loop over regions
      {
        if ($want_to_add_menu || ($want_to_add_menu_option)) # need to add stuff
        {
          if ($in_menu_for_confirm || ($in_menu_for_add))
          {
            #make sure there are the correct number of options in menu
            if ($in_menu_for_add && ($line =~ m|\(([0-9]+)\)|)) 
            {
              $1 eq ($add_option_to_menu[1] - 1) && ($correct_no_options = 1);
              $1 >= $add_option_to_menu[1] && ($correct_no_options = 0);
            } 
            if ($line =~ m|^\s*([0-9]+)|)
            {
              my $option_number = $1;
              if ($in_menu_for_confirm)
              {
                if ($option_number eq 
                        $options_to_confirm[$options_confirmed][1])
                {
                  $options_confirmed++;
                  $in_menu_for_confirm = 0;
                }
                else # confirm failed
                {
                  $want_to_add_menu = 0; # don't want to add menu 
                  $in_menu_for_confirm = 0;
                }
              }
              elsif ($in_menu_for_add)
              {
                if ($correct_no_options)
                {
                  my $spaces = "  ";
                  ($add_option_to_menu[1] > 9) || ($spaces = "   ");
              
                  $line = $spaces . "(" . $add_option_to_menu[1] . ") " . 
                    $add_option_to_menu[2] . "\n" .$line;
              
                  $in_menu_for_add = 0;
                  $want_to_add_menu_option = 0;
                  $correct_no_options = 0;
                }
                else
                {
                  #weren't correct number of options
                  $in_menu_for_add = 0;
                  $want_to_add_menu_option = 0;
                } 
              }
            }
          } 
          elsif ($want_to_add_menu &&  # not in menu to confirm or add
            ($options_confirmed > $#options_to_confirm))
          {
            if ($line =~ m|$insert_above_this_line|)
            {
              # make sure menu hasn't already been added 
              if ($prev_line2 =~ m|Grid-based Finite Element| 
                    || ($prev_line3 =~ m|Grid-based Finite Element|))
              { 
                # menu already added
                $want_to_add_menu = 0
              }
              else
              {
                $line =~ s|^($insert_above_this_line)|$lines_to_insert\n$1| 
                    && ($want_to_add_menu = 0);
              }
            }
          }
          else # not in menu to confirm or add
          {
            if ($want_to_add_menu 
                  && ($options_confirmed <= $#options_to_confirm))
            {
              $menu_name = $options_to_confirm[$options_confirmed][0];
              $line =~ m|^\s*\Q$menu_name\E\s*\[.*| 
                          && ($in_menu_for_confirm = 1);
            }
            if ($want_to_add_menu_option)
            {
              $menu_name = $add_option_to_menu[0];
              $line =~ m|^\s*\Q$menu_name\E\s*\[.*| && ($in_menu_for_add = 1);
            }
          }
        }
        # store previous lines
        $prev_line3 = $prev_line2;
        $prev_line2 = $prev_line1;
        $prev_line1 = $line;
        print OUTPUT $line;
      }
	  }
	  close (INPUT);
	  close (OUTPUT);
  }

#==========================================================
# Normally you should not need to change things below here.
#==========================================================

my $cvs_start;
my $start_directory;
my $usage = "Usage: $0 [-h] [-c examples/path_in_cvs]\n";
my %options;
if (!(getopts('hc:', \%options)))
  {
    print $usage;
    exit;
  }

for (keys %options)
  {
    if (m%^h$%)
      {
		  print $usage;
		  exit;
      }
    elsif (m%^c$%)
      {
		  #Directory in cvs below which all files will be
		  #checked out and tested.
		  $cvs_start = $options{$_}; 
      }
  }

if ($cvs_start)
  {
	 #We are going to work with a cvs repository!
	 if (! ($cvs_start =~ m"^examples"))
		{
		  die "cvs path must begin with examples.";
		}
	 if (-d "examples")
		{
		  die "examples directory already exists in current working directory.";
		}
	 if (! $ENV{CVSROOT})
		{
		  die "you must have CVSROOT defined.";
		}
	 if (! -d "$ENV{CVSROOT}/$cvs_start")
		{
		  die "$cvs_start is not found in repository $ENV{CVSROOT}.";
		}
	 #Change the search string
	 $FILESTOCHANGE =~ s%\$$%,v\$%;
	 $start_directory = $ENV{PWD};
  }

# Find the files (excluding links) needing to be changed
#
if ($cvs_start)
  {
	 find(\&foreachfile, "$ENV{CVSROOT}/$cvs_start");

	 print "\n\n: CVS final update\n";
	 !system("cvs update -P $cvs_start") || die "Unable to prune directories from $cvs_start";
  }
else
  {
	 find(\&foreachfile, ".");
  }

sub foreachfile 
  {
	 # This routine is executed for every file in the tree.  It manages the
	 # temporary and old files and log file.
	 if ($_ =~ $FILESTOCHANGE)
		{
		  my $date = localtime(time());
		  my $login;
		  my $keep_ = $_; #Ensure we don't mess this up
		  my $filename;
	  
		  if (!($login = getlogin()))
			 {
				$login = "unknown";
			 }
		  
		  print "Converting file $File::Find::name ....";
		  if ($cvs_start)
			 {
				chdir $start_directory;
				$filename = $File::Find::name;
				$filename =~ s%,v$%%;
				$filename =~ s%^$ENV{CVSROOT}/%%;
				!system("cvs -Q checkout $filename") || die "Unable to checkout $filename";
			 }
		  else
			 {
				$filename = $_;
			 }
		  if (-s "$filename")
			 {
				if (update_file("$filename", "$filename.tmp"))
				  {
					 if (-s "$filename.tmp")
						{
						  if (compare("$filename", "$filename.tmp") == 1)
							 {
								#Files are different
								rename ("$filename", "$filename.old");
								rename ("$filename.tmp", "$filename");
								print " done.\n";
								if (-w $logfile)
								  {
									 if (open(LOGFILE, ">>$logfile"))
										{
										  print LOGFILE "$login: Updated $File::Find::name $date\n";
										  close(LOGFILE);
										}
								  }
							 }
						  else
							 {
								#No changes made
								unlink "$filename.tmp";
								if ($cvs_start)
								  {
									 !system("cvs -Q update -r0 $filename") || die "Unable to release $filename";
								  }
								print " no changes required.\n";
								if (-w $logfile)
								  {
									 if (open(LOGFILE, ">>$logfile"))
										{
										  print LOGFILE "$login: Parsed $File::Find::name $date\n";
										  close(LOGFILE);
										}
								  }
							 }
						}
					 else
						{
						  print "\n";
						  die "ERROR : Corrected file is zero size.. reverting.";
						}
				  }
				else
				  {
					 print "\n";
					 die "ERROR : Update failed.. reverting.";
				  }
			 }
		  else
			 {
				print " empty - no change.\n";
				if ($cvs_start)
				  {
					 !system("cvs -Q update -r0 $filename") || die "Unable to release $filename";
				  }
			 }
		  if ($cvs_start)
			 {
				chdir $File::Find::dir;
			 }
		  $_ = $keep_; #Ensure we don't mess this up
		}
	 elsif ($_ =~ m/^(CVS|Attic)$/)
		{
		  $File::Find::prune = 1;
		}
  }

