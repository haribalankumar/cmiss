#!/usr/bin/perl -w
#
# This is a perl script to update comfiles with changing commands
#
# Usage:
#   ${CMISS_ROOT}/cmiss_utils/filefixes/
#
# Created:
#   Shane Blackett 2 October 2000
# Modified:
#   Richard Christie 27 October 2000
#   (added fixes for composite and fibre_sheet_axes computed fields)
#   Richard Christie 2 November 2000
#   (added fix for gradient computed field)
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

my $FILESTOCHANGE = q/(\.com|_com\.cmiss|cmiss_test.out)$/;

#
# Location of the log file recording script usage
#
my $logfile = "$ENV{CMISS_ROOT}/filefixes/20000825/fixcom.log";
#
# The name of this script as written into the comfiles when changes
# are made.
#
my $scriptname = "fixcom.pl";

sub make_changes
  {
	 my $line = shift;
	 my $change = 0;

	 #Change the given line into the new line you want.
	 #The new $line variable may contain multiple lines
	 my $gfx_match_string = q/(gfx[^\#]*)/;
	 if ($line =~ m/$gfx_match_string/)
		{
		  my $gfx_command = $1;

      # "modify texture evaluate width_texels height_texels depth_texels format"
		#  -> "gfx modify texture specify_width specify_height specify_depth specify_format evaluate"
 		  my $mod_texture =
				q/((mod|cre)\w*\s+tex\w*)(.*)/;
			if ($gfx_command =~ m/$mod_texture/i)
			{
			   my $command = $1;
				my $parameters = $3;

				$parameters =~ s/(eval\w+\s+.*)(\swid\w+\s+)([\$\w]+)/specify_width $3 $1/;
				$parameters =~ s/(eval\w+\s+.*)(\shei\w+\s+)([\$\w]+)/specify_height $3 $1/;
				$parameters =~ s/(eval\w+\s+.*)(\sdep\w+\s+)([\$\w]+)/specify_depth $3 $1/;
				$parameters =~ s/(eval\w+\s+.*)(\sfor\w+\s+)([\$\w]+)/specify_format $3 $1/;

				$gfx_command =~ s/$mod_texture/$command$parameters/i;
			}


		  $line =~ s/$gfx_match_string/$gfx_command/;
		}

	 return ($line);
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

sub update_file
  {
	 # Given the input filename and output filename this routine parses
	 # the given file and makes any required changes.
	 my $input = shift;
	 my $output = shift;
	 my $line;
	 my $new_line;
	 my $comment_line;
	 my $uncomment_line;
	 my $date = localtime(time());

	 open(INPUT, "<$input") || die "Could not open file $input: $!\n";
	 open(OUTPUT, ">$output") || die "Could not open file $output: $!\n";

	 while (defined($line = <INPUT>))
		{
		  #Don't process comments
		  if ($line =~ m/^([^#]*)(#.*)$/)
			 {
				$comment_line = $2;
			   if ($uncomment_line = $1)
			    {
					$new_line = make_changes($uncomment_line) . $comment_line ."\n";
			    }
			  else
			    {
					$new_line = $line;
				 }
			 }
		  else
			 {
				$new_line = make_changes($line);
			 }
		  if ($new_line ne $line)
			 {
				print OUTPUT "\n# Command updated by $scriptname on $date\n" .
				  "# Old command: $line";
				print OUTPUT $new_line;
				print OUTPUT "\n";
			 }
		  else
			 {
				print OUTPUT $line;
			 }
		}
	 close (INPUT);
	 close (OUTPUT);
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

