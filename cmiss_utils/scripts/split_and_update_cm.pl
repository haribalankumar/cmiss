
use strict;
#For the preprocess/split_nad_update_cm_split script
$ENV{ABI} = "32";
$ENV{DEBUG} = "true";

if (! (-d "source"))
  {
	 die "This script should be run in your cm directory at the same level as your Makefile";
  }

#Keep a log of who is running the script
my $log = "$ENV{CMISS_ROOT}/cmiss_utils/scripts/split_and_update_cm.log";
if (-f $log)
  {
	 system("echo $ENV{USER} `date` $ENV{PWD} >> $log");
  }

#Check that we can get up to date presplit (do it twice so we get the conflicts as a return code)
system ("cvs update -r pre-cm-split && cvs update -r pre-cm-split")
  and die "Unable to update changes made to cm prior to splitting";

#Update to split time in preparation of split
system ("cvs update -r cm-split && cvs update -r cm-split")
  and die "Unable to update to split point";

#For each fe*.f check out the equivalent fe*/*.f directory
#and split over the top of this, removing changes as we process them.
opendir(DIR, "source") or die "can't opendir source directory: $!\n";
#Make sure that there is a split directory for each source file we are going to work on
my @source_files = grep {!system("cvs log source/$_ >/dev/null")}
  map {$_ =~ s/\.f$//; $_} grep { /^fe.*\.f$/ } readdir(DIR);
closedir DIR;
my $source_file;
for $source_file (@source_files)
  {
	 #Check out the entire directory for each fe*.f
	 system ("cd source ; cvs update -r cm-split -d $source_file")
		and die "Unable to update source/$source_file directory";
	 #Split the file
	 system ("perl $ENV{CMISS_ROOT}/cmiss_utils/scripts/split_and_update_cm_split.pl source/$source_file.f source")
		and die "Unable to split changes in $source_file.f";
	 #Find out which ones have changed
	 opendir(DIR, "source/$source_file") or die "can't opendir source/$source_file directory: $!\n";
	 my @split_files = grep { /.f$/ } readdir(DIR);
	 closedir DIR;
	 my $split_file;
	 my @changed_split_files;
	 for $split_file (@split_files)
	 {
		my $result = system "cvs diff source/$source_file/$split_file > /dev/null";
		if ($result & 0xff)
		  {
			 die "Unable to cvs diff file source/$source_file/$split_file";
		  }
		if ($result << 8)
		  {
			 push (@changed_split_files, $split_file);
		  }
		else
		  {
			 #There is nothing to keep so just remove it
			 unlink "source/$source_file/$split_file"
				or die "Unable to unlink source/$source_file/$split_file";
		  }
	 }
	 #Clean up the entries file
	 system "rm source/$source_file/CVS/Entries && echo 'D' > source/$source_file/CVS/Entries"
		and die "Unable to clean source/$source_file/CVS/Entries";
	 #Make the entries static so that it only checks out the files we want
	 system "touch source/$source_file/CVS/Entries.Static"
		and die "Unable to make source/$source_file/CVS/Entries.Static";
	 #Re check out only those that have been modified and copy the modified split files over
	 for $split_file (@changed_split_files)
	 {
		#Move the file with the changes out of the way
		rename "source/$source_file/$split_file", "source/$source_file/$split_file.keep"
		  or die "Unable to move source/$source_file/$split_file out of the way";
		#Update the file so CVS recognises it again
		system ("cd source/$source_file ; cvs update -r cm-split $split_file")
		  and die "Unable to update source/$source_file/$split_file";
		#Move the changes back
		rename "source/$source_file/$split_file.keep", "source/$source_file/$split_file"
		  or die "Unable to move source/$source_file/$split_file back";
	 }
	 #Move the old changed unsplit file out of the way otherwise we will get conflicts when we update further
	 rename "source/$source_file.f", "source/$source_file.f.presplit"
		or die "Unable to rename presplit file source/$source_file.f";
	 #Update to the final unsplit version so we do not generate errors on the next update
	 system ("cvs update -r cm-split source/$source_file.f")
		and die "Unable to update source/$source_file directory";
  }

#Update to post split
system ("cvs update -r post-cm-split && cvs update -r post-cm-split")
  and die "Unable to remove old unsplit files and update to post split point";

#Update to current state removing all sticky tags
system ("cvs update -A && cvs update")
  and die "Conflict or error after splitting";
