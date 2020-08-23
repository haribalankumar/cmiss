#!/usr/bin/perl -w

use strict;

my $executable;
my @executable_list;
my $source_directory;
my $target_directory;
my $mirror_remote_directory;

my $uname_node = `uname -n`;
chomp $uname_node;
if ($uname_node =~ /(esu[0-9]*|hpc2)/)
  {
	 @executable_list = ("cmgui", "cmgui-debug", "cmgui64");
	 $source_directory = "$ENV{CMISS_ROOT}/cmgui/bin/mips-irix";
	 $target_directory = "$ENV{CMISS_ROOT}/bin/mips-irix";
  }
elsif ($uname_node =~ /bioeng[0-9]*/)
  {
	 @executable_list = ("cmgui", "cmgui-debug", "cmgui-gtk");
	 $source_directory = "$ENV{CMISS_ROOT}/cmgui/bin/i686-linux";
	 $target_directory = "$ENV{CMISS_ROOT}/bin/i686-linux";
  }
elsif ($uname_node =~ /hpc/)
  {
	 @executable_list = ("cmgui", "cmgui-debug", "cmgui64", "cmgui64-debug");
	 $source_directory = "$ENV{CMISS_ROOT}/cmgui/bin/rs6000-aix";
	 $target_directory = "$ENV{CMISS_ROOT}/bin/rs6000-aix";
  }
else
  {
	 die ("Unknown machine $uname_node.");
  }

print ("Updating cmgui executables on $uname_node...\n");

for $executable (@executable_list)
  {
	 if (! -d "$target_directory")
	 {
		system ("mkdir -p $target_directory");
	 }
	 if ((! -e "$target_directory/$executable") ||
		  (-M "$source_directory/$executable" < -M "$target_directory/$executable"))
		{
		  rename ("$target_directory/$executable", "$target_directory/$executable.save");
		  print ("copying $source_directory/$executable to $target_directory/$executable\n");
		  system ("cp $source_directory/$executable $target_directory/$executable");

		  if (defined $mirror_remote_directory)
			 {
				print ("copying $source_directory/$executable to $mirror_remote_directory/$executable\n");
				system ("scp $target_directory/$executable $mirror_remote_directory/$executable;");
			 }
		}
  }
