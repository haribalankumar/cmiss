#!/usr/bin/perl -w

use strict;

my $usage = "Usage: $0 hpc2 [example_subtree]\n";

sub dirname ($) {
    # return the directory part of a filename
    $_[0] =~ m'^(.*)/[^/]*$' ? $1 : '.';
}
my $examples_root;
BEGIN { # begin for use lib below
  $examples_root = $ENV{CMISS_EXAMPLES} or
    die "Need to set environment variable CMISS_EXAMPLES\n";
}

# for sh_quote
use lib dirname __FILE__, "$examples_root/common/cmiss_type1";
use test_cmiss;

if (@ARGV < 1)
  {
    print $usage;
    exit;
  }

my $test_output;
my $example_pattern;
my $machine_address;

if ($ARGV[0] eq "hpc2")
  {
	 $example_pattern = '(hpc_)?cmo?(64)?_(irix|n32)(_.*)?';
	 $machine_address = "hpc2";
  }
else
  {
	 die ("ERROR: Invalid machine name " . $ARGV[0] . "\n$usage");
  }

my @subdirectory_option = ();
if (@ARGV == 2)
  {
	 if (! -d "$examples_root/$ARGV[1]")
		{
		  die ("Directory $examples_root/$ARGV[1] not found.  \n$usage");
		}
	 @subdirectory_option = ('-e',$ARGV[1]);
  }

0 == system( qw(ssh -n), $machine_address,
	     '$CMISS_EXAMPLES/common/cmiss_type1/send_output.pl',
	     test_cmiss::sh_quote( @subdirectory_option,
				   test_cmiss::hostname() . ":$examples_root/",
				   '-v', $example_pattern ))
  or die "send_output.pl failed on $machine_address";

# !!! this will actually merge other examples and versions from previous
# unmerged send_output runs.

0 == system( "$examples_root/common/cmiss_type1/merge_histories.pl",
	     $machine_address )
  or die "merge_histories failed for $machine_address";
