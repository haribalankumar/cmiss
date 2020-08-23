#!/usr/bin/perl -w

use strict;

my $usage = "Usage: $0 cvs_checkedout_dir since_timestamp\n";

if (@ARGV < 2)
  {
    print $usage;
    exit;
  }

my $root_dir = $ARGV[0];
my $since_timestamp = $ARGV[1];
my $filename;
my $date;
my $author;
my $revision;
my $change_string;
# I put everything into this hash so we can sort it by time.
my $changes_log;

open (CVS_LOG, "cd $root_dir ; cvs log -d \"$since_timestamp<\" 2>&1 |");

while (<CVS_LOG>)
  {
    if (m%Working file: (.*)%)
      {
	$filename = $1;
      }
    if (m%revision (.*)%)
      {
	$revision = $1;
	$_ = <CVS_LOG>;
	if (m%date: (.*);  author: (.*);  state: Exp;%)
	  {
	    $date = $1;
	    $author = $2;
	    $change_string = "Author $author changed $root_dir/$filename to version $revision at $date\n";
	    $_ = <CVS_LOG>;
	    while (! m%-----------% && ! m%==========%)
	      {
		$change_string .= "      $_";
		$_ = <CVS_LOG>;
	      }
	    push (@{$changes_log->{$date}}, $change_string);
	  }
      }
  }

close CVS_LOG;

for $date (sort {$a cmp $b} keys %{$changes_log})
  {
    for (@{$changes_log->{$date}})
      {
	print "$_\n";
      }
  }
