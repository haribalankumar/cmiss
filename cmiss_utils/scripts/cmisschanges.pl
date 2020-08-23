eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
#
# Usage : cmisschanges
#						-tree					(directory to search)				[]
#						-branch 			(auckland|oxford)						[all]
#						-leaf					(fe21.f|cmgui.c)						[all]
#
#						-sincelast		(day|month|year|all)				[day]
#						-user					(all|cheng|blackett)				[all]
#
#						-outputformat	(html|text)									[text]
#
#
use strict;
use Getopt::Long;
use Pod::Usage;
use Time::Local;
use POSIX qw(strftime);

$SIG{__WARN__} = sub { warn "$0: ",@_ };

# Get CVSROOT first to determine whether we need remote access.

my $cvsroot = $ENV{CVSROOT};

{
	my $parser = new Getopt::Long::Parser config => ['pass_through'];
	$parser->getoptions( "cvsroot=s" => \$cvsroot ) or pod2usage(2);
}

if ( ! $cvsroot )
{
	die "Must set environment variable CVSROOT or use --cvsroot option\n";
}

if( $cvsroot =~ m|^(?:ext:)?([a-z0-9_@.-]+):(/.*)|i )
	{ # remote access
		my $remote_host = $1;
		$cvsroot = $2;

		sub sh_quote (@) {
			# quotes and joins arguments for input to shell.
			join ' ', map {
				my $arg = $_;
				# use quotes if arg is zero length or if there are meta characters
				if( $arg eq '' || $arg =~ m'[^-_.:/=A-Z0-9]'i ) {
					# KAT I can't find a consistent method for quoting
					# newlines that works in sh and csh.  We don't know what
					# the (login) shell (that ssh uses) is so abort on
					# newlines.
					$arg =~ m/\n/ and
						die "don't know how to quote newlines for shell";
					# In csh, the history substitution character and newlines
					# need a backslash quote.
					# Assume the history substitution character hasn't been
					# changed through histchars.
					$arg =~ s/['!]/'\\$1'/g;	# quote quotes
					$arg = "'$arg'";			# quote result
				}
				$arg;
			} @_;
		}

		sub exec_command (@)
			{
				my (@command) = @_;
				exec @command or die "can't exec $command[0]: $!\n";
			}

		my $rsh = $ENV{CVS_RSH} || 'ssh';

		my $remote_command = $ENV{CMISSCHANGES_SERVER} ||
			'/people/cmiss/cmiss_utils/scripts/cmisschanges.pl';

		exec_command $rsh, $remote_host,
			sh_quote $remote_command, "--cvsroot=$cvsroot", @ARGV;
	}

#
# Setup the option flags
#
my $searchtree;
my $searchbranch;
my $searchuser;
my $sinceperiod = 'default';
my $after_stamp; # string
my $after_time; # in seconds
my $before_stamp; # string
my $before_time; # in seconds
my $searchleaf;
my $usehtml = 0;

sub normalize_timestamp ($)
	# Put timestamp into the form in the logfile.
	{
    my ($indate) = @_;
    my @matches =
			$indate =~ m%^(\d{4})(?:[-./](\d+)(?:[-./](\d+)(?:[. ;,](\d+)(?:[.:](\d+)(?:[.:](\d+))?)?)?)?)?$%;
    @matches or die "can't parse date $indate\n";
		# Fill in unspecified less significant components
		for ( my $i = 1; $i < 6; $i++ )
			{
				# Remember months and days start at 1; year is always specified
				defined $matches[$i] or
					$matches[$i] = (undef, 1, 1, 0, 0, 0)[$i];
			}
    return sprintf '%.4d-%.2d-%.2d,%.2d:%.2d:%.2d', @matches;
	}

sub stamp_from_time ($)
	{
    my @gmtime = gmtime $_[0];
    my $sec = $gmtime[0];
    my $min = $gmtime[1];
    my $hour = $gmtime[2];
    my $day = $gmtime[3];
    my $month = $gmtime[4] + 1;
    my $year = $gmtime[5] + 1900;
    sprintf "%.4d-%.2d-%.2d,%.2d:%.2d:%.2d",$year,$month,$day,$hour,$min,$sec;
	}

sub time_from_stamp ($)
	{
    my ($time_stamp) = @_;
    $time_stamp =~ m/(\d+)-(\d+)-(\d+),(\d+):(\d+):(\d+)/
      or die "unrecognized time stamp `$time_stamp'\n";
    my $sec = $6;
    my $min = $5;
    my $hour = $4;
    my $mday = $3;
    my $mon = $2 - 1;
    my $year = $1;
    return Time::Local::timegm($sec,$min,$hour,$mday,$mon,$year);
	}

GetOptions
	(
		"tree=s" => \$searchtree,
		"leaf=s" => \$searchleaf,
		"branch=s" => sub
			{
				my $value = $_[1];

				if ( $value eq 'auckland' )
					{
						# Swap auckland branch to a null - it is blank in the commit.log file
						$searchbranch = '';
					}
				elsif ( $value eq 'all' )
					{
						undef $searchbranch
					}
				else
					{
						$searchbranch = $value;
						grep { $searchbranch eq $_ } qw(oxford physiome release)
							or warn "unknown cvs branch $searchbranch\n";
					}
			},
		"sincelast=s" => \$sinceperiod,
		'after=s' => sub
			{
				$after_stamp = normalize_timestamp $_[1];
				$after_time = time_from_stamp $after_stamp;
			},
		'before=s' => sub
			{
				$before_stamp = normalize_timestamp $_[1];
				$before_time = time_from_stamp $before_stamp;
			},
		"user=s" => \$searchuser,
		"outputformat=s" => sub
			{
				my $format = $_[1];

				if($format eq "text")
					{
						$usehtml=0;
					}
				elsif ($format eq "html")
					{
						$usehtml=1;
					}
				else
					{
						die "Unknown output format: try text or html\n\n";
					}
			},
		"help" => sub { pod2usage 0 }
	 ) or pod2usage(2);

if(@ARGV != 0)
{
  print STDERR "Unknown options were specified ...\n  @ARGV \n\n";
	pod2usage(2);
}

#
# Must specify the tree
#
defined $searchtree or do
{
	print "Must specify the tree to search \n\n";
	pod2usage(2);
};

# Check that there are not two starting times specified.
if ( defined $after_stamp )
	{
		$sinceperiod eq 'default' or
			die "can't specify both -after and -sincelast options\n";
		undef $sinceperiod;
	}
elsif ( $sinceperiod eq 'default' )
	{
		# Apply default to $sinceperiod.
		$sinceperiod = 'day'
	}
elsif ( $sinceperiod eq 'all' )
	{
		# There is no interval for the changes to list.
		undef $sinceperiod;
	}

if ( defined $sinceperiod )
	{
		my $period = do
			{
				if( $sinceperiod eq "day")
					{
						#LKC - doesn't work when first day of month
						# change to a simple previous 24 hours
						#$sincetime=timegm(0,0,12,$mday-2,$mon-1,$year-1900);
						24 * 60 * 60;
					}
				elsif( $sinceperiod eq "week")
					{
						7 * 24 * 60 * 60;
					}
				elsif( $sinceperiod eq "month")
					{
						30 * 24 * 60 * 60;
					}
				elsif( $sinceperiod eq "year")
					{
						365 * 24 * 60 * 60;
					}
				else # unknown
					{
						print STDERR
							"Unknown Time Since format - try day, week, month, year, all\n\n";
						pod2usage(2);
					}
			};

		$after_time = $^T - $period;
		$after_stamp = stamp_from_time $after_time;
	}

$cvsroot =~ s%^:local:%%;

#
# The logfile which records all our commits
#
my $logfile = "$cvsroot/CVSROOT/logs/commit.log";

if($usehtml==0)
{
	print   "CMISS CVS Changes for ....\n";
	print   "  Tree             : $searchtree\n";
	defined $searchbranch and
		print "  Branch           : $searchbranch\n";
	defined $searchuser and
		print "  User             : $searchuser\n";
	defined $sinceperiod and
		print "  Since Last       : $sinceperiod\n";
	defined $after_stamp and
		print "  After            : ",(join ' ', split ',', $after_stamp),
			' (',(strftime "%a %b %e %T %Z %Y", localtime $after_time),")\n";
	defined $before_stamp and
		print "  Before           : ",(join ' ', split ',', $before_stamp),
			' (',(strftime "%a %b %e %T %Z %Y", localtime $before_time),")\n";
	print "\n";
}

my @allfiles;
my $filename=$logfile;
open (LOGFILE,"<$filename") || die "Can't Open $filename: $!\n";

while (<LOGFILE>)
{
	my $file;		
	my ($time_stamp,$directory,$branch,$user,@allfiles) = split(/\|/);
	chomp @allfiles;	


	# tree
  if( $searchtree ne "all" )
	{
		next unless ($directory =~ /^$searchtree/xo );
	}
	

	#since date
	if ( defined $after_stamp )
	{
		next unless $time_stamp ge $after_stamp;
	}

	
	#before date
	if ( defined $before_stamp )
	{
		next unless $before_stamp gt $time_stamp;
	}

	
	# branch
	if( defined $searchbranch )
	{
		next if $branch ne "$searchbranch";
	}

	# user
	if( defined $searchuser )
	{
		next if $user ne "$searchuser";
	}


	#
	# Outputs the modified files ...
	#
	my $found=0;
	my $version=0;
	for($file=0;$file<=$#allfiles;$file++)
	{
		#split $file into actual $file, $action and version
		($filename,$version)=split(",",$allfiles[$file]);

		if ( ! defined $searchleaf or $filename eq $searchleaf )
		{
#			print "****** $file\n";

			#split $time_stamp into components	
			my ($date,$time)= split(",",$time_stamp);

			if ($found==0)
			{
				$found=1;
				if ($usehtml==0)
				{
					print "\nAt $date $time, $user \n";	
				}
				else
				{
					print "<BR><B>At $date $time, $user</B><BR> \n<UL> \n";
				} #usehtml
			} #first file

			if ($usehtml==1)
			{
				print "  <LI>";
			}
			if ( $version =~ s/^R//x )
			{
				print "     removed  $directory/$filename -- v$version  \n";
			}
			elsif ( $version =~ s/^A//x )
			{
				print "     added    $directory/$filename  -- v$version  \n";
			}
			else
			{
				print "     modified $directory/$filename -- v$version  \n";
			}

		} #filename valid
		
		if($file==$#allfiles && $found==1) #last file
		{
			if($usehtml==1)
			{
				print "</UL>\n";
			} #usehtml
		}
	} #end for



	#
	# Outputs the comments from the rcs file ...
	#
	if ($found== 1)
	{
		my $rcsfilename; #an rcs filename which we know still exists
		
		FILE: for($file=0;$file<=$#allfiles;$file++)
		{
			#split $file into actual $file, $action and version
			($filename,$version)=split(",",$allfiles[$file]);

#	print "*** totalfiles $#allfiles $filename $version\n";

			foreach my $rcssubdir ('','Attic/') {
				my $rcsfilenametest = "$cvsroot/$directory/$rcssubdir$filename,v";
				if( -f $rcsfilenametest ) {
					$rcsfilename = $rcsfilenametest;
					last FILE;
				}
			}
		}

		unless( defined $rcsfilename ) {
			print "         *** RCS files no longer exist ***\n";
		}
		else { 
			#
			# Remove the Add (A) and Remove (R) tags 
			#
			$version =~ s/^[AR]//x;

			unless( open (RCSFILE,"<$rcsfilename") ) {
				warn "can't open $rcsfilename: $!\n";
			} else {
				while (<RCSFILE>)	
				{
					$found = grep(/^$version/,$_); #look for the current version number

#	print "** parsing log message -- found == $found version == $version\n";

					if( $found == 1 )
					{
						$_=<RCSFILE>; # the word log (?)
						$found = grep(/^log/,$_);
						if( $found==1 )
						{
							last;
						}
					}
				}

				unless( $found ) {
					warn "version $version not found in $rcsfilename\n";
				} else {
					if($usehtml==1)
					{
					  print "<FONT COLOR=\"blue\"><PRE>\n";
					}					# while printing out the comment paragraph
					while ($found==1) 
					{
						$_=<RCSFILE>; # the start of the comment
						$found = ! grep(/^\@$/,$_);
						chomp $_;
						$_=~s/^@//x;    #remove the @
						print "$_\n";
					}
					if($usehtml==1)
					{
					  print "</PRE></FONT>\n";
					}
				} # version found
				
				close RCSFILE;
			} # rcs file opened
		} # rcs file found

		if($usehtml==1)
		{
			print "<HR>";
		}
	} #if found
} #while
close LOGFILE;


#
# Documentation Stuff
#

=head1 NAME

cmisschanges.pl - shows the cvs transactions that have occured.

=head1 SYNOPSIS

=over 12

=item B<cmisschanges.pl>
B<-tree treename>
[B<-leaf leafname>]
[B<-branch branchname>]
[B<-sincelast period>|B<-after timestamp>]
[B<-before timestamp>]
[B<-user username>]
[B<-outputformat formattype>]
[B<-help>]

=back

=head1 OPTIONS AND ARGUMENTS

=over 8

=item B<-tree>

The cvs directory. eg cm/source, example/a, cmgui/source
(Required).

=item B<-leaf>

The cvs file to exmine. eg fe07.f, cmgui.c
Default is all.

=item B<-branch>

The branch of the the tree. eg. auckland, oxford.
Default is all.

=item B<-sincelast>

Restrict changes listed to those that occurred in the previous time interval
of the specified period.
Options are day, week, month (30 days), year (365 days), and all.
Default is "day" unless B<-after> is specified.

=item B<-after>

Restrict changes listed to those that occurred after the specified date UTC.
Date formats recognized include "2003-05-21 05:53:23", "2003-05-21,05:53:23",
and "2003.05.21.05.53.23".  The less significant components need not be
specified.

=item B<-before>

Restrict changes listed to those that occurred before the specified date UTC.
Date formats recognized are the same as those for B<-after>.

=item B<-user>

The user who modified the files. eg cheng, blackett.
Default is all.

=item B<-outputformat>

Outputs the results in "text" or "html" format.
Default is "text".

=item B<-help>

Outputs this screen

=cut

### Local variables:
### tab-width: 2
### perl-indent-level: 2
### End: 
