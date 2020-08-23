#! /usr/paterson/local/bin/perl
#
# Perl script to parse the cvs log and issue the appropriate cvs 
# commands (with sensible log messages) back to the Auckland CVS repository.
#
# Usage:
#   updateaklcvs.pl
#
# Created:
#   Chris Bradley, 9 July 2002.
#
# Updates:
#

use Time::Local;

#$debug = 0;
$debug = 1;

$cvslibrary = $ENV{CMISS_ROOT}."/cvs/cm/source/";
$cvslogfile = $ENV{CMISS_ROOT}."/cvs/CVSROOT/logs/commit.log";
$crondirectory = $ENV{CMISS_ROOT}."/cmiss_utils/cronjobs";
$cvsexe = "/usr/paterson/local/bin/cvs";
$aklcvsexe = "/usr/paterson/local/bin/scvs";
$aklcvsdirectory = $ENV{CMISS_ROOT}."/aklcmiss/cm/source";
$aklcvslocalportnum = "2400";
$aklcvsservername = "localhost";
$cvsroot = "/usr/paterson/cmiss/cvs";
$aklcvsroot = "/product/cmiss/cvs";

%emailhash = ( 
	       "cpb"=>"bradley", 
	       "mn"=>"nash",
	       "nps"=>"nsmith",
	       "gry"=>"garny",
	       "dpn"=>"nickerso",
	       "pjh"=>"hunter",
	       "ajp2"=>"pullan"
	       );

$nummissed = 0;
#Process any cvs commands that didn't succeed last time.
if( -e "$crondirectory/missed_updateaklcvs.dat" ) {
    print "Found missed update commands. Processing missed commands...\n";
    $lasttime = 108000;
    rename("$crondirectory/missed_updateaklcvs.dat","$crondirectory/missed_updateaklcvs_work.dat");
    open(MISSED,">$crondirectory/missed_updateaklcvs.dat");
    open(OLDMISSED,"$crondirectory/missed_updateaklcvs_work.dat");
    while( $logline = <OLDMISSED> ) {
	chomp($logline);
	&updateaklcvs($logline);
    }
    close(MISSED);
    if( $nummissed == 0 ) {
	unlink("$crondirectory/missed_updateaklcvs.dat");
    }
    close(OLDMISSED);
    unlink("$crondirectory/missed_updateaklcvs_work.dat");
}

#Determine the last time the script was run and process since that time. If
#we can't determine the last time use the current time less one day.
if( -e "$crondirectory/last_updateaklcvs.dat" ) {
    printf "IN HERE\n";
    open(LASTTIME,"$crondirectory/last_updateaklcvs.dat");
    if( $line = <LASTTIME> ) {
	chomp($line);
	$lasttime = $line;
    } else {
	$lasttime = time - 86400;
    }
    close(LASTTIME);
} else {
    $lasttime = time - 86400;
}

#$lasttime = 1;

if( $lasttime < 108000 ) { 
    $lasttime = 108000; 
}

#Open up the missed files
if( -e "$crondirectory/missed_updateaklcvs.dat" ) {
    open(MISSED,">>$crondirectory/missed_updateaklcvs.dat");
} else {
    open(MISSED,">$crondirectory/missed_updateaklcvs.dat");
}

#Process the cvs log file
print "Processing cvs log...\n";
open(LOG,$cvslogfile);
while( $logline = <LOG> ) {
    chomp($logline);
    &updateaklcvs($logline);
}
close(LOG);

close(MISSED);
#Delete the missed file if necessary
if( $nummissed == 0 ) {
    unlink("$crondirectory/missed_updateaklcvs.dat");
} else {
    print "\nERROR: cvs commands were missed.\n";
}

#Set the last time the script was run
$lasttime = timestamp2time( time2timestamp ( time ) );
system("echo $lasttime > $crondirectory/last_updateaklcvs.dat");

sub updateaklcvs { #Takes a cvs log line and puts that change into Auckland's cvs repository
    my $cvslogline = $_[0];
    my @allfiles;
    my ($timestamp,$directory,$branch,$user,@allfiles) = split(/\|/,$cvslogline);
    chomp @allfiles;

    if( $debug ) {
	printf "\n*** updateaklcvs\n";
	printf "*** cvslogline = $cvslogline\n";
	printf "*** timestamp = $timestamp\n";
	printf "*** directory = $directory\n";
	printf "*** branch = $branch\n";
	printf "*** user = $user\n";
	printf "*** allfiles = @allfiles\n";
    }

    if( $directory eq "cm/source" ) {

	$changetime = timestamp2time($timestamp);

	if( $debug ) {
	    printf "\n*** changetime = $changetime\n";
	    printf "*** lasttime = $lasttime\n";
	}

	if( timestamp2time($timestamp) > $lasttime ) {

	    my $filename = "";
	    my $version = "";

	    if( $debug ) {
		printf "\n*** #allfiles = $#allfiles\n";
		}

	    for($file=0;$file<=$#allfiles;$file++) {
		#split $file into actual $file, $action and version
		($filename,$version)=split(",",$allfiles[$file]);
		
		if( $debug ) {
		    printf "\n*** allfiles = $allfiles[$file]\n";
		    printf "*** filename = $filename\n";
		    printf "*** version = $version\n";
		}

		$rcsfile = $cvslibrary.$filename.",v";

		if ( $version =~ s/^R//x ) {
		    #Removed a file from the cvs repository

		    if( $debug ) {
			printf "\n*** File was removed from the cvs repository\n";
		    }

		    $rcslogcommand="rlog -N -r$version $rcsfile |";

		    if( $debug ) {
			printf "\nrcslogcommand = $rcslogcommand\n";
		    }

		    if( open(RCSLOG,$rcslogcommand ) ) {
			&getcomments;
			close(RCSLOG);

			if( $revisionnumber > 0 ) {
			    if( open(COMMENTFILE,">/tmp/cvscomments.txt" ) ) {
				&writecomments;
				close(COMMENTFILE);
				$checkoutcommand="co -f -r$revisionnumber $rcsfile";
				if( $debug ) {
				    print "\ncheckoutcommand = $checkoutcommand\n";
				}

				if( ! system("cd $aklcvsdirectory; $checkoutcommand") ) {
				    if( $filename =~ "cmiss_archive" ) {
					if ( system("cd $aklcvsdirectory; mv -f $filename archive/.") ){
					    error("Failed to move $filename to archive directory");
					    return;
					}
					$filename="archive/".$filename;
				    }
				    if( exists $emailhash{$user} ) {
					$cvscommand="$aklcvsexe -d :pserver:$emailhash{$user}"."@"."$aklcvsservername:$aklcvslocalportnum$aklcvsroot remove $filename";
					if( $debug ) {
					    print "\ncvscommand = $cvscommand\n";
					}
					if( system("cd $aklcvsdirectory; $cvscommand") ) {
					    error("Failed to complete cvs command: $cvscommand");
					} else {
					    $cvscommand="$aklcvsexe -d :pserver:$emailhash{$user}"."@"."$aklcvsservername:$aklcvslocalportnum$aklcvsroot commit -F /tmp/cvscomments.txt -roxford $filename";
					    if( $debug ) {
						print "\ncvscommand = $cvscommand\n";
					    }
					    if( system("cd $aklcvsdirectory; $cvscommand") ) {
						error("Failed to complete cvs command: $cvscommand");
					    }
 					}
				    } else {
					error("$who does not exist in the email hash");
				    }
				} else {
				    error("Could not check out revision $revisionnumber of element $element");
				}
				close(COMMENTFILE);
				#unlink(COMMENTFILE);
			    } else {
				error("Could not open comment file");
			    }
			} else {
			    error("Could not obtain revision number");
			}
		    } else {
			error("Could not obtain rcs log");
		    }

		} elsif ( $version =~ s/^A//x ) {
		    #Added a file into the cvs repository

		    if( $debug ) {
			printf "\n*** File was added into the cvs repository\n";
		    }

		    $rcslogcommand="rlog -N -r$version $rcsfile |";

		    if( $debug ) {
			printf "\nrcslogcommand = $rcslogcommand\n";
		    }

		    if( open(RCSLOG,$rcslogcommand ) ) {
			&getcomments;
			close(RCSLOG);

			if( $revisionnumber > 0 ) {
			    if( open(COMMENTFILE,">/tmp/cvscomments.txt" ) ) {
				&writecomments;
				close(COMMENTFILE);
				$checkoutcommand="co -f -r$revisionnumber $rcsfile";
				if( $debug ) {
				    print "\ncheckoutcommand = $checkoutcommand\n";
				}

				if( ! system("cd $aklcvsdirectory; $checkoutcommand") ) {
				    if( $filename =~ "cmiss_archive" ) {
					if ( system("cd $aklcvsdirectory; mv -f $filename archive/.") ){
					    error("Failed to move $filename to archive directory");
					    return;
					}
					$filename="archive/".$filename;
				    }
				    if( exists $emailhash{$user} ) {
					$cvscommand="$aklcvsexe -d :pserver:$emailhash{$user}"."@"."$aklcvsservername:$aklcvslocalportnum$aklcvsroot add $filename";
					if( $debug ) {
					    print "\ncvscommand = $cvscommand\n";
					}
					if( system("cd $aklcvsdirectory; $cvscommand") ) {
					    error("Failed to complete cvs command: $cvscommand");
					} else {
					    $cvscommand="$aklcvsexe -d :pserver:$emailhash{$user}"."@"."$aklcvsservername:$aklcvslocalportnum$aklcvsroot commit -F /tmp/cvscomments.txt -roxford $filename";
					    if( $debug ) {
						print "\ncvscommand = $cvscommand\n";
					    }
					    if( system("cd $aklcvsdirectory; $cvscommand") ) {
						error("Failed to complete cvs command: $cvscommand");
					    }
 					}
				    } else {
					error("$who does not exist in the email hash");
				    }
				} else {
				    error("Could not check out revision $revisionnumber of element $element");
				}
				close(COMMENTFILE);
				#unlink(COMMENTFILE);
			    } else {
				error("Could not open comment file");
			    }
			} else {
			    error("Could not obtain revision number");
			}
		    } else {
			error("Could not obtain rcs log");
		    }


		} else {
		    #Modified a file in the cvs repository

		    if( $debug ) {
			printf "\n*** File was modified in the cvs repository\n";
		    }

		    $rcslogcommand="rlog -N -r$version $rcsfile |";

		    if( $debug ) {
			printf "\nrcslogcommand = $rcslogcommand\n";
		    }

		    if( open(RCSLOG,$rcslogcommand ) ) {
			&getcomments;
			close(RCSLOG);

			if( $revisionnumber > 0 ) {
			    if( open(COMMENTFILE,">/tmp/cvscomments.txt" ) ) {
				&writecomments;
				close(COMMENTFILE);
				$checkoutcommand="co -f -r$revisionnumber $rcsfile";
				if( $debug ) {
				    print "\ncheckoutcommand = $checkoutcommand\n";
				}
				if( ! system("cd $aklcvsdirectory; $checkoutcommand") ) {
				    if( $filename =~ "cmiss_archive" ) {
					if ( system("cd $aklcvsdirectory; mv -f $filename archive/.") ){
					    error("Failed to move $filename to archive directory");
					    return;
					}
					$filename="archive/".$filename;
				    }
				    if( exists $emailhash{$user} ) {
					$cvscommand="$aklcvsexe -d :pserver:$emailhash{$user}"."@"."$aklcvsservername:$aklcvslocalportnum$aklcvsroot commit -F /tmp/cvscomments.txt -roxford $filename";
					if( $debug ) {
					    print "\ncvscommand = $cvscommand\n";
					}
					if( system("cd $aklcvsdirectory; $cvscommand") ) {
					    error("Failed to complete cvs command: $cvscommand");
					}
				    } else {
					error("$who does not exist in the email hash");
				    }
				} else {
				    error("Could not check out revision $revisionnumber of element $element");
				}
				close(COMMENTFILE);
				#unlink(COMMENTFILE);
			    } else {
				error("Could not open comment file");
			    }
			} else {
			    error("Could not obtain revision number");
			}
		    } else {
			error("Could not obtain rcs log");
		    }
		}


	    }
	}
    }

}
sub getcomments { #Obtains comments from a RCS log file
    $foundrevision = 0;
    $revisionnumber = -1;
    $numcomments = -2;
    $processingcomment = 0;
    @commentlines = ();
    while ( $rcsline = <RCSLOG> ) { 
	chomp($rcsline);
	if( $debug ) {
	    print "rcsline = $rcsline\n";
	}
	if( $foundrevision ) {
	    @revisionwords = split(/\s+/,$rcsline);
	    if( $debug ) {
		print "revision words:\n";
		foreach $revisionword ( @revisionwords ) {
		    print "   $revisionword\n"
		    }
	    }
	    $revisionnumber = $revisionwords[1];
	    $foundrevision = 0;
	    $processingcomment = 1;
	}
	if( $rcsline =~ /----/ ) {
	    $foundrevision = 1;
	} elsif( $rcsline =~ /====/ ) {
	    $processingcomment = 0;
	}
	if( $processingcomment ) {
	    $numcomments++;
	    if( $numcomments > 0 ) {
		push(@commentlines,$rcsline);
	    }
	}
    }
    if( $debug ) {
	print "\nrevisionwords = @revisionwords\n";
	print "revisionnumber = $revisionnumber\n";
	print "numcomments = $numcomments\n";
	print "commentlines = @commentlines\n";
    }
}

sub writecomments { #Writes the comments lines to a comments file
    foreach $commentline ( @commentlines ) {
	print COMMENTFILE "$commentline\n";
    }
}

sub time2timestamp { #Converts a time to a time stamp i.e., into a date and time
    my @mygmtime = gmtime(@_[0]);
    my $sec = $mygmtime[0];
    my $min = $mygmtime[1];
    my $hour = $mygmtime[2];
    my $day = $mygmtime[3];
    my $month = $mygmtime[4] + 1;
    my $year = $mygmtime[5] + 1900;
    $equivtimestamp=sprintf("%.4d-%.2d-%.2d,%.2d:%.2d:%.2d",$year,$month,$day,$hour,$min,$sec);

    if( $debug ) {
	printf "*** time2timestamp\n";
	printf "*** time = @_[0]\n";
	printf "*** sec = $sec\n";
	printf "*** min = $min\n";
	printf "*** hour = $hour\n";
	printf "*** day = $day\n";
	printf "*** month = $month\n";
	printf "*** year = $year\n";
	printf "*** equivtimestamp = $equivtimestamp\n";
    }

    return $equivtimestamp;
}

sub timestamp2time { #Converts a time to a time stamp i.e., into a date and time
    my ($datepart,$timepart) = split(",",@_[0]);
    my ($year,$month,$day) = split("-",$datepart);
    my ($hour,$minute,$second) = split(":",$timepart);

    if( $debug ) {
	printf "*** timestamp2time\n";
	printf "*** timestamp = @_[0]\n";
	printf "*** datepart = $datepart\n";
	printf "*** timepart = $timepart\n";
	printf "*** year = $year\n";
	printf "*** month = $month\n";
	printf "*** day = $day\n";
	printf "*** hour = $hour\n";
	printf "*** minute = $minute\n";
	printf "*** second = $second\n";
	printf "*** equivtime = $equivtime\n";
    }

    $month = $month - 1;
    $year = $year - 1900;
    $equivtime = timegm($second,$minute,$hour,$day,$month,$year);

    if( $debug ) {
	printf "*** timestamp2time\n";
	printf "*** timestamp = @_[0]\n";
	printf "*** datepart = $datepart\n";
	printf "*** timepart = $timepart\n";
	printf "*** year = $year\n";
	printf "*** month = $month\n";
	printf "*** day = $day\n";
	printf "*** hour = $hour\n";
	printf "*** minute = $minute\n";
	printf "*** second = $second\n";
	printf "*** equivtime = $equivtime\n";
    }

    return $equivtime;
}

sub error { #Error condition
    my $errorstr=shift;
    print "ERROR: ".$errorstr."\n";
    print MISSED "$logline\n";
    $nummissed++;
}


