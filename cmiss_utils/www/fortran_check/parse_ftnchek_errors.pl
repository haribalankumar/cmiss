eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
# Parses ftnchek output, splitting information into various files for
# conversion to html.  Unrecognized lines are sent to STDOUT.

use strict;

sub dirname ($) {
    # return the directory part of a filename
    $_[0] =~ m'^(.*)/[^/]*$' ? $1 : '.';
}

my $line = &getline();

#-----------------------------------------------------------------------------
# Subroutines

sub shortname ($)
  # return a shorted version of the filename if it is a full path name
  {
    return $_[0] =~ m{^/.*/cm/(?:source|generated/[^/]+)/(.*)$} ? $1 : $_[0];
  }

#-----------------------------------------------------------------------------

# Open the output files.

my $output_dir =
    exists $ENV{CM_PARSE_FTNCHK} ? $ENV{CM_PARSE_FTNCHK}
    : dirname __FILE__;

&open_or_die( \*GENERAL_ERROR, ">$output_dir/general_error.dat" );
print GENERAL_ERROR "2 files lines\n";
print GENERAL_ERROR "Fortran Errors\n";

&open_or_die( \*SUBPROGRAM_USAGE,
	      ">$output_dir/argument_usage_mismatch.dat" );
print SUBPROGRAM_USAGE "2 modules arguments\n";
print SUBPROGRAM_USAGE "Argument Mismatches\n";

# parse ftnchek output

if( ! grep { $_ eq $line } ('FTNCHEK Version 3 September 2000',
			    'FTNCHEK Version 3.1 May 2001',
			    'FTNCHEK Version 3.2 November 2002',
			    'FTNCHEK Version 3.3 November 2004') ) {
    print STDERR "$0: ftnchek version not recognized\n";
} else {
    $line = &getline();
}

my %used_blocks;
my %block_files;
my %module_files;

while( defined $line ) {
    if( $line =~ m/^(?:File \S+:| +\d+ (?:syntax errors? detected|warnings? issued) in file \S+)$/ ) {
	# File header: do nothing
	    $line = &getline();
    } elsif( $line =~ m/^"(\S+)", line \d+: Warning in module (\S+): (.*)$/ ) {
	# Module based warning
	my ($file,$module,$message) = ($1,$2,$3);
	if( $message eq 'Variables declared but never referenced:' ) {
	    # variables not referenced
	    $line = &getline();
	    while( defined $line and
		   $line =~
m/^"\S+", line (\d+):     (\S+) (declared(?: \(dummy argument\))?)$/
		   ) {
		my ($line_number,$variable,$message) = ($1,$2,$3);
		print GENERAL_ERROR
"$file $module variable $variable $message line $line_number but never referenced\n";
		$line = &getline();
		if( $line =~ m/^"\S+", line \d+:     \(where included\)$/ ) {
		    $line = &getline();
		}
	    }
	} elsif( $message eq 'Labels defined but not used:' ) {
	    # Labels defined but not used
	    $line = &getline();
	    while( defined $line and
		   $line =~ m/^"\S+", line (\d+):     (\S+) defined$/ ) {
		my ($line_number,$label) = ($1,$2);
		print GENERAL_ERROR
"$file $module label $label defined line $line_number but not used\n";
		$line = &getline();
		if( $line =~ m/^"\S+", line \d+:     \(where included\)$/ ) {
		    $line = &getline();
		}
	    }
	} elsif( $message eq 'Variables set but never used:' ) {
	    # Variables set but not used
	    $line = &getline();
	    while( defined $line and
		   $line =~ m/^"\S+", line (\d+):     (\S+) set$/ ) {
		my ($line_number,$variable) = ($1,$2);
		print GENERAL_ERROR
"$file $module variable $variable set line $line_number but never used\n";
		$line = &getline();
		if( $line =~ m/^"\S+", line \d+:     \(where included\)$/ ) {
		    $line = &getline();
		}
	    }
	} elsif( grep { $message eq $_ }
		 ('Variables used before set',
		  'Variables may be used before set:') ) {
	    # Used before set.
	    # For 'may be' only print a line for used; not set also.
	    $line = &getline();
	    while( defined $line and
		   $line =~
		   m/^"\S+", line (\d+):     (\S+) (used|used; never set|set|may be set)$/
		   ) {
		my ($line_number,$variable,$operation) = ($1,$2,$3);
		unless( $operation =~ m/set$/ ) {
		    if( $operation eq 'used' ) { $variable .= ' may be' }
		    print GENERAL_ERROR
			"$file $module variable $variable used line $line_number before set\n";
		}
		$line = &getline();
		if( $line =~ m/^"\S+", line \d+:     \(where included\)(; never set)?$/ ) {
		    $line = &getline();
		}
	    }
	} elsif( $message eq 'Identifiers of undeclared type' ) {
	    # Undeclared indentifiers
	    $line = &getline();
	    while( defined $line and
		   $line =~ m/^"\S+", line (\d+):     (\S+) (referenced|first occurrence)$/ ) {
		my ($line_number,$identifier,$use) = ($1,$2,$3);
		unless( grep { $identifier eq $_ } qw(%VAL %REF %LOC ACHAR) ) {
		    print GENERAL_ERROR
"$file $module identifier $identifier is of undeclared type; $use line $line_number\n";
		}
		$line = &getline();
		if( $line =~ m/^"\S+", line \d+:     \(where included\)$/ ) {
		    $line = &getline();
		}
	    }
	} else {
	    # Unrecognized line: send to STDOUT
	    print "$line\n";
	    $line = &getline();
	}

    } elsif( $line =~ m/^ *(\d+) / ) {
	# line is an error line with line number
	my $line_number_1 = $1;
	my $error_line = $line;
	my $error_found = 0;
	my $error_pointer;
	$line = &getline();
	while( defined $line and $line =~ m/^ *\d+ # \d+( +".*")? *$/ )
	  { # These preprocessor directives seem to turn up if the error was
            # on the last line of a file, but are not helpful: leave them out.
	    $line = &getline();
	  }
	while( defined $line ) {
	    my $message_regexp;
	    if( $line =~ m/^ +\^$/ ) {
		$error_pointer = $line;
		$line = &getline();
		$message_regexp = qr/^"(\S+)", line (\d+) col (?:\d+): (.*)$/;
	    } else {
		$message_regexp = qr/^"(\S+)", near line (\d+): (.*)$/;
	    }
	    unless( $line =~ $message_regexp ) { last }
	    my ($file,$line_number,$message) =
		(shortname $1, $2, $3);
	    unless( $line_number eq $line_number_1 ) {
		print STDERR
"$0: line number $line_number does not match $line_number_1\n";
	    }
	    unless( grep { $message =~ $_ } (
# %LOC is interpreted as real*2.
qr/^Warning: real\*2 expr %LOC\(.*\) promoted to intg\*4 \w+: may not give desired precision$/,
# ACHAR is interpreted as real
qr/^Error: type mismatch: real expr ACHAR\(\d+\) assigned to char\*1 \w+$/,
# ignore truncation of character expressions
qr/^Warning: char\*\d+ expr .* truncated to char\*\d+ \w+$/ )
		    ) { 
		$message .= defined $error_pointer ?
		    " <PRE> $error_line <BR> $error_pointer </PRE>":
			" <PRE> $error_line </PRE>";
		print GENERAL_ERROR "$file $message\n";
	    }
	    $error_found = 1;
	    undef $error_pointer;
	    $line = &getline();
	    # Extra newline after Zero-length string not allowed.
	    defined $line and $line eq '' and $line = &getline();
	}
	# if no error found for error line; send error line to STDOUT
	if( defined $error_pointer ) { print "$error_pointer\n" }
	unless( $error_found ) { print "$error_line\n" }

    } elsif( $line =~ m/^"(\S+)", (line \d+ col \d+: .*)$/ ) {
	# errors without the line of code supplied
	my ($file,$message) = (shortname $1, $2);
	print GENERAL_ERROR "$file $message\n";
	$line = &getline();

    } elsif( $line =~
m/^"(\S+)", line \d+: Warning: Subprogram (\S+) (.*)$/ ) {
	# Subprogram argument mismatch errors
	my ($file,$module,$warning) = ($1,$2,$3);
	my $position = '';
	if( $warning =~ m/ at position \d+:$/ ) {
	    $warning = $`; #`
	    $position = $&;
	}
	my $mismatch =
	    grep { $warning eq $_ } ('argument data type mismatch',
				     'argument arrayness mismatch',
				     'argument usage mismatch');
	my $usage_mismatch = $warning eq 'argument usage mismatch';
	$line = &getline();
	while( 1 ) {
	    my @messages = ();
	    my @callers = ();
	    my $args = 0;
	    my $may_be = 1;
	    my $uncertain_alias = 1;
	    my $start_save;
	    while( defined $line and
		   $line =~ m/^("\S+", line \d+:) (?:\(location of error\)|(\(where included\))?   (.*))$/ ) {
		my ($start,$where_included,$rest) = ($1,$2,$3);
		if( $rest ) {
		    if( $where_included ) { $start = $start_save }
		    if( $rest =~
m{^(?:Actual arg (\S+|\'.*)|Invoked|Declared) in module (\S+)(?: (.*))?$} ) {
			my ($arg,$module,$rest) = ($1,$2,$3);
			unless( grep { $module eq $_ } @callers ) {
			    push @callers, $module;
			}
			# Filter out mismatches with %VAL() or %REF()
			# in the actual arg.
			if( $mismatch ) {
			    $arg =~ m{^\%(VAL|REF)\(} or $args++;
			}
# In arg-common-alias warnings, ftnchek doesn't actually check that
# the common block variable may be used in the subprogram.  There are
# too many of this sort of error to look through them all for possible
# problems.  But: I think our code should be designed so that, if a
# common variable is passed into a routine, that routine should not
# use that common block.  If it does need that common block then the
# common block file probably needs splitting up.  So arg-common-alias
# errors should be reported.
# For now filter out arg-common-alias reports unless they state
# something about the common variable (e.g. `which may be modified').
			( $usage_mismatch
			  and $rest =~ m{^is in common block \S+$} )
			    or $uncertain_alias = 0;
		    } elsif( $rest =~
m/^(?:Dummy arg \S+|Defined) in module (\S+)(?: (.*))?$/ ) {
			my ($module2,$rest) = ($1,$2);
			if( $module2 ne $module ) {
			    print STDERR
				"$0: module $module2 does not match $module\n";
			}
			$args++;
			if( $usage_mismatch ) {
			    if( $rest =~
m/^is aliased to common var \d+: \S+ in block \S+(?: (.*))?$/ ) {
				my $rest = $1;
				if( $rest ) {
				    $uncertain_alias = 0;
				    $rest ne 'which may be modified'
					and $may_be = 0;
				}
			    } else {
# Filter out `may be modified's because there are too many.  There may
# be many warnings associated with each position.  If we find any
# warning that is not `may be' then make this report.
				$rest ne 'may be modified' and $may_be = 0;
			    }
			}
		    } else {
			# Don't recognize this.  Report anyway.
			warn "`$line' not expected";
		    }
		    $start =~ m/"(\S+)"(, line \d+:)/ or
			print STDERR "$0: regexp match failed";
		    push @messages, (shortname $1) . "$2 $rest";
		} else {
		    # record location of error
		    $start_save = $start;
		}

		$line = &getline();
	    }

	    # Filter out mismatches with %VAL() in the actual arg,
	    # may be modifieds, and vague common aliasing.
	    unless( ($mismatch && $args < 2) ||
		    ($usage_mismatch && ($may_be || $uncertain_alias)) ) {
		if( @callers ) {
		    print SUBPROGRAM_USAGE
			join('-',@callers,$module),
			join(' <BR> '," $module $warning$position",@messages),
			"\n";
		} else {
		    print GENERAL_ERROR
			"$file $module",
			join(' <BR> '," $warning$position",@messages),
			"\n";
		}
	    }

	    if( defined $line ) {
		if( $line =~ m/^"\S+", line \d+:  and( at position \d+:)?$/ ) {
		    if( defined $1 ) {
			$position = $1;
		    } else {
			$position = '';
		    }
		    $line = &getline();
		} elsif( $line eq ' etc...' ) {
		    $position = $line;
		} else {
		    last;
		}
	    } else {
		last;
	    }
	}

    } elsif( $line =~
m/^"(\S+)", line \d+: Warning: Common block (\S+) (.*)$/ ) {
	# Common block warnings
	my ($file,$block,$warning) = (shortname $1, $2, $3);
	my $position = '';
	if( $warning =~ m/ at position \d+:$/ ) {
	    $warning = $`; #`
	    $position = $&;
	}
	$line = &getline();

	# Filter out unused messages for modules.
	# What we really want to know is whether the .cmn file should
	# be included.  This is not implemented here.
	if( $warning eq 'unused in the following modules:' ) {
	    while( defined $line
		   and $line =~ 
		   m/^"(\S+)", line \d+: (?:\(location of error\)|(\(where included\))?   Unused in module (\S+))$/ ) {
		my ($file,$where_included,$module) = (shortname $1, $2, $3);
		# store the block file against the block key
		unless( defined $where_included ) {
		    ${$block_files{$block}}{$file} = 1;
		}
		# store the unused block against the module key.
		if( defined $module ) {
		    push @{$used_blocks{$module}}, $block;
		}
		$line = &getline();
	    }
	} else { while( 1 ) {
	    my @files = ($file);
	    my @messages = ();

	    while( defined $line ) {
		if( $line =~ 
m/^"(\S+)", line \d+: (?:\(location of error\)$|(\(where included\))?   )/
		    ) {
		    # seems to be part of the same error
		    my ($file,$where_included) = (shortname $1, $2);
		    # if the error is occuring across more than one file
		    # then include all filenames in the level heading.
		    unless( defined $where_included
			    or grep { $file eq $_ } @files ) {
			    push @files, $file;
			}
		} elsif( $line !~ m/^   ( [A-Z_]\w*)+$/ ) {
		    # line not part of same error (not list of elements).
		    last;
		}
		push @messages, $line;
		$line = &getline();
	    }

	    print GENERAL_ERROR
		join('-',@files),
		join(' <BR> ',
		     " Common block $block $warning$position",@messages),
		"\n";

	    if( defined $line ) {
		if( $line =~ m/^"\S+", line \d+:  and( at position \d+:)?$/ ) {
		    if( defined $1 ) {
			$position = $1;
		    } else {
			$position = '';
		    }
		    $line = &getline();
		} elsif( $line eq ' etc...' ) {
		    $position = $line;
		} else {
		    last;
		}
	    } else {
		last;
	    }
	}}

    } elsif( $line =~
	     m/^"(\S+)": ([^:]+): (?:File|Included file (\S+)) (.*)$/ ) {
	# Nonportable usage: contains tabs
	my ($file,$message) =
	    (shortname(defined $3 ? $3:$1),"$2: $4");
	print GENERAL_ERROR "$file $message\n";
	$line = &getline();

    } else {
	# Unrecognized line: send to STDOUT
	print "$line\n";
	$line = &getline();

    }

}

# close output files

close GENERAL_ERROR or die "error closing general_error.dat: $!\n";
close SUBPROGRAM_USAGE or
  die "error closing argument_usage_mismatch.dat: $!\n";

#-----------------------------------------------------------------------------
# Subroutines

sub open_or_die {
    my ($fh,$filename) = @_;

    open $fh, $filename or
	die "$0: can't open $filename: $!\n"; 
}

sub getline {
    my $line = <>;
    chomp $line if defined $line;
    return $line;
}

sub tail {
    my ($filename) = @_;
    return substr $filename, 1 + rindex $filename, '/';
}
