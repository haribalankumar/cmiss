eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

# Generates from a CMISS fortran file, a file of fortran source in
# <output-directory> ready for compilation.  If the extension is .f, then a
# dependency file <input-file-base>.d is produced in <output-directory>
# including all the object dependencies.

# Created: Karl Tomlinson 2/10/00

use strict;

my $usage = "usage: $0 <input-file> <output-directory>\n";

# remove path from warning and error messages
$SIG{__WARN__} = sub { $0 =~ m'([^/]*)$'; warn "$1: warning: ",@_ };
$SIG{__DIE__} = sub { $0 =~ m'([^/]*)$'; die "$1: ",@_ };

# collect arguments
if( @ARGV != 2 ) { die $usage }

my ($input_file, $output_dir) = @ARGV;
unless( $input_file ) {
    die "empty input file name";
}

if( substr( $output_dir, -1 ) ne '/' ) { $output_dir .= '/' }

# environment variables exported by Makefile
unless( exists $ENV{ABI} ) {
    die "environment variable ABI not set\n";
}
my $abi = $ENV{"ABI"};
unless( exists $ENV{DEBUG} ) {
    die "environment variable DEBUG not set\n";
}
my $opt = $ENV{DEBUG} eq 'false';
# unless( exists $ENV{FC} ) {
#     die "environment variable FC not set\n";
# }
my $fc = $ENV{FC};
my $sysname = do {
    if( exists $ENV{SYSNAME} ) {
	$ENV{SYSNAME};
    } else {
	my $output = `uname`; $? and die; chomp $output;
	$output;
    }
};
my $f90 =
  defined $fc ? $fc =~ /^(pathf90$|xlf9)/
  : $sysname eq 'AIX' ; # for backward compatibility (FC may not be set).

# make local output directory
my ($file_base,$file_ext) = (&fileparse($input_file))[1,2];
my $output_file = $output_dir.$file_base.$file_ext;

# Derive read permissions of output file from input file.
# Executable permission should not be required.
# Remove write permission to resist accidental modification of generated files.
# We should really make the groups the same or remove group privileges
# but cvs doesn't do this so there doesn't seem much point.

stat $input_file or die "can't stat $input_file: $!\n";
# my $group = (stat _)[5];
my $permissions = 00444 & (stat _)[2];

# open input
open INPUT_FILE, $input_file or
    die "can't open $input_file: $!\n";

# open output file
if( -f $output_file )
  {
    chmod $permissions | 00200, $output_file or
      die "can't chmod $permissions $output_file: $!\n";
  }

open OUTPUT_FILE, ">$output_file" or
  die "can't open $output_file: $!\n";
#     chown -1, $group, $output_file or
# 	die "can't chgrp $group $output_file, $!\n";

# write the output file and collect any dependencies

my @include_files;

while( defined(my $line = <INPUT_FILE>) )
  {
    # Search for dependencies and modify code
    &preprocess_line( \$line, \@include_files );
    print OUTPUT_FILE $line;
  }


close INPUT_FILE;

close OUTPUT_FILE or die "error closing $output_file: $!\n";

chmod $permissions, $output_file or
	die "can't chmod $permissions $output_file: $!\n";

# if( $sysname eq 'AIX' ) {
#   # checking for null bytes in output files
#   open NEW_FILE, "$output_file" or
#     die "can't open $output_file: $!\n";
#   grep {/\0/} <NEW_FILE> and die "$output_file now contains null bytes";
#   close NEW_FILE;
# }

# The target in the Makefile should be made last so that if there
# is an error with or interrupt during anything else the Makefile
# will try to remake the target.

if( $file_ext eq '.f' )
  {

    # need to write dependencies

    # write dependency file for module

    my $depend_file = $output_dir.$file_base.'.d';
    open DEPEND_FILE, ">$depend_file" or
	die "can't open $depend_file: $!\n";

    # unit .o dependencies.

    print DEPEND_FILE join( " ",
			    "$output_dir$file_base.o:",
			    $output_file,
			    @include_files ), "\n";

    close DEPEND_FILE or die "error closing $depend_file: $!\n";

  }
elsif( @include_files )
  {
    die "dependencies @include_files of $input_file were not expected.\n";
  }

#-----------------------------------------------------------------------------
# subroutines

sub preprocess_line {
# preprocess and collect any dependencies
    my ($line,$include_files) = @_;
    # INTEGER*4 declarations are pointers.  Convert to INTEGER*8 if 64-bit
    if( $$line =~ /^ {6}POINTER\b/ ) {
      if ($abi eq '64'){
        $$line =~ s/POINTER/INTEGER*8/;
      }
      else
      {
        $$line =~ s/POINTER/INTEGER*4/;
      }
    }
    elsif( $abi eq '64' && $$line =~ /^ {6}INTEGER\*4\b/ ) {
	$$line =~ s/4/8/;
    } elsif( $$line =~
	     /^(?:[Cc]\$|  ) {4}INCLUDE +'([^']+)'/i ) {
	my $file = $1;
	push @$include_files, $file;
    } elsif( $opt and
	     $$line =~ /^([\d ]{6})CALL (ENTERS|EXITS)\([^\(]+\) *(!|\n)/ ) {
	# strip calls to ENTERS and EXITS from optimized version
	if( $1 eq '      ' ) {
	    $$line = "\n";
	} else {
	    $$line = "${1}CONTINUE\n";
	}
    } elsif( $$line =~ /^ {5}[^ ] +'cmiss\$examples:'/ ) {
	if( exists $ENV{CMISS_EXAMPLES} ) {
	    $$line =~ s%cmiss\$examples:%$ENV{CMISS_EXAMPLES}%;
	} else {
	    warn "CMISS_EXAMPLES not defined\n";
	    $$line =~ s%cmiss\$examples:%.%;
	}
    #F90 uses LOC() rather than %LOC(), change for AIX machine
    } elsif( $f90 && $$line =~ m'\%LOC\(' ) {
	$$line =~ s/\%LOC/LOC/;
#       } elsif( $fc eq 'pathf90' && $$line =~ m'\%REF\(' ) {
# 	# handles only one nested set of brackets
# 	$$line =~ s/\%REF\(([^(]*(\([^(]*\))*)\)/$1/g;
    }
}

# sub update_unit {
#     my ($prog_unit,$filename) = @_;

#     my $tmpfile = "$filename.tmp";

#     # write generated fortran file
#     open OUTPUT_FILE, ">$tmpfile" or
# 	die "can't open >$tmpfile, $!\n";
#     print OUTPUT_FILE @$prog_unit;
#     close OUTPUT_FILE or die "error closing $tmpfile, $!";

#     rename $tmpfile, $filename or die "can't rename $tmpfile, $filename, $!\n";

# #     chown -1, $group, $filename or
# # 	die "can't chgrp $group $filename, $!\n";
#     chmod $permissions, $filename or
# 	die "can't chmod $permissions $filename, $!\n";
# }

sub fileparse {
    my ($filename) = @_;
    my $lastslash = rindex( $filename, '/' ) + 1;
    my ($dirname,$tail);
    if( $lastslash ) {
	$dirname = substr( $filename, 0, $lastslash );
	$tail = substr( $filename, $lastslash );
    }
    else {
	$dirname = '';
	$tail = $filename;
    }
    my $lastdot = rindex( $tail, '.' );
    my ($basename,$suffix);
    if( $lastdot != -1 ) {
	$basename = substr( $tail, 0, $lastdot );
	$suffix = substr( $tail, $lastdot );
    }
    else {
	$basename = $tail;
	$suffix = '';
    }
    return $dirname, $basename, $suffix;
}
