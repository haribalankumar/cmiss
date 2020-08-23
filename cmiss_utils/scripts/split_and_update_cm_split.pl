eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

# Generates from a cmiss fortran file, a directory
# <input-file-base>.split in <output-directory> of fortran source
# files ready for compilation.  If some of the temp files already
# exist and there is no need to change them, no change is made.  A
# dependency file <input-file-base>.d is produced in
# <output-directory> including all the object dependencies.  An object
# list file <input-file-base>.olist is produced containing a list of
# all the object file names.

# If the file extension is not .f or .F then the file is simply
# preprocessed with the output file placed in <output-directory>.

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
my $sysname = do {
    if( exists $ENV{SYSNAME} ) {
	$ENV{SYSNAME};
    } else {
	my $output = `uname`; $? and die; chomp $output;
	$output;
    }
};
my $curdir = do {
    if( exists $ENV{CURDIR} ) {
	$ENV{CURDIR};
    } else {
	my $output = `pwd`; $? and die; chomp $output;
	$output;
    }
} . '/';
my $f90 = $sysname eq 'AIX';

# make local output directory
my ($file_base,$file_ext) = (&fileparse($input_file))[1,2];
my $output_file = $output_dir.$file_base.$file_ext;

# Derive read permissions of output file from input file.
# Executable permission should not be required.
# Remove write permission to resist accidental modification of generated files.
# We should really make the groups the same or remove group privileges
# but cvs doesn't do this so there doesn't seem much point.

stat $input_file or die "can't stat $input_file, $!\n";
# my $group = (stat _)[5];
my $permissions = 00644 & (stat _)[2];

# open input
open INPUT_FILE, $input_file or
    die "can't open $input_file, $!\n";

#-----------------------------------------------------------------------------

if( $file_ext eq '.f' || $file_ext eq '.F' ) {

    my $default_unit_name = '-unknown-';

    my $split_dir = "$file_base/";
    my $local_dir = $output_dir.$split_dir;
    unless( -d $local_dir ) {
    # let the umask and parent directory determine the mode
	mkdir $local_dir, 0777 or
	    die "can't make directory $local_dir, $!\n";
    }

    # make files in the local output directory and get dependencies

    my %dependencies;

    # Set up default unit
    my $unit_name = $default_unit_name;
    my @prog_unit = ();
    my $no_code = 1; # don't create a file unless code is found
    my $output_file = $local_dir.$unit_name.'.f';
    my $same = open( OLD_FILE, $output_file );

    print "Updated: ";

    while( defined(my $line = <INPUT_FILE>) ) {
	# Check for start of a unit
	if( $line =~ /^ {6}(?:SUBROUTINE|BLOCK ?DATA|PROGRAM|(?:(?:INTEGER|REAL|LOGICAL|COMPLEX)(?:\*\d*)?|CHARACTER(?:\*\([*0-9]*\))?) +FUNCTION) +(\w*)/i ) {
	    $unit_name = $1;
	    # Clean up last unit
	    if( $same ) {
		$same = ! defined <OLD_FILE>;
		close OLD_FILE;
	    } 
	    unless( $same || $no_code ) {
		&update_unit( \@prog_unit, $output_file );
	    }
	    # Set up this unit
	    @prog_unit = ();
	    @{$dependencies{$unit_name}} = ();
	    $no_code = 0; # code found
	    $output_file = $local_dir.$unit_name.'.f';
	    $same = open( OLD_FILE, "$output_file" );
	}
	my $dependency;
	# Search for dependencies and modify code
	&preprocess_line( \$line, \$dependency );
	if( $dependency ) { # an INCLUDE statement was found
	    $no_code = 0;
	    push @{$dependencies{$unit_name}}, $dependency;
	}
	# Save current line
	push @prog_unit, $line;
	# Compare unit with last version
	if ( $same && $line ne <OLD_FILE> ) {
	    close( OLD_FILE );
	    $same = 0;
	}
	# Test of our parsing / file validity
	if( $no_code && $line !~ /^([Cc*]| *($|!))/ ) {
	    warn "code found before start of unit in $input_file:\n $line";
	}
    }
    # Clean up last unit
    if( $same ) {
	$same = ! defined <OLD_FILE>;
	close OLD_FILE;
    } 
    unless( $same || $no_code ) {
	&update_unit( \@prog_unit, $output_file );
    }
    print "\n";

    # tidy up local output directory

    opendir( LOCAL_DIR, $local_dir ) or
	die "can't open $local_dir, $!\n";
    while( defined( my $file = readdir( LOCAL_DIR ) ) ) {
	if( $file =~ /\.[fo]$/ && ! exists $dependencies{$`} ) { #`}){
	    my $filename = "$local_dir$file";
	    print "rm $filename\n";
	    unlink $filename;
	}
    }
    closedir( LOCAL_DIR );

}

close INPUT_FILE;

#-----------------------------------------------------------------------------
# subroutines

sub preprocess_line {
# preprocess and collect any dependencies
    my ($line,$dependency) = @_;
	 if( $$line =~
	     /^(?:[Cc]\$|  ) {4}INCLUDE +'((?:cmiss|gx)\$[^:]*:)?(.*)'/i ) {
	my ($tag,$file) = ($1,$2);
	if( $tag ) { $$line =~ s/'(?:cmiss|gx)\$[^:]*:.*'/'$file'/i }
	$$dependency = $file;
    }
}

sub update_unit {
    my ($prog_unit,$filename) = @_;

    my $tmpfile = "$filename.tmp";

    # write generated fortran file
    open OUTPUT_FILE, ">$tmpfile" or
	die "can't open >$tmpfile, $!\n";
    print OUTPUT_FILE @$prog_unit;
    close OUTPUT_FILE or die "error closing $tmpfile, $!";

    if( $sysname eq 'AIX' ) {
	# checking for null bytes in output files
	open NEW_FILE, "$tmpfile" or
	    die "can't open >$tmpfile, $!\n";
	grep {/\0/} <NEW_FILE> and die "$tmpfile now contains null bytes";
	close NEW_FILE;
    }

    rename $tmpfile, $filename or die "can't rename $tmpfile, $filename, $!\n";

#     chown -1, $group, $filename or
# 	die "can't chgrp $group $filename, $!\n";
    chmod $permissions, $filename or
	die "can't chmod $permissions $filename, $!\n";
    # progress message
    print ' ',(&fileparse($filename))[1];
}

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
