eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;

# executes ftnchek on cm source
# optional argument specifies the source directory to use

use strict;

$0 =~ s%^.*/%%; # remove path from warning messages
$SIG{__DIE__} = sub { die "$0: ",@_ }; # die misbehaves if set in BEGIN

exists $ENV{CMISS_ROOT} or die "$0: CMISS_ROOT not defined\n";
my $cmiss_root = $ENV{CMISS_ROOT};
my $local_root = @ARGV ? $ARGV[0] : "$cmiss_root/cm";
substr($local_root,-1) eq '/' and chop $local_root;

my $cmiss_locale =
	exists $ENV{CMISS_LOCALE} ? $ENV{CMISS_LOCALE} : "DEFAULT";

my $global_source_dir = "$cmiss_root/cm/generated/i686-linux-debug";
my $source_dir = "$local_root/generated/i686-linux-debug";
#OLD IRIX files
#my $global_source_dir = "$cmiss_root/cm/generated/mips-n32-irix-debug";
#my $source_dir = "$local_root/generated/mips-n32-irix-debug";

# There are too many files to supply as arguments to ftnchek so we
# supply the source on STDIN.  This is done by forking a child to
# receive the source.

unless( my $pid = open FORTRAN, '|-' ) {
    defined $pid or die "can't fork: $!";

    # Child
    # Source is supplied on STDIN by the parent.
    # Run ftnchek on this source.

    my $command = '/hpc/cmiss/cmiss_utils/www/fortran_check/ftnchek-3.3.1/ftnchek';
    #my $command = 'ftnchek';

    my @ftnchek_options =
	( "-include=$source_dir",
	  "-include=$local_root/source",
	  "-include=$global_source_dir",
	  "-include=$cmiss_root/cm/source",
	  "-include=$cmiss_root/linear_solvers/solver",
	  "-include=$cmiss_root/gx/source",
	  '-nonovice',
	  '-quiet',
	  '-nowrap',
	  '-errors=0',
	  '-array=0',
	  '-nopure',
	  '-identifier-chars=_%',
	  '-pretty=no-continuation',
	  '-portability=real-do,tabs',
	  '-style=do-construct,program-stmt',
# 	  '-cross-ref=calls',
#  	  '-list',

# With arg-array-alias, ftnchek doesn't check the subscripts so too
# many warnings are given.
	  '-usage=no-arg-array-alias,no-ext-undefined'
	  );

#     my $output_file = "$local_root/source/checkcmiss.out";
#     open STDOUT, ">$output_file" or die "$0: can't open $output_file: $!\n";

    exec $command,@ftnchek_options;
    die "can't exec $command: $!\n";
}

#??? Should we now close STDOUT which is written to by ftnchek?

# Fork so that we can receive the list of source files
# from STDOUT of the system command find.

unless( my $pid = open FILE_LIST, '-|' ) {
    defined $pid or die "can't fork: $!";

    # Child
    # Supply a list of source files on STDOUT to the parent.

    my $command = 'find';
    exec $command, $source_dir, '-name', '*.f';
    die "can't exec $command: $!\n";
}

# Read the list of source files from FILE_LIST;

my @file_list;
while( my $filename = <FILE_LIST> )
  {
    chomp $filename;
    push @file_list, $filename;
  }

unless( close FILE_LIST ) {
    $! and die "error closing pipe from find: $!\n";
    $? & 0xFF and die "find received signal: $?\n";
    die "find failed: ", $? >> 8, "\n";
}

my $file_count = 0;
foreach my $filename ( sort @file_list ) {

    my $shortname =
      $filename =~ m{^\Q$source_dir\E/(.*)} ? $1 : $filename;

    # strip out generated fz modules which would produce many errors
    if( $shortname !~ m{^fz.*.f$} ) {

	$file_count++;

# With this simpler approach vars not ref'd and set not used are ignored
#  	print FORTRAN "      INCLUDE \"$shortname\"\n";

	# Use the #line directive to tell ftnchek the reference filename
	# and reset the line counter to one.
	# If ftnchek only takes the first filename, then
	# modify get_cpp_directive like this:
#  Index: forlex.c
#  ===================================================================
#  --- forlex.c    3 Nov 2001 00:55:53 -0000       1.1.1.1
#  +++ forlex.c    11 Mar 2002 05:02:38 -0000      1.2
#  @@ -406,18 +406,22 @@
 
#                  /* A #line directive on first line of toplevel source file
#                     gives name of real original file.  Replace our idea
#  -                  of top_filename with that. But ignore an initial
#  -                  # 1 "" since that means cpp was working with stdin,
#  -                  probably from ftnpp.  Likewise ignore # 1 "stdin"
#  +                  of top_filename with that. But a
#  +                  # 1 "" means cpp was working with stdin,
#  +                  probably from ftnpp.  Likewise # 1 "stdin"
#                     which is a variant form, e.g. from fpp.
#                  */
#  -    if( cpp_start_of_file ) {
#  -      if( next_filename[0] != '\0' &&
#  -         strcmp(next_filename,"stdin") != 0 ) {
#  +    if( next_filename[0] == '\0' ||
#  +       strcmp(next_filename,"stdin") == 0 ) {
#          top_filename = next_filename;
#          current_filename = next_filename;
#  -       cpp_start_of_file = FALSE;
#  +      cpp_inc_depth = 0;
#  +      cpp_start_of_file = TRUE;
#         }
#  +    else if( cpp_start_of_file ) {
#  +      top_filename = next_filename;
#  +      current_filename = next_filename;
#  +      cpp_start_of_file = FALSE;
#       }
#       else {
#         if( cpp_inc_depth > 0 &&

	# Tell ftnchek we are stdin so that it knows to treat
	# the next #line directive as the top file.
	print FORTRAN "# $file_count \"stdin\"\n";
	# Give ftnchek the filename and next line is line 1.
  	print FORTRAN "# 1 \"$shortname\"\n";

	# Send the source file to STDOUT for ftnchek

	open FILE, $filename or die "can't open $filename: $!\n";
	my $have_newline = 1;
	while( defined( my $line = <FILE> ) ) {
	  $have_newline = $line =~ /\n$/;
	  print FORTRAN $line;
	}
	# Make sure there is a newline before the next #line directive.
	$have_newline or print FORTRAN "\n";
    }
}

# Remind ftnchek of the input source so that it doesn't attribute summary
# errors misleadingly to the last file input.
print FORTRAN "# $file_count \"stdin\"\n";

# If there is no main program in the fortran source, then ftnchek can't
# determine which routines are never invoked.  So supply a main
# program that accesses the root subroutines.
print FORTRAN <<FORTRAN;
      program %main

      external cmdestroy,cminitialise,
     &  getnextcomfileline,getprevcomfileline,
     &  maincmloop,setcmisscommandparams

      call %dummy(cmdestroy,cminitialise,
     &  getnextcomfileline,getprevcomfileline,
     &  maincmloop,setcmisscommandparams)

      end
FORTRAN

unless( close FORTRAN ) {
    $! and die "error closing pipe to ftnchek: $!\n";
    $? & 0xFF and die "ftnchek received signal: $?\n";
    die "ftnchek failed: ", $? >> 8, "\n";
}
