eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

my $usage = "usage: $0 <basename>\n";

# Reads a list of ftnchek errors from <basename>.dat, compares them to a
# record <basename>.date of dates of occurence of previous errors and updates
# the record.

use strict;

BEGIN
  {
    $SIG{__WARN__} = sub { warn "$0: ",@_ };
  }

@ARGV == 1 or die $usage;

$SIG{__DIE__} = sub { die "$0: ",@_ }; # die misbehaves if set in BEGIN

my($basename) = @ARGV;

sub dirname ($)
  {
    # return the directory part of a filename
    $_[0] =~ m'^(.*)/[^/]*$' ? $1 : '.';
  }

sub generalize_error ($)
  # Remove line numbers from error messages and ignore indentation changes.
  {
    my $error = $_[0];

    $error =~ s/ line \d+//g;
    $error =~ s/  +/ /g;
    $error =~ s/ <PRE> \d+ / /;
    return $error;
  }

sub time_stamp ($)
  {
    my @gmtime = gmtime $_[0];
#      my $sec = $gmtime[0];
    my $min = $gmtime[1];
    my $hour = $gmtime[2];
    my $mday = $gmtime[3];
    my $mon = $gmtime[4] + 1;
    my $year = $gmtime[5] + 1900;
    sprintf "%.4d-%.2d-%.2d,%.2d:%.2d",$year,$mon,$mday,$hour,$min;
  }

my $output_dir =
    exists $ENV{CM_PARSE_FTNCHK} ? $ENV{CM_PARSE_FTNCHK}
    : dirname __FILE__;

my $record_file = "$output_dir/$basename.date";
my $errors_file = "$output_dir/$basename.dat";
my $temp_file = "/tmp/$basename.date$$";

my %date;
# If there is an old record file then make a hash of dates keyed by the errors.
if( -e $record_file )
  {
    open RECORD, $record_file or die "can't open $record_file: $!\n";

    # Skip the header.  We could check that it matches $errors_file.
    defined <RECORD> and defined <RECORD>
      or die "error in format of $record_file\n";

    while( defined( my $line = <RECORD> ) )
      {
	chomp $line;
	my ($date,$error) = split / /, $line, 2;
	$date{generalize_error $error} = $date;
      }

    close RECORD;
  }

open ERRORS, $errors_file or die "can't open $errors_file: $!\n";
stat ERRORS or die "can't stat $errors_file: $!\n";
my $time_stamp = time_stamp( (stat _)[9] );

open RECORD, ">$temp_file" or die "can't open >$temp_file: $!\n";

# Copy the header to the new record file.
for( my $i = 2; $i > 0; $i-- )
  {
    my $line = <ERRORS>;
    defined $line or die "error in format of $errors_file\n";

    print RECORD $line;
  }

while( defined( my $error = <ERRORS> ) )
  {
    chomp $error;
    printf RECORD "%s %s\n",
      $date{generalize_error $error} || $time_stamp, $error;
  }

close ERRORS;
close RECORD or die "error closing $temp_file: $!\n";

# Use mv rather than rename in case moving to another file system.
0 == system 'mv', $temp_file, $record_file
  or die "can't mv $temp_file to $record_file: $!\n";
