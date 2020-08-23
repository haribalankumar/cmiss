eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;

my $usage = "usage: $0 <basename> <input_file>\n";

# Builds html files from lists of errors in the format of
# parse_ftnchek_errors.pl.

use strict;

BEGIN
  {
    $SIG{__WARN__} = sub { warn "$0: ",@_ };
  }

# collect arguments
if( @ARGV != 2 ) { die $usage }

$SIG{__DIE__} = sub { die "$0: ",@_ }; # die misbehaves if set in BEGIN

my ($basename, $input_file) = @ARGV;
unless( $input_file ) {
    die "empty input file name\n";
}

# open input
open INPUT_FILE, $input_file or
    die "can't open $input_file, $!\n";

# read the header

my ($line,$num_levels,@names);

unless( defined($line = <INPUT_FILE>)
	and ($num_levels, @names) = split ' ', $line
	and chomp $line
	and $num_levels and @names == $num_levels
	and defined($line = <INPUT_FILE>)
	and chomp $line ) {
    die "error in input file format:\n  $line\n";
}

my $title = $line;

# read the entries

my %entries = ();

while( defined($line = <INPUT_FILE>) ) {

    chomp $line;
    my ($date,@levels) = split ' ',$line, $num_levels + 1;

    unless( @levels == $num_levels ) {
	die "not enough entries in input line:\n  $line\n";
    }

    my $level_entries = \%{$entries{join ':', @levels[0..@levels-2]}};
#     my $level_entries = \%entries;
#     foreach my $level ( @levels[0..@levels-2] ) {
# 	$level_entries = \%{$level_entries->{$level}};
#     }

    my $level_array = \@{$level_entries->{''}};

    if( ! @$level_array or $date gt $level_array->[0] )
      {
	$level_array->[0] = $date;
# 	print scalar @$level_array,"\n";
# 	print "$date\n";
      }

    push @$level_array, $levels[-1];
}

close INPUT_FILE;

# write out the files
my $summary_file = "summary.txt";
my $temp_file = "/tmp/$basename.summary$$";

open SUMMARY, ">$temp_file" or die "can't open >$temp_file: $!\n";

my $lookup_url = 'http://cmiss.bioeng.auckland.ac.nz/cgi-bin/cmiss-lookup.cgi';

&output_level( \%entries, 0, $basename );

sub mangle_filename (@)
  {
    my $filename = (join ',', @_) . '.html';
    $filename =~ s|/|:|g;
    $filename =~ s|%|percent|g;
    return $filename;
  }

sub output_level {
    my ($level_entries, $level_num, @levels) = @_;

    my $filename = mangle_filename @levels;
    # This works with perl v5.6.0:
#      open my $output, "> $filename" or die "can't open $filename: $!\n";
    # Perl 5.005 needs something like this.
    local *HTML;
    open HTML, "> $filename" or die "can't open $filename: $!\n";

    my $subtitle = join( ' - ', $title, @levels[1..@levels-1]);
    $subtitle =~ s/:/: /g;
    print HTML "<HTML>\n";
    print HTML "<HEAD>\n";
    print HTML "<TITLE> CMISS Fortran Check: $subtitle </TITLE>\n";
    print HTML "</HEAD>\n";
    print HTML "<BODY>\n";
    # Provide links to routine documentation
    my $heading = $subtitle;
    $heading =~ s%^(Fortran Errors - fe\d+/)(\w+)(.f)$%$1<A HREF="$lookup_url?name=$2">$2</A>$3%;
    $heading =~ s%^(Argument Mismatches - )(\w+)(-)(\w+)$%$1<A HREF="$lookup_url?name=$2">$2</A>$3<A HREF="$lookup_url?name=$4">$4</A>%;
    print HTML "<H1> $heading </H1>\n";
#     print HTML "<P>\n";
#     print HTML "The following $names[$level_num] may contain errors:\n";
#     print HTML "</P>\n";
    print HTML "<UL>\n";

    foreach my $level ( sort keys %$level_entries ) {
	if( $level eq '' )
	  {
	    my @level_array = @{$level_entries->{''}};
	    my $date = shift @level_array;
	    printf SUMMARY "%s %s %s\n", $date, $filename, $subtitle;
	    foreach my $entry ( @level_array )
	      {
		print HTML "<LI> $entry\n";
	      }
	  }
	else
	  {
	    my $spaced_level = $level;
	    $spaced_level =~ s/:/: /g;
	    print HTML "<LI> <A HREF = \"" .
	      (mangle_filename @levels, $level) .
		"\">$spaced_level</A>\n";
	    &output_level( $level_entries->{$level}, $level_num + 1, @levels, $level );
	  }
    }
    print HTML "</UL>\n";
    close HTML;
}

close SUMMARY or die "error closing $temp_file: $!\n";

# Use mv rather than rename in case moving to another file system.
0 == system 'mv', $temp_file, $summary_file
  or die "can't mv $temp_file to $summary_file: $!\n";
