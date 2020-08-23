package CmUtils::IpFileFix;
require 5.006;

use strict;
use warnings;
use Carp;
use File::Compare;
use File::Find;
use Getopt::Std;

=head1 CmUtils::IpFileFix

The F<CmUtils::IpFileFix> package provides a simple mechanism for writing
scripts to update CMISS ipfiles.

=cut

BEGIN {
  use Exporter;  
  our ($VERSION, @ISA, @EXPORT);

  $VERSION  = 0.2;
  @ISA      = qw(Exporter);
  @EXPORT   = qw(&fixfiles);
}

=head2 VERSION

0.3 (18 June 2002)

=head2 CHANGES

0.3

GBS: Now gives sensible error message if run not using perl >= 5.6.0

0.2 

GBS: Added documentation.

GBS: Made pattern-checking optional for C<add_before> and C<add_after>, and
     removed it from the C<replace>.

=head2 Script files

The current set of scripts (in F<$CMISS_ROOT/cmiss_utils/filefixes/20010215/>)
have been rewritten to use this module, and subsequent fixes are most easily
copied from these files.

Each script should be named according to the files to be changed (e.g.
F<fixipopti.pl>) because the list of files to be changed is automatically
generated from the script name.  The file suffix is determined by the letters
between "fixip" and the next "." or "_", allowing multiple scripts for one file
type by including a description that starts with an underscore following the
filetype (e.g. F<fixipopti_10dec99.pl>).

=head2 Loading the module

The following header should be used to load the module.

  #!/usr/local/bin/perl -w
  use strict;
  use CmUtils::IpFileFix qw(fixfiles);

The standard interface is through the I<fixfiles()> routine.  See below in the
L</EXPERTS> section for lower-level interaction.

=cut

#
# Code for finding files.  Modified from old scripts by Shane Blackett.
#
# "common" variables for the find routine "foreachfile"
my ($FILESTOCHANGE, $logfile, $cvs_start, $start_directory, $login, $filter);

=head1 B<fixfiles()>

This routine takes one named argument called I<filter> which should specify a
code reference containing the commands required to "fix" each ipfile.  This is
most easily defined as an anonymous subroutine.  The subroutine is called with
one argument, which is a pointer to the file being processed.

  fixfiles(
    filter => sub {
      my $file = shift;
      ...
    }
  );

Earlier scripts processed the file line-by-line, but this led to confusion and
inefficiencies.  The C<fixfiles> routine reads the entire file before calling
the subroutine defined by I<filter> for processing.

=cut

sub fixfiles {
  my ($filetype) = ($0 =~ /fixip([^_.]*).*pl$/);
  my %params = @_;
  croak 'Must define an filter routine' unless ref($params{filter}) =~ /CODE/;
  $filter = $params{filter};

  my $usage = "Usage: $0 [-h] [-c examples/path_in_cvs]\n";
  my %options;
  if (!(getopts('hc:', \%options))) {
    print $usage;
    exit;
  }
  for (keys %options) {
    if (m/^h$/) {
      print $usage;
      exit;
    } elsif (m/^c$/) {
      # Directory in cvs below which all files will be checked out and tested
      $cvs_start = $options{$_}; 
    }
  }

  $FILESTOCHANGE = qr/(\.i[pr]$filetype|_i[pr]$filetype\.cmiss)$/;
  my $log = "/tmp/fix_$filetype.log";
  if (open LOGFILE, ">>$log") {
    $logfile = $log;
  }
  if (!($login = getlogin())) {
    $login = "unknown";
  }

  if ($cvs_start) {
    #We are going to work with a cvs repository!
    croak 'cvs path must begin with "examples/"'
      unless ($cvs_start =~ m/^examples/);
    croak "examples directory already exists in this directory" if -d "examples";
    croak "you must have CVSROOT defined" unless exists $ENV{CVSROOT};
    croak "$cvs_start not found in repository $ENV{CVSROOT}"
      unless -d "$ENV{CVSROOT}/$cvs_start";
    #Change the search string
    $FILESTOCHANGE =~ s/\$\)$/,v\$\)/;
    $start_directory = $ENV{PWD};
  }
  
  # Find the files (excluding links) needing to be changed
  if ($cvs_start) {
    find(\&foreachfile, "$ENV{CVSROOT}/$cvs_start");
    print "\n\n: CVS final update\n";
    !system("cvs update -P $cvs_start") || die "Unable to prune directories from $cvs_start";
  } else {
    find(\&foreachfile, ".");
  }
  if ($logfile) {
    close LOGFILE;
  }
}

=head1 Modification functions

One or more file modification functions describe the actions to take to update
the ipfiles, and are specified in the I<filter> subroutine in the region C<...>
above.

All functions have a named interface, which means that the arguments may occur
in any order.  Arguments are typically given as strings, although regular
expressions may be used for greater flexibility.  B<In all cases where answers
will be included in the match, regular expressions should be used, as they
allow for various answers/formats/spaces.>  Arguments given as strings will be
quoted in any pattern matches that are done, and are therefore treated as
literals.

Because the file is read into a single multi-line variable, if you want to tie
your regular expression to the start/end of any line, you will need to use the
multi-line C<(?m)> option:

  qr/(?m)^  Enter/

otherwise C<^> and C<$> will only match the start and end of the file.

Alternatively, if you wish to treat the file as one single line (for matching
across more than one line), use the single-line C<(?s)> option:

  qr/(?s)this.*that/

in which case the C<.> will also match newlines.  The C<\n> character may
always be used to specify a known newline.

=head1 B<add_before()>

  $file->add_before(target => "", insert => "" [, unless => ""]);

Adds a string I<insert> before each occurrence string/regexp I<target>.  Most
often, a third string I<unless> needs to be checked so that the addition is
not applied several times.

  $file->add_before(
    target => " USE_MINOS     (0 or 1)",
    insert => " USE_MAPS      (0 or 1)[1]: 0\n",
    unless => " USE_MAPS      (0 or 1)"
  );

  $file->add_before(
    target => qr/(?m)^ USE_MINOS/,
    insert => " USE_MAPS      (0 or 1)[1]: 0\n",
    unless => qr/USE_MAPS/
  );

=over 4

=item unless

String/regexp to test for.  The string insertion is made I<unless> this
matches.

=item target

String/regexp before which the new string is inserted.  As always, a string is
treated as a literal.

=item insert

New string to be inserted.  Makes no sense for this to be a regexp.  Multi-line
strings are most easily defined using C<<<>.

  $file->add_before(
    target => " Do you want to use parallel element stiffness computations",
    insert => <<EOF,
 Specify the convergence criteria [1]:
   (1) Ratio of unconstrained to constrained residuals
   (2) Ratio of unconstrained residuals to maximum Gauss point value
    1
EOF
    unless => "Specify the convergence criteria"
  );

=back

=cut

sub add_before {
  my $file = shift;
  croak 'File not read yet' unless $file->{text};

  my %args = @_;
  croak 'usage: $file->add_before(target => "", insert => "" [, unless => ""])' 
    unless exists $args{target} && exists $args{insert};
  $args{unless} = qr/\Q$args{unless}\E/ unless ref($args{unless}) =~ /regexp/i;
  $args{target} = qr/\Q$args{target}\E/ unless ref($args{target}) =~ /regexp/i;

  unless ($file->{text} =~ /$args{unless}/) {
    $file->{text} =~ s/($args{target})/$args{insert}$1/g;
  }
  return 1;
}

=head1 B<add_after()>

  $file->add_after(target => "", insert => "" [, unless => ""]);

Adds a string I<insert> after each occurrence of another string/regexp
I<target>, unless a third string I<unless> is matched.  As for C<add_before>.

=cut

sub add_after {
  my $file = shift;
  croak 'File not read yet' unless $file->{text};

  my %args = @_;
  croak 'usage: $file->add_after(target => "", insert => "" [, unless => ""])' 
    unless exists $args{target} && exists $args{insert};
  $args{unless} = qr/\Q$args{unless}\E/ unless ref($args{unless}) =~ /regexp/i;
  $args{target} = qr/\Q$args{target}\E/ unless ref($args{target}) =~ /regexp/i;

  unless ($file->{text} =~ /$args{unless}/) {
    $file->{text} =~ s/($args{target})/$1$args{insert}/g;
  }
  return 1;
}

=head1 B<match()>

Returns the result of attempting to match a given I<pattern> in the file.  Most
often, a regular expression should be used, so that varying formatting will
give correct answers.  Results may be captured (using C<()>) for later use if
necessary.

  my ($nodenum) = $file->match(pattern => qr/number of nodes:\s*(\d+)/);
  
  my @values = $file->match(
    pattern => qr/values for C2.*:\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)/
  );

  my @values = map {split} $file->match(pattern => qr/for C2.*:(.*)\n/);

If the argument is a string, it will do a literal string match:

  if ( $file->match(pattern => "number of nodes:\n") ) {
    ...
  }

=cut

sub match {
  my $file = shift;
  croak 'File not read yet' unless $file->{text};

  my %args = @_;
  croak 'usage: $file->match(pattern => ...)' unless exists $args{pattern};
  $args{pattern} = qr/\Q$args{pattern}\E/ unless ref($args{pattern}) =~ /regexp/i;

  return $file->{text} =~ /$args{pattern}/;  
}

=head1 B<delete()>

Deletes the first occurence of the given I<string> from the file.

  $file->delete(string => " (default 1)");
  $file->delete(string => qr/(?m)^ Enter the step limit.*$/);

Global deletion (if ever needed) could be accomplished with:

  my $pattern = "something";
  while ($file->match(pattern => $pattern)) {
    $file->delete(string => $pattern);
  }

=cut

sub delete {
  my $file = shift;
  croak 'File not read yet' unless $file->{text};

  my %args = @_;
  croak 'usage: $file->delete(string => ...)' unless exists $args{string};
  $args{string} = qr/\Q$args{string}\E/ unless ref($args{string}) =~ /regexp/i;

  $file->{text} =~ s/$args{string}//;
  return 1;
}

=head1 B<replace()>

Replaces each occurrence of one given string or regexp with another string.

  $file->replace(target => "nubmer", replace => "number");

=over 4

=item target

String/regexp searched for in the file.  As always, a string is treated as a
literal.

=item replace

Replacement string.  Makes no sense for this to be a regexp.

=back

=cut

sub replace {
  my $file = shift;
  croak 'File not read yet' unless $file->{text};

  my %args = @_;
  croak 'usage: $file->replace(target => ..., replace => ...)'
    unless exists $args{target} && exists $args{replace};
  $args{target} = qr/\Q$args{target}\E/ unless ref($args{target}) =~ /regexp/i;

  $file->{text} =~ s/$args{target}/$args{replace}/g;
  return 1;
}

=head1 EXPERTS

The variable holding the file is actually a blessed hash, with a single entry
C<{text}> containing the entire file as a single string.  Any manipulations
that cannot be done with the provided routines may be performed directly on the
C<$file-E<gt>{text}> entry e.g.

  $file->{text} =~ tr[A-Z][a-z];

Additionally, files may be read and written directly.

=head2 B<new>

Creates a CmUtils::IpFileFix hash.  An optional argument specifies the filename.

  $file = new CmUtils::IpFileFix;
  $file = new CmUtils::IpFileFix ("file.ippara");

=cut

sub new {
  my $class = shift;
  my $file = {};
  bless $file, $class;
  if ( @_ ) {
    $file->read( @_ );
  }
  return $file;
}

=head2 B<read>

Stores a file in the C<{text}> entry of the file hash.

  $file->read($filename);

=cut

sub read {
  my $file = shift;
  croak "usage: file->read(filename)" unless @_ == 1;
  my $filename = shift;
  
  local $/ = undef;
  open INPUT, $filename or croak "Could not open file $filename: $!\n";
  $file->{text} = <INPUT>;
  close INPUT;
}

=head2 B<write>

Writes the C<{text}> entry in the file hash to an output file.

  $file->write("output.ipfile");

=cut

sub write {
  my $file = shift;
  croak 'File not read yet' unless $file->{text};
  croak "usage: file->write(filename)" unless @_ == 1;
  my $filename = shift;

  open OUTPUT, ">$filename" or croak "Could not open file $filename: $!\n";
  print OUTPUT $file->{text};
  close OUTPUT;
}

sub foreachfile {
  # This routine is executed for every file in the tree.  It manages the
  # temporary and old files and log file.
  if ($_ =~ $FILESTOCHANGE) {
    my $date = localtime(time());
    my $keep_ = $_; #Ensure we don't mess this up
    my $filename = $_;
  
    printf ("Converting file %30s .... ",$File::Find::name);
    if ($cvs_start) {
      chdir $start_directory;
      $filename = $File::Find::name;
      $filename =~ s/,v$//;
      $filename =~ s/^$ENV{CVSROOT}\///;
      !system("cvs -Q checkout $filename") || die "Unable to checkout $filename";
    }

# LKC Ensure it is only working on plain files, not soft links
#    if (-s $filename) {
    if (-s $filename) {   
       if(! -l $filename) {
        my $tmpfile = "$filename.tmp";
        my $file = new CmUtils::IpFileFix;
        $file->read($filename);
        if ($filter->($file)) {
          $file->write($tmpfile);
          if (-s $tmpfile) {
            if (compare($filename, $tmpfile) == 1) {
              #Files are different
              rename ($filename, "$filename.old");
              rename ($tmpfile, $filename);
              print " done.\n";
              if ($logfile) {
                print LOGFILE "$login: Updated $File::Find::name $date\n";
               }
            } else {
              #No changes made
              unlink $tmpfile;
              if ($cvs_start) {
                !system("cvs -Q update -r0 $filename") || die "Unable to release $filename";
               }
              print " no changes required.\n";
               if ($logfile) {
                print LOGFILE "$login: Parsed $File::Find::name $date\n";
              }
            }
          } else {
            print "\n";
            die "ERROR : Corrected file is zero size.. reverting.";
          }
        } else {
          print "\n";
          die "ERROR : Update failed.. reverting.";
        }
      }  else {
        print " symbolic link - no change.\n";
      }      
    } else {
      print " empty - no change.\n";
      if ($cvs_start) {
        !system("cvs -Q update -r0 $filename") || die "Unable to release $filename";
      }
    }

    if ($cvs_start) {
    chdir $File::Find::dir;
    }
    $_ = $keep_; #Ensure we don't mess this up
  } elsif ($_ =~ m/^(CVS|Attic)$/) {
     $File::Find::prune = 1;
  }
}

1;
