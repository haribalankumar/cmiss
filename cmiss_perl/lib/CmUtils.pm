package CmUtils;
require 5.006;

=head1 CmUtils

Some basic utility functions for the CMISS Perl interpreter. Includes functions
for directories, creating lists, and printing defined variables.

=cut

use strict qw(vars);
if (eval {require warnings}) {
  eval 'use warnings';
} else {
  print STDERR <<EOF;
WARNING: You really should be using modules for perl 5.6.0.
         Try "setenv CMISS_PERL /usr/local/perl5.6/bin/perl" or
         whatever your path to perl5.6 is.

         We'll try to make things work anyway, but it won't be as nice!

EOF
}
use CmUtils::Pod::Usage qw(pod2usage);
my $Dumper = 1;
if (eval 'use Data::Dumper',$@) {
  $Dumper = 0;
  print STDERR <<EOF;
WARNING: Module Data::Dumper would be helpful.
         Please install, or upgrade to perl 5.005 or later.

EOF
}

BEGIN {
  use CmUtils::Exporter ();
  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION = 1.10;

=head1 VERSION

1.10 (7 September 2001)

=head1 CHANGES

=item 1.10

GBS:

Added B<max>, B<min> and B<constrain> as utility routines.

Added B<choose> routine to prompt the user to choose between several options.
Enabled this routine and B<ask> to use Gtk2 if running under cmgui.

=item 1.09

Added B<pushdir> and B<popdir> routines for keeping track of directory changes.

=item 1.08

Added functionality when external modules (Data::Dumper etc) are not
available.

Test for valid variables in pdefined commands is now on the definition,
rather than whether a variable is true or not.  

=item 1.07  

Changed package name from "CMISS::Utils" to "CmUtils"

Moved help routines back to this module.

=item 1.06 

Added "pother" to print other refs (usually objects)

=item 1.05  

Moved "h" and "man" routines to the Exporter module.

=item 1.04

Changed to use "CmUtils::Exporter" as the ISA - gives automatic
definition of subroutines and help routines

Added "man" function for help on Perl modules

=head1 SUBROUTINES

=cut

  # Base class for module
  @ISA         = qw(CmUtils::Exporter);

  # Default subroutines and variables to export to caller
  @EXPORT      = qw(
    &max &min &constrain
    &h &man
    &ask &choose
    &createwin 
    &getcurrentpath &getcurrentdir &cd &pushdir &popdir
    &step &range
    &proutine &phash &parray &pscalar &pother &pdefined &penv
  );

  # Routines and variables exportable on request
  @EXPORT_OK   = qw();

  # Import a subset of routines/variables
  # e.g. use CmUtils qw(:listvar)
  #      would only import the routines which list out defined variables etc  
  %EXPORT_TAGS = (
    listvar => [qw(&proutine &phash &parray &pscalar &pother &pdefined &penv)],
  );

}

=head2 B<max(list)>

Returns the maximum of a list of numbers.

=cut
# GBS
sub max (@) {
  my ($max);
  foreach (@_) {
    !defined $max || $max < $_ and $max = $_;
  }
  $max;
}

=head2 B<min(list)>

Returns the minimum of a list of numbers.

=cut
# GBS
sub min (@) {
  my ($min);
  foreach (@_) {
    !defined $min || $min > $_ and $min = $_;
  }
  $min;
}

=head2 B<constrain(value, minimum, maximum)>

Constrains a numeric variable to lie between two other values.
Ensures min < max.

=cut
# GBS
sub constrain {
  my ($val, $min, $max) = @_;
  ($min, $max) = ($max, $min) if $min > $max;
  return $min if $val < $min;
  return $max if $val > $max;
  return $val;
}

=head2 B<h([name][,topic])>

Searches for, and prints, documentation for the NAMEd subroutine/module/object
from the CMISS perl modules.

"name" can be either a subroutine name (e.g. C<'pdefined'>), a module name
(e.g. C<'CmUtils::Utils'>) or an object (e.g. C<$nodegroup>).

Additional help available is listed, and can be accessed by specifying a topic
along with the name e.g.

  h 'CmUtils::Objects::NodeGroup';
  h 'CmUtils::Objects::NodeGroup','list';

=cut
# GBS
sub h {
  my $name = shift || 'h';
  my ($topic, $msg);
  if (ref $name) {$name = ref $name};
  if ($CmUtils::Exporter::SUBS{$name}) {
    $msg = qq/Subroutine "$name" defined in $CmUtils::Exporter::SUBS{$name}/;
    $topic = '\b' . $name . '\b';
    pod2usage( {
      -verbose => -1,
      -exitval => -1,
      -input   => $CmUtils::Exporter::SUBS{$name},
      -topic   => ".*$topic.*",
      -msg     => $msg,
    } );
  } else {
    # See if a module is defined with this name
    if (@_) {
      $topic = shift;
      $msg = qq/Looking for help on "$topic" in $name/;
      pod2usage( {
        -verbose => -2,
        -exitval => -1,
        -input   => $name,
        -topic   => ".*$topic.*",
        -msg     => $msg,
      } );
    } else {
      $msg = "\nModule:";
      my $pod = pod2usage( {
        -verbose => -2,
        -exitval => -1,
        -input   => $name,
        -topic   => ".*$name.*",
        -msg     => $msg,
      } );
      $msg = "Help available on:";
      $pod = pod2usage( {
        -verbose       => -2,
        -exitval       => -1,
        -input         => $name,
        -topic         => ".*",
        -headings_only => 2,
        -msg           => $msg,
      } ) if $pod;
#       $pod = pod2usage( {
#         -verbose       => -2,
#         -exitval       => -1,
#         -input         => $name,
#         -topic         => "[BF]<.*",
#         -headings_only => 2,
#         -msg           => $msg,
#       } ) if $pod;
      if ($pod) {
        my @parents = @{"${name}::ISA"};
        my $i=0;
        while (my $parent = $parents[$i]) {
          push @parents, @{"${parent}::ISA"};
          $i++;
        }
        if (@parents) {
          print join "\n  ", "$name is derived from:", @parents, "\n";
        }
      }
    }
  }
}

=head2 B<man(name)>

Prints the full manpage for the NAMEd module.

=cut
# GBS
sub man {
  my $name = shift || 'Overview.pod';
  pod2usage( {
    -verbose => 2,
    -exitval => -1,
    -input => $name,
  } );
  print "\n";
}


###
### Functions related to commands
###   ask
###

=head2 B<ask(question, [default])>

Prompts for user-input with QUESTION.  Returns the entered string, or DEFAULT
if nothing is entered.  Now uses Gtk2 if running cmgui, otherwise prompts in
the terminal window.

=cut
# GBS
sub ask {
  use CmUtils::Utils qw(runningCmgui);
  my ($question, $default) = @_;
  my $prompt = join "", "$question ", defined $default ? "[$default] " : "", "? ";
  my $answer;
  if (runningCmgui()) {
    use CmGtk2;
    unless ($CmGtk2::CmGtk2started) {
      CmGtk2->init();
    }
            # A dialog box
    my $window = new Gtk2::Dialog();
    $window->set_title("Question:");
    $window->set_position(mouse);
            # size width,height 300,100 
    $window->set_default_size(300,100);
            # with a prompt string at the top
    my $label = new Gtk2::Label($prompt);
    $window->vbox->add($label);
            # and a text entry area, with the default string
    my $text = new Gtk2::Entry();
    $text->set_text($default);
    $text->select_region(0, length $default);
    $text->signal_connect(activate => 
      sub { 
        $answer = $text->get_text();
        $window->destroy();
        undef $window;
      });
    $window->action_area->add($text);
    $window->show_all();

    do {
      CmGtk2->main();
      $text->grab_focus();
    } while $window;
#     print "Answer: $answer\n";
  } else {
    print STDOUT $prompt;
    $answer = <STDIN>;
    print "\n";
    chomp $answer;
  }
  return $answer || $default;
}

=head2 B<choose(question, default, choices)>

Prompts for user to choose between several options.  CHOICES is an array
of options, with DEFAULT being the default option in that array (first is
zero).  Uses Gtk2 if running cmgui, otherwise prompts in terminal window.

=cut
# GBS

# AJC 
sub choose {
  use CmUtils::Utils qw(runningCmgui);
  my ($question, $default, @choices) = @_;
  my $answer;
  if (runningCmgui()) {
    use CmGtk2;
    CmGtk2->init() unless $CmGtk2::CmGtk2started;
            # A dialog box
    my $window = new Gtk2::Dialog();
    $window->set_title("Choose:");
    $window->set_position(mouse);
            # size width,height 300,100 
    $window->set_default_size(300,100);
    
            # with a prompt at the top
    my $label = new Gtk2::Label($question);
    $window->vbox->add($label);
    my $row;
    use Gtk2::SimpleList;
    my $list = Gtk2::SimpleList->new('Choice' => 'scalar');
    $list->get_selection->set_mode('browse');
#    $list->signal_connect(select_row => sub { shift; $row = shift; });
    $list->set_data_array(\@choices);
    $list->select($default);

    $window->vbox->add($list);
            # and a button to apply the choice
    my $button = new Gtk2::Button("OK");
    $window->action_area->add($button);
    $button->signal_connect(clicked =>
      sub {
        my @answers;
        @answers = $list->get_selected_indices;
	$answer = $choices[$answers[0]];
        $window->destroy();
        undef $window;
      });


    $window->show_all();
    do {
      CmGtk2->main();
    } while $window;

  # not running Cmgui  
  } else {
    do {
      print STDOUT "$question (* is default)\n";
      my $num = 0;
      foreach (@choices) {
        print $num == $default ? "*" : " ";
        printf "[%2d]:\t$_\n", $num;
        $num++;
      }
      print "? ";
      $answer = <STDIN>;
      print "\n";
      chomp $answer;
    } until ($answer >=0 and $answer < @choices);
    $answer = $choices[$answer eq '' ? $default : $answer];
  }
  return $answer;
}


sub choose1 {
  use CmUtils::Utils qw(runningCmgui);
  my ($question, $default, @choices) = @_;
  my $answer;
  if (runningCmgui()) {
    use CmGtk2;
    CmGtk2->init() unless $CmGtk2::CmGtk2started;
            # A dialog box
    my $window = new Gtk2::Dialog();
    $window->set_title("Choose:");
    $window->set_position(mouse);
            # size width,height 300,100 
    $window->set_default_size(300,100);
#    $window->default_width(300);
    
            # with a prompt at the top
    my $label = new Gtk2::Label($question);
    $window->vbox->add($label);
            # and a list of choices
    my $row;
    my $list = new Gtk2::CList(1);
    $list->set_selection_mode(browse);
    $list->signal_connect(select_row => sub { shift; $row = shift; });
    $list->append($_) foreach @choices;
    $list->select_row($default, 0);
    $window->vbox->add($list);
            # and a button to apply the choice
    my $button = new Gtk2::Button("OK");
    $window->action_area->add($button);
    $button->signal_connect(clicked =>
      sub {
        $answer = $list->get_text($row, 0);
        $window->destroy();
        undef $window;
      });
    
    $window->show_all();
    do {
      CmGtk2->main();
    } while $window;
#     print "Answer: $answer\n";
  } else {
    do {
      print STDOUT "$question (* is default)\n";
      my $num = 0;
      foreach (@choices) {
        print $num == $default ? "*" : " ";
        printf "[%2d]:\t$_\n", $num;
        $num++;
      }
      print "? ";
      $answer = <STDIN>;
      print "\n";
      chomp $answer;
    } until ($answer >=0 and $answer < @choices);
    $answer = $choices[$answer eq '' ? $default : $answer];
  }
  return $answer;
}

###
### Graphical functions
###   createwin
###

=head2 B<createwin(name)>

Creates a NAMEd 3d window in Cmgui, default is window 1.  Useful in com files
as it doesn't return an warning (and die) if the window already exists.

=cut
# GBS
# doesn't crash on creation error (if window already exists)
sub createwin {
  my $name = shift;
  $name = 1 unless (defined($name));
  if (cmiss::cmiss("gfx create window \"$name\"")) {
    print "Creating window $name\n";
  }
}

###
### Directory functions
###   getcurrentpath getcurrentdir cd
###

=head2 B<getcurrentpath()>

Returns the current directory path (no trailing slash).

=cut
# GBS
sub getcurrentpath {
  use Cwd;
  return getcwd();
}

=head2 B<getcurrentdir()>

Returns the name of the current directory.

=cut
# GBS
sub getcurrentdir {
  use Cwd;
  my ($dir)=(getcwd() =~ m(/([^/]*)$));
  return $dir;
}

# list of directories for pushing and popping
my @dirlist = ();

=head2 B<pushdir(dirname)>

Stores current directory and changes current directory to DIRNAME.

=cut
# GBS
sub pushdir {
  my $dir = shift;
  if (defined $dir && -d $dir) {
    push @dirlist, getcurrentpath();
    return cd($dir);
  } else {
    print STDERR "Directory $dir does not exist\n";
  }
}

=head2 B<popdir()>

Changes to the previous directory.

=cut
# GBS
sub popdir {
  return unless @dirlist;
  return cd(pop @dirlist);
}

=head2 B<cd(dirname)>

Changes current directory to DIRNAME for both perl and CMISS.  Cmgui handles
this alright already.

=cut
# GBS
# CMGUI handles directory changing OK, but CMISS doesn't (yet!) => use "cd"
sub cd {
  my $dir = shift;
  if (defined $dir && -d $dir) {
    chdir $dir;
    $dir = getcurrentpath();
    cmiss::cmiss("set dir $dir");
  } else {
    print STDERR "Directory $dir does not exist\n";
  }
  return getcurrentpath();
}

###
### List creation functions
###   step range
###

=head2 B<step(start,end,step)>

Returns a list of numbers from START to END in increments of STEP (default 1).

  foreach $value (step(4,20,4)) {  # 4,8,12,16,20

=cut
# GBS
# step(1,10,2) returns values from 1 to 10 in steps of 2 (1,3,5,7,9)
# stepsize defaults to 1
sub step {
  my ($start, $end, $step) = @_;
  $step = 1 unless defined($step);
  return $start if $step == 0;
  return map {$_*$step+$start} ( 0 .. ($end-$start)/$step );
}

=head2 B<range(start,end,number)>

Returns a list of NUMBER values evenly spaced between START and END.

  foreach $value (range(4,20,3)) {  # 4,12,20

=cut
# GBS
# range(0,10,6) returns 6 values from 0 to 10 (0,2,4,6,8,10)
sub range {
  my ($start, $end, $num) = @_;
  $num -= 1;
  return map {$_*($end-$start)/$num+$start} 0..$num;
}

###
### Functions for examining defined objects
###   proutine phash parray pscalar pother pdefined penv
###

=head2 B<proutine(name)>

Prints a list of NAMEd or all defined subroutines, along with their stored
help-text if available.  Specify name as a string.

=cut
# GBS
sub proutine {
  package cmiss;
  my @name = @_;

  # either use passed names, or grab all keys from cmiss:: that are subroutines and
  # create an array containing info for printing
  unless (@name) { @name = keys %cmiss::};
  my @subs = sort grep {exists &$_} @name;
  print "\nSubroutines:\n" if @subs;

  # find maximum length of first elements of the array
  my $maxlen = CmUtils::max( map{ length } @subs);
  foreach my $key (@subs) {
    if (exists $CmUtils::Exporter::SUBS{$key}) {
      print "$key:  "," " x ($maxlen-length($key)), "$CmUtils::Exporter::SUBS{$key}\n";
    } else {
      print "$key:\n";
    }
  }
  print "\n" if @subs;
}

=head2 B<phash(name)>

Prints a list of NAMEd or all defined hashes and hashrefs and their values.
Specify name as a string.

=cut
# GBS
sub phash {
  package cmiss;
  my @name = @_;
#   use Data::Dumper;
  local $Data::Dumper::Indent=1;

  # either use passed names or grab all keys from cmiss:: 
  # that are hashes or hashrefs
  unless (@name) { @name = keys %cmiss::};

  my @hashkeys = sort grep {defined %$_} @name;
  print "\nHashes:\n" if @hashkeys;
  foreach my $key (@hashkeys) {
    if ($Dumper) {
      print Data::Dumper->Dump([\%$key],["*$key"]);
    } else {
      print "\%$key:\n";
      foreach my $k (keys %$key) {
        print "  $k : $$key{$k}\n";
      }
    }
  }

  my @hashrefs = sort grep {defined $$_ && ref($$_)=~/HASH/} @name;
  print "Hashrefs:\n" if @hashrefs;
  foreach my $key (@hashrefs) {
    if ($Dumper) {
      print Data::Dumper->Dump([$$key],["$key"]);
    } else {
      print "\$$key:\n";
      foreach my $k (keys %{$$key}) {
        print "  $k : $$key->{$k}\n";
      }
    }
  }

  print "\n" if (@hashkeys + @hashrefs);
}

=head2 B<parray(name)>

Prints a list of NAMEd or all defined arrays and arrayrefs and their values.
Specify name as a string.

=cut
# GBS
sub parray {
  package cmiss;
  my @name = @_;
#   use Data::Dumper;
  local $Data::Dumper::Indent=1;

  # either use passed names or grab all keys from cmiss:: 
  # that are arrays or arrayrefs
  unless (@name) { @name = keys %cmiss::};

  my @arraykeys = sort grep {defined @$_} @name;
  print "\nArrays:\n" if @arraykeys;
  foreach my $key (@arraykeys) {
    if ($Dumper) {
      print Data::Dumper->Dump([\@$key],["*$key"]);
    } else {
      print "\@$key:\n  ";
      print join "\n  ", @$key;
    }
  }

  my @arrayrefs = sort grep {defined $$_ && ref($$_)=~/ARRAY/} @name;
  print "\nArrayrefs:\n" if @arrayrefs;
  foreach my $key (@arrayrefs) {
    if ($Dumper) {
      print Data::Dumper->Dump([$$key],["$key"]);
    } else {
      print "\$$key:\n  ";
      print join "\n  ", @$$key;
    }
  }

  print "\n" if (@arraykeys + @arrayrefs);
}

=head2 B<pscalar(name)>

Prints a list of NAMEd or all defined scalars and their values, expect those
which are references.  Specify name as a string.

=cut
# GBS
sub pscalar {
  package cmiss;
  my @name = @_;
  
  # either use passed names or grab all keys from cmiss:: that are scalars
  unless (@name) { @name = keys %cmiss::};
  
  my @scalarkeys = sort grep {defined $$_ && !ref($$_)} @name;
  print "\nScalars:\n" if @scalarkeys;
  # find maximum length of scalar variable names
  my $maxlen = CmUtils::max(map {length} @scalarkeys);
  foreach my $key (@scalarkeys) {
    print "\$$key ", " " x ($maxlen-length($key)), "= $$key\n";
  }
  print "\n" if @scalarkeys;
}

=head2 B<pother(name)>

Prints a list of NAMEd or all variables other than those printed above.
Specify name as a string.

=cut
# GBS
sub pother {
  package cmiss;
  my @name = @_;
  
  # either use passed names or grab all keys from cmiss:: that are refs, 
  # but not Hash, Array or Scalar refs
  unless (@name) { @name = keys %cmiss::};

  my @keys = sort grep { 
    defined $$_ && 
    ref($$_) && 
    ref($$_) !~ /HASH|ARRAY|SCALAR/ 
    } @name;
  print "\nOther refs and objects:\n" if @keys;
  # find maximum length of variable names
  my $maxlen = CmUtils::max(map {length} @keys);
  foreach my $key (@keys) {
    print "\$$key ", " " x ($maxlen-length($key)), ": ",ref $$key, "\n";
  }
  print "\n" if @keys;
}

=head2 B<pdefined(name)>

Prints a list of NAMEd or all defined variables and subroutines and their values.
Specify name as a string.

  pdefined('INC');

=cut
# GBS
sub pdefined {
  &proutine;
  &phash;
  &parray;
  &pscalar;
  &pother;
}

=head2 B<penv(name)>

Prints a list of NAMEd or all defined environment variables.

=cut
# GBS
sub penv {
  my @name = @_;

  unless (@name) { @name = keys %ENV};
  my @scalarkeys = sort grep {$ENV{$_}} @name;

  my $maxlen = max(map {length} @scalarkeys);
  foreach my $key (@scalarkeys) {
    print "$key:  ", " " x ($maxlen-length($key)), "$ENV{$key}\n";
  }
}

1;
