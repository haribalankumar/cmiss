package CmGtk2;
require 5.008;

=head1 CmGtk2

Core routines to initialise and start Gtk2 running under CMGUI.  Usage:

  use CmGtk2;
  CmGtk2->init();    # or   init CmGtk2;
  CmGtk2->main();    # or   main CmGtk2;

=cut

BEGIN {
  use Exporter ();
  our ($VERSION, @ISA, @EXPORT);
  $VERSION = 1.00;
  @ISA     = qw( Exporter );
  @EXPORT  = qw( );
}

eval "use Gtk2";

our $CmGtk2started = 0;

sub _step {
  while (Gtk2->events_pending()) {
    Gtk2->main_iteration();
  }
  return ("OK");
}

=head1 B<init>

Initialise Gtk2 to run under CMISS.

=cut

sub init {
  my $i;
  my $index;
  my $original_env_display;
  my @original_display_args;
  cmiss::cmiss("attach gtk start_detection perl_action &CmGtk2::_step()");
  #Change the X form of the argument for the Gtk2 one
  for ($i = 0 ; $i < scalar @ARGV - 1 ; $i++)
  {
	 my $arg = $ARGV[$i];
	 if ($arg =~ m/^-display$/)
	 {
	   $index = $i;

	   #In older versions of Gtk2 the --display arg doesn't seem
	   #to work correctly.
	   $original_env_display = $ENV{DISPLAY};
	   $ENV{DISPLAY} = $ARGV[$i+1];

	   @original_display_args = ($ARGV[$i], $ARGV[$i+1]);
	   $ARGV[$i] = "--display=$ARGV[$i+1]";
	   splice (@ARGV, $i+1, 1);
	 }
  }
  Gtk2->init();
  #Put the arguments back how we found them
  if (defined @original_display_args)
  {
    $ENV{DISPLAY} = $original_env_display;
    splice (@ARGV, $index, 1, @original_display_args);
  }
  cmiss::cmiss("attach gtk end_detection");
  Gtk2->set_locale();
}

=head1 B<main>

Process all Gtk2 events to get loop started.

=cut

sub main {
  CmGtk2::_step();
  $CmGtk2started = 1;
#   die("return to main cmgui mainloop");
}

1;
