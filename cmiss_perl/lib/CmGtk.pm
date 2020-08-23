package CmGtk;
require 5.006;

=head1 CmGtk

Core routines to initialise and start Gtk running under CMGUI.  Usage:

  use CmGtk;
  CmGtk->init();    # or   init CmGtk;
  CmGtk->main();    # or   main CmGtk;

=cut

BEGIN {
  use Exporter ();
  our ($VERSION, @ISA, @EXPORT);
  $VERSION = 1.00;
  @ISA     = qw( Exporter );
  @EXPORT  = qw( );
}

eval "use Gtk";

our $CmGtkstarted = 0;

sub _step {
  while (Gtk->events_pending()) {
    Gtk->main_iteration();
  }
  return ("OK");
}

=head1 B<init>

Initialise Gtk to run under CMISS.

=cut

sub init {
  my $i;
  my $index;
  my $original_env_display;
  my @original_display_args;
  cmiss::cmiss("attach gtk start_detection perl_action &CmGtk::_step()");
  #Change the X form of the argument for the Gtk one
  for ($i = 0 ; $i < scalar @ARGV - 1 ; $i++)
  {
	 my $arg = $ARGV[$i];
	 if ($arg =~ m/^-display$/)
	 {
	   $index = $i;

	   #In older versions of Gtk the --display arg doesn't seem
	   #to work correctly.
	   $original_env_display = $ENV{DISPLAY};
	   $ENV{DISPLAY} = $ARGV[$i+1];

	   @original_display_args = ($ARGV[$i], $ARGV[$i+1]);
	   $ARGV[$i] = "--display=$ARGV[$i+1]";
	   splice (@ARGV, $i+1, 1);
	 }
  }
  Gtk->init();
  #Put the arguments back how we found them
  if (defined @original_display_args)
  {
    $ENV{DISPLAY} = $original_env_display;
    splice (@ARGV, $index, 1, @original_display_args);
  }
  cmiss::cmiss("attach gtk end_detection");
  Gtk->set_locale();
}

=head1 B<main>

Process all Gtk events to get loop started.

=cut

sub main {
  CmGtk::_step();
  $CmGtkstarted = 1;
#   die("return to main cmgui mainloop");
}

1;
