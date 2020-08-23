eval 'exec perl -wS $0 ${1+"$@"}'
    if 0;

# Finds the list of source files, determining which are local and which are
# global and prints them to <depend-file>, together with dependencies of the
# DSOs, in a format to be included in the cm Makefile.  Also produces object
# list files in a format suitable for the linker containing lists of the
# object file names associated with each DSO and with the static executable.

# Object files are removed if they are no longer required.

my $usage = "usage: $0 <source-directory> <source-patterns> <working-directory> <global-root> <depend-file>\n";

use strict;

use File::Path;
use Fcntl;

# remove path from warning and error messages
BEGIN {
    $^W = 1; # turn on warnings
    $SIG{__WARN__} = sub { $0 =~ m'([^/]*)$'; warn "$1: warning: ",@_ };
}
# die misbehaves if set in BEGIN
$SIG{__DIE__} = sub { $0 =~ m'([^/]*)$'; die "$1: ",@_ };

@ARGV == 5 or die $usage;
my ($source_dir,$source_patterns,$working_dir,$global_root,$depend_file)
  = @ARGV;

# This should be an argument but not changing argument list yet for backward
# compatibility with Makefile.
my $object_list = "$working_dir/object.list";

my $sysname = do {
    if( exists $ENV{SYSNAME} ) {
	$ENV{SYSNAME};
    } else {
	my $output = `uname`; $? and die; chomp $output;
	$output;
    }
};

# my $machname = do {
#     if( exists $ENV{MACHNAME} ) {
# 	$ENV{MACHNAME};
#     } else {
# 	my $output = `uname -m`; $? and die; chomp $output;
# 	$output;
#     }
# };

# printf STDERR "%d: %s", time - $^T, 'searching sources...';

# Make sure we have a pattern so that glob returns nothing if nothing matches.
my @source_patterns =
  map {$_ =~ s%\.(\w)$%\.[$1]%; $_} split " ",$source_patterns ;

my (@local_sources,@global_sources);

foreach my $source_pattern (@source_patterns)
  {
    my @pattern_local =
      map {$_ =~ s%$source_dir/%%; $_} glob "$source_dir/$source_pattern";
    my @pattern_global =
      grep { my $match = $_;
	     ! scalar grep {m/^$match$/} @pattern_local }
	map {$_ =~ s%$global_root/$source_dir/%%; $_}
	  glob "$global_root/$source_dir/$source_pattern";

    if (!(scalar @pattern_local + scalar @pattern_global))
      {
	if ($source_pattern =~ s/fz(.*.f)/fe$1/)
	  {
	    @pattern_local =
	      map {$_ =~ s%$source_dir/%%; $_}
		glob "$source_dir/$source_pattern";
	    @pattern_global =
	      grep { my $match = $_;
		     ! scalar grep {m/^$match$/} @pattern_local }
		map {$_ =~ s%$global_root/$source_dir/%%; $_}
		  glob "$global_root/$source_dir/$source_pattern";
	    @pattern_local =
	      map {$_ =~ s%fe(.*\.f)%fz$1%; $_}
		(@pattern_local,@pattern_global);
	    @pattern_global = ();
	  }
      }
    if (!(scalar @pattern_local + scalar @pattern_global))
      {
	die ("Unable to find match for pattern $source_pattern");
      }
    push(@local_sources, @pattern_local);
    push(@global_sources, @pattern_global);
  }

my @source_list = sort (@local_sources,@global_sources);
my @local_units = map {$_ =~ s%\.[^/]*$%%; $_} map {$_} @local_sources;
my @global_units = map {$_ =~ s%\.[^/]*$%%; $_} map {$_} @global_sources;
my @units = sort (@local_units,@global_units);

my %objects;

foreach my $unit (@global_units)
  {
    $objects{$unit} = "$global_root/$working_dir/$unit.o";
  }
foreach my $unit (@local_units)
  {
    $objects{$unit} = "$working_dir/$unit.o";
  }

my %dso_unit_lists;

# Units containing common symbol definitions (initializations) should be in
# the first DSO to ensure that other DSOs referencing the common symbol do not
# provide a definition (as g77 and GNU ld do with symbols from a common
# blocks).
my $definition_units_regexp = '(00|_INIT)$|/(BLK|BD)';
my $definition_dso_unit = 'fedef';

foreach my $unit (@units)
  {
    my $dso =
      do { if( $unit =~ /$definition_units_regexp/o )
	     { $definition_dso_unit }
	   else
	     { $unit =~ m|^([^/]*)|; $1 }
	 };
    push @{$dso_unit_lists{$dso}}, $unit;
  }

my @dso_units =
  sort {
    if( $a eq $definition_dso_unit )
      { -1 }
    elsif( $b eq $definition_dso_unit )
      { 1 }
    else
      { $a cmp $b }
  } keys %dso_unit_lists;

# printf STDERR "done\n%d: %s", time - $^T, 'updating files...';

# Find a temporary directory if it is used for building objects on a local
# filesystem.  Normally it is faster to build on the destination filesystem
# and rename but GNU ld appears to interact poorly with nfs when writing so
# build somewhere on a local file system first then copy and delete.

# TODO:
# Should use stat -f to check whether the destination filesystem is local.

my $mytmpdir;
if( $sysname eq 'Linux' )
  {
    # Assuming the current directory is the root of the tree.  It is possible
    # that the user has made symlinks from this tree to the local filesystem,
    # but, if so, they are most likely using the same filesystem as TMPDIR so
    # there should be no significant penalty in moving the file between
    # directories.
    my $filesystem_type = `stat --file-system --format="%T" .`;
    $? and warn "failed to determine filesystem type from command stat\n";
    chomp $filesystem_type;

    if( $filesystem_type eq 'nfs' )
      {
	my $tmpdir = $ENV{TMPDIR} || '/tmp';

	# If the temporary directory is only writable by the user then it makes it
	# easier to ensure that no one else creates a file of the same name as
	# intended in the build.

	if( (stat $tmpdir)[2] & 0022 )
	  {		     # Writable by others.  Make our own subdirectory.
	    # Assume that the user has a TMPDIR having the sticky bit set.
	    my $username = getpwuid $>;

	    $mytmpdir = "$tmpdir/$username";
	    # We're fine if this is already a directory and not a link and not
	    # writable by others.
	    unless( lstat $mytmpdir && -d _ && !( (stat _)[2] & 0022 ) )
	      {
		print "mkdir $mytmpdir\n";
		mkdir $mytmpdir, 0777 or die "can't mkdir $mytmpdir: $!\n";
	      }
	  }
	else
	  {
	    $mytmpdir = $tmpdir;
	  }
      }
  }

#===========================================================================
# Subroutines for writing to files only if they need updating.

sub update_open (\$$)
  {
    my ($handle,$filename) = @_;

    local *FH;
    my $changed = 0;
    # When using open with +< it doesn't create the file when it doesn't
    # exist, and with +>> the file pointer seems to get moved to the end of
    # the file with perl 5.004_05 for irix at least.  sysopen can give what we
    # want.
    sysopen FH, $filename, O_CREAT|O_RDWR
      or die "can't open +>> $filename: $!\n";
#     seek FH, 0, 0 or die "seek failed on $filename: $!\n";

    $$handle = { filename => $filename, fh => *FH, changed => $changed };
  }

sub update_print ($@)
  {
    my ($handle,@list) = @_;

    my $fh = $handle->{fh};
    my $changed = \$handle->{changed};

    my $new_data = join defined($,) ? $, : '', @list;

    if( ! $$changed )
      {
	my $data_beg = tell $fh;

	defined( read $fh, my ($old_data), length $new_data )
	  or die "error reading $handle->{filename}: $!\n";

	if( $new_data ne $old_data )
	  {
	    # SEEK_SET is not available in perl 5.004_05
	    seek $fh, $data_beg, 0
	      or die "seek failed on $handle->{filename}: $!\n";
	    $$changed = 1;
	  }
      }

    if( $$changed )
      {
	print { $fh } $new_data
	  or die "error printing to $handle->{filename}: $!\n";
      }
  }

sub update_close ($)
  {
    my ($handle) = @_;

    my $fh = $handle->{fh};
    my $changed = \$handle->{changed};

    my $new_end = tell $fh;

    if( ! $$changed )
      {
	if( $new_end ne -s $fh )
	  {
	    $$changed = 1;
	  }
      }

    if( $$changed )
      {
	truncate $fh, $new_end
	  or die "can't truncate $handle->{filename}: $!\n";
      }

    close $fh or die "error closing $handle->{filename}: $!\n";
  }

#===========================================================================

sub update_object_list ($@)
  {
    my ($file,@objects) = @_;

    update_open my ($handle), $file;
#    my $is_linker_script = $machname ne 'ia64' && $sysname eq 'Linux' || $sysname eq 'win32';
    my $is_linker_script = $sysname eq 'Linux' || $sysname eq 'win32';

    if( $is_linker_script )
      {
	update_print $handle, "INPUT(\n";
      }

    foreach my $object (@objects)
      {
	update_print $handle, "$object\n"
      }

    if( $is_linker_script )
      {
	update_print $handle, ")\n";
      }

    update_close $handle;
  }

#===========================================================================

open DEPEND, "> $depend_file" or die "can't open > $depend_file: $!\n";

print DEPEND "SOURCES_LOCAL := @local_sources\n";
print DEPEND "SOURCES_GLOBAL := @global_sources\n";
print DEPEND "SOURCES := @source_list\n";
print DEPEND "UNITS := @units\n";

print DEPEND "DSO_UNITS := " . join(" ", @dso_units) . "\n";
print DEPEND "DSOS := " .
  join(" ", map {$_ = "\$(DSO_DIR)/lib$_.so"} map {$_} @dso_units) . "\n";

if( defined $mytmpdir )
  {
    print DEPEND "get_build_filename = $mytmpdir/\$(subst /,!,\$(CURDIR)!\$(dir \$1))\$(notdir \$1)\n";
  }
else
  {
    print DEPEND "get_build_filename = \$1.build\n";
  }

print DEPEND "ifeq (\$(TASK),objects)\n";

foreach my $dso (@dso_units)
  {
    my @objects = map { $objects{$_} } @{$dso_unit_lists{$dso}};
    print DEPEND "\$(call get_build_filename,\$(DSO_DIR)/lib$dso.so): @objects\n";

    update_object_list "$working_dir/$dso.olist", @objects;
  }

print DEPEND "endif\n";

close DEPEND or die "error closing $depend_file: $!\n";

update_object_list( $object_list, map { $objects{$_} } @units );

#===========================================================================

# printf STDERR "done\n%d: %s", time - $^T, 'tidying working directory...';

my %local_dirs_required;
@local_dirs_required{ map { $_ =~ m{^(.*)/} ? $1 : () } @local_units } = undef;

my %local_units;
@local_units{ @local_units } = undef;

sub possibly_remove_working_file ($$)
  # Removes a specified file in the generated directory if no longer required.
  # Returns 1 if the file still exists or 0 if removed.
  {
    my ($subname,$fullname) = @_;

    my ($subroot,$ext) = $subname =~ m{^(.*)\.(.*)$};

    my $keep = do
      {
	if( ! defined $ext )
	  { # No extension; not our file so don't remove it
	    1
	  }
	elsif( $ext =~ /^[fod]$/ )
	  { # Seems to have been generated from a source file.
	    # Keep if the unit is still local.
	    exists $local_units{$subroot};
	  }
	elsif( $ext eq 'cmn' )
	  { # Looks like a preprocessed Fortran common file.
	    # Keep if the file is still local.
	    -f "$source_dir/$subname";
	  }
	elsif( $ext eq 'olist' )
	  { # List of objects in dso; keep if dso still should exist.
	    exists $dso_unit_lists{$subroot};
	  }
	else
	  { # Not recognized; may not be ours; leave it.
	    1;
	  }
      };

    unless( $keep )
      { print "unlink $fullname\n";
	unlink $fullname or do
	  { warn "can't unlink $fullname: $!\n";
	    $keep = 1;
	  };
	# TODO: include something to force remaking of
	# anything created from these files.
      }

    return $keep;
  }

sub prune_working_dir ($$); # separate prototype because sub is recursive
sub prune_working_dir ($$)
  # Removes files in the specified subdir of the working directory if no
  # longer required.  Returns 1 if the directory still contains files or 0
  # otherwise.
  {
    my ($subdir,$fulldir) = @_;

    opendir DIR, $fulldir or do
      { warn "can't opendir $fulldir: $!\n";
	return 1;
      };

    my @contents = readdir DIR;

    closedir DIR;

    my $kept; # no contents kept yet

    foreach my $item ( @contents )
      {
	if( $item =~ m{^\.\.?$} )
	  { #ignore
	    next;
	  }

	my $subname = $subdir eq '.' ? $item : "$subdir/$item";
	my $fullname = "$fulldir/$item";

	if( -d $fullname )
	  { if( prune_working_dir $subname, $fullname )
	      { $kept = 1;
		delete $local_dirs_required{subname};
	      }
	    else # empty directory
	      { rmdir $fullname or do
		  { warn "can't rmdir $fullname: $!\n";
		    $kept = 1;
		  };
	      }
	  }
	else # file (not directory)
	  { possibly_remove_working_file $subname, $fullname and $kept = 1;
	  }
      }

    return $kept;
  }

prune_working_dir '.',$working_dir;

foreach my $dir ( keys %local_dirs_required )
  {
    my $fullname = "$working_dir/$dir";

    -d $fullname or do
      { print "mkdir $fullname\n";
#	mkdir -p $fullname, 0777 or die "can't mkdir $fullname: $!\n";
	&File::Path::mkpath($fullname) or die "can't mkdir $fullname: $!\n";
      };
  }

# printf STDERR "done %d\n", time - $^T;
