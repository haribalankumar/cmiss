# CmUtils::Pod::Text::Mono -- Overrides for printing POD in Mono
#
# Modified by Greg Sands 2001. from Pod::Text::Color
#
# Copyright 1999 by Russ Allbery <rra@stanford.edu>
#
# This program is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.

############################################################################
# Modules and declarations
############################################################################

package CmUtils::Pod::Text::Mono;

require 5.004;

use Pod::Text ();

use strict;
use vars qw(@ISA $VERSION);

@ISA = qw(Pod::Text);
$VERSION = 1.00;

############################################################################
# Overrides
############################################################################

# Fix the various interior sequences.
sub seq_b { return '*' . $_[1] .'*'  }
sub seq_c { return $_[1] }
sub seq_f { return '"' . $_[1] .'"'   }
sub seq_i { return '*' . $_[1] .'*' }
sub seq_l { return '_' . $_[1] .'_' }

############################################################################
# Module return value and documentation
############################################################################

1;
__END__

=head1 NAME

CmUtils::Pod::Text::Mono - Overrides for printing POD in Mono

=head1 SYNOPSIS

    use CmUtils::Pod::Text::Mono;
    my $parser = CmUtils::Pod::Text::Mono->new (sentence => 0, width => 78);

    # Read POD from STDIN and write to STDOUT.
    $parser->parse_from_filehandle;

    # Read POD from file.pod and write to file.txt.
    $parser->parse_from_file ('file.pod', 'file.txt');

=head1 DESCRIPTION

CmUtils::Pod::Text::Mono is a simple subclass of Pod::Text that overrides a few
routines so as to be more useful in CMISS.  See L<Pod::Text> for details and
available options.

=head1 SEE ALSO

L<Pod::Text|Pod::Text>, L<Pod::Parser|Pod::Parser>

=head1 AUTHOR

Greg Sands

=cut
