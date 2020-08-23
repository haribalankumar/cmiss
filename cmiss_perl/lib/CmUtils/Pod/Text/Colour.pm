# CmUtils::Pod::Text::Colour -- Convert POD data to formatted colour ASCII text
#
# Modified by Greg Sands 2001. from Pod::Text::Colour
# Copyright 1999 by Russ Allbery <rra@stanford.edu>
#
# This program is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.

############################################################################
# Modules and declarations
############################################################################

package CmUtils::Pod::Text::Colour;

require 5.004;

use Pod::Text ();
use Term::ANSIColor qw(colored);

use strict;
use vars qw(@ISA $VERSION);

@ISA = qw(Pod::Text);
$VERSION = 1.10;

############################################################################
# Overrides
############################################################################

# Make level one headings bold.
sub cmd_head1 {
    my $self = shift;
    local $_ = shift;
    s/\s+$//;
    $self->SUPER::cmd_head1 (
      colored('*', 'black on_yellow') . ' ' . colored ($_, 'bold'));
}

# Make level two headings bold.
# sub cmd_head2 {
#     my $self = shift;
#     local $_ = shift;
#     s/\s+$//;
#     s/^/  /;
#     $self->SUPER::cmd_head2 (colored ($_, 'bold'));
#     $_ = $self->interpolate ($_, shift);
#     $self->output (' ' x ($$self{indent} / 2) . $_ . "\n");
# }

# Fix the various interior sequences.
sub seq_b { return colored ($_[1], 'bold')   }
sub seq_c { return colored ($_[1], 'cyan') }
sub seq_f { return colored ($_[1], 'green')   }
sub seq_i { return colored ($_[1], 'underline') }
sub seq_l { return colored ($_[1], 'cyan underline') }

# We unfortunately have to override the wrapping code here, since the normal
# wrapping code gets really confused by all the escape sequences.
sub wrap {
    my $self = shift;
    local $_ = shift;
    my $output = '';
    my $spaces = ' ' x $$self{MARGIN};
    my $width = $$self{width} - $$self{MARGIN};
    while (length > $width) {
        if (s/^((?:(?:\e\[[\d;]+m)?[^\n]){0,$width})\s+//
            || s/^((?:(?:\e\[[\d;]+m)?[^\n]){$width})//) {
            $output .= $spaces . $1 . "\n";
        } else {
            last;
        }
    }
    $output .= $spaces . $_;
    $output =~ s/\s+$/\n\n/;
    $output;
}

############################################################################
# Module return value and documentation
############################################################################

1;
__END__

=head1 NAME

CmUtils::Pod::Text::Colour - Convert POD data to formatted colour ASCII text

=head1 SYNOPSIS

    use CmUtils::Pod::Text::Colour;
    my $parser = CmUtils::Pod::Text::Colour->new (sentence => 0, width => 78);

    # Read POD from STDIN and write to STDOUT.
    $parser->parse_from_filehandle;

    # Read POD from file.pod and write to file.txt.
    $parser->parse_from_file ('file.pod', 'file.txt');

=head1 DESCRIPTION

CmUtils::Pod::Text::Colour is a simple subclass of Pod::Text that highlights
output text using ANSI colour escape sequences.  Apart from the colour, it in all
ways functions like Pod::Text.  See L<Pod::Text> for details and available
options.

Term::ANSIColor is used to get colours and therefore must be installed to use
this module.

=head1 SEE ALSO

L<Pod::Text|Pod::Text>, L<Pod::Text::Color|Pod::Text::Color>,
L<Pod::Parser|Pod::Parser>

=head1 AUTHORS

Greg Sands E<lt>g.sands@auckland.ac.nzE<gt>.

Russ Allbery E<lt>rra@stanford.eduE<gt>.

=cut
