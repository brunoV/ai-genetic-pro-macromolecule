package AI::Genetic::Pro::Macromolecule::Types;

# ABSTRACT: Specific types for AI::Genetic::Pro::Macromolecule

use Moose;
use MooseX::Types::Moose qw(Str);
use MooseX::Types -declare => [qw(Probability AIGeneticPro)];
use namespace::autoclean;

class_type AIGeneticPro, { class => 'AI::Genetic::Pro' };

subtype Probability, as Str, where { $_ < 1 and $_ > 0 };


__PACKAGE__->meta->make_immutable;

__END__

=head1 DESCRIPTION

This module defines specific types and type coercions to be used by
AI::Genetic::Pro::Macromolecule.

=cut
