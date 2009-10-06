package AI::Genetic::Pro::Macromolecule::Types;
use Moose;
use MooseX::Types::Moose qw(Str);
use MooseX::Types -declare => [qw(Probability AIGeneticPro)];
use AI::Genetic::Pro;
use namespace::autoclean;

class_type AIGeneticPro, { class => 'AI::Genetic::Pro' };

subtype Probability, as Str, where { $_ < 1 and $_ > 0 };

1;
