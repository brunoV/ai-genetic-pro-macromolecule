package AI::Genetic::Pro::Macromolecule;
use Moose;
use MooseX::Types::Moose qw(Str Bool Int Num ArrayRef CodeRef);
use AI::Genetic::Pro::Macromolecule::Types qw(AIGeneticPro Probability);
use AI::Genetic::Pro;
use Moose::Util::TypeConstraints;
use List::Util 'max';
use Modern::Perl;
use namespace::autoclean;

# AI::Genetic::Pro::Macromolecule - Genetic Algorithms to evolve protein, DNA and RNA sequences

=attr fitness

Accepts a CodeRef that should assign a numeric score to each string
sequence that it's passed to it as an argument. Required.

    sub fitness {
        my $seq = shift;

        # Do something with $seq and return a score
        my $score = f($seq);

        return $score;
    }

=cut

has fitness => (
   is  => 'ro',
   isa => CodeRef,
   required => 1,
);

=method evolve

    $m->evolve($n);

Evolve the sequence population for the specified number of generations.
If $n is 0 or undef, it will evolve undefinitely or C<terminate> returns
true.

=method generation

Returns the current generation number.

=cut

has _ga => (
    is  => 'ro',
    isa => AIGeneticPro,
    handles => [qw(evolve generation)],
    lazy_build => 1,
);

sub _build__ga {
    my $self = shift;

    my $ga = AI::Genetic::Pro->new(

        -type            => 'listvector',
        -population      => $self->population_size,
        -crossover       => $self->crossover,
        -mutation        => $self->mutation,
        -parents         => $self->parents,
        -selection       => $self->selection,
        -strategy        => $self->strategy,
        -cache           => $self->cache,
        -history         => 1,
        -preserve        => $self->preserve,
        -variable_length => $self->variable_length,
        -fitness         => $self->_actual_fitness,
    );

    if (
        !$self->variable_length         and
         $self->_has_initial_population and
         $self->_seq_lengths_are_different
    ) { die "initial population lengths cannot be different when variable_length is set to 0\n"; }

    my $initial_population_size = scalar @{$self->initial_population} // 0;
    if ( $initial_population_size > $self->population_size ) {
        warn "initial_population has more sequences than population_size allows\n"
    }

    if ($self->_has_terminate) { $ga->terminate($self->_actual_terminate) };

    $ga->init(
        [ map { $self->_alphabet } (1 .. $self->length) ]
    );

    $ga->inject([ map { [ split '', $_ ] } @{$self->initial_population} ]);

    return $ga;
}

sub _seq_lengths_are_different {
    # returns true if lengths of the inserted sequences are equal
    my $self = shift;

    my $initial_length = length($self->initial_population->[0]);

    return grep { length $_ != $initial_length } @{$self->initial_population};
}

=method fittest

Returns a hash reference with the desired number of top scoring
sequences. The hash reference has two keys, 'seq' which points to the
sequence string, and 'score' which points to the sequence's score.

    my @top_2 = $m->fittest(2);
    # [
    #     { seq => 'VIKP', score => 10 },
    #     { seq => 'VLKP', score => 9  },
    # ]

When called with no arguments, it returns the top scoring sequence.

    my $fittest = $m->fittest;
    # { seq => 'VIKP', score => 10 }

=cut

sub fittest {
    my ($self, $n) = @_;
    $n //= 1;

    my @fittest;
    my @chromosomes = $self->_ga->getFittest($n, 1);

    foreach my $chrom (@chromosomes) {
        my $seq = $self->_ga->as_string($chrom);
        $seq =~ s/_//g;

        push @fittest, {
            seq   => $seq,
            score => $self->_ga->as_value ($chrom),
        };
    }

    return ( $n == 1 ) ? $fittest[0] : @fittest;
}

=method history

Returns a hash reference with the minimum, maximum and mean score for
each generation.

    my $history = $m->history;
    # {
    #     min  => [ 0, 0, 0, 1, 2, ... ],
    #     max  => [ 1, 2, 2, 3, 4, ... ],
    #     mean => [ 0.2, 0.3, 0.5, 1.5, 3, ... ],
    # }

To access the mean score for the C<$n>-th generation, for instance:

    $m->history->{mean}->[$n - 1];

=cut

sub history {
    my $self = shift;

    my $history = $self->_ga->getHistory;

    return {
        max  => $history->[0],
        mean => $history->[1],
        min  => $history->[2],
    };
}

=method current_stats

Returns a hash reference with the minimum, maximum and mean score fore
the current generation.

    $m->current_stats;
    # { min => 2, max => 10, mean => 3.5 }

=cut

sub current_stats {
    my $self = shift;

    my ($max, $mean, $min) = $self->_ga->getAvgFitness;

    return { max => $max, mean => $mean, min => $min };
}

=method current_population

Returns a list with all the sequences of the current generation and
their scores, in no particular order.

    my @seqs = $m->current_population;
    # [
    #     { seq => 'VIKP', score => 10 },
    #     { seq => 'VLKP', score => 9  },
    #     ...
    # ]

=cut

sub current_population {
    my $self = shift;

    my @population;

    my $chromosomes = $self->_ga->people;
    foreach my $chrom (@$chromosomes) {

        my $seq = $self->_ga->as_string( $chrom );
        $seq =~ s/_//g;

        my $score = $self->_ga->as_value($chrom);

        push @population, { seq => $seq, score => $score };

    }

    return @population;
}

has _alphabet => (
    is  => 'ro',
    isa => ArrayRef,
    lazy_build => 1,
);

our %alphabet_for = (
    protein => [qw(A C D E F G H I K L M N P Q R S T V W Y)],
    dna     => [qw(A C G T)],
    rna     => [qw(A C G U)],
);

sub _build__alphabet {
    my $self = shift;

    return $alphabet_for{ lc($self->type) };
}

=attr length

Manually set the allowed maximum length of the sequences.

This attribute is required unless an initial population is provided. In
that case, C<length> will be set as equal to the length of the longest
sequence provided if it's not explicity specified.

=cut

has length => (
    is  => 'ro',
    isa => Num,
    lazy_build => 1,
);

sub _build_length {
    my $self = shift;

    unless ( $self->_has_initial_population ) {
        die "Either length or initial_population should be defined\n";
    }

    my $max_length = max( map { length } @{$self->initial_population} );

    return $max_length;
}

=attr type

Macromolecule type: protein, dna, or rna. Required.

=cut

has type => (
    is  => 'ro',
    isa => enum([qw(protein Protein dna DNA rna RNA)]),
    required => 1,
);

=attr initial_population

Sequences to add to the initial pool before evolving. Accepts an array
reference of strings.

    my $m = AI::Genetic::Pro::Macromolecule->new(
        initial_population => ['ACGT', 'CAAC', 'GTTT'],
        ...
    );

=cut

has initial_population => (
    is  => 'ro',
    isa => ArrayRef[Str],
    predicate => '_has_initial_population',
);

=attr cache

Accepts a Bool value. When true, score results for each sequence will be
stored, to avoid costly and unnecesary recomputations. Set to 1 by
default.

=cut

has cache => (
   is  => 'ro',
   isa => Bool,
   default => 1,
);

=attr mutation

Mutation rate, a number between 0 and 1. Default is 0.05.

=cut

has mutation => (
   is  => 'ro',
   isa => Probability,
   default => 0.05,
);

=attr crossover

Crossover rate, a number between 0 and 1. Default is 0.95.

=cut

has crossover => (
   is  => 'ro',
   isa => Probability,
   default => 0.95,
);

=attr population_size

Number of sequences per generation. Default is 300.

=cut

has population_size => (
   is  => 'ro',
   isa => Int,
   default => 300,
);

=attr parents

Number of parents sequences in recombinations. Default is 2.

=cut

has parents => (
   is  => 'ro',
   isa => Int,
   default => 2,
);

=attr selection

Defines how sequences are selected to crossover. It expects an ArrayRef:

    selection => [ $type, @params ]

See docs in L<AI::Genetic::Pro> for details on available selection
strategies, parameters, and their meanings. Default is Roulette, in
which at first the best individuals/chromosomes are selected. From this
collection parents are selected with probability poportionaly to its
fitness.

=cut

has selection => (
   is  => 'ro',
   isa => ArrayRef,
   default => sub { ['Roulette'] },
);

=attr strategy

Defines strategy of crossover operation. It expects an ArrayRef:

    strategy => [ $strategy, @params ]

See docs in L<AI::Genetic::Pro> for details on available crossover
strategies, parameters, and their meanings. Default is [ Points, 2 ], in
which parents are crossed at 2 points and the best child is moved to the
next generation.

=cut

has strategy => (
   is  => 'ro',
   isa => ArrayRef,
   default => sub { [ 'Points', 2 ] },
);

=attr preserve

Whether to inject the best sequences for next generation, and if so, how
many. Defaults to 5.

=cut

has preserve => (
   is  => 'ro',
   isa => Int,
   default => '5',
);

=attr terminate

Accepts a CodeRef. It will be applied once at the end of each
generation. If returns true, evolution will stop, disregarding the
generation steps passed to the C<evolve> method.

It accepts an C<AI::Genetic::Pro::Macromolecule> object as argument, and
should return either true or false.

    sub terminate {
        my $m = shift;  # an AI::G::P::Macromolecule object

        my $highest_score = $m->fittest->{score};

        if ( $highest_score > 9000 ) {
            warn "It's over 9000!";
            return 1;
        }
    }

In the above example, evolution will stop the moment the top score in
any generation exceeds the value 9000.

=cut

has terminate => (
   is  => 'ro',
   isa => CodeRef,
   predicate => '_has_terminate',
);

has '_actual_' . $_ => (
    is => 'ro',
    isa => CodeRef,
    lazy_build => 1,
) for qw(fitness terminate);

sub _build__actual_fitness {
    my $self = shift;

    return sub {
        my ($ga, $chromosome) = @_;
        my $seq = $ga->as_string($chromosome);
        $seq =~ s/_//g;

        return $self->fitness->($seq);
    }
}

sub _build__actual_terminate {
    my $self = shift;

    return sub { return $self->terminate };
}

=attr variable_length

Decide whether the sequences can have different lengths. Accepts a Bool
value. Defaults to 1.

=cut

has variable_length => (
    is  => 'ro',
    isa => Bool,
    default => 1,
);

__PACKAGE__->meta->make_immutable;
1;

__END__

=head1 SYNOPSIS

    use AI::Genetic::Pro::Macromolecule;

    my @proteins = ($seq1, $seq2, $seq3, ... );

    my $m = AI::Genetic::Pro::Macromolecule->new(
        type    => 'protein',
        fitness => \&hydrophobicity,
        initial_population => \@proteins,
    );

    sub hydrophobicity {
        my $seq = shift;
        my $score = score($seq)

        return $score;
    }

    $m->evolve(10) # evolve for 10 generations;

    my $most_hydrophobic = $m->fittest->{seq};   # get the best sequence
    my $highest_score    = $m->fittest->{score}; # get top score

    # Want the score stats throughout generations?
    my $history = $m->history;

    my $mean_history = $history->{mean}; # [ mean1, mean2, mean3, ... ]
    my $min_history  = $history->{min};  # [ min1,  min2,  min3,  ... ]
    my $max_history  = $history->{max};  # [ max1,  max2,  max3,  ... ]

=head1 DESCRIPTION

AI::Genetic::Pro::Macromolecule is a wrapper over L<AI::Genetic::Pro>,
aimed at easily evolving protein, DNA or RNA sequences using arbitrary
fitness functions.

Its purpose it to allow optimization of macromolecule sequences using
Genetic Algorithms, with as little set up time as possible.

Standing atop L<AI::Genetic::Pro>, it is reasonably fast and memory
efficient. It is also highly customizable, although I've chosen what I
think are sensible defaults for every parameter, so that you don't have
to worry about them if you don't know what they mean.

=cut


