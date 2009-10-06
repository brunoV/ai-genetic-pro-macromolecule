use Test::More;
use Test::Exception;
use Test::Warn;
use Modern::Perl;

use_ok( 'AI::Genetic::Pro::Macromolecule' );

sub fitness {
    my $seq = shift;

    if ($seq eq 'VIKP') { return 10 }
    else { return 0 }
}

my $gm = AI::Genetic::Pro::Macromolecule->new(
    fitness            => \&fitness,
    initial_population => ['VIKP', 'LEP'],
    type               => 'protein',
);

isa_ok( $gm, 'AI::Genetic::Pro::Macromolecule' );

my %defaults = (
    population_size => 300,              # population
    crossover       => 0.95,             # probab. of crossover
    mutation        => 0.05,             # probab. of mutation
    parents         => 2,                # number  of parents
    selection       => [ 'Roulette' ],   # selection strategy
    strategy        => [ 'Points', 2 ],  # crossover strategy
    cache           => 1,                # cache results
    preserve        => 5,                # remember the bests
    variable_length => 1,                # turn variable length ON
);

# Check for default attrs
foreach my $attr (keys %defaults) {
    is_deeply( $gm->$attr, $defaults{$attr}, "Default for $attr ok" );
}

lives_ok { $gm->evolve(1) };

is( $gm->length, 4, 'Chromosome length was guessed correctly' );

is( $gm->fittest->{ seq }, 'VIKP', 'Asking for the top seq returns the top seq' );
is( $gm->fittest->{ score },   10, 'Asking for the top score returns the top score' );

is( $gm->fittest(2), 2, 'Asking for top two returns a list with 2 things' );

is(
    $gm->current_population,
    $gm->population_size + 2,
    'current_population returns list of correct size'
);

# methods history and current_stats return a hashref with correct keys
foreach my $cmd ("history", "current_stats") {

    my $stats = $gm->$cmd;

    ok( defined $stats->{$_}, "$cmd has $_" ) for qw(max mean min);

}

{
    # Testing for failure upon using wrong alphabet
    $gm = AI::Genetic::Pro::Macromolecule->new(
        type               => 'dna',
        fitness            => sub { 1 },
        initial_population => [qw(VIKP ALEP)],
    );

    dies_ok { $gm->evolve(1) };
}

{
    # Either length or initial_population should be defined
    $gm = AI::Genetic::Pro::Macromolecule->new(
        type    => 'dna',
        fitness => sub { 1 },
    );

    dies_ok { $gm->evolve(1) };
}

{
    # terminate works
    $gm = AI::Genetic::Pro::Macromolecule->new(
        type => 'dna',
        fitness   => sub { 1 },
        terminate => sub { (shift)->fittest->{seq} eq 'ACGT' },
        initial_population => [qw(ACGT)],
    );

    # Try to evolve ten times...
    $gm->evolve(10);

    # But since the terminate function tells us that we have what we
    # want, it doesn't evolve any generations
    is( $gm->generation, 0, 'terminate works' );
}

# What happens if the user inputs more sequences than the population
# size?

{
    $gm = AI::Genetic::Pro::Macromolecule->new(
        type => 'dna', fitness => sub { 1 },
        population_size => 2,
        preserve        => 0,
        initial_population => [qw(C A G T)],
    );

    warning_like { $gm->evolve(1) } qr/has more sequences than/;

}

# What happens when we call after-evolution methods before evolving?

{
    $gm = AI::Genetic::Pro::Macromolecule->new(
        type => 'dna', fitness => sub { 1 }, length => 10,
        population_size => 10,
    );

    is( $gm->generation, 0 );

    is( ref $gm->current_stats, 'HASH' );
    is_deeply( [ sort keys %{$gm->current_stats} ], [ qw(max mean min) ] );

    my @pop = $gm->current_population;

    is( @pop, 10 );
    is_deeply( [ sort keys %{$pop[0]} ], [ qw(score seq) ] );

}

done_testing;
