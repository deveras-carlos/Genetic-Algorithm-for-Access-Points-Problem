#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"

void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

unsigned short randgen(  ){
	return ( unsigned short )( rand(  ) % AMT_ACCESS_POINTS );
}

double fitness( unsigned short* genes, int n, Problem* problem ){
    register int i;
    double sum = 0;
    double x_diff, y_diff;
    for ( i = 0; i < n; i++ ){
        x_diff = ( double )( problem->clients[ i ].x - problem->access_points[ genes[ i ] ].x );
        y_diff = ( double )( problem->clients[ i ].y - problem->access_points[ genes[ i ] ].y );
        sum += sqrt( pow( x_diff, 2 ) + pow( y_diff, 2 ) );
    }

    return sum;
}

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( unsigned short*, int n, Problem* problem ), Problem* problem
){
    unsigned int i, j;
    unsigned short gene;
    double fitness, sum = 0;

    for ( i = 0; i < population_size; i++ ){
        population->reading_individuals[ i ].genes = ( unsigned short* ) malloc( individual_size * sizeof( unsigned short ) );
        population->writing_individuals[ i ].genes = ( unsigned short* ) malloc( individual_size * sizeof( unsigned short ) );

        if ( i == 0 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = MAX_POPULATION - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = MAX_POPULATION - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = 2;
        } else if ( i == 1 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = MAX_POPULATION - 1;
            population->writing_individuals[ i ].neighbors[ 1 ] = 0;
            population->writing_individuals[ i ].neighbors[ 2 ] = i + 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = i + 2;
        } else if ( i == MAX_POPULATION - 1 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = i - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = i - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = 0;
            population->writing_individuals[ i ].neighbors[ 3 ] = 1;
        } else if ( i == MAX_POPULATION - 2 ){
            population->writing_individuals[ i ].neighbors[ 0 ] = i - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = i - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = MAX_POPULATION - 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = 0;
        } else {
            population->writing_individuals[ i ].neighbors[ 0 ] = i - 2;
            population->writing_individuals[ i ].neighbors[ 1 ] = i - 1;
            population->writing_individuals[ i ].neighbors[ 2 ] = i + 1;
            population->writing_individuals[ i ].neighbors[ 3 ] = i + 2;
        }

        for ( j = 0; j < AMT_NEIGHBORS; j++ )
            population->reading_individuals[ i ].neighbors[ j ] = population->writing_individuals[ i ].neighbors[ j ];

        for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
            population->reading_individuals[ i ].amt_allocated_clients[ j ] = 0;
            population->writing_individuals[ i ].amt_allocated_clients[ j ] = 0;
        }

        for ( j = 0; j < individual_size; j++ ){
            gene = randgen(  );
            population->reading_individuals[ i ].genes[ j ] = gene;
            population->writing_individuals[ i ].genes[ j ] = gene;

            population->reading_individuals[ i ].amt_allocated_clients[ gene ] += 1;
            population->writing_individuals[ i ].amt_allocated_clients[ gene ] += 1;
        }

        fitness = fitness_function( population->reading_individuals[ i ].genes, individual_size, problem );
        population->writing_individuals[ i ].fitness = population->reading_individuals[ i ].fitness = fitness;

        if ( i == 0 || i == MAX_POPULATION / 4 || i == MAX_POPULATION / 2 || i == MAX_POPULATION - MAX_POPULATION / 4 ){
            population->reading_individuals[ i ].random = TRUE;
            population->writing_individuals[ i ].random = TRUE;
        } else {
            population->reading_individuals[ i ].random = FALSE;
            population->writing_individuals[ i ].random = FALSE;
        }

        population->reading_individuals[ i ].sels = 0;
        population->writing_individuals[ i ].sels = 0;

        sum += fitness;
    }

    population->population_size = population_size;
    population->individual_size = individual_size;
    population->sum_fitness = sum;
    population->best = 0;
    population->worst = 0;
    population->mutation_amount = 0;
    population->equals = 0;
}

int fix_unfeasible( unsigned short* xr, unsigned int* amt_allocated_clients, Problem* problem ){
    if ( amt_allocated_clients[ *xr ] < problem->access_points_capacity[ *xr ] )
        return 1;
    
    unsigned short current_ap = *xr;

    for (int i = 0; i < AMT_ACCESS_POINTS; i++) {
        if (i == current_ap) {
            continue;
        }

        if (amt_allocated_clients[i] < problem->access_points_capacity[i]) {
            *xr = i;

            amt_allocated_clients[i]++;

            amt_allocated_clients[current_ap]--;

            return 1;
        }
    }
    return 0;
}

void crossover( Population* population, int father, int mother, int son, float alpha, Problem* problem ){
    int aux, i;
    int point1 = rand() % population->individual_size;
    int point2 = rand() % population->individual_size;
    if (point1 > point2){
        aux = point1;
        point1 = point2;
        point2 = aux;
    }

    for ( i = point1; i < point2; i++ ) {
        population->writing_individuals[ son ].genes[i] = population->reading_individuals[ father ].genes[i];
    }


    for ( i = 0; i < point1; i++ ) {
        population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ mother ].genes[i];
    }

    for ( i = point2; i < population->individual_size; i++ ) {
        population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ mother ].genes[i];
    }
}

void mutate(Chromossome* individual, int individual_size, float mutation_rate, Problem* problem) {
    for (int i = 0; i < individual_size; i++) {
        if ( (float)rand() / RAND_MAX < mutation_rate) {
            unsigned short original_gene = individual->genes[i];
            individual->genes[i] = randgen();
            fix_unfeasible(&(individual->genes[i]), individual->amt_allocated_clients, problem);

            if (individual->genes[i] != original_gene) {
                individual->amt_allocated_clients[original_gene]--;
                individual->amt_allocated_clients[individual->genes[i]]++;
            }
        }
    }
}

void downhill_local_search(
    Chromossome* individual, int individual_size, int step,
    double ( *fitness_function )( unsigned short*, int n, Problem* problem ), Problem* problem
){
    int i;
    double step_size = 1 + ( double ) ( ( rand(  ) % 101 / 100 ) * step );
    double new_fitness, original_fitness = individual->fitness;
    double original_gene;

    for (i = 0; i < individual_size; i++) {
        original_gene = ( double )individual->genes[i];

        individual->genes[i] = ( ( unsigned short )( original_gene + step_size ) ) % AMT_ACCESS_POINTS;
        fix_unfeasible(&(individual->genes[i]), individual->amt_allocated_clients, problem);
        new_fitness = fitness_function(individual->genes, individual_size, problem);

        if (new_fitness < original_fitness) {
            original_fitness = new_fitness;
        } else {
            individual->genes[i] = ( ( unsigned short )( original_gene - step_size ) ) % AMT_ACCESS_POINTS;
            fix_unfeasible(&(individual->genes[i]), individual->amt_allocated_clients, problem);
            new_fitness = fitness_function(individual->genes, individual_size, problem);

            if (new_fitness < original_fitness) {
                original_fitness = new_fitness;
            } else {
                individual->genes[i] = ( unsigned short ) original_gene;
            }
        }
    }

    individual->fitness = original_fitness;
}

void genetic_algorithm( Problem* problem ){
    clock_t start, end;

    Population population;
    int pa1, pa2;
    double fit, dvp, med;
    int cfo = 0, generation = 0;
    int i, j;

    unsigned short gene;
    unsigned short *auxiliar_genes;
    unsigned int* auxiliar_allocated_clients;

    int num_threads = 8;

    srand((unsigned) time(0));

    start_population( &population, ( unsigned int ) MAX_POPULATION, ( unsigned int ) MAX_GENES, fitness, problem );

    for ( int i = 0; i < MAX_POPULATION; i++ ){
        printf( "Individual[ %d ] = {\n", i );
        printf( "\tFitness: %lf\n", population.reading_individuals[ i ].fitness );
        printf( "\tGenes : [ " );
        for ( int j = 0; j < MAX_GENES - 1; j++ )
            printf( "%d, ", population.reading_individuals[ i ].genes[ j ]);
        printf( "%d ]\n", population.reading_individuals[ i ].genes[ MAX_GENES - 1 ] );
        printf( "\tVizinhos: [ %d, %d, %d, %d ]\n",
            population.reading_individuals[ i ].neighbors[ 0 ], population.reading_individuals[ i ].neighbors[ 1 ],
            population.reading_individuals[ i ].neighbors[ 2], population.reading_individuals[ i ].neighbors[ 3 ]
        );
        printf( "}\n" );
    }

    for ( i = 0; i < population.population_size; i++ ){
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    printf( "\n\tAverage distance: %lf\n", ( population.sum_fitness / ( MAX_POPULATION * AMT_CLIENTS ) ) );
    
    start = clock();

    do {
        #pragma omp parallel for num_threads( num_threads ) \
        private( pa1, pa2, fit, gene, i, j ) shared( population, problem ) reduction( +:cfo )
        for ( i = 0; i < population.population_size; i++ ){
            if ( population.writing_individuals[ i ].random == TRUE ){
                for ( j = 0; j < AMT_ACCESS_POINTS; j++ )
                    population.writing_individuals[ i ].amt_allocated_clients[ j ] = 0;
                
                for ( j = 0; j < population.individual_size; j++ ){
                    gene = (unsigned short) randgen(  );
                    population.writing_individuals[ i ].genes[ j ] = gene;
                    population.writing_individuals[ i ].amt_allocated_clients[ gene ] += 1;
                }

                fit = fitness( population.writing_individuals[ i ].genes, population.individual_size, problem );
                population.writing_individuals[ i ].fitness = fit;
            } else {
                pa1 = population.reading_individuals[ i ].neighbors[ 0 + (rand(  ) % AMT_NEIGHBORS ) ];
                pa2 = population.reading_individuals[ i ].neighbors[ 0 + (rand(  ) % AMT_NEIGHBORS ) ];

                if ( pa1 != pa2 ){
                    crossover( &population, pa1, pa2, i, XALPHA, problem );
                    for ( j = 0; j < population.individual_size; j++ ){
                        fix_unfeasible( &( population.writing_individuals[ i ].genes[ j ] ), population.writing_individuals[ i ].amt_allocated_clients, problem );
                    }
                    population.writing_individuals[ i ].fitness = fitness( population.writing_individuals[ i ].genes, population.individual_size, problem);
                    if ( population.writing_individuals[ i ].fitness < population.reading_individuals[ i ].fitness )
                        population.writing_individuals[ i ].sels++;
                    
                    cfo++;
                }
            }
        }

        #pragma omp parallel for num_threads( num_threads ) shared(population, problem)
        for (int i = 0; i < population.population_size; i++) {
            mutate(&population.writing_individuals[i], population.individual_size, MUTATION_RATE, problem);
        }

        #pragma omp parallel for num_threads( num_threads ) \
        private( auxiliar_genes, auxiliar_allocated_clients, i, j ) shared( population, problem )
        for ( i = 0; i < population.population_size; i++ ){
            if (population.writing_individuals[ i ].sels > population.reading_individuals[ i ].sels ){
                population.reading_individuals[ i ].fitness = population.writing_individuals[ i ].fitness;

                auxiliar_genes = population.reading_individuals[ i ].genes;
                population.reading_individuals[ i ].genes = population.writing_individuals[ i ].genes;
                population.writing_individuals[ i ].genes = auxiliar_genes;

                for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
                    population.reading_individuals[ i ].amt_allocated_clients[ j ] = population.writing_individuals[ i ].amt_allocated_clients[ j ];
                }

                population.reading_individuals[ i ].sels++;
            }
        }
        // printf( "\nFez as trocas lá\n" );

        #pragma omp parallel for num_threads(num_threads) private(i, j) shared( population, problem )
        for (i = 0; i < population.population_size; i++) {
            for ( j = 0; j < MAX_ITER_LOCAL_SEARCH; j++ )
                downhill_local_search(
                    &(population.reading_individuals[i]), 
                    population.individual_size,
                    generation,
                    fitness, 
                    problem
                );
        }

        generation++;
        printf( "CFO: %d\n", cfo );
    } while ( cfo < MAX_GEN );

    end = clock();

    printf("\n\tExecutado em tempo = %.4f,CFO=%d ", (double) ( (end - start) / (CLOCKS_PER_SEC*num_threads) ), cfo);

    population.sum_fitness = 0.0F;

    for ( i = 0; i < population.population_size; i++ ){
        population.sum_fitness += population.reading_individuals[ i ].fitness;
        printf( "\n\tIndividual (%d) = %.10f | %d", i, population.reading_individuals[ i ].fitness, population.reading_individuals[ i ].sels );
        if ( population.reading_individuals[ i ].fitness < population.reading_individuals[ population.best ].fitness )
            population.best = i;
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    printf( "\n\tAverage distance: %lf\n", ( population.sum_fitness / ( double )MAX_POPULATION ) / ( double )AMT_CLIENTS );
}

#endif
