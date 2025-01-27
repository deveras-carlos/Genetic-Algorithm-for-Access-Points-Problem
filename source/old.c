#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"

unsigned short randgen(  ){
	return rand(  ) % AMT_ACCESS_POINTS;
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
    int i, j;
    unsigned short gene;
    double fit, sum = 0;

    for ( i = 0; i < population_size; i++ ){
        population->reading_individuals[ i ].genes = ( unsigned short* ) malloc( MAX_GENES * sizeof( unsigned short ) );
        population->writing_individuals[ i ].genes = ( unsigned short* ) malloc( MAX_GENES * sizeof( unsigned short ) );

        population->reading_individuals[ i ].neighbors[ 0 ] = i % COLS == 0 ? i + ( COLS - 1 ) : i - 1;         // Left
        population->reading_individuals[ i ].neighbors[ 1 ] = i % COLS == COLS - 1 ? i - ( COLS - 1 ) : i + 1;  // Right

        population->reading_individuals[ i ].neighbors[ 2 ] = ( i / COLS ) > ROWS - 3 ? i - COLS * 2 : i + COLS * 2;         // Top or Bottom
        population->reading_individuals[ i ].neighbors[ 3 ] = ( i / COLS ) == ROWS - 2 || ( i / COLS ) == 0 ? i + COLS : i - COLS * 2;  // Bottom or Top

        // Correction if it's on line 1 or last line
        population->reading_individuals[ i ].neighbors[ 3 ] += ( i / COLS ) == 1 || ( i / COLS ) > ROWS - 2 ? COLS : 0;
        population->reading_individuals[ i ].neighbors[ 2 ] += ( i / COLS ) == 1 || ( i / COLS ) > ROWS - 2 ? COLS : 0;    
        // printf( "%d : %d\n\n", i, ( i / COLS ) );

        for ( j = 0; j < 4; j++ )
            population->writing_individuals[ i ].neighbors[ j ] = population->reading_individuals[ i ].neighbors[ j ];

        for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
            population->reading_individuals[ i ].amt_allocated_clients[ j ] = 0;
            population->writing_individuals[ i ].amt_allocated_clients[ j ] = 0;
        }

        for ( j = 0; j < individual_size; j++ ){
            gene = randgen(  );
            population->reading_individuals[ i ].genes[ j ] = gene;
            population->writing_individuals[ i ].genes[ j ] = gene;

            population->reading_individuals[ i ].amt_allocated_clients[ j ] += 1;
            population->writing_individuals[ i ].amt_allocated_clients[ j ] += 1;
        }

        fit = fitness_function( population->reading_individuals[ i ].genes, individual_size, problem );
        population->writing_individuals[ i ].fitness = population->reading_individuals[ i ].fitness = fit;

        if ( i == 0 || i == COLS / 4 || i == COLS / 2 || i == COLS - COLS / 4 ){
            population->reading_individuals[ i ].random = TRUE;
            population->writing_individuals[ i ].random = TRUE;
        } else {
            population->reading_individuals[ i ].random = FALSE;
            population->writing_individuals[ i ].random = FALSE;
        }

        population->reading_individuals[ i ].sels = 0;
        population->writing_individuals[ i ].sels = 0;

        sum += ( fit ); // fabs
    }
    
    population->population_size = population_size;
    population->individual_size = individual_size;
    population->sum_fitness = ( sum / ( double )population_size ) / ( double )individual_size ;
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

void blend_crossover( Population* population, int father, int mother, int son, float alpha, Problem* problem ){
    double a, b, f, m, s;
    int i;
    int r;

    for ( i = 0; i < AMT_ACCESS_POINTS; i++ ){
        population->writing_individuals[ son ].amt_allocated_clients[ i ] = 0;
    }

    for ( i = 0; i < population->individual_size; i++ ){
        r = rand(  ) % 3;

        printf( "Gene %d, son: %d - father: %d - mother: %d\n", i, son, father, mother );
        switch ( r ){
            case 0:
                population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ father ].genes[ i ];
                break;
            case 1:
                population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ mother ].genes[ i ];
                break;
            default:
                a = -alpha;
                b = 1 + alpha;
                f = a + ( rand(  ) / 101 / 100 ) * ( b - a );
                population->writing_individuals[ son ].genes[ i ] = ( ( unsigned short )(
                    ( double )population->reading_individuals[ father ].genes[ i ] +
                    r * ( double )( population->reading_individuals[ mother ].genes[ i ] - population->reading_individuals[ father ].genes[ i ] )
                ) ) % AMT_ACCESS_POINTS;
                break;
        }

        // printf( "Access point: %d\n", population->writing_individuals[ son ].genes[ i ] );

        population->writing_individuals[ son ].amt_allocated_clients[ population->writing_individuals[ son ].genes[ i ] ] += 1;
        fix_unfeasible( &( population->writing_individuals[ son ].genes[ i ] ), population->writing_individuals[ son ].amt_allocated_clients, problem );
    }
}

void downhill_local_search(
    Chromossome* individual, int individual_size, int step,
    double ( *fitness_function )( unsigned short*, int n, Problem* problem ), Problem* problem
){
    int i;
    double step_size = 1 - 0.0001 * ( double ) step;
    double new_fitness, original_fitness = individual->fitness;

    for (i = 0; i < individual_size; i++) {
        double original_gene = ( double )individual->genes[i];

        individual->genes[i] = ( unsigned short )( original_gene + step_size ) % AMT_ACCESS_POINTS;
        fix_unfeasible(&(individual->genes[i]), individual->amt_allocated_clients, problem);
        new_fitness = fitness_function(individual->genes, individual_size, problem);

        if (new_fitness < original_fitness) {
            original_fitness = new_fitness;
        } else {
            individual->genes[i] = ( unsigned short )( original_gene - step_size ) % AMT_ACCESS_POINTS;
            fix_unfeasible(&(individual->genes[i]), individual->amt_allocated_clients, problem);
            new_fitness = fitness_function(individual->genes, individual_size, problem);

            if (new_fitness < original_fitness) {
                original_fitness = new_fitness;
            } else {
                individual->genes[i] = original_gene;
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

    int num_threads = 4;

    srand((unsigned) time(0));

    start_population( &population, MAX_POPULATION, MAX_GENES, fitness, problem );

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
            if ( population.reading_individuals[ i ].random == TRUE ){
                for ( j = 0; j < AMT_ACCESS_POINTS; j++ )
                    population.writing_individuals->amt_allocated_clients[ j ] = 0;
                
                for ( j = 0; j < population.individual_size; j++ ){
                    gene = (unsigned short) randgen(  );
                    printf( "Alelo %d, indivíduo %d - gene: %d", j, i, gene );
                    population.writing_individuals[ i ].genes[ j ] = gene;
                    population.writing_individuals->amt_allocated_clients[ j ] += 1;
                }

                fit = fitness( population.writing_individuals[ i ].genes, population.individual_size, problem );
                population.writing_individuals[ i ].fitness = fit;
                // printf( "\nGerou indivíduos aleatórios.\n" );
            } else {
                pa1 = population.reading_individuals[ i ].neighbors[ 0 + (rand(  ) % AMT_NEIGHBORS ) ];
                pa2 = population.reading_individuals[ i ].neighbors[ 0 + (rand(  ) % AMT_NEIGHBORS ) ];

                if ( pa1 != pa2 ){
                    blend_crossover( &population, pa1, pa2, i, XALPHA, problem );
                    population.writing_individuals[ i ].fitness = fitness( population.writing_individuals[ i ].genes, population.individual_size, problem);
                    if ( population.writing_individuals[ i ].fitness < population.reading_individuals[ i ].fitness )
                        population.reading_individuals[ i ].sels++;
                    
                    cfo++;
                }
                // printf( "\nTentou cruzar.\n" );
            }
        }
        // printf( "\nFez cruzamento\n" );

        #pragma omp parallel for num_threads( num_threads ) \
        private( auxiliar_genes, auxiliar_allocated_clients, i, j ) shared( population, problem )
        for ( i = 0; i < population.population_size; i++ ){
            if (population.reading_individuals[ i ].sels > population.writing_individuals[ i ].sels ){
                population.reading_individuals[ i ].fitness = population.writing_individuals[ i ].fitness;

                auxiliar_genes = population.reading_individuals[ i ].genes;
                population.reading_individuals[ i ].genes = population.writing_individuals[ i ].genes;
                population.writing_individuals[ i ].genes = auxiliar_genes;

                for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
                    population.reading_individuals[ i ].amt_allocated_clients[ j ] = population.writing_individuals[ i ].amt_allocated_clients[ j ];
                }

                population.writing_individuals[ i ].sels++;
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
        // printf( "\nFez a busca local\n" );

        generation++;
        printf( "Evaluated individuals %d\n", cfo );
        // printf( "Generation number %d\n", generation );
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
