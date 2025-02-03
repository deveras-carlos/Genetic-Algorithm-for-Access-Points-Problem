#ifndef GENETIC_ALGORITHM_C
#define GENETIC_ALGORITHM_C
#include "genetic_algorithm.h"

unsigned short randgen(  ){
	return ( unsigned short )( rand(  ) % AMT_ACCESS_POINTS );
}

double distance(Point a, Point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double fitness( Chromossome* individual, Problem* problem ){
    unsigned short* genes = individual->genes;
    register int i, j;
    double sum = 0;

    // for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
    //     individual->amt_allocated_clients[ j ] = 0;
    // }

    for ( i = 0; i < AMT_CLIENTS; i++ ){
        sum += distance( problem->clients[ genes[ i ] ], problem->access_points[ genes[ i ] ] );
        // individual->amt_allocated_clients[ genes[ i ] ]++;
    }

    for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
        if ( individual->amt_allocated_clients[ j ] > problem->access_points_capacity[ j ] )
            sum += 1e9;
    }

    return sum;
}

void fix_unfeasible( unsigned short gene, Chromossome* individual, Problem* problem) {
    if ( individual->amt_allocated_clients[ gene ] < problem->access_points_capacity[ gene ] )
        return;
    
    unsigned short current_ap = individual->genes[ gene ];

    for (int i = 0; i < AMT_ACCESS_POINTS; i++) {
        if (i == current_ap) {
            continue;
        }

        if (individual->amt_allocated_clients[i] < problem->access_points_capacity[i]) {
            individual->genes[ gene ] = i;
            individual->amt_allocated_clients[i]++;
            individual->amt_allocated_clients[current_ap]--;

            return;
        }
    }
    return;
}

void start_population( Population* population, Problem* problem ){
    unsigned int i, j;
    unsigned short gene;
    double fit, sum = 0;

    for ( i = 0; i < MAX_POPULATION; i++ ){
        population->reading_individuals[ i ].genes = ( unsigned short* ) malloc( AMT_CLIENTS * sizeof( unsigned short ) );
        population->writing_individuals[ i ].genes = ( unsigned short* ) malloc( AMT_CLIENTS * sizeof( unsigned short ) );

        for ( j = 0; j < AMT_ACCESS_POINTS; j++ ){
            population->reading_individuals[ i ].amt_allocated_clients[ j ] = 0;
            population->writing_individuals[ i ].amt_allocated_clients[ j ] = 0;
        }

        for ( j = 0; j < AMT_CLIENTS; j++ ){
            gene = randgen(  );
            population->reading_individuals[ i ].genes[ j ] = gene;
            population->writing_individuals[ i ].genes[ j ] = gene;

            population->reading_individuals[ i ].amt_allocated_clients[ gene ]++;
            population->writing_individuals[ i ].amt_allocated_clients[ gene ]++;

            fix_unfeasible( j, &population->reading_individuals[ i ], problem );
        }

        fit = fitness( &population->reading_individuals[ i ], problem );
        population->writing_individuals[ i ].fitness = population->reading_individuals[ i ].fitness = fit;

        if ( i == 0 || i == MAX_POPULATION / 4 || i == MAX_POPULATION / 2 || i == MAX_POPULATION - MAX_POPULATION / 4 ){
            population->reading_individuals[ i ].random = TRUE;
            population->writing_individuals[ i ].random = TRUE;
        } else {
            population->reading_individuals[ i ].random = FALSE;
            population->writing_individuals[ i ].random = FALSE;
        }

        population->reading_individuals[ i ].sels = 0;
        population->writing_individuals[ i ].sels = 0;

        sum += fit;
    }

    population->population_size = MAX_POPULATION;
    population->sum_fitness = sum;
    population->best = 0;
    population->worst = 0;
    population->mutation_amount = 0;
    population->equals = 0;
}

int tournament_selection(Population* population) {
    int best = rand() % MAX_POPULATION;
    for (int i = 1; i < TOURNAMENT_SIZE; i++) {
        int challenger = rand() % MAX_POPULATION;
        if ( population->reading_individuals[ i ].fitness  < population->reading_individuals[ best ].fitness) {
            best = challenger;
        }
    }
    return best;
}

void crossover( Population* population, int father, int mother, int son, float alpha, Problem* problem ){
    int aux, i;
    int point1 = rand() % AMT_CLIENTS;
    int point2 = rand() % AMT_CLIENTS;
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

    for ( i = point2; i < AMT_CLIENTS; i++ ) {
        population->writing_individuals[ son ].genes[ i ] = population->reading_individuals[ mother ].genes[i];
    }
}

void mutate(Chromossome* individual, Problem* problem) {
    for (int i = 0; i < AMT_CLIENTS; i++) {
        if ( ( rand() / (double) RAND_MAX ) < MUTATION_RATE ) {
            unsigned short original_gene = individual->genes[i];
            individual->genes[i] = randgen();
            individual->amt_allocated_clients[ original_gene ]--;
            individual->amt_allocated_clients[ individual->genes[ i ] ]++;
            fix_unfeasible( i, individual, problem );
        }
    }
}

void downhill_local_search( Chromossome* individual, Problem* problem ){
    for (int i = 0; i < AMT_CLIENTS; i++) {
        int best_ap = individual->genes[ i ];
        double best_dist = distance( problem->clients[ i ], problem->access_points[ individual->genes[ i ] ] );
        for ( int j = 0; j < AMT_ACCESS_POINTS; j++ ){
            if ( j != best_ap ){
                if ( individual->amt_allocated_clients[ j ] < problem->access_points_capacity[ j ] ){
                    double new_dist = distance( problem->clients[ i ], problem->access_points[ j ] );
                    if ( new_dist < best_dist ){
                        individual->genes[ i ] = j;
                        best_dist = new_dist;
                    }
                }
            }
        }
        if ( individual->genes[ i ] != best_ap ){
            individual->amt_allocated_clients[ best_ap ]--;
            individual->amt_allocated_clients[ individual->genes[ i ] ]++;
        }
    }
    individual->fitness = fitness( individual, problem );
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

    int num_threads = 40;

    srand( time( NULL) );

    start_population( &population, problem );

    for ( int i = 0; i < MAX_POPULATION; i++ ){
        printf( "Individual[ %d ] = {\n", i );
        printf( "\tFitness: %lf\n", population.reading_individuals[ i ].fitness );
        printf( "\tGenes : [ " );
        for ( int j = 0; j < AMT_CLIENTS - 1; j++ )
            printf( "%d, ", population.reading_individuals[ i ].genes[ j ]);
        printf( "%d ]\n", population.reading_individuals[ i ].genes[ AMT_CLIENTS - 1 ] );
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
                for ( j = 0; j < AMT_CLIENTS; j++ ){
                    gene = (unsigned short) randgen(  );
                    population.writing_individuals[ i ].genes[ j ] = gene;
                }
                population.writing_individuals[ i ].fitness = fit;
            } else {
                pa1 = tournament_selection( &population );
                pa2 = tournament_selection( &population );

                if ( pa1 != pa2 ){
                    crossover( &population, pa1, pa2, i, XALPHA, problem );

                    for ( j = 0; j < AMT_CLIENTS; j++ ){
                        fix_unfeasible( j, &population.reading_individuals[ i ], problem );
                    }

                    mutate(&population.writing_individuals[i], problem);
                    population.writing_individuals[ i ].fitness = fitness( &population.writing_individuals[ i ], problem);

                    if ( population.writing_individuals[ i ].fitness < population.reading_individuals[ i ].fitness )
                        population.writing_individuals[ i ].sels++;
                    
                    cfo++;
                }
            }
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

        #pragma omp parallel for num_threads(num_threads) private(i, j) shared( population, problem )
        for (i = 0; i < population.population_size; i++) {
            for ( j = 0; j < MAX_ITER_LOCAL_SEARCH; j++ )
                downhill_local_search( &population.reading_individuals[i], problem );
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
        
        printf("\n\tTotal clients per access point:\n");
        for (int j = 0; j < AMT_ACCESS_POINTS; j++) {
            printf("\t\tAP %d: %d clients\n", j, population.reading_individuals[ i ].amt_allocated_clients[ j ]);
        }
    }

    printf( "\n\nThe best is %d with fitness equals to %.4f\n", population.best, population.reading_individuals[ population.best ].fitness );

    printf( "\n\tAverage distance: %lf\n", ( population.sum_fitness / ( double )MAX_POPULATION ) / ( double )AMT_CLIENTS );
}

#endif
