#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#define FALSE 0
#define TRUE 1

#define MAX_GENES 350
#define ROWS 10
#define COLS 4
#define MAX_POPULATION ROWS * COLS
#define AMT_NEIGHBORS 4
#define MAX_GEN 10000

#define XALPHA  0.35

#define TOURNAMENT_SIZE 3
#define MAX_ITER_LOCAL_SEARCH 1

#define AMT_ACCESS_POINTS 4
#define AMT_CLIENTS 350

typedef struct _point_{
    unsigned int x;
    unsigned int y;
} Point;

typedef struct _problem_{
    Point access_points[ AMT_ACCESS_POINTS ];
    unsigned int access_points_capacity[ AMT_ACCESS_POINTS ];
    Point clients[ AMT_CLIENTS ];
} Problem;

#ifdef GENETIC_ALGORITHM_C
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>

typedef struct _chromossome_ {
    unsigned short* genes;
    unsigned int amt_allocated_clients[ AMT_ACCESS_POINTS ];
    unsigned int neighbors[ AMT_NEIGHBORS ];
    double fitness;
    char random; // 1 if it's a random chromossome
    char sels;
} Chromossome;

typedef struct _population_ {
        Chromossome     reading_individuals[ MAX_POPULATION ];
        Chromossome     writing_individuals[ MAX_POPULATION ];
        double          sum_fitness;
        int             population_size;
        int             individual_size;
        int             best;
        int             worst;
        int             mutation_amount;
        int             equals;
        int             best_generation;
} Population;

void start_population(
    Population* population, unsigned int population_size, unsigned int individual_size,
    double ( *fitness_function )( unsigned short*, int n, Problem* problem ), Problem* problem
);

int fix_unfeasible( unsigned short* xr, unsigned int* amt_allocated_clients, Problem* problem );

void blend_crossover( Population* population, int father, int mother, int son, float alpha, Problem* problem );

void genetic_algorithm( Problem* problem );

#else

extern void genetic_algorithm( Problem* problem );

#endif

#endif