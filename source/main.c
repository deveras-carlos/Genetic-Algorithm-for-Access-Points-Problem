#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "genetic_algorithm.h"

Problem* read_csv(const char* filename) {
    // Allocate memory for the Problem structure
    Problem* problem = (Problem*)malloc(sizeof(Problem));
    if (!problem) {
        perror("Failed to allocate memory for Problem");
        return NULL;
    }

    // Initialize access points and their capacities
    problem->access_points[0] = (Point){0, 0};
    problem->access_points[1] = (Point){80, 0};
    problem->access_points[2] = (Point){0, 80};
    problem->access_points[3] = (Point){80, 80};

    problem->access_points_capacity[0] = 64;
    problem->access_points_capacity[1] = 64;
    problem->access_points_capacity[2] = 128;
    problem->access_points_capacity[3] = 128;

    // Open the CSV file for reading
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        free(problem);
        return NULL;
    }

    // Read the file line by line and parse client data
    char line[128];
    int client_index = 0;

    // Skip the header line
    fgets(line, sizeof(line), file);

    while (fgets(line, sizeof(line), file)) {
        if (client_index >= AMT_CLIENTS) {
            fprintf(stderr, "Exceeded maximum number of clients\n");
            break;
        }

        // Parse the CSV line
        unsigned int id, x, y;
        if (sscanf(line, "%u;%u;%u", &id, &x, &y) == 3) {
            problem->clients[client_index] = (Point){x, y};
            client_index++;
        } else {
            fprintf(stderr, "Failed to parse line: %s\n", line);
        }
    }

    fclose(file);

    // Check if the number of clients matches the expected value
    if (client_index != AMT_CLIENTS) {
        fprintf(stderr, "Warning: Number of clients read (%d) does not match expected (%d)\n", client_index, AMT_CLIENTS);
    }

    return problem;
}

Problem* generate_random_problem() {
    // Allocate memory for the Problem structure
    Problem* problem = (Problem*)malloc(sizeof(Problem));
    if (!problem) {
        perror("Failed to allocate memory for Problem");
        return NULL;
    }

    // Initialize access points and their capacities
    problem->access_points[0] = (Point){0, 0};
    problem->access_points[1] = (Point){80, 0};
    problem->access_points[2] = (Point){0, 80};
    problem->access_points[3] = (Point){80, 80};

    problem->access_points_capacity[0] = 64;
    problem->access_points_capacity[1] = 64;
    problem->access_points_capacity[2] = 128;
    problem->access_points_capacity[3] = 128;

    // Seed the random number generator
    srand(time(NULL));

    // Generate random client data
    for (int i = 0; i < AMT_CLIENTS; i++) {
        problem->clients[i].x = rand() % 81; // Generate x in range [0, 80]
        problem->clients[i].y = rand() % 81; // Generate y in range [0, 80]
    }

    return problem;
}

int main( int argc, char* argv[  ] ){
    const char* filename = "ag_data.csv";
    Problem* problem = read_csv(filename);
    // Problem* problem = generate_random_problem(  );

    if (problem) {
        printf("Access Points:\n");
        for (int i = 0; i < AMT_ACCESS_POINTS; i++) {
            printf("AP%d: (%u, %u) Capacity: %u\n",
                   i + 1, problem->access_points[i].x, problem->access_points[i].y,
                   problem->access_points_capacity[i]);
        }

        printf("\nClients:\n");
        for (int i = 0; i < AMT_CLIENTS; i++) {
            printf("Client %d: (%u, %u)\n", i + 1, problem->clients[i].x, problem->clients[i].y);
        }

        genetic_algorithm( problem );

        // Free the allocated memory
        free(problem);
    }


    return 0;
}