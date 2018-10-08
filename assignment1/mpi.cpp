#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "functions.h"

using namespace std;

/**
 *
 * GLOBAL VAR
 *
 */
points_container **point_containers;
int **final_paths;
double *results;
int *grid_degrees;

bool degreeAllowed(int g_i, int count) {
    int degree = grid_degrees[g_i];

    for (int i = 0; i < count; i++) {
        if (grid_degrees[i] < degree) {
            return false;
        }
    }

    return true;
}

void updateDegree(int g_i) {
    grid_degrees[g_i] += 2;
}

int degreesLeft(int count) {
    int counter = 0;

    for (int i = 0; i < count; i++) {
//        printf("%i ", grid_degrees[i]);
        if (grid_degrees[i] < 2) {
            counter += 1;
        }
    }

    //    printf("   left=%i\n", counter);
    return counter;
}

bool degreesAllEqual2(int count) {
    return degreesLeft(count) == 0;
}

void add_sub_graph(points_container **full_path, int *full_path_index, int g_i, int previous_nn_p, int *p_to_search,
                   int p_i) {
    if (point_containers[g_i]->count > 1) {
        // printf("{%i} {%i}", p_to_search[1], p_to_search[0]);
        int index_of_initial = find_index_in_path_for_point(previous_nn_p, point_containers[g_i],
                                                            final_paths[g_i]);

        int value_modifier = p_to_search[p_i] == final_paths[g_i][index_of_initial - 1] ? 1 : -1;

        int count_for_this_grid = point_containers[g_i]->count;
        for (int i = 1; i < count_for_this_grid; i++) {
            int calculated_path_index = (
                    ((i * value_modifier) + index_of_initial + count_for_this_grid) %
                    count_for_this_grid);
            int calculated_point_index = final_paths[g_i][calculated_path_index];
            (*full_path)->points[(*full_path_index)++] = point_containers[g_i]->points[
                    calculated_point_index - 1];
        }
    }
}

/**
 *
 * Traveling Salesman Algorithm
 *
 */

#define INF DBL_MAX;

/**
 *
 * MAIN
 *
 */
int main(int argc, char *argv[]) {
    /**
     *
     * MPI Code
     *
     */

    int rank, total_tasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &total_tasks); // This will need to be a power of 4 right now.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    int number_of_points;
    int tag = 1;
    int path_tag = 2;
    points_container *current_process_point_container;
    int *current_processes_final_path;

    points_container *points;

    struct timespec start, end;

    if (argc != 2) {
        cout << "MPI-Usage: ./main tmp.txt\n" << endl;
        exit(0);
    }

    int NUM_THREADS = total_tasks;

    if (rank == 0) {
        // load data, then broadcast to other nodes it is ready
        char *filename = argv[1];

        clock_gettime(CLOCK_MONOTONIC_RAW, &start);

        // create the point containers that will be used per thread.
        point_containers = (points_container **) malloc(sizeof(points_container *) * NUM_THREADS);
        final_paths = (int **) malloc(NUM_THREADS * sizeof(int *));

        points = get_the_points(filename);


        double min_x = points->points[0].x;
        double max_x = points->points[0].x;
        double min_y = points->points[0].y;
        double max_y = points->points[0].y;
        for (int i = 0; i < points->count; i++) {
            min_x = min(min_x, floor(points->points[i].x));
            max_x = max(max_x, ceil(points->points[i].x));
            min_y = min(min_x, floor(points->points[i].y));
            max_y = max(max_x, ceil(points->points[i].y));
        }
        double range_x = max_x - min_x;
        double block_range_x = ceil((range_x + 1) / sqrt(NUM_THREADS));
        double range_y = max_y - min_y;
        double block_range_y = ceil((range_y + 1) / sqrt(NUM_THREADS));

        // TODO: count how many elements are stored in each grid block so that we can malloc it
        int *block_sizes = (int *) calloc((size_t) NUM_THREADS, sizeof(int));

        // count how many points each block should have
        for (int i = 0; i < points->count; i++) {
            point p = points->points[i];
            int x_block = (int) floor((p.x - min_x) / block_range_x);
            int y_block = (int) floor((p.y - min_y) / block_range_y);

            block_sizes[(x_block * (int) sqrt(NUM_THREADS)) + y_block]++;
        }

        // initialize data structure for each grid block
        for (int j = 0; j < NUM_THREADS; j++) {
            point_containers[j] = (points_container *) malloc(sizeof(points_container));
            point_containers[j]->points = (point *) malloc(block_sizes[j] * sizeof(point));
            point_containers[j]->count = block_sizes[j];
        }

        // populate each grid block
        int *block_pointer_indices = (int *) calloc((size_t) NUM_THREADS, sizeof(int));
        for (int i = 0; i < points->count; i++) {
            point p = points->points[i];
            int x_block = (int) floor((p.x - min_x) / block_range_x);
            int y_block = (int) floor((p.y - min_y) / block_range_y);

            int block_i = (x_block * (int) sqrt(NUM_THREADS)) + y_block;

            int index = block_pointer_indices[block_i];

            memcpy(&point_containers[block_i]->points[index], &p, sizeof(point));

            block_pointer_indices[block_i]++;
        }

        /*
         * Send each block to each individual process (except rank = 0)
         */
        for (int k = 0; k < NUM_THREADS; k++) {
            final_paths[k] = (int *) malloc(point_containers[k]->count * sizeof(int));

            if (k != 0) {
                // send each point (2 doubles each)
                MPI_Send(&point_containers[k]->count, 1, MPI_INT, k, tag, MPI_COMM_WORLD);

                double *points_extrapolated = (double *) malloc(point_containers[k]->count * 2 * sizeof(double));

                for (int p = 0; p < point_containers[k]->count; p++) {
                    points_extrapolated[p * 2] = point_containers[k]->points[p].x;
                    points_extrapolated[(p * 2) + 1] = point_containers[k]->points[p].y;
                }

                MPI_Send(points_extrapolated, point_containers[k]->count * 2, MPI_DOUBLE, k, tag, MPI_COMM_WORLD);
            } else {
//                current_process_point_container = point_containers[0];

                number_of_points = point_containers[0]->count;

                current_process_point_container = (points_container *) malloc(sizeof(points_container));
                current_process_point_container->count = number_of_points;
                current_process_point_container->points = point_containers[0]->points;

            }
        }
    } else {
        MPI_Recv(&number_of_points, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &stat);

        double *points_extrapolated = (double *) malloc(number_of_points * 2 * sizeof(double));

        MPI_Recv(points_extrapolated, number_of_points * 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &stat);

        current_process_point_container = (points_container *) malloc(sizeof(points_container));
        current_process_point_container->count = number_of_points;
        current_process_point_container->points = (point *) malloc(number_of_points * sizeof(point));

        for (int i = 0; i < number_of_points; i++) {
            point *p = &(current_process_point_container->points[i]);
            *p = point();
//
            (*p).x = points_extrapolated[i * 2];
            (*p).y = points_extrapolated[(i * 2) + 1];
        }
    }

    current_processes_final_path = (int *) malloc(number_of_points * sizeof(int));

    run(current_process_point_container, current_processes_final_path);


    if (rank == 0) {
        final_paths[0] = (int *) malloc(point_containers[0]->count * sizeof(int));
        final_paths[0] = current_processes_final_path;

        for (int i = 1; i < NUM_THREADS; i++) {
            final_paths[i] = (int *) malloc(point_containers[i]->count * sizeof(int));
            MPI_Recv(&final_paths[i][0], point_containers[i]->count, MPI_INT, i, path_tag, MPI_COMM_WORLD, &stat);
        }
    } else {
        // send data to master
        MPI_Send(current_processes_final_path, current_process_point_container->count, MPI_INT, 0, path_tag,
                 MPI_COMM_WORLD);

        MPI_Finalize();
    }

    if (rank == 0) {
        grid_degrees = (int *) calloc(NUM_THREADS, sizeof(int));

        points_container *full_path = (points_container *) malloc(sizeof(points_container));
        full_path->points = (point *) malloc(points->count * sizeof(point));
        full_path->count = points->count;
        int full_path_index = 0;

        // First, mark all "already completed" grids as completed
        for (int g_i = 0; g_i < NUM_THREADS; g_i++) {
            points_container *g = point_containers[g_i];
            if (g->count == 0) {
                grid_degrees[g_i] = 2;
            }
        }

        // for all g \in G // this will run in parallel
        for (int g_i = 0; g_i < 1; g_i++) { // TODO: run for all g_i < NUM_THREADS
            points_container *g = point_containers[g_i];
            grid_degrees[g_i] = 1;
            // p-to-search = all p \in g
            int count_p_to_search = g->count;
            int *p_to_search = (int *) malloc((count_p_to_search + 1) * sizeof(int));
            p_to_search = final_paths[g_i];
            bool done = false;
            // while true {
            while (!done) {
                // for all p \in p-to-search
                for (int p_i = 0; p_i < count_p_to_search; p_i++) {
                    /// Get Nearest neighbor
                    int nn_grid_index = -1;
                    int nn_p;
                    int previous_nn_p;
                    double nn_d = INF;
                    // for all g' \in G, g' != g and degreeAllowed(g')
                    for (int g1_i = 0; g1_i < NUM_THREADS; g1_i++) {
                        if (g1_i != g_i && degreeAllowed(g1_i, NUM_THREADS)) {
                            // for all p' \in g'
//                            printf("LISTED %i\n", grid_degrees[g1_i]);
                            for (int p1_i = 1; p1_i <= point_containers[g1_i]->count; p1_i++) {
                                if (grid_degrees[g1_i] == 0 || point_are_neighbors(full_path->points[0],
                                                                                   point_containers[g1_i]->points[p1_i -
                                                                                                                  1],
                                                                                   point_containers[g1_i],
                                                                                   final_paths[g1_i])) {
                                    // d = distance(p, p')
                                    double d = distance(&point_containers[g_i]->points[p_to_search[p_i] - 1],
                                                        &point_containers[g1_i]->points[p1_i - 1]);
                                    // if d < nn_d
                                    if (d < nn_d) {
                                        // nn_grid_index = g'-index
                                        // nn_p = p'
                                        // nn_d = d
                                        nn_grid_index = g1_i;
                                        nn_p = p1_i;
                                        nn_d = d;
                                    }
                                }
                            }
                        }
                    }

                    /// We've got the Nearest neighbor
//                printf("we got a nearest neighbor \n");
//                    point *a = &point_containers[g_i]->points[p_to_search[p_i] - 1];
//                    point *b = &point_containers[nn_grid_index]->points[nn_p - 1];
//                    printf("connector.put(%i, [%lf %lf;%lf %lf]);\n", NUM_THREADS - degreesLeft(NUM_THREADS), a->x,
//                           a->y,
//                           b->x, b->y);
//
//                    // if best_selected_point
//
//                    for (int i = 0; i < point_containers[g_i]->count; i++) {
//                        printf("%i -> ", final_paths[g_i][i]);
//                    }

                    point *point_we_are_about_to_add = &(point_containers[g_i]->points[p_to_search[p_i] - 1]);

                    if (full_path_index != 0) { // don't add the previous graph if we haven't added anything as of yet
                        // push the values for the grid we are leaving behind
                        // if best neighbor was previous neighbor
                        add_sub_graph(&full_path, &full_path_index, g_i, previous_nn_p, p_to_search, p_i);
                    }

                    full_path->points[full_path_index++] = *point_we_are_about_to_add;
                    full_path->points[full_path_index++] = point_containers[nn_grid_index]->points[nn_p - 1];

                    updateDegree(nn_grid_index);

                    if (degreesAllEqual2(NUM_THREADS)) {
                        int candidate_city_identifier = get_city_identifier(full_path->points[0],
                                                                            point_containers[nn_grid_index]);
//                    int candidate_point_path_index = find_index_in_path_for_point(candidate_city_identifier, point_containers[nn_grid_index], final_path);
                        int *fake_p_to_search = (int *) malloc(sizeof(int));
                        fake_p_to_search[0] = candidate_city_identifier;
                        add_sub_graph(&full_path, &full_path_index, nn_grid_index, nn_p, fake_p_to_search, 0);

                        // break
                        done = true;
                        break;
                    }

                    previous_nn_p = nn_p;
                    g_i = nn_grid_index;

                    // count_p_to_search = get_count_p_to_search(G[nn_grid_index])
                    count_p_to_search = get_count_p_to_search(point_containers[nn_grid_index]);
                    // p-to-search = get_p_to_search(nn_p, G[nn_grid_index])
                    p_to_search = get_p_to_search(nn_p, point_containers[nn_grid_index], final_paths[nn_grid_index]);
                }
            }
        }


        // handle inversions
        int point_count = full_path_index - 1;
        full_path->count = point_count;
        for (int i = 0; i < 100; i++) {
//        printf("INVERSIONS?\n");
            int intersection_count = 0;
            for (int x0 = 0; x0 < point_count - 2; x0++) {
                int x1 = x0 + 1;

                for (int y_i = 0; y_i < point_count - 1 - 2; y_i++) {
                    int y0 = (y_i + x1 + 1) % full_path->count;
                    int y1 = (y0 + 1) % full_path->count;

//            printf("max[%i] %i,%i  -> %i,%i\n", points->count, x0, x1, y0, y1);
                    if (pointsIntersect(full_path->points[x0], full_path->points[x1], full_path->points[y0],
                                        full_path->points[y1])) {
                        handleInversion(&full_path, x1, y0);
                        intersection_count++;
                    }
                }
            }
        }

//        printf("\nconnected_connectors = [");
//        for (int k = 0; k < full_path_index - 1; k++) {
//            printf("%lf %lf;", full_path->points[k].x, full_path->points[k].y);
//        }
//        printf("]\n");

        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;

        printf("time (ms): %llu\n", diff);
    }

    if (rank == 0) {

        // check if things are correct around here.
        MPI_Finalize();
    }

    /**
     *
     * End MPI Specific Code
     *
     */
}