#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

/**
 * Point as created in `gen_cities.sh`.
 */
struct point {
    double x;
    double y;
};

/**
 * Points container is used to keep track of the
 * array of points and most importantly the total
 * number of points in the array.
 */
struct points_container {
    int count;
    point *points;
};

/**
 *
 * GLOBAL VAR
 *
 */
points_container **point_containers;
int **final_paths;
double *results;
int *grid_degrees;

/**
 *
 * Load data from file
 *
 */
int get_file_line_count(char *filename) {
    int count = 0;

    ifstream file(filename);

    for (string line; getline(file, line);) {
        count++;
    }

    return count;
}

points_container *get_the_points(char *filename) {
    int line_count = get_file_line_count(filename);

    point *points = (point *) malloc(sizeof(point) * line_count);

    ifstream file(filename);

    int i = 0;

    for (string line; getline(file, line);) {
        points[i] = point();

        double *x = 0, y = 0;

        char *line_chars = const_cast<char *>(line.c_str());

        sscanf(line_chars, "%lf %lf", &points[i].x, &points[i].y);

        i++;
    }

    points_container *container = (points_container *) malloc(sizeof(points_container));

    container->count = line_count;
    container->points = points;

    return container;
}

double distance(point *a, point *b) {
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
}

/**
 * Calculate Distances
 */

double *get_distance_matrix(points_container *container) {
    double *distances = (double *) malloc((container->count * container->count) * sizeof(double));

    /*
     * Understanding we are repeating calculations that we don't need to,
     * we start from left to right on each row.
     *
     * It is understood that this matrix could always contains zeros along the diagonal (no need to compute this)
     * and that the values are mirrors across the diagonal.
     */
    for (int i = 0; i < container->count; i++) {
        for (int j = 0; j < container->count; j++) {
            distances[(j * container->count) + i] = distance(&container->points[i], &container->points[j]);
//            distances[(j * container->count) + i] = sqrt(pow(container->points[i].x - container->points[j].x, 2) +
//                                                         pow(container->points[i].y - container->points[j].y, 2));
//            printf("distance between %i and %i is %lf\n", i, j, distances[(j * container->count) + i]);
        }
    }

    return distances;
}

/**
 *
 * Traveling Salesman Algorithm
 *
 */

#define INF DBL_MAX;

int combinations_calculation(int n, int k) {
    int result = n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

int count_of_subsets_of_size(int size, int distances_count) {
    return combinations_calculation(distances_count, size);
}

void combinations(int *arr, int size, int distances_count, int *index, int bit_position, int num_bits_set, int data) {
    // if bit_position less than distance_count
    if (bit_position < distances_count) {
        // "place `0` at bit position
        // recurse.
        combinations(arr, size, distances_count, index, bit_position + 1, num_bits_set, data);

        // place `1` bit at bit_position
        data = (1 << bit_position) | data;

        // if num_bits_set + 1 == size
        if (num_bits_set + 1 == size) {
            // save in arr[*index]
            arr[*index] = data;
            (*index)++;
        } else {
            combinations(arr, size, distances_count, index, bit_position + 1, num_bits_set + 1, data);
        }
    }
}

int *get_subsets_of_size(int size, int distances_count) {
    int *subsets = (int *) calloc(count_of_subsets_of_size(size, distances_count), sizeof(int *));
    int index = 0;
    combinations(subsets, size, distances_count, &index, 0, 0, 0);
    return subsets;
}


/*
 * One of the hardest parts of this is the difference between 0-based physical
 * and 1-based logical city indexing.
 *
 * This function gets the array-index from distances[][] as _logical 1-based numbering_
 */
int get_combined_x_y_from_logical(int x, int y, int width) {
    return (y - 1) * width + (x - 1);
}

/*
 * This function gets the array-index from distances[][] as _physical 0-based numbering_
 */
int get_combined_x_y_from_physical(int x, int y, int width) {
    return y * width + x;
}

double traveling_salesman(double *distances, int count, int *final_path) {
    // bugfix -> handle when count == 1
    if (count == 0) { return 0; }
    // bugfix -> handle when count == 2
//    if (count == 2) { return 2 * distances[1]; }

//    printf("Making %li * sizeof(double)\n", static_cast<long>((1 << count) * count));
    double **c = (double **) malloc(static_cast<long>((1 << count)) * sizeof(double *));

    for (int i = 0; i < (1 << count); i++) {
        c[i] = (double *) calloc(count + 1, sizeof(double));
    }

    // C({1},1) = 0
    c[1][1] = 0;
    // for s = 2 to n:
    for (int s = 2; s <= count; s++) {
        // for all subsets S ⊆ {1,2,...,n} of size s and containing 1:
        int *subsets = get_subsets_of_size(s, count);
        for (uint s_i = 0; s_i < count_of_subsets_of_size(s, count); s_i++) {
            if ((subsets[s_i] & 1) != 0) { // containing 1

                // C(S,1) = ∞
                c[subsets[s_i]][1] = INF;
                // for all j∈S, j≠1:
                for (int j = 2; j <= count; j++) { // check if j≠1 too
                    int j_bitmask = 1 << (j - 1);
                    if (subsets[s_i] & j_bitmask && j != 1) { // check if j∈S
                        // C(S, j) = min{C(S−{j},i)+dij: i∈S, i≠j}
                        double min = INF;
                        int min_i = 0;
                        int min_j_bitmask = 0;

                        for (int i = 1; i <= count; i++) { // for each C(S−{j},i)+dij
                            int i_bitmask = 1 << (i - 1);
                            if ((i_bitmask & subsets[s_i]) != 0 && i != j) { // only if i∈S and i≠j
                                int array_position = get_combined_x_y_from_logical(i, j, count);
                                double distance = c[subsets[s_i] & (~j_bitmask)][i] + distances[array_position];
                                if (min > distance) {
                                    min_j_bitmask = j_bitmask;
                                    min_i = i;
                                    min = distance;
                                }
                            }
                        }

                        c[subsets[s_i]][j] = min;

//                        string previous_bitmask_string = std::bitset<12>(
//                                static_cast<unsigned long long int>(subsets[s_i] & ~min_j_bitmask)).to_string();
//                        string x = std::bitset<12>(static_cast<unsigned long long int>((subsets[s_i]))).to_string();
//                        printf("\n\n    dl[%i][%i] + c[%s][%i]", j, min_i, previous_bitmask_string.c_str(), min_i);
//                        printf("\n   %lf + %lf = %lf stored at c[%s][%i]",
//                               distances[get_combined_x_y_from_logical(min_i, j, count)],
//                               c[subsets[s_i] & ~min_j_bitmask][min_i],
//                               min, x.c_str(), j);
                    }
                }
            }
        }
//        free(subsets);
    }

    // return minjC({1,...,n},j)+dj1
    double min = INF;
    int min_i = 0;

    uint full_bitmask = 0;
    for (int i = 0; i < count; i++) {
        full_bitmask = full_bitmask | (1 << i);
    }

    for (int i = 2; i <= count; i++) {
        double distance = c[full_bitmask][i] + distances[get_combined_x_y_from_logical(1, i, count)];

        if (min > distance) {
            min = distance;
            min_i = i;
        }
    }

    /*
     * Go backwards and determine the path
     */
//    int *final_path_cities = (int *) malloc((count + 1) * sizeof(int));

    final_path[0] = 1;
    // we know the next city because we calculated that just now:
    final_path[1] = min_i;

//    printf("\n hey at least we know [0] -> 0 and [1] -> %i", min_i);

    // find each city travelled to
    double current_sum =
            min - distances[get_combined_x_y_from_logical(final_path[0], final_path[1], count)];
    uint current_bitmask = full_bitmask;
    // for each city-slot available
    for (int i = 2;
         i <= count; i++) { // NOTE: skip 1 because we've already added a final_path_cities[0] and final_path_cities[1]
//        printf("\n\n\n");
//    for (int i = 1; i < count + 1; i++) { // NOTE: skip 1 because we've already added a final_path_cities[0]
        // for each candidate city path taken
        for (int j = 1; j <= count; j++) {
            int x = final_path[i - 1]; // last city travelled to (used to figure out the distances to the current
            int y = j;
            if (x != y) {
                uint tmp_bitmask = (uint) current_bitmask ^(uint) (1 << (x - 1));
//                string a_str = std::bitset<12>(static_cast<unsigned long long int>(tmp_bitmask)).to_string();
//                string sum_str = std::bitset<12>(static_cast<unsigned long long int>(current_bitmask)).to_string();
//
                double calculated_prev = c[tmp_bitmask][j];
                double calculated_distance = distances[get_combined_x_y_from_logical(x, y, count)];
                double calculated_for_current_iteration = calculated_prev + calculated_distance;

//                printf("\nc[%s][%i] + d[%i][%i] = c[%s][%i]\n", a_str.c_str(), j, x, y, sum_str.c_str(), j);
//                printf("  %lf + %lf", calculated_prev, calculated_distance);
//                printf(" = %lf", calculated_for_current_iteration);
//                printf(" ==? %lf ???", current_sum);

                bool is_this_it = abs(calculated_for_current_iteration - current_sum) < 0.000001;

//                if (is_this_it) {
//                    printf(" YES!!!  <---");
//                }

                if (is_this_it) {
                    final_path[i] = j;
//                int x = final_path_cities[i - 1];
//                int y = final_path_cities[i];
                    current_sum = current_sum - distances[get_combined_x_y_from_logical(x, y, count)];
                    current_bitmask = tmp_bitmask;
                    break;
                }
            }
        }
    }

    return min;
}

void *run(void *ptr) {
    int thread_identifier = *(int *) ptr;

    int city_count = point_containers[thread_identifier]->count;

    double min = 0;
    double *distances = get_distance_matrix(point_containers[thread_identifier]);

    final_paths[thread_identifier] = (int *) malloc(city_count * sizeof(int));

    min = traveling_salesman(distances, city_count, final_paths[thread_identifier]);

    results[thread_identifier] = min;
}


int find_index_in_path_for_point(int p_index, points_container *container, int *final_path) {
    for (int i = 0; i < container->count; i++) {
        if ((p_index) == final_path[i]) {
            return i;
        }
    }

    return -1; // this should never occur
}

/*
 * p_index is actually the city-identifier! (which is 1-based)
 */
// func next_neighbor_point(p, g):
int next_neighbor_point(int p_index, points_container *container, int *final_path) {
//    for (int i = 0; i < container->count; i++) {
//        printf("%i -> ", final_path[i]);
//    }

    // g[&p + 1 % g->points->count]
    int index_for_current_city = find_index_in_path_for_point(p_index, container, final_path);
//    printf(" for [%i] which is (%i) next is [%i] which is (%i) \n", index_for_current_city, p_index, (index_for_current_city + 1) % container->count, final_path[(index_for_current_city + 1) % container->count]);
    return final_path[(index_for_current_city + 1) % container->count];
}

/**
 *
 * @param p_index this is the index for the point from the container
 * @param container
 * @param final_path
 * @return
 */
// func prev_neighbor_point(p, g):
int prev_neighbor_point(int p_index, points_container *container, int *final_path) {
//    for (int i = 0; i < container->count; i++) {
//        printf("%i -> ", final_path[i]);
//    }


    // g[&p - 1 % g->points->count]
    int index_for_current_city = find_index_in_path_for_point(p_index, container, final_path);
//    int index = find_index_in_path_for_point(p_index, container, final_path);
//    printf(" for [%i] which is (%i) prev is [%i] which is (%i) \n", index_for_current_city, p_index, (index_for_current_city + (container->count - 1)) % container->count, final_path[(index_for_current_city + (container->count - 1)) % container->count]);
    return final_path[(index_for_current_city + (container->count - 1)) % container->count];
}


int get_count_p_to_search(points_container *container) {
    if (container->count == 1 || container->count == 2) {
        return 1;
    }

    return 2;
}

/**
 *
 * @param p_index this is the index for the point from the container
 * @param container
 * @param final_path
 * @return
 */
// nn_p, &point_containers[nn_grid_index], final_paths[nn_grid_index]
int *get_p_to_search(int p_index, points_container *container, int *final_path) {
    int *result = (int *) malloc(get_count_p_to_search(container) * sizeof(int));

    // if g->points->count == 1:
    if (container->count == 1) {
        // return [p]
        result[0] = p_index;
    } else if (container->count == 2) {
        // return [next_neighbor_point(p)];
        result[0] = next_neighbor_point(p_index, container, final_path);
    } else {
        result[0] = next_neighbor_point(p_index, container, final_path);
        result[1] = prev_neighbor_point(p_index, container, final_path);
    }

    return result;
}

bool are_these_points_the_same(point a, point b) {
    return a.x == b.x && a.y == b.y;
}

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

int get_city_identifier(point p, points_container *container) {
    int city_identifier = -1;

    for (int i = 0; i < container->count; i++) {
        if (p.x == container->points[i].x && p.y == container->points[i].y) {
            city_identifier = i + 1;
        }
    }

    return city_identifier;
}

bool point_are_neighbors(point key_point, point candidate_point, points_container *container, int *path) {
    int city_identifier = get_city_identifier(key_point, container);
    int candidate_city_identifier = get_city_identifier(candidate_point, container);

    int city_identifier_0 = prev_neighbor_point(city_identifier, container, path);
    int city_identifier_1 = next_neighbor_point(city_identifier, container, path);

    int candidate_point_path_index = find_index_in_path_for_point(candidate_city_identifier, container, path);

    return candidate_point_path_index == city_identifier_0 || candidate_point_path_index == city_identifier_1;
}

void add_sub_graph(points_container **full_path, int *full_path_index, int g_i, int previous_nn_p, int *p_to_search, int p_i) {
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

int triangleOrientation(point a, point b, point c) {
    double val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);

    // colinear
    if (val == 0) {
        return 0;
    }

    // clock or counterclock wise
    return val > 0 ? 1 : 2;
}

bool pointsIntersect(point a1, point a2, point b1, point b2) {
    int o1 = triangleOrientation(a1, a2, b1);
    int o2 = triangleOrientation(a1, a2, b2);
    int o3 = triangleOrientation(b1, b2, a1);
    int o4 = triangleOrientation(b1, b2, a2);

    return o1 != o2 && o3 != o4;

//     General case
//    if (o1 != o2 && o3 != o4) {
//        return true;
//    }

    return false;
}

void swap_points(points_container **full_path, int a, int b) {
//    if (a == b) {
//        printf("WOW They equal\n");
//    }
//    printf("swap %i <-> %i\n", a, b);
    point tmp = (*full_path)->points[a];

    (*full_path)->points[a] = (*full_path)->points[b];
    (*full_path)->points[b] = tmp;
}

void handleInversion(points_container **full_path, int x1, int y0) {
//    printf("handle inversion for %i %i\n", x1, y0);
    if (x1 < y0) {
        int range = x1 < y0 ? y0 - x1 : ((*full_path)->count) - x1 + y0;
//    printf("range2: %i, %i\n", range - (2*floor(range/2)), range);

        int path_length = ((*full_path)->count);

        for (int i = 0; i <= ceil(range / 2); ++i) {
            swap_points(full_path, (x1 + i) % path_length, (y0 - i + path_length) % path_length);
        }
    }
}

/**
 *
 * MAIN
 *
 */
int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Threaded-Usage: ./main tmp.txt\n" << endl;
        exit(0);
    }
    char *filename = argv[1];

    struct timespec start, end;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    int NUM_THREADS = 16;

    // create the point containers that will be used per thread.
    point_containers = (points_container **) malloc(sizeof(points_container *) * NUM_THREADS);
    final_paths = (int **) malloc(NUM_THREADS * sizeof(int *));

    points_container *points = get_the_points(filename);


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

    // create the threads to get the grid of tsp sub-solutions
    pthread_t *threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));
    int *thread_identifier = (int *) malloc(NUM_THREADS * sizeof(int));
    results = (double *) malloc(NUM_THREADS * sizeof(double));

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_identifier[i] = i;
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_create(&threads[i], NULL, &run, (void *) &thread_identifier[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        void *result = 0;
        pthread_join(threads[i], NULL);

//        printf("returned to main thread: %lf\n", results[i]);

        int *final_path = final_paths[i];
        points_container *p = point_containers[i];

//        /**
//         * Print the paths for display in matlab
//         */
//        printf("\nseq_paths.put(%i, [", i);
//        for (int j = 0; j < p->count; j++) {
//            printf("%lf %lf;", p->points[final_path[j] - 1].x, p->points[final_path[j] - 1].y);
//        }
//        printf("]);");
    }


    grid_degrees = (int *) calloc(NUM_THREADS, sizeof(int));

    points_container *full_path = (points_container *) malloc(sizeof(points_container));
    full_path->points = (point *) malloc(points->count * sizeof(point));
    full_path->count = points->count;
    int full_path_index = 0;

    /// First, mark all "already completed" grids as completed
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
                        for (int p1_i = 1; p1_i <= point_containers[g1_i]->count; p1_i++) {
                            if (grid_degrees[g1_i] == 0 || point_are_neighbors(full_path->points[0],
                                                                               point_containers[g1_i]->points[p1_i - 1],
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
                point *a = &point_containers[g_i]->points[p_to_search[p_i] - 1];
                point *b = &point_containers[nn_grid_index]->points[nn_p - 1];

                point *point_we_are_about_to_add = &point_containers[g_i]->points[p_to_search[p_i] - 1];

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

                if (pointsIntersect(full_path->points[x0], full_path->points[x1], full_path->points[y0],
                                    full_path->points[y1])) {
                    handleInversion(&full_path, x1, y0);
                    intersection_count++;
                }
            }
        }
    }

//    printf("\nconnected_connectors = [");
//    for (int k = 0; k < full_path_index - 1; k++) {
//        printf("%lf %lf;", full_path->points[k].x, full_path->points[k].y);
//    }
//    printf("]\n");

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;

//    printf("time (ms): %llu\n", diff);

    double min = 0;

    for (int k = 0; k < full_path_index - 1 - 1; k++) {
        min += distance(&full_path->points[k], &full_path->points[k + 1]);
//        printf("%lf %lf;", full_path->points[k].x, full_path->points[k].y);
    }

    printf("%lf", min);
}