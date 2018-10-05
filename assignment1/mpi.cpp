#include <stdio.h>
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
            distances[(j * container->count) + i] = sqrt(pow(container->points[i].x - container->points[j].x, 2) +
                                                         pow(container->points[i].y - container->points[j].y, 2));
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
    if (count == 1) { return 0; }
    // bugfix -> handle when count == 2
    if (count == 2) { return 2 * distances[1]; }

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
    for (int i = 2; i <= count; i++) { // NOTE: skip 1 because we've already added a final_path_cities[0] and final_path_cities[1]
//        printf("\n\n\n");
//    for (int i = 1; i < count + 1; i++) { // NOTE: skip 1 because we've already added a final_path_cities[0]
        // for each candidate city path taken
        for (int j = 1; j <= count; j++) {
            int x = final_path[i - 1]; // last city travelled to (used to figure out the distances to the current
            int y = j;
            if (x != y) {
                uint tmp_bitmask = (uint) current_bitmask ^(uint) (1 << (x - 1));
                string a_str = std::bitset<12>(static_cast<unsigned long long int>(tmp_bitmask)).to_string();
                string sum_str = std::bitset<12>(static_cast<unsigned long long int>(current_bitmask)).to_string();

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

/**
 *
 * MAIN
 *
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Threaded-Usage: ./main tmp.txt\n" << endl;
        exit(0);
    }
    char *filename = argv[1];

//    printf("%s\n", argv[1]);

    points_container *points = get_the_points(filename);


    struct timespec start, end;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    double *distances = get_distance_matrix(points);
    int *final_path = (int *) malloc(points->count * sizeof(int));
    double min = traveling_salesman(distances, points->count, final_path);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    uint64_t diff = (1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec) / 1e6;

//    for (int i = 0; i < points->count; i++) {
//        if (i != 0) {
//            printf(" -> ");
//        }
//        printf("%i", final_path[i]);
//    }

//    printf("\n\nTSP-min is: %lf\n", min);

    printf("time (ms): %llu\n", diff);
}