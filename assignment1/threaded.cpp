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

int factorial(int n) {
    int total = 1;

    for (int i = n; i > 0; i--) {
        total *= i;
    }

    return total;
}

int count_of_subsets_of_size(int size, int distances_count) {
    return factorial(distances_count) / (factorial(size) * factorial(distances_count - size));
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


double traveling_salesman(double *distances, int count) {
    // bugfix -> handle when count == 1
    if (count == 1) { return 0; }
    // bugfix -> handle when count == 2
    if (count == 2) { return 2 * distances[1]; }

    printf("Making %d * sizeof(double)\n", (1 << count) * count);
    double **c = (double **) malloc((1 << count) * sizeof(double *));

    for (int i = 0; i < (1 << count); i++) {
        c[i] = (double *) calloc(count, sizeof(double));
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
                for (int j = 2; j <= count; j++) {
                    if (subsets[s_i] & (1 << (j - 1))) { // check if j∈S
                        // C(S, j) = min{C(S−{j},i)+dij: i∈S, i≠j}
                        double min = INF;
                        for (int i = 1; i <= count; i++) {
                            int bitmask = 1 << (i - 1);
                            if ((bitmask & subsets[s_i]) != 0 && i != j) { // only if i∈S and i≠j
                                if (s != 2) {
                                    // pass
                                }
                                if (s == 2 || ~bitmask & subsets[s_i] & 1) { // do not check bit 1
                                    double distance =
                                            distances[((j - 1) * count) + (i - 1)] + c[~bitmask & subsets[s_i]][j];
                                    if (min > distance) {
                                        min = distance;
                                    }
                                }
                            }
                        }
                        c[subsets[s_i]][j] = min;

//                        std::cout << "\n    stored at c[" << std::bitset<16>(static_cast<unsigned long long int>((subsets[s_i])))
//                        << "][" << j << "]" << std::endl;

//                        std::cout << "C(" << j << ",{"
//                                << std::bitset<4>(static_cast<unsigned long long int>((subsets[s_i] & ~(1 << (j-1)) & ~1)))
//                                << "}) = " << d_a << " + " << "C(" << d_b_1 << ",{"
//                                << std::bitset<4>(static_cast<unsigned long long int>(d_b_2))
//                                << "}) = " << min_a << " + " << min_b << " = " << min << std::endl;

//                        std::cout << std::bitset<16>(static_cast<unsigned long long int>(subsets[s_i]))
//                                << " " << j << ": " << min << std::endl;
                    }
                }
            }
        }
//        free(subsets);
    }

    // return minjC({1,...,n},j)+dj1
    double min = INF;

    uint full_bitmask = 0;
    for (int i = 0; i < count; i++) {
        full_bitmask = full_bitmask | (1 << i);
    }

    for (int i = 1; i < count; i++) {
        double distance = c[full_bitmask][i];

        if (min > distance) {
            min = distance;
        }
    }

//    for (int i = 0; i < (1 << count); i++) {
//        free(c[i]);
//    }
//    free(c);

    return min;
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

    printf("%s\n", argv[1]);

    int NUM_THREADS = 16;

    // create the point containers that will be used per thread.
    points_container **point_containers = (points_container **) malloc(sizeof(points_container *) * NUM_THREADS);


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
    double block_range_x = ceil(range_x / sqrt(NUM_THREADS));
    double range_y = max_y - min_y;
    double block_range_y = ceil(range_y / sqrt(NUM_THREADS));

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

//    // validate splits TODO: remove
//    for (int j = 0; j < NUM_THREADS; j++) {
//        for (int i = 0; i < point_containers[j]->count; i++) {
//            printf("block %i has a point %lf,%lf (%i out of %i)\n", j, point_containers[j]->points[i].x, point_containers[j]->points[i].y, i, point_containers[j]->count);
//        }
//    }

    // TODO: create threads
    for (int j = 0; j < NUM_THREADS; j++) {
        double min = 0;
        if (point_containers[j]->count > 1) {
            double *distances = get_distance_matrix(point_containers[j]);
            min = traveling_salesman(distances, point_containers[j]->count);
        }

        printf("TSP-min is: %lf\n", min);
        printf(". . . ");
    }
}