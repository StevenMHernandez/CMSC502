#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

#define INF DBL_MAX

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

double distance(point *a, point *b) {
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2));
}

uint combinations_calculation(int n, int k) {
    uint result = (uint) n;
    for (int i = 2; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

uint count_of_subsets_of_size(int size, int distances_count) {
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
    int *subsets = (int *) calloc(count_of_subsets_of_size(size, distances_count), sizeof(int));
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

double traveling_salesman(double *distances, int count, int *final_path) {
    // bugfix -> handle when count == 1
    if (count == 0) { return 0; }

    double **c = (double **) malloc(static_cast<long>((1 << count)) * sizeof(double *));

    for (int i = 0; i < (1 << count); i++) {
        c[i] = (double *) malloc((count + 1) * sizeof(double));
        memset(c[i], 999, (count + 1) * sizeof(double));
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

                        for (int i = 1; i <= count; i++) { // for each C(S−{j},i)+dij
                            int i_bitmask = 1 << (i - 1);
                            if ((i_bitmask & subsets[s_i]) != 0 && i != j) { // only if i∈S and i≠j
                                int array_position = get_combined_x_y_from_logical(i, j, count);
                                double distance = c[subsets[s_i] & (~j_bitmask)][i] + distances[array_position];
                                if (min > distance) {
                                    min = distance;
                                }
                            }
                        }

                        c[subsets[s_i]][j] = min;
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
    double current_sum = min - distances[get_combined_x_y_from_logical(final_path[0], final_path[1], count)];
    uint current_bitmask = full_bitmask;
    // for each city-slot available
    for (int i = 2;
         i <= count; i++) { // NOTE: skip 1 because we've already added a final_path_cities[0] and final_path_cities[1]
        // for each candidate city path taken
        for (int j = 1; j <= count; j++) {
            int x = final_path[i - 1]; // last city travelled to (used to figure out the distances to the current
            int y = j;
            if (x != y) {
                uint tmp_bitmask = (uint) current_bitmask ^(uint) (1 << (x - 1));
                double calculated_prev = c[tmp_bitmask][j];
                double calculated_distance = distances[get_combined_x_y_from_logical(x, y, count)];
                double calculated_for_current_iteration = calculated_prev + calculated_distance;
                bool is_this_it = calculated_for_current_iteration != INF
                                  && abs(calculated_for_current_iteration - current_sum) > -0.000001
                                  && abs(calculated_for_current_iteration - current_sum) <  0.000001; // TODO: it may be required to change this threshold

                if (is_this_it) {
                    final_path[i] = j;
                    current_sum = calculated_prev;
                    current_bitmask = tmp_bitmask;
                    break;
                }
            }
        }
    }

    return min;
}

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
        }
    }

    return distances;
}

void run_tsp(points_container *container, int *final_path) {
    int city_count = container->count;

    double min = 0;
    double *distances = get_distance_matrix(container);

    traveling_salesman(distances, city_count, final_path);
}

int sign(double x) {
    return static_cast<int>(x / abs(x));
}

/*
 * This link explains triangle orientation a bit https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
 */
int triangleOrientation(point a, point b, point c) {
    return sign((b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y));
}

bool pointsIntersect(point a1, point a2, point b1, point b2) {
    return triangleOrientation(a1, a2, b1) != triangleOrientation(a1, a2, b2) &&
           triangleOrientation(b1, b2, a1) != triangleOrientation(b1, b2, a2);
}

void swap_points(points_container **full_path, int a, int b) {
    point tmp = (*full_path)->points[a];

    (*full_path)->points[a] = (*full_path)->points[b];
    (*full_path)->points[b] = tmp;
}

void handleInversion(points_container **full_path, int x1, int y0) {
    if (x1 < y0) {
        int range = x1 < y0 ? y0 - x1 : ((*full_path)->count) - x1 + y0;

        int path_length = ((*full_path)->count);

        for (int i = 0; i <= ceil(range / 2); ++i) {
            swap_points(full_path, (x1 + i) % path_length, (y0 - i + path_length) % path_length);
        }
    }
}