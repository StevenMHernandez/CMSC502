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
    int count = 0;
    point *points;
};

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

        char* line_chars = const_cast<char*>(line.c_str());

        sscanf(line_chars, "%lf %lf", &points[i].x, &points[i].y);

        i++;
    }

    points_container *container = (points_container *) malloc(sizeof(points_container));

    container->count = line_count;
    container->points = points;

    return container;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage: ./main tmp.txt\n" << endl;
        exit(0);
    }

    printf("%s\n", argv[1]);

    points_container *points = get_the_points(argv[1]);

    for (int i = 0; i < points->count; i++) {
        printf(" || I got: x -> %f, y -> %f\n", points->points[i].x, points->points[i].y);
    }
}