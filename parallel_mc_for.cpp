
#include <omp.h>
#include <stdio.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include "tables.hpp"
using namespace std;

struct Point {
    double x, y, z;
};

struct Triangle {
    Point p[3];
};

struct GridCell {
    Point p[8];
    double values[8];
};

Point interpolate(Point p1, Point p2, double val1, double val2, double iso_level = 0)  // para saber si esta fuera o dentro de la superficie
{
    // comparacion de los dos vértices de ccada arista
    if (fabs(iso_level - val1) < 0.00001)
        return p1;
    if (fabs(iso_level - val2) < 0.00001)
        return p2;

    // las anteriores significan q la superf. pasa justo por p1 o p2
    if (fabs(val1 - val2) < 0.00001)
        return p1;  // no se interseca la superf.
    // esto es como el else
    double mu = (iso_level - val1) / (val2 - val1);
    return {// 3 coord x,y,z
            p1.x + mu * (p2.x - p1.x),
            p1.y + mu * (p2.y - p1.y),
            p1.z + mu * (p2.z - p1.z)};
    // comb. lineal de los vértices y sus valores(mu)
}

GridCell create_grid_cell(double (*func)(double, double, double),
                          double x_min, double y_min, double z_min,
                          double x_max, double y_max, double z_max) {
    GridCell grid;
    // vertic. del cubo
    grid.p[0] = {x_min, y_min, z_min};  // 0
    grid.p[1] = {x_max, y_min, z_min};  // 1
    grid.p[2] = {x_max, y_max, z_min};  // 2
    grid.p[3] = {x_min, y_max, z_min};  // 3
    grid.p[4] = {x_min, y_min, z_max};  // 4
    grid.p[5] = {x_max, y_min, z_max};  // 5
    grid.p[6] = {x_max, y_max, z_max};  // 6
    grid.p[7] = {x_min, y_max, z_max};  // 7 sentido antihorario

    for (int i = 0; i < 8; i++) {
        grid.values[i] = func(grid.p[i].x, grid.p[i].y, grid.p[i].z);
    }

    return grid;
}

int make_polygon(GridCell grid, Triangle *triangles, double iso_level = 0) {
    int n_triangles = 0;  // to count triangles
    int cube_index = 0;   // which edges are intersected
    Point vert_list[12];  // each edge can have a vertex (interpolated)
    // manejo por bits
    if (grid.values[0] < iso_level)
        cube_index |= 1;
    if (grid.values[1] < iso_level)
        cube_index |= 2;
    if (grid.values[2] < iso_level)
        cube_index |= 4;
    if (grid.values[3] < iso_level)
        cube_index |= 8;
    if (grid.values[4] < iso_level)
        cube_index |= 16;
    if (grid.values[5] < iso_level)
        cube_index |= 32;
    if (grid.values[6] < iso_level)
        cube_index |= 64;
    if (grid.values[7] < iso_level)
        cube_index |= 128;

    if (edge_table[cube_index] == 0)
        return 0;
    // 12 as the edges
    if (edge_table[cube_index] & 1)
        vert_list[0] = interpolate(grid.p[0], grid.p[1], grid.values[0], grid.values[1], iso_level);
    if (edge_table[cube_index] & 2)
        vert_list[1] = interpolate(grid.p[1], grid.p[2], grid.values[1], grid.values[2], iso_level);
    if (edge_table[cube_index] & 4)
        vert_list[2] = interpolate(grid.p[2], grid.p[3], grid.values[2], grid.values[3], iso_level);
    if (edge_table[cube_index] & 8)
        vert_list[3] = interpolate(grid.p[3], grid.p[0], grid.values[3], grid.values[0], iso_level);
    if (edge_table[cube_index] & 16)
        vert_list[4] = interpolate(grid.p[4], grid.p[5], grid.values[4], grid.values[5], iso_level);
    if (edge_table[cube_index] & 32)
        vert_list[5] = interpolate(grid.p[5], grid.p[6], grid.values[5], grid.values[6], iso_level);
    if (edge_table[cube_index] & 64)
        vert_list[6] = interpolate(grid.p[6], grid.p[7], grid.values[6], grid.values[7], iso_level);
    if (edge_table[cube_index] & 128)
        vert_list[7] = interpolate(grid.p[7], grid.p[4], grid.values[7], grid.values[4], iso_level);
    if (edge_table[cube_index] & 256)
        vert_list[8] = interpolate(grid.p[0], grid.p[4], grid.values[0], grid.values[4], iso_level);
    if (edge_table[cube_index] & 512)
        vert_list[9] = interpolate(grid.p[1], grid.p[5], grid.values[1], grid.values[5], iso_level);
    if (edge_table[cube_index] & 1024)
        vert_list[10] = interpolate(grid.p[2], grid.p[6], grid.values[2], grid.values[6], iso_level);
    if (edge_table[cube_index] & 2048)
        vert_list[11] = interpolate(grid.p[3], grid.p[7], grid.values[3], grid.values[7], iso_level);

    for (int i = 0; tri_table[cube_index][i] != -1; i += 3) {
        triangles[n_triangles].p[0] = vert_list[tri_table[cube_index][i]];
        triangles[n_triangles].p[1] = vert_list[tri_table[cube_index][i + 1]];
        triangles[n_triangles].p[2] = vert_list[tri_table[cube_index][i + 2]];
        n_triangles++;
    }
    return n_triangles;
}

// ----------

double f(double x, double y, double z) {
    return x * x + y * y + z * z - 1.0;
}

// Number of additional sample points to check when all corners have same sign
const int N_SAMPLES_PER_DIM = 16;

unsigned long long next_power_of_two(unsigned long long n) {
    if (n == 0) return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    n++;
    return n;
}

// Iterative Marching Cubes paralelizado
void iterative_marching_cubes(double (*func)(double, double, double),
                              vector<Triangle> &triangles,
                              double x_min, double y_min, double z_min,
                              double x_max, double y_max, double z_max,
                              double precision, double iso_level = 0.0) {
    double required_res_x = (x_max - x_min) / precision;
    double required_res_y = (y_max - y_min) / precision;
    double required_res_z = (z_max - z_min) / precision;

    double max_required_res = std::max({required_res_x, required_res_y, required_res_z});

    const unsigned long long grid_dim = next_power_of_two(static_cast<unsigned long long>(ceil(max_required_res)));
    const long long nx = grid_dim;
    const long long ny = grid_dim;
    const long long nz = grid_dim;

    const double cell_size_x = (x_max - x_min) / nx;
    const double cell_size_y = (y_max - y_min) / ny;
    const double cell_size_z = (z_max - z_min) / nz;

    vector<vector<Triangle>> local_triangles_per_thread;

#pragma omp parallel
    {
        vector<Triangle> local_triangles;

#pragma omp for collapse(3) schedule(dynamic)
        for (long long k = 0; k < nz; ++k) {
            for (long long j = 0; j < ny; ++j) {
                for (long long i = 0; i < nx; ++i) {
                    double x = x_min + i * cell_size_x;
                    double y = y_min + j * cell_size_y;
                    double z = z_min + k * cell_size_z;

                    GridCell grid = create_grid_cell(func, x, y, z,
                                                     x + cell_size_x, y + cell_size_y, z + cell_size_z);

                    Triangle tris[5];
                    int n = make_polygon(grid, tris, iso_level);

                    if (n > 0) {
                        local_triangles.insert(local_triangles.end(), tris, tris + n);
                    }
                }
            }
        }

#pragma omp critical
        {
            local_triangles_per_thread.push_back(std::move(local_triangles));
        }
    }

    for (const auto &local_triangles : local_triangles_per_thread) {
        triangles.insert(triangles.end(), local_triangles.begin(), local_triangles.end());
    }
}

void draw_surface(double (*func)(double, double, double), const string &output_filename,
                  double x_min, double y_min, double z_min,
                  double x_max, double y_max, double z_max, double precision) {
    vector<Triangle> triangles;
    double iso_level = 0.0;  // Level set for the surface (0 is border)

    // Use iterative
    iterative_marching_cubes(func, triangles, x_min, y_min, z_min,
                             x_max, y_max, z_max, precision, iso_level);

    // Write triangles to PLY file
    ofstream file(output_filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << output_filename << endl;
        return;
    }

    // PLY header
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "element vertex " << (triangles.size() * 3) << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";
    file << "element face " << triangles.size() << "\n";
    file << "property list uchar int vertex_indices\n";
    file << "end_header\n";

    // Write vertices
    for (const auto &triangle : triangles) {
        for (int i = 0; i < 3; i++) {
            file << triangle.p[i].x << " " << triangle.p[i].y << " " << triangle.p[i].z << "\n";
        }
    }

    // Write faces
    for (size_t i = 0; i < triangles.size(); i++) {
        file << "3 " << (i * 3) << " " << (i * 3 + 1) << " " << (i * 3 + 2) << "\n";
    }

    file.close();
    cout << "Generated " << triangles.size() << " triangles" << endl;
    cout << "Surface saved to " << output_filename << endl;
}

double noise(double x, double y) {
    return sin(x * 3) * 0.5 * cos(y * 3);
}

double height_map_2d(double x, double y, double z) {
    return noise(x, y) - z;
}

int main(int argc, char* argv[]) {
    cout << "Iterative Marching Cubes Implementation" << endl;
    
    int n = 512; // default value
    
    // parse n value
    if (argc > 1) {
        n = atoi(argv[1]);
        if (n <= 0) {
            cerr << "Error: n must be a positive integer" << endl;
            cerr << "Usage: " << argv[0] << " [n]" << endl;
            cerr << "  n: grid resolution (default: 512)" << endl;
            return 1;
        }
    }
    
    double precision = 10.0 / n;
    cout << "Using n = " << n << endl;
    cout << "Con precision: " << precision << endl;
    int threads[] = {1, 2, 4, 8, 16, 32, 64, 128};
    for (int t = 0; t < 8; t++) {
        omp_set_num_threads(threads[t]);
        double t_f;
        double t_0 = omp_get_wtime();
        vector<Triangle> triangles;
        iterative_marching_cubes(height_map_2d, triangles, -5, -5, -5, 5, 5, 5, precision);
        t_f = omp_get_wtime() - t_0;
        cout << "Con " << threads[t] << " threads, tiempo: " << t_f << endl;
    }

    return 0;
}