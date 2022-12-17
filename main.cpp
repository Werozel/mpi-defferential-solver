#include <iostream>
#include "vector"
#include "cmath"
#include "mpi.h"

#define MAX_T 20

class Point;

template<typename T>
using Grid = std::vector<std::vector<std::vector<T> > >;

// region Point
class Point {
public:
    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    double x;
    double y;
    double z;
};

std::vector<double> divide_axis(long n, double h) {
    std::vector<double> result;
    for (int i = 0; i < n; i++) {
        result.push_back(h * (double) (i + 1));
    }
    return result;
}

Grid<Point> make_points_grid(
        const std::vector<double> &axis_points_x,
        const std::vector<double> &axis_points_y,
        const std::vector<double> &axis_points_z
) {
    std::vector<std::vector<std::vector<Point> > > result;
    for (double x_axi: axis_points_x) {
        std::vector<std::vector<Point> > yz_plane;
        for (double y_axi: axis_points_y) {
            std::vector<Point> z_line;
            for (double z_axi: axis_points_z) {
                z_line.emplace_back(x_axi, y_axi, z_axi);
            }
            yz_plane.push_back(z_line);
        }
        result.push_back(yz_plane);
    }
    return result;
}
//endregion Point

// region helpers
void print_values(const Grid<double> &grid) {
    for (const auto &yz_plane: grid) {
        for (const auto &z_line: yz_plane) {
            for (auto value: z_line) {
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void print_points(const Grid<Point> &points_grid) {
    for (const auto &yz_plane: points_grid) {
        for (const auto &z_line: yz_plane) {
            for (auto point: z_line) {
                std::cout << point.x << " " << point.y << " " << point.z << std::endl;
            }
        }
    }
}
// endregion helpers

// region math functions
double u_analytical(
        const Point &point,
        int t,
        double lx, double ly, double lz,
        double a_t
) {
    return std::sin(3 * M_PI * point.x / lx) * sin(2 * M_PI * point.y / ly) * sin(2 * M_PI * point.z / lz)
           * cos(a_t * (double) t + 4 * M_PI);
}

double phi(
        const Point &point,
        double lx, double ly, double lz,
        double a_t
) {
    return u_analytical(point, 0, lx, ly, lz, a_t);
}

double laplassian(
        int i, int j, int k,
        const Grid<double> &prev_v,
        double hx, double hy, double hz
) {
    double result = 0;

    auto n = prev_v.size();
    result += (i == 0 || i == n - 1)
              ? 0
              : (prev_v[i - 1][j][k] - 2 * prev_v[i][j][k] + prev_v[i + 1][j][k]) / (hx * hx);
    auto m = prev_v[i].size();
    result += (j == 0 || j == m - 1)
              ? 0
              : (prev_v[i][j - 1][k] - 2 * prev_v[i][j][k] + prev_v[i][j + 1][k]) / (hy * hy);
    auto t = prev_v[i][j].size();
    result += (k == 0 || k == t - 1)
              ? 0
              : (prev_v[i][j][k - 1] - 2 * prev_v[i][j][k] + prev_v[i][j][k + 1]) / (hz * hz);

    return result;
}
// endregion math functions

// region build values
Grid<double> build_initial_prev_values_1(
        const Grid<Point> &points_grid,
        double lx,
        double ly,
        double lz,
        double a_t
) {
    Grid<double> result;
    for (const auto &yz_plane: points_grid) {
        std::vector<std::vector<double> > yz_plane_values;
        for (const auto &z_line: yz_plane) {
            std::vector<double> z_line_values;
            for (auto point: z_line) {
                z_line_values.push_back(phi(point, lx, ly, lz, a_t));
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}

Grid<double> build_initial_prev_values_2(
        const Grid<Point> &points_grid,
        const Grid<double> &previous_values_1,
        double hx, double hy, double hz
) {
    Grid<double> result;
    for (int i = 0; i < points_grid.size(); i++) {
        const auto &yz_plane = points_grid[i];
        std::vector<std::vector<double> > yz_plane_values;
        for (int j = 0; j < yz_plane.size(); j++) {
            const auto &z_line = yz_plane[j];
            std::vector<double> z_line_values;
            for (int k = 0; k < z_line.size(); k++) {
                z_line_values.push_back(previous_values_1[i][j][k] +
                                        0.5 * laplassian(i, j, k, previous_values_1, hx, hy, hz));
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}

Grid<double> build_current_values(
        const Grid<Point> &points_grid,
        const Grid<double> &previous_values_1,
        const Grid<double> &previous_values_2,
        double hx, double hy, double hz,
        double t
) {
    Grid<double> result;
    for (int i = 0; i < points_grid.size(); i++) {
        const auto &yz_plane = points_grid[i];
        std::vector<std::vector<double> > yz_plane_values;
        for (int j = 0; j < yz_plane.size(); j++) {
            const auto &z_line = yz_plane[j];
            std::vector<double> z_line_values;
            for (int k = 0; k < z_line.size(); k++) {
                z_line_values.push_back(
                        (t * t) * laplassian(i, j, k, previous_values_2, hx, hy, hz)
                        + 2 * previous_values_2[i][j][k]
                        - previous_values_1[i][j][k]
                );
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}
// endregion build values

int main(int argc, char *argv[]) {
    // region init
    double lx, ly, lz;
    long n;
    if (argc == 3) {
        lx = ly = lz = std::stod(argv[1]);
        n = strtol(argv[2], nullptr, 10);
    } else if (argc == 5) {
        lx = std::stod(argv[1]);
        ly = std::stod(argv[2]);
        lz = std::stod(argv[3]);
        n = strtol(argv[4], nullptr, 10);
    } else {
        std::cout << "Usage:" << std::endl << "\t./main Lx Ly Lz n" << std::endl << "\t./main L n" << std::endl;
        return 0;
    }

    double hx = lx / (double) n;
    double hy = ly / (double) n;
    double hz = lz / (double) n;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a_t = M_PI * std::sqrt(9 / (lx * lx) + 4 / (ly * ly) + 4 / (lz * lz));

    std::vector<double> axis_points_x = divide_axis(n, hx);
    std::vector<double> axis_points_y = divide_axis(n, hy);
    std::vector<double> axis_points_z = divide_axis(n, hz);
    Grid<Point> points_grid = make_points_grid(axis_points_x, axis_points_y, axis_points_z);
    // endregion init

    MPI_Barrier(MPI_COMM_WORLD);

    double start_time = MPI_Wtime();

    Grid<double> previous_values_1;
    Grid<double> previous_values_2;
    Grid<double> current_values;
    for (int t = 0; t <= MAX_T; t++) {
        if (t == 0) {
            previous_values_1 = build_initial_prev_values_1(points_grid, lx, ly, lz, a_t);
            continue;
        }
        if (t == 1) {
            previous_values_2 = build_initial_prev_values_2(points_grid, previous_values_1, hx, hy, hz);
            continue;
        }

        current_values = build_current_values(points_grid, previous_values_1, previous_values_2, hx, hy, hz, t);

        previous_values_1 = previous_values_2;
        previous_values_2 = current_values;
    }

    double end_time = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Time: " << end_time - start_time << std::endl;
    }

    if (!rank) {
        print_values(current_values);
    }

    MPI_Finalize();

    return 0;
}
