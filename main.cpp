#include <iostream>
#include "vector"
#include "cmath"

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

Grid<Point> make_points_grid(const std::vector<double> &axis_points) {
    std::vector<std::vector<std::vector<Point> > > result;
    for (double x_axi: axis_points) {
        std::vector<std::vector<Point> > yz_plane;
        for (double y_axi: axis_points) {
            std::vector<Point> z_line;
            for (double z_axi: axis_points) {
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
double u_analytical(const Point &point, int t, double l, double a_t) {
    return std::sin(3 * M_PI * point.x / l) * sin(2 * M_PI * point.y / l) * sin(2 * M_PI * point.z / l)
           * cos(a_t * (double) t + 4 * M_PI);
}

double phi(const Point &point, double l, double a_t) {
    return u_analytical(point, 0, l, a_t);
}

// TODO index out of range
double laplassian(int i, int j, int k, const Grid<double> &prev_v, double h) {
    double result = 0;

    result += prev_v[i - 1][j][k] - 2 * prev_v[i][j][k] + prev_v[i + 1][j][k];
    result += prev_v[i][j - 1][k] - 2 * prev_v[i][j][k] + prev_v[i][j + 1][k];
    result += prev_v[i][j][k - 1] - 2 * prev_v[i][j][k] + prev_v[i][j][k + 1];

    return result / (h * h);
}
// endregion math functions

// region build values
Grid<double> build_initial_prev_values_1(const Grid<Point> &points_grid, double l, double a_t) {
    Grid<double> result;
    for (const auto &yz_plane: points_grid) {
        std::vector<std::vector<double> > yz_plane_values;
        for (const auto &z_line: yz_plane) {
            std::vector<double> z_line_values;
            for (auto point: z_line) {
                z_line_values.push_back(phi(point, l, a_t));
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}

Grid<double>
build_initial_prev_values_2(const Grid<Point> &points_grid, const Grid<double> &previous_values_1, double h) {
    Grid<double> result;
    for (int i = 0; i < points_grid.size(); i++) {
        const auto &yz_plane = points_grid[i];
        std::vector<std::vector<double> > yz_plane_values;
        for (int j = 0; j < yz_plane.size(); j++) {
            const auto &z_line = yz_plane[j];
            std::vector<double> z_line_values;
            for (int k = 0; k < z_line.size(); k++) {
                z_line_values.push_back(previous_values_1[i][j][k] +
                                        0.5 * laplassian(i, j, k, previous_values_1, h));
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
        double h,
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
                        (t * t) * laplassian(i, j, k, previous_values_2, h)
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
    if (argc < 3) {
        std::cout << "Usage: ./main L n" << std::endl;
        return 0;
    }

    double l = std::stod(argv[1]);
    long n = strtol(argv[2], nullptr, 10);
    double h = l / (double) n;

    double l_2 = l * l;
    double a_t = M_PI * std::sqrt(9 / l_2 + 4 / l_2 + 4 / l_2);

    std::vector<double> axis_points = divide_axis(n, h);
    Grid<Point> points_grid = make_points_grid(axis_points);

    Grid<double> previous_values_1;
    Grid<double> previous_values_2;
    Grid<double> current_values;
    for (int t = 0; t <= MAX_T; t++) {
        if (t == 0) {
            previous_values_1 = build_initial_prev_values_1(points_grid, l, a_t);
            continue;
        }
        if (t == 1) {
            previous_values_2 = build_initial_prev_values_2(points_grid, previous_values_1, h);
            continue;
        }

        current_values = build_current_values(points_grid, previous_values_1, previous_values_2, h, t);

        previous_values_1 = previous_values_2;
        previous_values_2 = current_values;
    }

    print_values(current_values);

    return 0;
}
