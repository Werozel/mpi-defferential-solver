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
        double tau,
        double lx, double ly, double lz,
        double a_t
) {
    return std::sin(3 * M_PI * point.x / lx) * std::sin(2 * M_PI * point.y / ly) * std::sin(2 * M_PI * point.z / lz)
           * std::cos(a_t * tau + 4 * M_PI);
}

double phi(
        const Point &point,
        double lx, double ly, double lz,
        double a_t
) {
    return u_analytical(point, 0, lx, ly, lz, a_t);
}

double plain_laplassian(
        int i, int j, int k,
        const Grid<double> &curr_v,
        const int curr_block_dims[3],
        double hx, double hy, double hz
) {
    double result = 0;

    result += (i == 0 || i == curr_block_dims[0] - 1) ? 0 : (curr_v[i + 1][j][k] - 2 * curr_v[i][j][k] + curr_v[i - 1][j][k]) / (hx * hx);
    result += (j == 0 || j == curr_block_dims[1] - 1) ? 0 : (curr_v[i][j + 1][k] - 2 * curr_v[i][j][k] + curr_v[i][j - 1][k]) / (hy * hy);
    result += (k == 0 || k == curr_block_dims[2] - 1) ? 0 : (curr_v[i][j][k + 1] - 2 * curr_v[i][j][k] + curr_v[i][j][k - 1]) / (hz * hz);

    return result;
}

double laplassian(
        int i, int j, int k,
        const Grid<double> &curr_v,
        const Grid<double> &prev_v,
        const Grid<double> &next_v,
        const int curr_block_dims[3],
        double hx, double hy, double hz
) {
    double result = 0;

    if (i == 0) {
        result += (curr_v[i + 1][j][k] - 2 * curr_v[i][j][k] + prev_v[0][j][k]) / (hx * hx);
    } else if (i == curr_block_dims[0] - 1) {
        result += (next_v[0][j][k] - 2 * curr_v[i][j][k] + curr_v[i - 1][j][k]) / (hx * hx);
    } else {
        result += (curr_v[i + 1][j][k] - 2 * curr_v[i][j][k] + curr_v[i - 1][j][k]) / (hx * hx);
    }

    if (j == 0) {
        result += (curr_v[i][j + 1][k] - 2 * curr_v[i][j][k] + prev_v[1][i][k]) / (hy * hy);
    } else if (j == curr_block_dims[1] - 1) {
        result += (next_v[1][i][k] - 2 * curr_v[i][j][k] + curr_v[i][j - 1][k]) / (hy * hy);
    } else {
        result += (curr_v[i][j + 1][k] - 2 * curr_v[i][j][k] + curr_v[i][j - 1][k]) / (hy * hy);
    }

    if (k == 0) {
        result += (curr_v[i][j][k + 1] - 2 * curr_v[i][j][k] + prev_v[2][i][j]) / (hz * hz);
    } else if (k == curr_block_dims[2] - 1) {
        result += (next_v[2][i][j] - 2 * curr_v[i][j][k] + curr_v[i][j][k - 1]) / (hz * hz);
    } else {
        result += (curr_v[i][j][k + 1] - 2 * curr_v[i][j][k] + curr_v[i][j][k - 1]) / (hz * hz);
    }

    return result;
}
// endregion math functions

// region build values
Grid<double> build_initial_prev_values_1(
        const Grid<Point> &points_grid,
        double lx, double ly, double lz,
        double a_t,
        const int curr_block_dims[3]
) {
    Grid<double> result;
    for (int i = 0; i < curr_block_dims[0]; i++) {
        const auto& yz_plane = points_grid[i];
        std::vector<std::vector<double> > yz_plane_values;
        for (int j = 0; j < curr_block_dims[1]; j++) {
            const auto &z_line = yz_plane[j];
            std::vector<double> z_line_values;
            for (int k = 0; k < curr_block_dims[2]; k++) {
                const auto &point = z_line[k];
                z_line_values.push_back(phi(point, lx, ly, lz, a_t));
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}

Grid<double> build_initial_prev_values_2(
        double tau,
        const Grid<double> &curr_v,
        double hx, double hy, double hz,
        const int curr_block_dims[3]
) {
    Grid<double> result;
    for (int i = 0; i < curr_block_dims[0]; i++) {
        std::vector<std::vector<double> > yz_plane_values;
        for (int j = 0; j < curr_block_dims[1]; j++) {
            std::vector<double> z_line_values;
            for (int k = 0; k < curr_block_dims[2]; k++) {
                z_line_values.push_back(curr_v[i][j][k] +
                                                (tau * tau) * 0.5 * plain_laplassian(i, j, k, curr_v, curr_block_dims, hx, hy, hz));
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}

Grid<double> build_analytical_values(
        const Grid<Point> &points_grid,
        double tau,
        double lx, double ly, double lz,
        double a_t,
        const int curr_block_dims[3]
) {
    Grid<double> result;
    for (int i = 0; i < curr_block_dims[0]; i++) {
        const auto &yz_plane = points_grid[i];
        std::vector<std::vector<double> > yz_plane_values;
        for (int j = 0; j < curr_block_dims[1]; j++) {
            const auto &z_line = yz_plane[j];
            std::vector<double> z_line_values;
            for (int k = 0; k < curr_block_dims[2]; k++) {
                const auto &point = z_line[k];
                z_line_values.push_back(u_analytical(point, tau, lx, ly, lz, a_t));
            }
            yz_plane_values.push_back(z_line_values);
        }
        result.push_back(yz_plane_values);
    }
    return result;
}

double get_diff(
        const Grid<double> &values_1,
        const Grid<double> &values_2
) {
    double result = 0;
    int count = 0;
    for (int i = 0; i < values_1.size(); i++) {
        const auto &yz_plane_1 = values_1[i];
        const auto &yz_plane_2 = values_2[i];
        for (int j = 0; j < yz_plane_1.size(); j++) {
            const auto &z_line_1 = yz_plane_1[j];
            const auto &z_line_2 = yz_plane_2[j];
            for (int k = 0; k < z_line_1.size(); k++) {
                result += std::abs(z_line_1[k] - z_line_2[k]);
                count++;
            }
        }
    }
    return result / count;
}
// endregion build values

int main(int argc, char *argv[]) {
    // region init
    double lx, ly, lz;
    double tau;
    long n;
    if (argc == 4) {
        lx = ly = lz = std::stod(argv[1]);
        tau = std::stod(argv[2]);
        n = strtol(argv[3], nullptr, 10);
    } else if (argc == 6) {
        lx = std::stod(argv[1]);
        ly = std::stod(argv[2]);
        lz = std::stod(argv[3]);
        tau = std::stod(argv[4]);
        n = strtol(argv[5], nullptr, 10);
    } else {
        std::cout << "Usage:" << std::endl << "\t./main Lx Ly Lz tau n" << std::endl << "\t./main L tau n" << std::endl;
        return 0;
    }

    double hx = lx / (double) n;
    double hy = ly / (double) n;
    double hz = lz / (double) n;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm;
    int dims[3] = {0, 0, 0}, periods[3] = {0, 0, 0}, coords[3];
    MPI_Dims_create(size, 3, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm);
    MPI_Cart_coords(comm, rank, 3, coords);

    if (rank == 0) {
        std::cout << "dims = " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    }

    int curr_block_dims[3];
    for (int i = 0; i < 3; ++i) {
        curr_block_dims[i] = n / dims[i];
        if (coords[i] == dims[i] - 1) {
            curr_block_dims[i] += n % dims[i];
        }
    }

    double a_t = M_PI * std::sqrt(9 / (lx * lx) + 4 / (ly * ly) + 4 / (lz * lz));

    std::vector<double> axis_points_x = divide_axis(n, hx);
    std::vector<double> axis_points_y = divide_axis(n, hy);
    std::vector<double> axis_points_z = divide_axis(n, hz);
    Grid<Point> points_grid = make_points_grid(axis_points_x, axis_points_y, axis_points_z);
    // endregion init

    double start_time = MPI_Wtime();

    std::vector<Grid<double> > values(MAX_T + 1);

    for (int t = 0; t <= MAX_T; t++) {
        if (t == 0) {
            values[0] = build_initial_prev_values_1(points_grid, lx, ly, lz, a_t, curr_block_dims);
            continue;
        }
        if (t == 1) {
            values[1] = build_initial_prev_values_2(
                    tau,
                    values[0],
                    hx, hy, hz,
                    curr_block_dims
            );
            continue;
        }

        // region send and receive x
        long long max_number_of_values_x = curr_block_dims[1] * curr_block_dims[2];
        auto *send_prev_values_x = (double *) calloc(max_number_of_values_x, sizeof(double));
        auto *send_next_values_x = (double *) calloc(max_number_of_values_x, sizeof(double));
        auto *rcv_prev_values_x = (double *) calloc(max_number_of_values_x, sizeof(double));
        auto *rcv_next_values_x = (double *) calloc(max_number_of_values_x, sizeof(double));

        for (int j = 0; j < curr_block_dims[1]; j++) {
            for (int k = 0; k < curr_block_dims[2]; k++) {
                send_prev_values_x[j * curr_block_dims[2] + k] = values[t - 1][0][j][k];
                send_next_values_x[j * curr_block_dims[2] + k] = values[t - 1][curr_block_dims[0] - 1][j][k];
            }
        }

        int src_prev_rank, dest_prev_rank;
        if (coords[0] != 0) {
            MPI_Cart_shift(comm, 0, -1, &src_prev_rank, &dest_prev_rank);
            MPI_Sendrecv(
                    send_prev_values_x, max_number_of_values_x, MPI_DOUBLE, dest_prev_rank, 0,
                    rcv_prev_values_x, max_number_of_values_x, MPI_DOUBLE, dest_prev_rank, 0,
                    comm, MPI_STATUS_IGNORE
            );
        }

        int src_next_rank, dest_next_rank;
        if (coords[0] != curr_block_dims[0] - 1) {
            MPI_Cart_shift(comm, 0, 1, &src_next_rank, &dest_next_rank);
            MPI_Sendrecv(
                    send_next_values_x, max_number_of_values_x, MPI_DOUBLE, dest_next_rank, 0,
                    rcv_next_values_x, max_number_of_values_x, MPI_DOUBLE, dest_next_rank, 0,
                    comm, MPI_STATUS_IGNORE
            );
        }

        // endregion send and receive x

        // region send and receive y
        long long max_number_of_values_y = curr_block_dims[0] * curr_block_dims[2];
        auto *send_prev_values_y = (double *) calloc(max_number_of_values_y, sizeof(double));
        auto *send_next_values_y = (double *) calloc(max_number_of_values_y, sizeof(double));
        auto *rcv_prev_values_y = (double *) calloc(max_number_of_values_y, sizeof(double));
        auto *rcv_next_values_y = (double *) calloc(max_number_of_values_y, sizeof(double));

        for (int i = 0; i < curr_block_dims[0]; i++) {
            for (int k = 0; k < curr_block_dims[2]; k++) {
                send_prev_values_y[i * curr_block_dims[2] + k] = values[t - 1][i][0][k];
                send_next_values_y[i * curr_block_dims[2] + k] = values[t - 1][i][curr_block_dims[1] - 1][k];
            }
        }

        if (coords[1] != 0) {
            MPI_Cart_shift(comm, 1, -1, &src_prev_rank, &dest_prev_rank);
            MPI_Sendrecv(
                    send_prev_values_y, max_number_of_values_x, MPI_DOUBLE, dest_prev_rank, 0,
                    rcv_prev_values_y, max_number_of_values_x, MPI_DOUBLE, dest_prev_rank, 0,
                    comm, MPI_STATUS_IGNORE
            );
        }

        if (coords[1] != curr_block_dims[1] - 1) {
            MPI_Cart_shift(comm, 1, 1, &src_next_rank, &dest_next_rank);
            MPI_Sendrecv(
                    send_next_values_y, max_number_of_values_x, MPI_DOUBLE, dest_next_rank, 0,
                    rcv_next_values_y, max_number_of_values_x, MPI_DOUBLE, dest_next_rank, 0,
                    comm, MPI_STATUS_IGNORE
            );
        }
        // endregion send and receive y

        // region send and receive z
        long long max_number_of_values_z = curr_block_dims[0] * curr_block_dims[1];
        auto *send_prev_values_z = (double *) calloc(max_number_of_values_z, sizeof(double));
        auto *send_next_values_z = (double *) calloc(max_number_of_values_z, sizeof(double));
        auto *rcv_prev_values_z = (double *) calloc(max_number_of_values_z, sizeof(double));
        auto *rcv_next_values_z = (double *) calloc(max_number_of_values_z, sizeof(double));

        for (int i = 0; i < curr_block_dims[0]; i++) {
            for (int j = 0; j < curr_block_dims[1]; j++) {
                send_prev_values_z[i * curr_block_dims[1] + j] = values[t - 1][i][j][0];
                send_next_values_z[i * curr_block_dims[1] + j] = values[t - 1][i][j][curr_block_dims[2] - 1];
            }
        }

        if (coords[2] != 0) {
            MPI_Cart_shift(comm, 2, -1, &src_prev_rank, &dest_prev_rank);
            MPI_Sendrecv(
                    send_prev_values_z, max_number_of_values_x, MPI_DOUBLE, dest_prev_rank, 0,
                    rcv_prev_values_z, max_number_of_values_x, MPI_DOUBLE, dest_prev_rank, 0,
                    comm, MPI_STATUS_IGNORE
            );
        }

        if (coords[2] != curr_block_dims[2] - 1) {
            MPI_Cart_shift(comm, 2, 1, &src_next_rank, &dest_next_rank);
            MPI_Sendrecv(
                    send_next_values_z, max_number_of_values_x, MPI_DOUBLE, dest_next_rank, 0,
                    rcv_next_values_z, max_number_of_values_x, MPI_DOUBLE, dest_next_rank, 0,
                    comm, MPI_STATUS_IGNORE
            );
        }
        // endregion send and receive z

        // region convert array to vector
        Grid<double> prev_values, next_values;

        std::vector<std::vector<double> > prev_x_values, next_x_values;
        for (int j = 0; j < curr_block_dims[1]; j++) {
            std::vector<double> prev_row;
            std::vector<double> next_row;
            for (int k = 0; k < curr_block_dims[2]; k++) {
                prev_row.push_back(rcv_prev_values_x[j * curr_block_dims[2] + k]);
                next_row.push_back(rcv_next_values_x[j * curr_block_dims[2] + k]);
            }
            prev_x_values.push_back(prev_row);
            next_x_values.push_back(next_row);
        }
        prev_values.push_back(prev_x_values);
        next_values.push_back(next_x_values);

        std::vector<std::vector<double> > prev_y_values, next_y_values;
        for (int i = 0; i < curr_block_dims[0]; i++) {
            std::vector<double> prev_row;
            std::vector<double> next_row;
            for (int k = 0; k < curr_block_dims[2]; k++) {
                prev_row.push_back(rcv_prev_values_y[i * curr_block_dims[2] + k]);
                next_row.push_back(rcv_next_values_y[i * curr_block_dims[2] + k]);
            }
            prev_y_values.push_back(prev_row);
            next_y_values.push_back(next_row);
        }
        prev_values.push_back(prev_y_values);
        next_values.push_back(next_y_values);

        std::vector<std::vector<double> > prev_z_values, next_z_values;
        for (int i = 0; i < curr_block_dims[0]; i++) {
            std::vector<double> prev_row;
            std::vector<double> next_row;
            for (int j = 0; j < curr_block_dims[1]; j++) {
                prev_row.push_back(rcv_prev_values_z[i * curr_block_dims[1] + j]);
                next_row.push_back(rcv_next_values_z[i * curr_block_dims[1] + j]);
            }
            prev_z_values.push_back(prev_row);
            next_z_values.push_back(next_row);
        }
        prev_values.push_back(prev_z_values);
        next_values.push_back(next_z_values);
        // endregion convert array to vector

        Grid<double> result;
        for (int i = 0; i < curr_block_dims[0]; i++) {
            const auto &yz_plane = points_grid[i];
            std::vector<std::vector<double> > yz_plane_values;
            for (int j = 0; j < curr_block_dims[1]; j++) {
                const auto &z_line = yz_plane[j];
                std::vector<double> z_line_values;
                for (int k = 0; k < curr_block_dims[2]; k++) {
                    auto laplas = laplassian(i, j, k, values[t - 1], prev_values, next_values, curr_block_dims, hx, hy, hz);
                    z_line_values.push_back(
                            (tau * tau) *
                            laplas
                            + 2 * values[t - 1][i][j][k]
                            - values[t - 2][i][j][k]
                    );
                }
                yz_plane_values.push_back(z_line_values);
            }
            result.push_back(yz_plane_values);
        }

        values[t] = result;
    }

    double end_time = MPI_Wtime() - start_time;
    double max_time;
    MPI_Reduce(&end_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Time: " << max_time << std::endl;
    }

    std::vector<Grid<double> > analytical_values_vector(MAX_T + 1);
    for (int t = 0; t <= MAX_T; t++) {
        analytical_values_vector[t] = build_analytical_values(points_grid, tau * t, lx, ly, lz, a_t, curr_block_dims);
    }

    double diff = get_diff(values[MAX_T], analytical_values_vector[MAX_T]);
    double max_diff;
    MPI_Reduce(&diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (!rank) {
        std::cout << "Diff: " << max_diff << std::endl;
    }

    MPI_Finalize();

    return 0;
}
