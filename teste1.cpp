#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <mpi.h>

using namespace std;

int prime_number(int n, int id, int p);
void timestamp();

int main(int argc, char *argv[]) {
    int id, p, ierr;
    int n_lo = 1, n_hi = 1048576, n_factor = 2;
    int primes, primes_part;
    double wtime, total_time = 0;
    vector<int> local_counts(12, 0);  // To store counts for each core

    ierr = MPI_Init(&argc, &argv);
    if (ierr != 0) {
        cout << "PRIME_MPI - Fatal error! MPI_Init returned nonzero ierr.\n";
        exit(1);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id == 0) {
        timestamp();
        cout << "PRIME_MPI - C++/MPI version\n";
        cout << "An MPI program to count the number of primes.\n";
        cout << "The number of processes is " << p << "\n\n";
        cout << "N\tPi\tTime\tDistribution\n\n";
    }

    for (int n = n_lo; n <= n_hi; n *= n_factor) {
        if (id == 0) wtime = MPI_Wtime();

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        primes_part = prime_number(n, id, p);

        MPI_Reduce(&primes_part, &primes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Gather(&primes_part, 1, MPI_INT, local_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (id == 0) {
            wtime = MPI_Wtime() - wtime;
            total_time += wtime;

            cout << n << "\t" << primes << "\t" << fixed << setprecision(6) << wtime << "\t";
            for (int i = 0; i < p; ++i) {
                cout << local_counts[i] << " ";
            }
            cout << "\n";
        }
    }

    if (id == 0) {
        cout << "\nTotal execution time: " << total_time << " seconds\n";
        cout << "Average time per iteration: " << total_time / log2((n_hi/n_lo) + 1) << " seconds\n";
        timestamp();
    }

    MPI_Finalize();
    return 0;
}

int prime_number(int n, int id, int p) {
    int total = 0;
    for (int i = 2 + id; i <= n; i = i + p) {
        bool is_prime = true;
        for (int j = 2; j <= sqrt(i); j++) {
            if (i % j == 0) {
                is_prime = false;
                break;
            }
        }
        if (is_prime) total++;
    }
    return total;
}

void timestamp() {
    time_t now = time(nullptr);
    cout << ctime(&now);
}