#include <omp.h>
#include <iostream>

int main() {
    int num_threads = omp_get_max_threads();
    std::cout << "Number of OpenMP threads: " << num_threads << std::endl;
    return 0;
}