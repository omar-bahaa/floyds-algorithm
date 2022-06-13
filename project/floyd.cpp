#include <omp.h>
#include <vector>
#include <cstdint>
#include <iostream>
#include <algorithm> 
#include <math.h>

#include "../parallelprogrammingbook/include/hpc_helpers.hpp"

#define inf UINT64_MAX

template <typename value_t,
          typename index_t>
void floyd(std::vector<value_t>& A,     // adjacency matrix n*n
           std::vector<value_t>& Dk,    // result
           index_t n){                  // size of A
    
    Dk = A;
    for (index_t k = 1; k < n; k++)
        for (index_t i = 0; i < n; i++)
            for (index_t j = 0; j < n; j++)
                Dk[i*n + j] = std::min(Dk[i*n + j], (Dk[i*n + k] + Dk[k*n + j]));
       
}
  
int main(int argc, char* argv[]){
    uint64_t n = 4;
    std::vector<uint64_t> A = {0, 3, inf, 5, 2, 0, inf, 4, inf, 1, 0, inf, inf, inf, 2, 0};
    std::vector<uint64_t> Dk(n*n);
    floyd(A, Dk, n);
    for (uint64_t i = 0; i < n; i++)
    {
        for (uint64_t j = 0; j < n; j++)
        {
            std::cout << Dk[i*n + j] << " , ";
        }
        std::cout << std::endl;
    }
    
}

