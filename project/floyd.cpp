#include <vector>
#include <cstdint>
#include <type_traits>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <cstring>


#include "../parallelprogrammingbook/include/hpc_helpers.hpp"


template <typename value_t,
          typename index_t>
void init(std::vector<value_t>& A,              // A[n][n]
          index_t n){
    
    for(index_t i=index_t(0); i < n; i++){
        for(index_t j=index_t(0); j < n; j++){
            if (i == j){
                A[i*n + j] = value_t(0);
            }else{
            A[i*n + j] = (std::rand() % 100) + 1;
            }
        }
    }
}


template <typename value_t,
          typename index_t>
void floyd_par_omp(std::vector<value_t>& A,      // adjacency matrix n*n
           std::vector<value_t>& Dk,             // result
           index_t n,                            // size of A
           bool parallel){
    
    Dk = A;
    // const int num = std::ceil(sqrt(n));
    index_t i, j;
    for (index_t k = 1; k < n; k++){
        #pragma omp parallel for collapse(2) if(parallel)
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                if ((Dk[i*n + k] * Dk[k*n + j] != 0) && (i != j))   /* See if you can't get a shorter path between i and j by interspacing
                                                                                                    k somewhere along the current path */
                    Dk[i*n + j] = std::min(Dk[i*n + j], (Dk[i*n + k] + Dk[k*n + j]));             
            }
        }
    }
}


int main(int argc, char* argv[]){
    std::srand(1);                              // initalizing random seed
    
    TIMERSTART(overall)
    
    TIMERSTART(alloc)
    uint64_t n = 1024;
    std::vector<uint64_t> A(n*n);
    std::vector<uint64_t> Dk_s(n*n);
    std::vector<uint64_t> Dk_p_omp(n*n);
    TIMERSTOP(alloc)

    TIMERSTART(init)
    init(A, n);
    TIMERSTOP(init)
    
    TIMERSTART(seq)
    floyd_par_omp(A, Dk_s, n, false);
    TIMERSTOP(seq)

    TIMERSTART(parallel_omp)
    floyd_par_omp(A, Dk_p_omp, n, true);
    TIMERSTOP(parallel_omp)

    TIMERSTOP(overall)
    for (uint64_t i = 0; i < n; i++){
        for (uint64_t j = 0; j < n; j++)
            if ((Dk_s[i*n + j] != Dk_p_omp[i*n + j]))
                std::cout << "error at position " << i << ", " << j << std::endl;
    } 
    return 0;
}
