#include "miscfunctions.h"

#include <algorithm>
#include <execution>
#include <iostream>
#include <numeric>
#include <random>
#include <omp.h>

namespace DSLA{

  void generateRandomMatrix(const size_t nrow, const size_t ncol, const double min, const double max, double*& mat, bool dropOff){
    const auto dim = nrow*ncol;
    if(mat == nullptr)
      mat = new double[dim];

    std::random_device rnd;
    std::default_random_engine eng(rnd());
  
    std::uniform_real_distribution<double> urd(min,max);
  
    std::generate(mat, (mat)+dim, [&]{ return urd(eng);});
    for(auto i=nrow - 10;i<nrow;++i)
        printf("i=%lu, buf[i*nrow+i] = %.3f \n", i, mat[i*nrow+i]);

    if(dropOff){
      for(auto i=0ul;i<nrow;++i){
        for(auto j=0ul;j<ncol;++j){
          if(i>j){
            mat[i*nrow + j] /= (7.50*pow((double)(i-j), 1.55));
            //printf("i=%lu, j=%lu, i*nrow+j = %lu, mat[i*nrow+j]=%.5f\n", i,j,i*nrow+j,mat[i*nrow+j]);
          }
          else if(i<j){
            mat[i*nrow + j] /= (7.50*pow((double)(j-i), 1.55));
            //printf("i=%lu, j=%lu, i*nrow+j = %lu, mat[i*nrow+j]=%.5f\n", i,j,i*nrow+j,mat[i*nrow+j]);
          }
        }
        mat[i*nrow+i] *= 5.0; 
      }
    }

  
  }
  
  double fNorm(const double* __restrict__ mat, const size_t dim){    
    double sum = 0.0;  
  #pragma omp parallel for reduction(+:sum)
    for(size_t i=0;i<dim;++i){
      sum += mat[i]*mat[i];
    }
    return sum;
  }

  void _dscale(double* __restrict__ mat, const size_t dim, const double fac){
    #pragma omp parallel for 
    for(size_t i=0;i<dim;++i){
      mat[i]*=fac;
    }
  }

  void _dadd(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim){
    #pragma omp parallel for
    for(size_t i=0;i<dim;++i){
      mat[i] += rhs[i];
    }
  }

  void _dsub(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim){
    #pragma omp parallel for
    for(size_t i=0;i<dim;++i){
      mat[i] -= rhs[i];
    }
  }

  double _daddS(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for(size_t i=0;i<dim;++i){
      mat[i] += rhs[i];
      sum += mat[i]*mat[i];
    }
    return sum;
  }

  double _dsubS(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for(size_t i=0;i<dim;++i){
      mat[i] -= rhs[i];
      sum += mat[i]*mat[i];
    }
    return sum;
  }

  void printMatrix(double* mat, const size_t nrow, const size_t ncol){
    printf("\n");
    for(auto i=0ul;i<nrow;++i){
      for(auto j=0ul;j<ncol;++j){
        printf("%.6f ", mat[i*nrow+j]);
      }
      printf("\n");
    }
    printf("\n");
  }

}
