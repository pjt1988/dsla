#ifndef MISC_FUNCTIONS_H_
#define MISC_FUNCTIONS_H_

#include <iostream>
#include <cstddef>

namespace DSLA{

  void generateRandomMatrix(const size_t nrow, const size_t ncol, const double min, const double max, double*& mat, bool dropOff = false);
  
  double fNorm(const double* __restrict__ mat, const size_t dim);

  void _dscale(double* __restrict__ mat, const size_t dim, const double fac);

  void _dadd(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim);
  void _dsub(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim);


  //basic linalg functions returning the Frobenius2 norm of the result
  double _daddS(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim);
  double _dsubS(double* __restrict__ mat, const double* __restrict__ rhs, const size_t dim);


  template<class T>
  void printMatrix(T* mat, const size_t nrow, const size_t ncol){
    printf("\n");
    for(auto i=0ul;i<nrow;++i){
      for(auto j=0ul;j<ncol;++j){
        std::cout << mat[i*nrow+j] << " ";
        //std::cout << setprecision(6) << setw(9) << setfill(' ') << mat[i*nrow+j];
        //printf("%.6f ", mat[i*nrow+j]);
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

#endif
