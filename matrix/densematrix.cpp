#include "dsla.hpp"

#include "densematrix.h"
#include "imatrix.h"

#include <cstring>

namespace DSLA{

  DenseMatrix::DenseMatrix(const size_t nrow, const size_t ncol, double* buf)
    : IMatrix(MatrixType::Dense, buf, nrow, ncol){
      if(buf == nullptr)
        _buffer = new double[nrow*ncol];
  };
  
  DenseMatrix::DenseMatrix(const DenseMatrix& rhs)
    : IMatrix(MatrixType::Dense, nullptr, rhs._nrow, rhs._ncol) {
    _nrow = rhs._nrow;
    _ncol = rhs._ncol;
    if(_buffer == nullptr){
      _buffer = new double[_nrow * _ncol];
    }
    memcpy(_buffer, rhs._buffer, sizeof(double) * _nrow*_ncol);
  }


  DenseMatrix::~DenseMatrix(){
    delete[] _buffer;
  }


  void DenseMatrix::print() const{
    for(size_t i=0;i<_nrow;++i){
      for(size_t j=0;j<_ncol;++j){
        printf(" %.3f ", _buffer[i*_nrow + j]);
      }
      printf("\n");
    }
  }


  void DenseMatrix::clear(){
    if(_buffer != nullptr){
      delete[] _buffer;
      _buffer = nullptr;
    }
  }


  void DenseMatrix::zero(){
    const auto dim = _nrow * _ncol;
    std::fill(_buffer, _buffer+dim, 0.0);
  }

  void DenseMatrix::write(const std::string& oFile) const {
    FILE* fd = fopen(oFile.c_str(),"wb");
    if (fd == NULL){
      printf("Filename: '%s'\n",oFile.c_str());
      die("Could not open file for writing!");
    }

    if(fwrite(&_ncol,sizeof(size_t),1,fd)      != 1 ||
       fwrite(&_nrow,sizeof(size_t),1,fd)      != 1)
       die("Could not write dense dims");

    if(fwrite(&_buffer[0],_nrow*_ncol*sizeof(double),1,fd) != 1)
      die("Could not save dense matrix");
  }

  void DenseMatrix::read(const std::string& iFile){
    clear(); 
    FILE* fd = fopen(iFile.c_str(),"r");
    if (fd == NULL){
      printf("Filename: '%s'\n",iFile.c_str());
      die("Could not open file for reading!");
    }
    if(fread(&_ncol,sizeof(size_t),1,fd)      != 1 ||
       fread(&_nrow,sizeof(size_t),1,fd)      != 1)
       die("Read error dense matrix dims");

    if(_buffer == nullptr){
      _buffer = new double[_nrow*_ncol];
    }else{
      die("Attempting to alloc existing matrix");
    }

    if(fread(&_buffer[0],_nrow*_ncol*sizeof(double),1,fd) != 1)
      die("Read error dense matrix");
  }

  template <typename T>
  void DenseMatrix::copy(const T& rhs){
    if(rhs.getType() == MatrixType::BCSR){

    }else if(rhs.getType() == MatrixType::Dense){
      //clear();


    }else if(rhs.getType() == MatrixType::MDBCSR){

    }else{
      //exception
    }

  }


}
