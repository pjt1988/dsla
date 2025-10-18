#include "dsla.hpp"

#include "densematrix.h"
#include "imatrix.h"

#include <cmath>
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

    _matrixNorm = rhs._matrixNorm;
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

  void DenseMatrix::clear(const bool clearBuffer){
    if(_buffer != nullptr){
      delete[] _buffer;
      _buffer = nullptr;
    }
  }

  void DenseMatrix::zero(){
    const auto dim = _nrow * _ncol;
    std::fill(_buffer, _buffer+dim, 0.0);
    _matrixNorm = 0.0;
  }

  double DenseMatrix::norm() const {
    if(_buffer == nullptr || _nrow == 0 || _ncol == 0) [[unlikely]] {
      return 0.0;
    }

    double sum{0.0};
    for(size_t i=0;i<_nrow;++i){
      for(size_t j=0;j<_ncol;++j){
        sum += _buffer[i*_nrow + j] * _buffer[i*_nrow + j];
      }
    }
    return sum;
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

  void DenseMatrix::savePixmap(const std::string& oFile) const {
        FILE* of = fopen(oFile.c_str(),"w");

    /*
    
      Legend:
      -------
    
      white:     zero
      lightgrey: zero, but part of submatrix
      darkgrey:  smaller than REM_SPARS_THRESH, but larger than sigthr
      black:     larger than REM_SPARS_THRESH
    
    */
      fprintf(of,"/* XPM */\nstatic char * matrix_xpm[] = {\n\"%i %i 13 1\",\n",(int)_ncol,(int)_nrow);
      fprintf(of,"\"a\tc #ff0000\",\n");   // > 1
      fprintf(of,"\"b\tc #ff3300\",\n"); // > e-01
      fprintf(of,"\"c\tc #ff3333\",\n"); // > e-02
      fprintf(of,"\"d\tc #ff6633\",\n"); // > e-03
      fprintf(of,"\"e\tc #ff6666\",\n"); // > e-04
      fprintf(of,"\"f\tc #ff9966\",\n"); // > e-05
      fprintf(of,"\"g\tc #ff9999\",\n"); // > e-06
      fprintf(of,"\"h\tc #ffbb99\",\n"); // > e-07
      fprintf(of,"\"i\tc #ffbbbb\",\n"); // > e-08
      fprintf(of,"\"j\tc #ffeebb\",\n"); // > e-09
      fprintf(of,"\"k\tc #ffeeee\",\n"); // > e-10
      fprintf(of,"\".\tc #ffffff\",\n"); // zero
      fprintf(of,"\"p\tc #bbbbbb\",\n"); // zero, bit part of block...

      for(size_t i=0;i<_nrow;++i){
          for(size_t j=0;j<_ncol;++j){
              double absval = fabs(_buffer[i*_ncol + j]);
              if (absval >= 1.e+00) fprintf(of,"a");     // the less tight one... SPARS_THRESH
              else if (absval >= 1.e-01) fprintf(of,"b");
              else if (absval >= 1.e-02) fprintf(of,"c");
              else if (absval >= 1.e-03) fprintf(of,"d");
              else if (absval >= 1.e-04) fprintf(of,"e");
              else if (absval >= 1.e-05) fprintf(of,"f");
              else if (absval >= 1.e-06) fprintf(of,"g");
              else if (absval >= 1.e-07) fprintf(of,"h");
              else if (absval >= 1.e-08) fprintf(of,"i");
              else if (absval >= 1.e-09) fprintf(of,"j");
              else if (absval >= 1.e-10) fprintf(of,"k");
              else fprintf(of,".");
          }
          fprintf(of,"\n");
      }
      fprintf(of,"};\n");
      fclose(of);     
  }


}
