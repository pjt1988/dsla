#include "dsla.hpp"

#include "bcsrmatrix.h"
#include "imatrix.h"
#include "miscfunctions.h"
#include "settings.h"

#include <cassert>
#include <cstring>
#include <execution>
#include <iostream>
#include <numeric>
#include <omp.h>

namespace DSLA{

  BCSRMatrix::BCSRMatrix(const size_t nrow, const size_t ncol, const size_t blockSize, double** buf, const size_t nbrX, const size_t nbrY)
   : IMatrix(MatrixType::BCSR, buf, nrow, ncol)
  , _rowPtr(nbrY+1,0)
  , _colIndx()
  , _bNorms()
  , _nbrX(nbrX)
  , _nbrY(nbrY)
  , _blockSize(blockSize == 0 ? Settings::Instance()->getDefaultBlockSize() : blockSize )
  , _id(0)
  {
    if(buf != nullptr && nbrX > 0 && nbrY > 0){
      _buffer = new double*[nbrX*nbrY];
      _colIndx.reserve(_nbrX*_nbrY);
      for(auto i = 0u; i<nbrX*nbrY;++i)
        _buffer[i] = nullptr;
    }
    _id = reinterpret_cast<size_t>(this);

  }


  BCSRMatrix::BCSRMatrix(const BCSRMatrix& rhs)
  : IMatrix(MatrixType::BCSR, nullptr, rhs._nrow, rhs._ncol){
    copyFromBCSR2(rhs);
  }


  BCSRMatrix::BCSRMatrix(BCSRMatrix&&){

  }

  BCSRMatrix::~BCSRMatrix(){
    freeBuf();
  }


  double BCSRMatrix::norm(const bool recalculate) const {
    if(!recalculate) [[likely]] {
      return std::reduce(std::execution::par, _bNorms.begin(), _bNorms.end());
    }else [[unlikely]]{
      double sum{0.0};
      const size_t bDim{_blockSize * _blockSize};
      const size_t numBlocks{_nbrY * _nbrX};
      const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};

    #pragma omp parallel for reduction(+:sum) num_threads(nt)
      for(size_t i=0;i<numBlocks;++i){
        if(_buffer[i] != nullptr)
          sum += fNorm(_buffer[i], bDim);
       }
       return sum;
    }
    return 0.0;
  }


  void BCSRMatrix::print() const {
    if(_buffer == nullptr)
      return;
    size_t bRow{0};
    size_t bCol{0};
    const size_t totalCol{_nbrX * _blockSize};
    const size_t totalRow{_nbrY * _blockSize};
    
    printf("\n  ");
    for(size_t i=0, rCount=0;i<totalRow;++i,++rCount){
      size_t bRowOld{bRow};
      bRow = i / _blockSize;
      if(bRow != bRowOld){
        rCount = 0;
        std::cout << "\n";
      }
      for(size_t j=0, cCount=0;j<totalCol;++j, ++cCount){
        size_t bColOld{bCol};
        bCol = j / _blockSize;
        if(bCol != bColOld){
          cCount = 0;
          std::cout << "  ";
        }
        if(_buffer[bRow * _nbrY + bCol] == nullptr){
          for(size_t k = 0; k<_blockSize;++k)
            std::cout << " 0.000 ";
          j+=_blockSize;
        }else{
          if(_buffer[bRow * _nbrY + bCol][rCount * _blockSize + cCount] >= 0)
            printf("  %.3f ", _buffer[bRow * _nbrY + bCol][rCount * _blockSize + cCount]);
          else
            printf(" %.3f ", _buffer[bRow * _nbrY + bCol][rCount * _blockSize + cCount]);
        }
      }
      std::cout << "\n";
    }
  }


  bool BCSRMatrix::equalDims(const BCSRMatrix& rhs) const {
    return this->_nbrY == rhs._nbrY &&
           this->_nbrX == rhs._nbrX &&
           this->_blockSize == rhs._blockSize &&
           this->_bNorms.size() == rhs._bNorms.size();
  }

  void BCSRMatrix::freeBuf(){
    if(_buffer != nullptr){
      for(auto i=0u;i<_nbrX*_nbrY;++i){
        if(_buffer[i] != nullptr){
          delete[] _buffer[i];
        }
      }
      delete[] _buffer;
      _buffer = nullptr;
    }
  }


  void BCSRMatrix::rebuildIndxFromNorms(){
    const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};
    std::vector<std::vector<size_t>> colIndxLoc(_nbrY);
    _colIndx.clear();
    _colIndx.reserve(_nbrX*_nbrY);
    _rowPtr[0] = 0;

    #pragma omp parallel for num_threads(nt)
    for(auto i=0ul;i<_nbrY;++i){
      colIndxLoc[i].reserve(_nbrX);
      for(auto j=0ul;j<_nbrX;++j){
        if(_bNorms[i*_nbrY+j] > 0.0){
          colIndxLoc[i].emplace_back(j);
        }
      }
    }

    for(size_t i=0;i<_nbrY;++i){
      _colIndx.insert(_colIndx.end(), colIndxLoc[i].begin(), colIndxLoc[i].end());
      _rowPtr[i+1] = _colIndx.size();
    }
  }


  void BCSRMatrix::generateBlocking(const size_t blockDim, const double* mat, const bool screen){
    _blockSize = blockDim;
    _nbrX = _ncol % blockDim == 0 ? _ncol / blockDim : _ncol / blockDim + 1;
    _nbrY = _nrow % blockDim == 0 ? _nrow / blockDim : _nrow / blockDim + 1;

    printf("nbrx,nbry %lu %lu \n", _nbrX, _nbrY);

    if(_buffer == nullptr){
      _buffer = new double*[_nbrX*_nbrY];
    }

    const auto bDim{_blockSize * _blockSize};
    const auto thresh{Settings::Instance()->getSparsityThresh()};
  
    _rowPtr.resize(_nbrY + 1);
    std::fill(_rowPtr.begin(), _rowPtr.end(), 0);
    _colIndx.clear();
    _colIndx.reserve(_nbrX*_nbrY);
    _bNorms.resize(_nbrX*_nbrY);
    std::fill(_bNorms.begin(), _bNorms.end(), 0.0);
    std::vector<std::vector<size_t> > colIndxLoc(_nbrY);
    const auto nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;

#pragma omp parallel for num_threads(nt)
    for(size_t i=0;i<_nbrY;++i){
      const size_t rIndx{i*_nbrY};
      for(size_t j=0;j<_nbrX;++j){
        _buffer[rIndx + j] = new double[bDim];
        
        memset(_buffer[rIndx + j], 0.0, sizeof(double)*bDim);
        const size_t iShift{i*_blockSize};
        const size_t jShift{j*_blockSize};

        for(size_t ii=0;ii<_blockSize && ii + iShift < _nrow; ++ii){
          for(size_t jj=0;jj<_blockSize && jj + jShift < _ncol; ++jj){
            _buffer[rIndx+j][ii*_blockSize + jj] = mat[(iShift + ii)*_nrow + jShift + jj];
          }
        }
        if(screen){
          const auto fn{fNorm(_buffer[rIndx + j], bDim)};
          //printf("i=%lu, j=%lu, rindx=%lu _buffer=%.5f, fn = %.5f \n", i,j, rIndx,_buffer[rIndx+j][0],fn);
          if(j==i){
            if(fn < thresh){
              printf("fail - i = %lu, j = %lu, fn = %.7f rIndx = %lu bDim = %lu iShift = %lu, jShift= %lu _nbrY=%lu _nbrX=%lu \n", i, j, fn, rIndx, bDim, iShift, jShift, _nbrY, _nbrX);
              printf("nrow=%lu ncol=%lu nbrY=%lu nbrX=%lu \n", _nrow, _ncol, _nbrY, _nbrX);
              //printMatrix(_buffer[rIndx+j], _blockSize, _blockSize);
            }
            //assert(fn > thresh);
          }
          if(fn > thresh){
            _bNorms[rIndx+j] = fn;
            colIndxLoc[i].emplace_back(j);
          }else{
            delete[] _buffer[i*_nbrY+j];
            _buffer[rIndx+j]=nullptr;
          }
        }else{
          colIndxLoc[i].emplace_back(j);
        }
      }
    }

    //printMatrix(_bNorms.data(),_nbrY,_nbrX);

    for(size_t i=0;i<_nbrY;++i){
      _colIndx.insert(_colIndx.end(), colIndxLoc[i].begin(), colIndxLoc[i].end());
      _rowPtr[i+1] = _colIndx.size();
    }
  }


  void BCSRMatrix::recompress(){
    const size_t bDim{_blockSize * _blockSize};
    const auto thresh{Settings::Instance()->getSparsityThresh()};
    const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};

    std::fill(_rowPtr.begin(), _rowPtr.end(), 0);
    _colIndx.clear();
    _colIndx.reserve(_nbrX*_nbrY);
    std::fill(_bNorms.begin(), _bNorms.end(), 0.0);
    std::vector<std::vector<size_t> > colIndxLoc(_nbrY);

  #pragma omp parallel for num_threads(nt)
    for(size_t i=0;i<_nbrY;++i){
      for(size_t j=0;j<_nbrX;++j){
        const auto fn{fNorm(_buffer[i*_nbrY + j], bDim)}; 
        if(fn > thresh) [[likely]] {
          _bNorms[i*_nbrY+j] = fn;
          colIndxLoc[i].emplace_back(j);
        }else [[unlikely]] {
          delete[] _buffer[i*_nbrY+j];
          _buffer[i*_nbrY+j]=nullptr;
        }
      }
    }
    for(size_t i=0;i<_nbrY;++i){
      _colIndx.insert(_colIndx.end(), colIndxLoc[i].begin(), colIndxLoc[i].end());
      _rowPtr[i+1] = _colIndx.size();
    }
  }


  void BCSRMatrix::clear(){
    std::fill(_rowPtr.begin(), _rowPtr.end(), 0);
    _colIndx.clear();
    std::fill(_bNorms.begin(), _bNorms.end(), 0.0);
    freeBuf();
  }


  void BCSRMatrix::zero(){
    clear();
  }


  void BCSRMatrix::occupancy() const {
    const auto thresh{Settings::Instance()->getSparsityThresh()};
    const auto dim{_bNorms.size()};
    auto count1{0};
    for(const auto& it : _bNorms){
      if(it > thresh)
        count1++;
    }

    auto count2{0};
    for(auto i=0ul;i<dim;++i){
      if(_buffer[i] != nullptr){
        count2++;
      }
    }

    std::cout << "Stored norms occ - " << count1 << "/" << dim << " = " << (double) count1/dim;
    std::cout << "  Buffer occ - " << count2 << "/" << dim << " = " << (double) count2/dim << std::endl;
  }


  void BCSRMatrix::read(const std::string& iFile){
    clear(); 
    FILE* fd{fopen(iFile.c_str(),"r")};
    if (fd == NULL){
      printf("Filename: '%s'\n",iFile.c_str());
      die("Could not open file for reading!");
    }
    if(fread(&_ncol,sizeof(size_t),1,fd)      != 1 ||
       fread(&_nrow,sizeof(size_t),1,fd)      != 1 ||
       fread(&_nbrX,sizeof(size_t),1,fd)      != 1 ||
       fread(&_nbrY,sizeof(size_t),1,fd)      != 1 ||
       fread(&_blockSize,sizeof(size_t),1,fd) != 1 ||
       fread(&_id,sizeof(size_t),1,fd)        != 1){
      die("read error dim data");
    }
    _rowPtr.resize(_nbrY+1);
    _bNorms.resize(_nbrY * _nbrY);
    if(_buffer == nullptr){
    _buffer = new double*[_nbrY*_nbrX];
     for(size_t i =0;i<_nbrX*_nbrY;++i)
      _buffer[i]=nullptr;
    }else{
      die("Attempting to alloc existing matrix");
    }

    //read key vectors
    if(fread(&(_rowPtr[0]), _rowPtr.size() * sizeof(size_t), 1,fd)   != 1)
      die("read error rowPtr");
    _colIndx.resize(_rowPtr.back(),0);
    if(_colIndx.size()){
      if(fread(&(_colIndx[0]), _rowPtr.back() * sizeof(size_t), 1,fd) != 1){
        die("read error key vectors");
      }
    }
    //read norms
    if(fread(&(_bNorms[0]), _bNorms.size() * sizeof(double),1,fd) != 1)
      die("write error norms");

    for(size_t i=0;i<_nbrY;++i){
      for(size_t j=_rowPtr[i];j<_rowPtr[i+1];++j){
        const size_t col = _colIndx[j];
        _buffer[i*_nbrY+col] = new double[_blockSize*_blockSize];
        if(fread(&(_buffer[i*_nbrY+col][0]), _blockSize*_blockSize*sizeof(double),1,fd) != 1)
          die("write error for block");
      }
    }
    fclose(fd);
  }


  void BCSRMatrix::write(const std::string& oFile) const {
    FILE* fd{fopen(oFile.c_str(),"wb")};
    if (fd == NULL){
      printf("Filename: '%s'\n",oFile.c_str());
      die("Could not open file for writing!");
    }
    //write basic dimensions
    if(fwrite(&_ncol,sizeof(size_t),1,fd)      != 1 ||
       fwrite(&_nrow,sizeof(size_t),1,fd)      != 1 ||
       fwrite(&_nbrX,sizeof(size_t),1,fd)      != 1 ||
       fwrite(&_nbrY,sizeof(size_t),1,fd)      != 1 ||
       fwrite(&_blockSize,sizeof(size_t),1,fd) != 1 ||
       fwrite(&_id,sizeof(size_t),1,fd)        != 1){
      die("write error dim data");
    }
    //write key vectors
    if(fwrite(&(_rowPtr[0]), _rowPtr.size() * sizeof(size_t), 1,fd)   != 1){
      die("write error rowPtr");
    }
    if(_colIndx.size()){
      if(fwrite(&(_colIndx[0]), _colIndx.size() * sizeof(size_t), 1,fd) != 1){
        die("write error colIndx");
      }
    }
    //write norms
    if(fwrite(&(_bNorms[0]), _bNorms.size() * sizeof(double),1,fd) != 1)
      die("write error norms");

    for(size_t i=0;i<_nbrY;++i){
      for(size_t j=_rowPtr[i];j<_rowPtr[i+1];++j){
        const size_t col{_colIndx[j]};
        if(fwrite(&(_buffer[i*_nbrY+col][0]), _blockSize*_blockSize*sizeof(double),1,fd) != 1)
          die("write error for block");
      }
    }
    fclose(fd);
  }


  template <typename T>
  void BCSRMatrix::copy(const T& rhs){
    if(rhs.getType() == MatrixType::BCSR){

    }else if(rhs.getType() == MatrixType::Dense){
      copyFromDense(rhs);



    }else if(rhs.getType() == MatrixType::MDBCSR){

    }else{
      //exception
    }

  }


  void BCSRMatrix::copyFromDense(const DenseMatrix& rhs){
    clear();
    generateBlocking(_blockSize, rhs.getBuffer(), true);
  }


  void BCSRMatrix::copyFromBCSR(const BCSRMatrix& rhs){
    if(rhs._nbrX != _nbrX || rhs._nbrY != _nbrY || rhs._blockSize != _blockSize)
      clear();

    _nrow = rhs._nrow;
    _ncol = rhs._ncol;
    _nbrX = rhs._nbrX;
    _nbrY = rhs._nbrY;
    _blockSize = rhs._blockSize;
    _rowPtr = rhs._rowPtr;
    _colIndx = rhs._colIndx;
    _bNorms = rhs._bNorms;

    const auto rhsBuffer = rhs.getBuffer();
    const auto dim = _nbrX * _nbrY;
    const auto bDim = _blockSize * _blockSize;
    const auto bDimSize = sizeof(double)*bDim;
    const int nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;

    if(_buffer == nullptr){
      _buffer = new double*[dim];
      for(auto i=0ul;i<dim;++i)
        _buffer[i] = nullptr;
    }    

#pragma omp parallel for num_threads(nt)
    for(auto i=0ul;i<dim;++i){
      if(rhsBuffer[i] != nullptr){
        if(_buffer[i] == nullptr)
          _buffer[i] = new double[bDim];
        memcpy(_buffer[i], rhsBuffer[i], bDimSize);
      }else if(rhsBuffer[i] == nullptr && _buffer[i] != nullptr){
        delete[] _buffer[i];
        _buffer[i] = nullptr;
      }
    }
  }

  //only faster when copying into a zero'd or near zero'd matrix
  void BCSRMatrix::copyFromBCSR2(const BCSRMatrix& rhs){
    freeBuf();

    _nrow = rhs._nrow;
    _ncol = rhs._ncol;
    _nbrX = rhs._nbrX;
    _nbrY = rhs._nbrY;
    _blockSize = rhs._blockSize;
    _rowPtr = rhs._rowPtr;
    _colIndx = rhs._colIndx;
    _bNorms = rhs._bNorms;

    const auto rhsBuffer = rhs.getBuffer();
    const auto dim = _nbrX * _nbrY;
    const auto bDim = _blockSize * _blockSize;
    const auto bDimSize = sizeof(double)*bDim;

    _buffer = new double*[dim];
    for(auto i =0ul;i<dim;++i)
      _buffer[i] = nullptr;

  const int nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;

#pragma omp parallel for num_threads(nt)
    for(auto i=0u;i<_nbrY;++i){
      for(auto j=_rowPtr[i];j<_rowPtr[i+1];++j){
        auto col = _colIndx[j];
        _buffer[i*_nbrY + col] = new double[bDim];
        memcpy(_buffer[i*_nbrY+col], rhsBuffer[i*_nbrY+col], bDimSize);
      }
    }
  }

}
