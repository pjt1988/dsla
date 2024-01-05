#include "dsla.hpp"

#include "bcsrmatrix.h"
#include "imatrix.h"
#include "miscfunctions.h"
#include "settings.h"

#include <assert.h>
#include <cstring>
#include <limits>
#include <math.h>
#include <omp.h>


namespace DSLA{

    void BCSRMatrix::scale(const double fac){
        const double fac2{fac*fac};
        const auto thresh{Settings::Instance()->getSparsityThresh()};
        const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};
        const size_t dim{_blockSize * _blockSize};
        
        if(fabs(fac2) < 1.0) [[unlikely]] {
            std::vector<std::vector<size_t> > colIndxLoc(_nbrY);        
    
        #pragma omp parallel for num_threads(nt)
            for(auto i=0ul;i<_nbrY;++i){
                const auto rIndx = i*_nbrY;
                for(auto j=_rowPtr[i];j<_rowPtr[i+1];++j){
                    const auto indx = rIndx + _colIndx[j];
                    if(_bNorms[indx]*fac2 < thresh ) [[unlikely]] {
                        _bNorms[indx] = 0.0;
                        delete[] _buffer[indx];
                        _buffer[indx] = nullptr;
                    }else [[likely]] {
                        _bNorms[indx] *= fac2;
                        //dscale _buffer[rIndx+col] * fac;
                        _dscale(_buffer[indx],dim,fac);
                        colIndxLoc[i].emplace_back(_colIndx[j]);
                    }
                }
            }
            _rowPtr[0]=0;
            _colIndx.clear();
            _colIndx.reserve(_nbrX*_nbrY);
            for(size_t i=0;i<_nbrY;++i){
                _colIndx.insert(_colIndx.end(), colIndxLoc[i].begin(), colIndxLoc[i].end());
                _rowPtr[i+1] = _colIndx.size();
            }
        }else [[likely]] {
        #pragma omp parallel for num_threads(nt)
            for(auto i=0ul;i<_nbrY;++i){
                const auto rIndx = i*_nbrY;
                for(auto j=_rowPtr[i];j<_rowPtr[i+1];++j){
                    const auto indx = rIndx + _colIndx[j];
                    _bNorms[indx] *= fac2;
                    _dscale(_buffer[indx], dim, fac);
                }
            }
        }
    }

    void BCSRMatrix::add(const BCSRMatrix& rhs){
        assert(equalDims(rhs));

        const auto thresh{Settings::Instance()->getSparsityThresh()};
        const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};
        const size_t dim{_blockSize * _blockSize};
        const size_t bDim{sizeof(double) * dim};

    #pragma omp parallel for num_threads(nt)
        for(auto i=0ul;i<_nbrY;++i){
            const auto rIndx = i*_nbrY;
            for(auto j=rhs._rowPtr[i];j<rhs._rowPtr[i+1];++j){
                const auto indx = rIndx + rhs._colIndx[j];
                if(this->_buffer[indx] == nullptr){
                    this->_buffer[indx] = new double[dim];
                    memcpy(this->_buffer[indx], rhs._buffer[indx], bDim);
                    this->_bNorms[indx] = rhs._bNorms[indx];
                }else{
                    this->_bNorms[indx] = _daddS(this->_buffer[indx], rhs._buffer[indx], dim);
                    if(this->_bNorms[indx] < thresh) [[unlikely]] {
                        delete[] this->_buffer[indx];
                        this->_buffer[indx] = nullptr;
                        this->_bNorms[indx] = 0.0;
                    }
                }
            }
        }
        rebuildIndxFromNorms();
    }

    void BCSRMatrix::sub(const BCSRMatrix& rhs){
        assert(equalDims(rhs));

        const auto thresh = Settings::Instance()->getSparsityThresh();
        const int nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;
        const size_t dim = _blockSize * _blockSize;
        const size_t bDim = sizeof(double) * dim;

    #pragma omp parallel for num_threads(nt)
        for(auto i=0ul;i<_nbrY;++i){
            const auto rIndx = i*_nbrY;
            for(auto j=rhs._rowPtr[i];j<rhs._rowPtr[i+1];++j){
                const auto indx = rIndx + rhs._colIndx[j];
                if(this->_buffer[indx] == nullptr){
                    this->_buffer[indx] = new double[dim];
                    memcpy(this->_buffer[indx], rhs._buffer[indx], bDim);
                    _dscale(this->_buffer[indx],dim,-1.0);
                    this->_bNorms[indx] = rhs._bNorms[indx];
                }else{
                    this->_bNorms[indx] = _dsubS(this->_buffer[indx], rhs._buffer[indx], dim);
                    if(this->_bNorms[indx] < thresh) [[unlikely]] {
                        delete[] this->_buffer[indx];
                        this->_buffer[indx] = nullptr;
                        this->_bNorms[indx] = 0.0;
                    }
                }
            }
        }
        rebuildIndxFromNorms();
    }

    std::pair<double,double> BCSRMatrix::gershgorinEstimate() const {
        //the matrix must be at least square - don't bother checking symmetry...
        assert(_nbrX == _nbrY && _nrow == _ncol);
        double epsn{0};
        double eps0{0};

        //const int nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;

        double min0{std::numeric_limits<double>::max()};
        double min1{std::numeric_limits<double>::max()};
        
        double max0{std::numeric_limits<double>::min()};
        double max1{std::numeric_limits<double>::min()};

        std::vector<double> diags(_blockSize);
        std::vector<double> sums(_blockSize);

        for(auto i=0ul;i<_nbrY;++i){
            std::fill(sums.begin(), sums.end(), 0.0);
            for(auto j=0ul;j<_blockSize;++j){
                assert(_buffer[i*_nbrY+i] != nullptr);
                for(auto k=0ul;k<j;++k) sums[k] += fabs(_buffer[i*_nbrY+i][j*_blockSize + k]);
                for(auto k=j+1;k<_blockSize;++k) sums[k] += fabs(_buffer[i*_nbrY+i][j*_blockSize + k]);

                diags[j] = _buffer[i*_nbrY+i][j*_blockSize+j];
            }


            for(auto j=_rowPtr[i];j<_rowPtr[i+1];++j){
                const auto col = _colIndx[j];
                if(col != i){
                    for(size_t k=0;k<_blockSize;++k){
                        for(size_t t=0;t<_blockSize;++t) sums[k] += fabs(_buffer[i*_nbrY+col][k*_blockSize+t]);
                    }
                }
            }

            for(auto j=0ul;j<_blockSize;++j){
                if ((max0 + max1) < (diags[j] + sums[j])){
                    max0 = diags[j];
                    max1 = sums[j];
                    epsn = max0 + max1;
                }
                if ((min0 - min1) > (diags[j] - sums[j])){
                    min0 = diags[j];
                    min1 = sums[j];
                    eps0 = min0 - min1;
                }
            }
        }
        return std::make_pair(eps0,epsn);
    }


std::pair<double,double> BCSRMatrix::gershgorinEstimate_old() const {
    double* diags = new double[_blockSize];
    double* sums  = new double[_blockSize];

    double eps0 = 0.0;
    double epsn = 0.0;
    double disc_min[2] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    double disc_max[2] = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
    for(size_t i=0;i<_nbrY;i++){
      size_t dim_row = _blockSize;
      // diags...
      double* Aact = _buffer[i+i*_nbrY];
      if (Aact == NULL){
        printf("Gershgorin..diag element of row i=%lu is zero!!\n", i);
        die("zero diag-block?!?");
      }    
      for(size_t j=0;j<dim_row;++j){
        sums[j]  = 0.0;
        diags[j] = Aact[j+j*dim_row];
      }    
      for(size_t j=0;j<dim_row;++j){ // should be symmetric, so...
        for(size_t k=0;k<j;++k)         sums[k] += fabs(Aact[j*dim_row+k]);
        for(size_t k=j+1;k<dim_row;++k) sums[k] += fabs(Aact[j*dim_row+k]);
      }    

      for(size_t t=_rowPtr[i]; t<_rowPtr[i+1]; ++t){
        size_t A_bcol = _colIndx[t];
        if(A_bcol != i){
          Aact = _buffer[A_bcol + i*_nbrY];
          for(size_t j=0;j<dim_row;++j){
            for(size_t k=0;k<dim_row;++k) sums[j] += fabs(Aact[k + j*dim_row]);
            //Aact += dim_row; //talk about pointer arithmetic..
          }
        }
      }    

      for(size_t j=0;j<dim_row;++j){ // should be symmetric, so...
        if ((disc_max[0] + disc_max[1]) < (diags[j] + sums[j])){
          disc_max[0] = diags[j];
          disc_max[1] = sums[j];
          epsn = disc_max[0] + disc_max[1];
        }
        if ((disc_min[0] - disc_min[1]) > (diags[j] - sums[j])){
          disc_min[0] = diags[j];
          disc_min[1] = sums[j];
          eps0 = disc_min[0] - disc_min[1];
        }
      }    
    }    

    delete[] sums;
    delete[] diags;

    return std::make_pair(eps0,epsn);




    }


}

