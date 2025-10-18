#include "dsla.hpp"

#include "bcsrmatrix.h"
#include "imatrix.h"
#include "miscfunctions.h"
#include "settings.h"

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <execution>
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


    void BCSRMatrix::transpose_in_place(){

        //std::swap is way too slow...
 
        std::vector<std::vector<size_t> > sigCol(_nbrY);
        std::vector<bool> doit(_nbrX*_nbrY,true);
        double tmp{0.0};
        
        for(size_t i=0;i<_nbrY;++i){
            for(size_t j=_rowPtr[i]; j<_rowPtr[i+1]; ++j){
                const size_t icol = _colIndx[j];
                sigCol[icol].emplace_back(i);
                auto& mat0 = _buffer[icol + i*_nbrY];
                if (i == icol){
                    for(size_t t=0;t<_blockSize;++t){
                        for(size_t jj=0;jj<t;++jj){
                            tmp=mat0[t+jj*_blockSize];
                            mat0[t+jj*_blockSize]=mat0[jj+t*_blockSize];
                            mat0[jj+t*_blockSize]=tmp;
                        }
                    }
                }else if(doit[icol+i*_nbrY] && doit[i+icol*_nbrY]){
                    // check if we did this one already...
                    auto& mat1 = _buffer[i + icol*_nbrY];
                    if (mat1 == nullptr){
                        mat1 = new double[_blockSize*_blockSize];
                        for(size_t t=0;t<_blockSize;++t){
                            for(size_t jj=0;jj<_blockSize;++jj){
                                tmp = mat1[t+jj*_blockSize];
                                mat1[t+jj*_blockSize] = mat0[jj+t*_blockSize];
                                mat0[jj+t*_blockSize] = tmp;
                            }
                        }
                        delete[] mat0;
                        mat0 = nullptr;
                    }else{
                        for(size_t t=0;t<_blockSize;++t){
                            for(size_t jj=0;jj<_blockSize;++jj){
                                tmp = mat0[t+jj*_blockSize];
                                mat0[t+jj*_blockSize] = mat1[jj+t*_blockSize];
                                mat1[jj+t*_blockSize] = tmp;
                            }
                        }
                    }
                    tmp = _bNorms[i+icol*_nbrY];
                    _bNorms[i+icol*_nbrY] = _bNorms[icol+i*_nbrY];
                    _bNorms[icol+i*_nbrY] = tmp;
                }
                doit[icol+i*_nbrY]=false;
                doit[i+icol*_nbrY]=false;
            }
        }
        _rowPtr[0] = 0;
        size_t count = 0;
        for(size_t i=0;i<_nbrY;++i){
            memcpy(&_colIndx[count], &sigCol[i][0], sizeof(size_t) * sigCol[i].size());
            count += sigCol[i].size();
            _rowPtr[i+1] = count;
        }

    }

    double BCSRMatrix::trace() const {
        assert(_nbrX == _nbrY && _nrow == _ncol);

        double sum = 0.0;        
        const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};

        #pragma omp parallel for reduction(+:sum) num_threads(nt)
        for(size_t i=0;i<_nbrY;++i){
            assert(_buffer[i+i*_nbrY]);
            const auto& mat = _buffer[i+i*_nbrY];
            for(size_t j=0;j<_blockSize;++j){
                sum += mat[j+j*_blockSize];
            }
        }
        return sum;
    }

    void BCSRMatrix::add_to_diag(const double v, const bool rebuildVecs){
        const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};

        const double fnIncr = _blockSize * v * v;
        
        #pragma omp parallel for num_threads(nt)
        for(size_t i=0;i<_nbrY;++i){
            assert(_buffer[i+i*_nbrY]);
            const auto& mat = _buffer[i+i*_nbrY];
            double fnSum = fnIncr;
            for(size_t j=0;j<_blockSize;++j){
                fnSum += 2*mat[j+j*_blockSize]*v;
                mat[j+j*_blockSize] += v;      
            }
            _bNorms[i+i*_nbrY] += fnSum;            
        }
        if(rebuildVecs){
            rebuildIndxFromNorms();
        }
    }

    void BCSRMatrix::add_to_diag(const std::vector<double>& v, const bool rebuildIndx){ 
        const int nt{omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY};
        
        #pragma omp parallel for num_threads(nt)
        for(size_t i=0;i<_nbrY;++i){
            assert(_buffer[i+i*_nbrY]);
            const auto& mat = _buffer[i+i*_nbrY];
            double fnSum = 0.0;
            for(size_t j=0;j<_blockSize && i*_nbrY +j < v.size();++j){
                fnSum += (2*mat[j+j*_blockSize] + v[i*_nbrY + j])*v[i*_nbrY + j];
                mat[j+j*_blockSize] += v[i*_nbrY + j];
            }
            _bNorms[i+i*_nbrY] += fnSum;            
        }
        if(rebuildIndx){
            rebuildIndxFromNorms();
        }
    }

    void BCSRMatrix::mult_old(BCSRMatrix& A, BCSRMatrix& B, const MatMultMode mode, const double alpha, const double beta){
        using enum MatMultMode;
        switch(mode){
            case NN:
                break;
            case NT:
                B.transpose_in_place();
                break;
            case TN:
                A.transpose_in_place();
                break;
            case TT:
                A.transpose_in_place();
                B.transpose_in_place();
                break;
        }

        const auto thresh = Settings::Instance()->getSparsityThresh();
        const auto threshInner = Settings::Instance()->getPreScreeningInnerThresh();
        const int nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;
        const double beta2{beta*beta};
        const double alpha2{alpha*alpha};
        const size_t dim = _blockSize * _blockSize;

        const auto& arp = A._rowPtr;
        const auto& aci = A._colIndx;
        const auto& afn = A._bNorms;

        const auto& brp = B._rowPtr;
        const auto& bci = B._colIndx;
        const auto& bfn = B._bNorms;

        std::set<size_t>* cols = new std::set<size_t>[_nbrY];
        if(beta == 0.0){
            this->clear();
        }else{
            for(size_t i=0;i<_nbrY;++i){
                for(size_t j=_rowPtr[i];j<_rowPtr[i+1];++j){
                    const size_t col = _colIndx[j];
                    cols[i].emplace_hint(cols[i].end(),col);                   
                    if(beta != 1.0){
                        _bNorms[col + i*_nbrY] *= beta2;
                        _dscale(_buffer[i*_nbrY+col],dim,beta);
                    }
                }
            }
        }

        //symbolic for outer screening
        for(size_t i=0;i<_nbrY;++i){
            for(size_t j=arp[i];j<arp[i+1];++j){
                size_t kk=aci[j];
                for(size_t k=brp[kk];k<brp[kk+1];++k){
                    size_t jj=bci[k];
                    const double prod = alpha * afn[i*_nbrY+kk] * bfn[kk*_nbrY+jj];
                    if(prod > threshInner){
                        if (_buffer[i*_nbrY+jj] == nullptr){ //making sure the block exists..
                            cols[i].insert(cols[i].end(), jj);
                            _buffer[i*_nbrY + jj] = new double[dim];
                            std::fill(_buffer[i*_nbrY+jj],_buffer[i*_nbrY+jj]+dim,0.0);
                        }
                        _bNorms[i*_nbrY+jj] += prod;
                    }
                }
            }
        }


        #pragma omp parallel for schedule(dynamic) num_threads(nt)
        for (size_t i=0; i<_nbrY ; ++i) { 
            for(size_t j=arp[i]; j<arp[i+1]; ++j){
                const size_t kk = aci[j];
                const double prescreen_loc = threshInner / afn[i*_nbrY+kk] / alpha2;

                for(size_t k = brp[kk]; k<brp[kk+1]; ++k){
                    const size_t jj = bci[k]; 
                    if(_bNorms[i*_nbrY+jj] > thresh){
                        if(bfn[jj+kk*_nbrY] > prescreen_loc){                    
                            _dgemm(A._buffer[i*_nbrY+kk],B._buffer[kk*_nbrY+jj],_buffer[i*_nbrY+jj],alpha,_blockSize);
                        }
                    }else{
                        if(_buffer[i*_nbrY+jj] != nullptr){
                            cols[i].erase(jj);
                            delete[] _buffer[i*_nbrY+jj];
                            _buffer[i*_nbrY+jj] = nullptr;
                            _bNorms[i*_nbrY+jj] = 0.0;
                        }
                    }
                }
            }

            for(auto it = cols[i].begin(); it != cols[i].end();){
                const size_t off1 = *it;
                if(_bNorms[i*_nbrY+off1] > thresh){
                    const double dot = fNorm(_buffer[i*_nbrY+off1],dim);
                    if (dot < thresh){
                        delete[] _buffer[i*_nbrY+off1];
                        _buffer[i*_nbrY+off1] = nullptr;
                        _bNorms[i*_nbrY+off1] = 0.0;
                        cols[i].erase(off1);
                    }else{
                        _bNorms[i*_nbrY+off1] = dot;
                    }
                }
                it++;
            }
        }

        _rowPtr[0] = 0;
        _colIndx.clear();
        for(size_t i=0;i<_nbrY;++i){
            _rowPtr[i+1] = _rowPtr[i] + cols[i].size();
            _colIndx.insert(_colIndx.end(), cols[i].begin(), cols[i].end());
        }
        delete[] cols;

        switch(mode){
            case NN:
                break;
            case NT:
                B.transpose_in_place();
                break;
            case TN:
                A.transpose_in_place();
                break;
            case TT:
                A.transpose_in_place();
                B.transpose_in_place();
                break;
        }



    }

    void BCSRMatrix::mult(BCSRMatrix& A, BCSRMatrix& B, const MatMultMode mode, const double alpha, const double beta){
        
        //TODO implement initial check if norms A * B < thresh to skip the whole multiplication
        
        
        
        using enum MatMultMode;
        switch(mode){
            case NN:
                break;
            case NT:
                B.transpose_in_place();
                break;
            case TN:
                A.transpose_in_place();
                break;
            case TT:
                A.transpose_in_place();
                B.transpose_in_place();
                break;
        }

        const auto thresh = Settings::Instance()->getSparsityThresh();
        const auto threshInner = Settings::Instance()->getPreScreeningInnerThresh();
        const int nt = omp_get_num_threads() > (int)_nbrY ? omp_get_num_threads() : (int)_nbrY;
        const double beta2{beta*beta};
        const double alpha2{alpha*alpha};
        const size_t dim = _blockSize * _blockSize;

        const auto& arp = A._rowPtr;
        const auto& aci = A._colIndx;
        const auto& afn = A._bNorms;

        const auto& brp = B._rowPtr;
        const auto& bci = B._colIndx;
        const auto& bfn = B._bNorms;

        if(beta == 0.0){
            this->clear();
        }else if(beta != 1.0){
            for(size_t i=0;i<_nbrY;++i){
                for(size_t j=_rowPtr[i];j<_rowPtr[i+1];++j){
                    const size_t col = _colIndx[j];                
                    _bNorms[col + i*_nbrY] *= beta2;
                    _dscale(_buffer[i*_nbrY+col],dim,beta);
                }
            }
        }

        //symbolic for outer screening
        for(size_t i=0;i<_nbrY;++i){
            for(size_t j=arp[i];j<arp[i+1];++j){
                auto kk=aci[j];
                for(size_t k=brp[kk];k<brp[kk+1];++k){
                    auto jj=bci[k];
                    const double prod = alpha * afn[i*_nbrY+kk] * bfn[kk*_nbrY+jj];
                    if(prod > threshInner){
                        if (_buffer[i*_nbrY+jj] == nullptr){ //making sure the block exists..
                            _buffer[i*_nbrY + jj] = new double[dim];
                            std::fill(_buffer[i*_nbrY+jj],_buffer[i*_nbrY+jj]+dim,0.0);
                        }
                        _bNorms[i*_nbrY+jj] += prod;
                    }
                }
            }
        }


        #pragma omp parallel for schedule(dynamic) num_threads(nt)
        for (size_t i=0; i<_nbrY ; ++i) { 
            for(size_t j=arp[i]; j<arp[i+1]; ++j){
                const auto kk = aci[j];
                const double prescreen_loc = threshInner / afn[i*_nbrY+kk] / alpha2;

                for(size_t k = brp[kk]; k<brp[kk+1]; ++k){
                    const auto jj = bci[k]; 
                    if(_bNorms[i*_nbrY+jj] > thresh){
                        if(bfn[jj+kk*_nbrY] > prescreen_loc){                    
                            _dgemm(A._buffer[i*_nbrY+kk],B._buffer[kk*_nbrY+jj],_buffer[i*_nbrY+jj],alpha,_blockSize);
                        }
                    }else{
                        if(_buffer[i*_nbrY+jj] != nullptr){
                            delete[] _buffer[i*_nbrY+jj];
                            _buffer[i*_nbrY+jj] = nullptr;
                            _bNorms[i*_nbrY+jj] = 0.0;
                        }
                    }
                }
            }

            for(auto off1=0ul;off1<_nbrX;++off1){
                if(_bNorms[i*_nbrY+off1] > thresh){
                    const auto dot = fNorm(_buffer[i*_nbrY+off1],dim);
                    if (dot < thresh){
                        delete[] _buffer[i*_nbrY+off1];
                        _buffer[i*_nbrY+off1] = nullptr;
                        _bNorms[i*_nbrY+off1] = 0.0;
                    }else{
                        _bNorms[i*_nbrY+off1] = dot;
                    }
                }
            }
        }

        rebuildIndxFromNorms();

        switch(mode){
            case NN:
                break;
            case NT:
                B.transpose_in_place();
                break;
            case TN:
                A.transpose_in_place();
                break;
            case TT:
                A.transpose_in_place();
                B.transpose_in_place();
                break;
        }



    }


}

