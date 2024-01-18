#include "dsla.hpp"

#include "densematrix.h"
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

    void DenseMatrix::add(const DenseMatrix& rhs){
        assert(this->_nrow == rhs._nrow && this->_ncol == rhs._ncol);

        #pragma omp parallel for
        for(size_t i=0;i<_nrow*_ncol;++i){
            _buffer[i] += rhs._buffer[i];
        }
    }

    void DenseMatrix::sub(const DenseMatrix& rhs){
        assert(this->_nrow == rhs._nrow && this->_ncol == rhs._ncol);

        #pragma omp parallel for
        for(size_t i=0;i<_nrow*_ncol;++i){
            _buffer[i] -= rhs._buffer[i];
        }
    }

    void DenseMatrix::scale(const double val){
        #pragma omp parallel for
        for(size_t i=0;i<_nrow*_ncol;++i){
            _buffer[i] *= val;
        }
    }

    double DenseMatrix::trace() const {
        assert(_nrow == _ncol);
        double sum = 0.0;

        #pragma omp parallel for reduction(+:sum)
        for(size_t i=0;i<_nrow;++i){
            sum += _buffer[i*_nrow + i];
        }
        return sum;
    }

}