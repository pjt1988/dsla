#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include "imatrix.h"

namespace DSLA{

  
  class DenseMatrix : public IMatrix<double*> {
  
  public:
    DenseMatrix() = delete;
    explicit DenseMatrix(const size_t nrow, const size_t ncol, double* buf = nullptr); 
    ~DenseMatrix();
    DenseMatrix(const DenseMatrix&);

    void print() const;
    void clear(); 
    void zero();

    void read(const std::string& iFile);
    void write(const std::string& oFile) const;
    void savePixmap(const std::string& oFile) const {};

    template <typename T>
    void copy(const T& rhs);

    //linalg
    void add(const DenseMatrix& rhs);
    void sub(const DenseMatrix& rhs);
    void scale(const double val);
    double trace() const;
    double norm() const;

  private:
    
  
  };
}

#endif
