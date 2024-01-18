#ifndef BCSRMATRIX_H_
#define BCSRMATRIX_H_

#include "imatrix.h"
#include "densematrix.h"
#include "enums.h"

#include <set>
#include <string>
#include <vector>

namespace DSLA{

  class BCSRMatrix : public IMatrix<double**> {
  
  public:
    BCSRMatrix() = delete;
    explicit BCSRMatrix(const size_t row, const size_t col, const size_t blockSize = 0, double** buf = nullptr, const size_t nbrX = 0, const size_t nbrY = 0); 
    BCSRMatrix(const BCSRMatrix&);
    BCSRMatrix(BCSRMatrix&&);
    ~BCSRMatrix();
  
    constexpr const std::vector<size_t>& getRowPtr() const {return _rowPtr;};
    constexpr std::vector<size_t>& getRowPtr() {return _rowPtr;};
    
    constexpr const std::vector<size_t>& getColIndx() const {return _colIndx;};
    constexpr std::vector<size_t>& getColIndx() {return _colIndx;};
  
    constexpr const size_t getNBRX() const {return _nbrX;};
    constexpr const size_t getNBRY() const {return _nbrY;};
    constexpr const size_t getBlockSize() const {return _blockSize;};

    void setNBRX(const size_t nbr){ _nbrX = nbr;};
    void setNBRY(const size_t nbr){ _nbrY = nbr;};
    void setBlocksize(const size_t bs){ _blockSize = bs;};

    double norm(const bool recalculate=false) const;

    //data managment
    template <typename T>
    void copy(const T& rhs);
    void copyFromDense(const DenseMatrix& rhs);
    void copyFromBCSR(const BCSRMatrix& rhs);
    void copyFromBCSR2(const BCSRMatrix& rhs);
    void generateBlocking(const size_t blockSize, const double* mat, const bool screen = true);
    void recompress();
    void clear();
    void zero();

    //IO stuff
    void print() const;
    void occupancy() const;
    void savePixmapNorms(const std::string& oFile, const size_t bs) const {};
    void savePixmap(const std::string& oFile) const {};
    void read(const std::string& iFile);
    void write(const std::string& oFile) const;

    //linalg
    void scale(const double fac);
    void add(const BCSRMatrix& rhs);
    void sub(const BCSRMatrix& rhs);
    void mult(BCSRMatrix& A, BCSRMatrix& B, const MatMultMode mode = MatMultMode::NN, const double alpha = 1.0, const double beta = 0.0) {} // this = alpha * A * B + beta * this
    void transpose_in_place();
    double trace() const;
    void add_to_diag(const double val, const bool rebuildVecs=false); //performs A = A + c*I
    void add_to_diag(const std::vector<double>& vals, const bool rebuildVecs=false); //performs A = A + vec*I
    std::pair<double,double> gershgorinEstimate() const;

  
  private:
    std::vector<size_t> _rowPtr;
    std::vector<size_t> _colIndx;
    std::vector<double> _bNorms;
    size_t              _nbrX;
    size_t              _nbrY;
    size_t              _blockSize;
    size_t              _id;

    bool equalDims(const BCSRMatrix& rhs) const;
    void freeBuf();
    void rebuildIndxFromNorms(); //rebuild the _rowPtr and _colIndx from _bNorms. Useful when additional blocks get added to an existing matrix
  
  
  
  };


}

#endif


