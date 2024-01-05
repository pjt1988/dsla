#include "matrix.h"


class MatrixFactory {
  public:

    static std::unique_ptr<Matrix> createUMatrix(size_t dim, const std::string& type){
      if(type == "matrix_dense")
        return std::unique_ptr<Matrix>(new DenseMatrix(dim));
      if(type == "matrix_trash")
        return std::unique_ptr<Matrix>(new TrashMatrix(dim));
      return nullptr;
    }
    static Matrix* createMatrix(size_t dim, const std::string& type){
      if(type == "matrix_dense")
        return new DenseMatrix(dim);
      if(type == "matrix_trash")
        return new TrashMatrix(dim);

      return nullptr;
    }

};


