#ifndef ENUMS_H_
#define ENUMS_H_

#include <string>

namespace DSLA{

  enum class MatrixType {Dense, BCSR, MDBCSR};

  const std::string enumAsString(const MatrixType& type);

  enum class Timings{
    BCSR_READ, BCSR_WRITE, BCSR_ADD, BCSR_SUB, BCSR_COPY, BCSR_SCALE, BCSR_TRANSPOSE, BCSR_TRACE, BCSR_MULT, /*...etc*/
    DENSE_ADD, DENSE_SUB, DENSE_COPY, DENSE_SCALE, DENSE_MULT,

    MISC

    };

  enum class MatMultMode{
    NN,NT,TN,TT
  };

}


#endif
