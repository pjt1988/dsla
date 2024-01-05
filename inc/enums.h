#ifndef ENUMS_H_
#define ENUMS_H_

#include <string>

namespace DSLA{

  enum class MatrixType {Dense, BCSR, MDBCSR};

  const std::string enumAsString(const MatrixType& type);

  enum class Timings{
    BCSR_ADD, BCSR_SUB, BCSR_COPY, BCSR_SCALE, /*...etc*/
    DENSE_ADD, DENSE_SUB, DENSE_COPY, DENSE_SCALE,

    MISC

    };

}


#endif
