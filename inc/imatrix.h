#ifndef IMATRIX_H_
#define IMATRIX_H_

#include "enums.h"

#include <vector>

namespace DSLA{

  template<typename Buf>
  class IMatrix{
    public:
      IMatrix(){};
      virtual ~IMatrix(){};
      IMatrix(const IMatrix&){};
      IMatrix(const MatrixType& mt, Buf bf, size_t nr, size_t nc) : _matrixType(mt), _buffer(bf), _nrow(nr), _ncol(nc) {};
  
      constexpr std::string getTypeAsString() const {return enumAsString(_matrixType);};
      constexpr MatrixType getType() const {return _matrixType;};
      constexpr size_t getNRow() const {return _nrow;};
      constexpr size_t getNCol() const {return _ncol;};
      constexpr const Buf&     getBuffer() const {return _buffer;};
      constexpr Buf&     getBuffer() {return _buffer;};


      virtual void print() const = 0;
      virtual void clear() = 0;
      virtual void zero() = 0;
      void setBuffer(Buf&& t){_buffer = std::move(t);};
      void setNBR(const size_t nbr){}
      void setBlocksize(const size_t bs){}

      virtual void copy(const IMatrix& rhs){};
      virtual void write(const std::string&) const = 0;
      virtual void read(const std::string&) = 0;
      virtual void savePixmap(const std::string&) const = 0;

      

    protected:
      MatrixType     _matrixType;
      Buf            _buffer;
      size_t         _nrow;
      size_t         _ncol;
  
  };

}


#endif
