#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "enums.h"

#include <memory>
#include <math.h>
#include <iostream>

//this is intended as a singleton

namespace DSLA{

  class Settings{
    public:
      Settings();
      Settings(const Settings&) = delete;
      Settings(Settings&&) = delete;
      ~Settings(){}
  
      constexpr void setSparsityThresh(const double st){_sparseThresh = st*st;};
      constexpr double getSparsityThresh() const {return sqrt(_sparseThresh);};
  
      constexpr void setPreScreeningInnerThresh(const double psit){_preScreenInner = psit*psit;};
      constexpr double getPreScreeningInnerThresh() const {return sqrt(_preScreenInner);};
  
      constexpr void setDefaultMatrixType(const MatrixType mt){_matrixType = mt;};
      constexpr MatrixType getDefaultMatrixType() const {return _matrixType;};
  
      constexpr void setDefaultBlockSize(const size_t bs){_blockSize=bs;};
      constexpr size_t getDefaultBlockSize() const {return _blockSize;};
  
      constexpr void printSettings();
      static std::shared_ptr<Settings> Instance();

    private:
        MatrixType _matrixType;
        double     _sparseThresh;
        double     _preScreenInner;
        size_t     _blockSize;
        static std::shared_ptr<Settings> _instance;
  
  };

}

#endif


