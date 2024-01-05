#include "settings.h"
#include "enums.h"

namespace DSLA{

  const std::string enumAsString(const MatrixType& type){
    using enum MatrixType;
    switch(type){
      case Dense:
        return "Dense";
      case BCSR:
        return "BCSR";
      case MDBCSR:
        return "Memory-Dense BCSR";
      default:
        return "undefined";
    }
  }
  Settings::Settings() : _matrixType(MatrixType::Dense), _sparseThresh(1e-6), _preScreenInner(1e-10), _blockSize(100){ }

  std::shared_ptr<Settings> Settings::Instance(){
    if(_instance == nullptr){
      _instance = std::make_shared<Settings>();
    }
    return _instance;
  }

  std::shared_ptr<Settings> Settings::_instance = nullptr;


}



