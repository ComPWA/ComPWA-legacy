// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_DATAIOINTERFACE_HPP_
#define COMPWA_DATA_DATAIOINTERFACE_HPP_

#include "Data/Data.hpp"

namespace ComPWA {
namespace Data {

class DataIOInterface {
public:
  virtual ~DataIOInterface() {}
  /// Write sample to file.
  /// Method is supposed to be implemented by derived classes
  virtual void writeData(std::shared_ptr<ComPWA::Data::Data> Data,
                         const std::string &OutputFilePath) const = 0;

  /// Read sample from file.
  /// Method is supposed to be implemented by derived classes
  virtual std::shared_ptr<ComPWA::Data::Data>
  readData(const std::string &InputFilePath) const = 0;

  // virtual Data *clone() const { return new Data(*this); };

  // virtual Data *emptyClone() const { return new Data(); };
};

} // namespace Data
} // namespace ComPWA

#endif
