// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Reader for data in Root-files.
///
#ifndef COMPWA_DATAREADER_ROOTREADER_HPP_
#define COMPWA_DATAREADER_ROOTREADER_HPP_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Data/Data.hpp"
#include "Data/DataCorrection.hpp"
#include "Data/DataIOInterface.hpp"
// PWA-Headers
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"

class TTree;

namespace ComPWA {
namespace Data {

///
/// \class RootReader
/// Data class for read/write of ROOT files.
/// This class reads event-based data from root-files. It implements the
/// interface of Data.hpp.
///
class RootReader : public DataIOInterface {
  std::string TreeName;
  int NumberEventsToProcess;

public:
  virtual ~RootReader();

  /// constructor.
  ///
  /// \param TreeName	Name of tree in input or output file
  /// \param size		Number of events to read or write
  RootReader(const std::string TreeName_ = "data",
             int NumberEventsToProcess_ = -1);

  // virtual RootReader *clone() const;

  /// Create empty clone
  // virtual RootReader *emptyClone() const;

  virtual void writeData(std::shared_ptr<ComPWA::Data::Data> Data,
                         const std::string &OutputFilePath) const;

  virtual std::shared_ptr<ComPWA::Data::Data>
  readData(const std::string &InputFilePath) const;
  // LOG(ERROR)
  //    << "Data::writeData() | Base class does not provide functionality"
  //       "to write data to file.";
  //};
};

} // namespace Data
} // namespace ComPWA

#endif
