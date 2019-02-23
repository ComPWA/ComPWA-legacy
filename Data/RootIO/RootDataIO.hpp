// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ROOTDATAIO_HPP_
#define COMPWA_DATA_ROOTDATAIO_HPP_

#include <memory>
#include <string>

class TTree;

namespace ComPWA {
namespace Data {

class DataSet;

///
/// \class RootDataIO
/// Class for reading/writing Physics Events from/to ROOT files.
///
class RootDataIO {
  std::string TreeName;
  int NumberEventsToProcess;

public:
  /// \param TreeName_	Name of tree in input or output file
  /// \param NumberEventsToProcess_	-1 processes all events
  RootDataIO(const std::string TreeName_ = "data",
             int NumberEventsToProcess_ = -1);

  std::shared_ptr<DataSet> readData(const std::string &InputFilePath) const;

  void writeData(std::shared_ptr<const DataSet> DataSample,
                 const std::string &OutputFilePath) const;
};

} // namespace Data
} // namespace ComPWA

#endif
