// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ROOTDATAIO_HPP_
#define COMPWA_DATA_ROOTDATAIO_HPP_

#include <memory>
#include <string>

#include "Core/Event.hpp"

class TTree;

namespace ComPWA {
namespace Data {
struct DataSet;
namespace Root {

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
             std::size_t NumberEventsToProcess_ = -1);

  /// @param InpufFilePath Input file(s); can take wildcards, because it uses [`TChain::Add`](https://root.cern.ch/doc/master/classTChain.html).
  std::vector<ComPWA::Event> readData(const std::string &InputFilePath) const;

  void writeData(const std::vector<ComPWA::Event> &Events,
                 const std::string &OutputFilePath) const;
};

} // namespace Root
} // namespace Data
} // namespace ComPWA

#endif
