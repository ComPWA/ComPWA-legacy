// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Reader for data in Root-files.
///
#ifndef _RootReader_HPP
#define _RootReader_HPP

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <utility>

// PWA-Headers
#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"
#include "DataReader/Data.hpp"
#include "DataReader/DataCorrection.hpp"

// Root-Headers
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"

namespace ComPWA {
namespace DataReader {

///
/// \class RootReader
/// Data class for read/write of ROOT files.
/// This class reads event-based data from root-files. It implements the
/// interface of Data.hpp.
///
class RootReader : public Data {

public:
  virtual ~RootReader();

  RootReader();

  /// Constructor read from TTree
  /// \param tr 		TTree to read
  /// \param size 		Number of events to read
  RootReader(TTree *tr, int size = -1);

  /// Read constructor.
  /// 
  /// \param inRootFile Input file name
  /// \param inTreeName Name of tree in input file
  /// \param size 		 Number of events to read in
  RootReader(const std::string inRootFile,
             const std::string inTreeName = "data", int size = -1);

  virtual RootReader *Clone() const;

  /// Create empty clone
  virtual RootReader *EmptyClone() const;

  virtual void WriteData(std::string file = "", std::string trName = "");

protected:
  /// Open ROOT file and set branch addresses to TTree
  void read(TTree* tr, double readSize = -1);

};

} // namespace DataReader
} // namespace ComPWA

#endif
