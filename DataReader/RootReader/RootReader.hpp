//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Reader for data in Root-files.
/*! \class RootReader
 * @file RootReader.hpp
 * This class reads event-based data from root-files. It implements the
 * interface of Data.hpp.
 */

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
namespace RootReader {

class RootReader : public Data {

public:
  //! Destructor
  virtual ~RootReader();

  //! Default constructor
  RootReader();

  /**! Constructor read from TTree
   *
   * @param tr 		TTree to read
   * @param size 		Number of events to read
   * @param binned 	Create binned/unbinned data set
   */
  RootReader(TTree *tr, int size = -1, const bool binned = false);

  /*! Read constructor.
   *
   * @param inRootFile Input file name
   * @param inTreeName Name of tree in input file
   * @param size 		 Number of events to read in
   * @param binned	 Create binning
   */
  RootReader(const std::string inRootFile,
             const std::string inTreeName = "data", int size = -1,
             const bool binned = false);

  //! Create clone
  virtual RootReader *Clone() const;

  //! Create empty clone
  virtual RootReader *EmptyClone() const;

  virtual void writeData(std::string file = "", std::string trName = "");

protected:
  // Open ROOT file and set branch addresses to TTree
  void read();

  // Read TTree and store Events in fEvents
  virtual void storeEvents();

  // Create binning (obsolete?)
  virtual void bin();

  // Number of events to read from TTree
  int readSize;

  // Input file name
  std::string fileName;

  // Input tree name
  std::string treeName;

  // Pointer to Tfile
  TFile *fFile;

  // Pointer to TTree
  TTree *fTree;

  // TTree branch variables
  TClonesArray *fParticles;
  double feventWeight;
  double feventEff;
  int fCharge;
  int fFlavour;
};

} /* namespace RootReader */
} /* namespace DataReader */
} /* namespace ComPWA */

#endif /* _RootReader_HPP */
