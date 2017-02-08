//-------------------------------------------------------------------------------
// Copyright (c) 2016 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Reader for data in Root-files.
/*! \class JakeReader
 * @file JakeReader.hpp
 * This class reads event-based data from BESIII root-files. It implements the
 * interface of Data.hpp.
 */

#ifndef _JakeReader_HPP
#define _JakeReader_HPP

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <utility>

// PWA-Headers
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"

// Root-Headers
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

namespace ComPWA {
namespace DataReader {
namespace JakeReader {

class JakeReader : public Data {

public:
  //! Default Constructor: Empty object, can be filled
  JakeReader();
  JakeReader(TTree *tr, int size = -1, const bool binned = false);
  /*! Read constructor.
   *
   * @param inRootFile Input file name
   * @param inTreeName Name of tree in input file
   * @param size 		Number of events to read in
   * @param binned	 Create binning
   */
  JakeReader(const std::string inRootFile, const std::string inTreeName = "kin",
             int size = -1, const bool binned = false);

  virtual const std::vector<std::string> &getVariableNames();

  //  virtual void writeToFile();
  virtual void pushEvent(const Event &evt) { fEvents.push_back(evt); };
  virtual Event &getEvent(const int);
  virtual const int getBin(const int, double &, double &);
  // virtual const int getEvent(const int, TLorentzVector& , TLorentzVector& ,
  // double&);
  virtual void writeData(std::string file = "", std::string trName = "");
  virtual void Clear();

  virtual const unsigned int getNEvents() const { return fEvents.size(); };
  virtual const unsigned int getNBins() const { return fmaxBins; };

  virtual std::vector<Event> getEvents() { return fEvents; }
  virtual std::vector<dataPoint> getDataPoints();
  virtual void Add(Data &otherSample) {
    std::vector<Event> otherEvents = otherSample.getEvents();
    fEvents.insert(fEvents.end(), otherEvents.begin(), otherEvents.end());
  }
  //! Destructor
  virtual ~JakeReader();

  //! Create clone
  virtual JakeReader *Clone() const;

  //! Create empty clone
  virtual JakeReader *EmptyClone() const;

  //! Set efficiency value for all stored events. Efficiency is taken from
  //! Efficiency object.
  void setEfficiency(std::shared_ptr<Efficiency> eff);
  //! Reset effciencies of all events
  void resetEfficiency(double e = 1.);
  //! Reset weights
  void resetWeights(double w = 1.);
  //! Weights set?
  bool hasWeights();

protected:
  void read();
  int readSize;
  bool _readFlag;

  // input tree
  std::string fileName;
  std::string treeName;
  TFile *fFile;
  TTree *fTree;
  double fe[3];
  double fpx[3];
  double fpy[3];
  double fpz[3];
  double feventWeight;

  // binning
  bool fBinned;
  unsigned int fmaxBins;
  std::map<int, std::pair<double, double>> fBins;

  // storage of events
  std::vector<std::string> fVarNames;
  std::vector<Event> fEvents;
  //	unsigned int fmaxEvents;
  //  unsigned int fEvent;

  virtual void storeEvents();
  virtual void bin();
};

} /* namespace JakeReader */
} /* namespace DataReader */
} /* namespace ComPWA */

#endif /* _JakeReader_HPP */
