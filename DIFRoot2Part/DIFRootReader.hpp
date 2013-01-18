//! Reader for data in Root-files.
/*! \class DIFRootReader
 * @file DIFRootReader.hpp
 * This class reads event-based data from root-files. It implements the
 * interface of Data.hpp.
*/

#ifndef _DIFRootReader_HPP
#define _DIFRootReader_HPP

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <utility>

//PWA-Headers
#include "Data.hpp"
#include "PWAEvent.hpp"

//Root-Headers
#include "TMath.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TRandom3.h"

class DIFRootReader : public Data {

public:
  /// Default Constructor (0x0)
  DIFRootReader(const std::string inConfigFile, const bool binned);

  virtual const int getEvent(const int, PWAEvent&);
  virtual const int getBin(const int, double&, double&);
  virtual const int getEvent(const int, TLorentzVector& , TLorentzVector& , double&);

  virtual const unsigned int getNEvents() const {return fmaxEvents;};
  virtual const unsigned int getNBins() const {return fmaxBins;};

  /** Destructor */
  virtual ~DIFRootReader();

protected:
  TFile* fFile;
  TTree* fTree;
  TClonesArray* fParticles;
  unsigned int fmaxEvents;
  unsigned int fEvent;
  bool fBinned;
  unsigned int fmaxBins;
  std::map<int, std::pair<double,double> > fBins;
  // vector<string> paramNames;

  virtual void bin();

};

#endif /* _DIFRootReader_HPP */
