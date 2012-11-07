//! Reader for data in Root-files.
/*! \class DIFRootReader
 * @file DIFRootReader.hpp
 * This class reads event-based data from root-files. It implements the
 * interface of DIFBase.hpp.
*/

#ifndef _DIFRootReader_HPP
#define _DIFRootReader_HPP

#include <vector>
#include <memory>
#include <string>

//PWA-Headers
#include "DIFBase.hpp"
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

class DIFRootReader : public DIFBase {

public:
  /// Default Constructor (0x0)
  DIFRootReader(std::string inConfigFile);

  virtual const int getEvent(const int, PWAEvent&);
  virtual const int getEvent(const int, TLorentzVector& , TLorentzVector& , double&);

  virtual const unsigned int getNEvents() const {return fmaxEvents;};

  /** Destructor */
  virtual ~DIFRootReader();

protected:

private:
  TFile* fFile;
  TTree* fTree;
  TClonesArray* fParticles;
  unsigned int fmaxEvents;
  unsigned int fEvent;
  // vector<string> paramNames;

};

#endif /* _DIFRootReader_HPP */
