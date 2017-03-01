//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------
#ifndef ROOTEFFICIENCY_HPP_
#define ROOTEFFICIENCY_HPP_

#include <vector>
#include <memory>

#include <TEfficiency.h>
#include <TH2.h>

#include "Core/Efficiency.hpp"

namespace ComPWA {
namespace DataReader {

/**
 *  \class RootEfficiency
 *  \brief Efficiency provided by a histogram
 */
class RootEfficiency : public Efficiency {
protected:
  std::shared_ptr<TEfficiency> effHist;

public:
  //! Construct RootEfficiency from TEfficiency object
  RootEfficiency(TEfficiency *eff);
  //! Construct RootEfficiency from two TH2 objects for passed and total
  //! events
  RootEfficiency(TH1 *passed, TH1 *total);
  RootEfficiency(const RootEfficiency &);
  ~RootEfficiency(){};

  //! returns efficiency for current datapoint
  virtual double Evaluate(const dataPoint &point) const;
};

/**
 *  \class RootAngleEfficiency
 *  \brief Uses also TEfficiency object, but the variables are one invariant
 *  mass and the corresponding helicity angle. This avoids binning effects near
 * the phsp boundaries.
 *  ATTENTION: We assume that the invariant mass of particle 2 and 3 and the
 * helicity angle
 *  between 1 and 2 are used!
 */
class RootAngleEfficiency : public RootEfficiency {
public:
  RootAngleEfficiency(TEfficiency *eff) : RootEfficiency(eff){};
  
  RootAngleEfficiency(TH1 *passed, TH1 *total)
      : RootEfficiency(passed, total){};
  
  RootAngleEfficiency(const RootAngleEfficiency &p)
      : RootEfficiency(p){};
  
  ~RootAngleEfficiency(){};

  virtual double Evaluate(const dataPoint &point) const;
};
  
} /* namespace DataReader */
} /* namespace ComPWA */

#endif /* ROOTEFFICIENCY_HPP_ */
