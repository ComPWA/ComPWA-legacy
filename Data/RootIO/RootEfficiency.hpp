// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ROOTEFFICIENCY_HPP_
#define COMPWA_DATA_ROOTEFFICIENCY_HPP_

#include <vector>
#include <memory>

#include <TEfficiency.h>

#include "Core/Efficiency.hpp"

class TH1;

namespace ComPWA {
namespace Data {

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
  virtual double evaluate(const DataPoint &point) const;
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

  virtual double evaluate(const DataPoint &point) const;
};
  
} /* namespace DataReader */
} /* namespace ComPWA */

#endif /* ROOTEFFICIENCY_HPP_ */
