// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ROOTEFFICIENCY_HPP_
#define COMPWA_DATA_ROOTEFFICIENCY_HPP_

#include <memory>
#include <vector>

#include "Core/Efficiency.hpp"

#include <TEfficiency.h>

class TH1;

namespace ComPWA {
namespace Data {
struct DataSet;

namespace Root {

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
  ~RootEfficiency() = default;

  virtual std::vector<double> evaluate(const DataSet &dataset) const;
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

  ~RootAngleEfficiency() = default;

  virtual std::vector<double> evaluate(const DataSet &dataset) const;
};

} // namespace Root
} // namespace Data
} // namespace ComPWA

#endif
