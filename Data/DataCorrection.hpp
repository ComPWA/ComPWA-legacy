// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_DATACORRECTION_HPP_
#define COMPWA_DATA_DATACORRECTION_HPP_

#include <cfloat>
#include <stdexcept>

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Data/CorrectionTable.hpp"

namespace ComPWA {
namespace Data {
class DataCorrection {
public:
  virtual ~DataCorrection() {}
  virtual double correction(Event &ev) const = 0;
};

class UnitCorrection : public DataCorrection {
public:
  virtual ~UnitCorrection() {}
  virtual double correction(Event &ev) const { return 1; }
};

class MomentumCorrection : public DataCorrection {
public:
  MomentumCorrection(ComPWA::ParticleList PartList_,
                     std::vector<ComPWA::Data::CorrectionTable> inCorr,
                     std::string t = "");
  virtual ~MomentumCorrection() {}

  virtual double correction(Event &ev) const;

  virtual double operator()(Event &ev) const { return correction(ev); };

  std::string title() { return Title; }

  void setTitle(std::string t) { Title = t; }

private:
  ParticleList PartList;
  std::vector<ComPWA::Data::CorrectionTable> Corrections;
  std::string Title;

  friend std::ostream &operator<<(std::ostream &out,
                                  const MomentumCorrection &b) {
    out << "MomentumCorrection::Print() | " << b.Title;
    for (const auto &c : b.Corrections)
      out << c;
    return out;
  }
};
} // namespace Data
} // namespace ComPWA
#endif
