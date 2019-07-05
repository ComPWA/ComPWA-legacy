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
  virtual double correction(Event &ev) = 0;

  virtual void print() const = 0;
};

class UnitCorrection : public DataCorrection {
public:
  virtual ~UnitCorrection() {}
  virtual double correction(Event &ev) { return 1; }

  virtual void print() const {
    LOG(INFO) << "UnitCorrection::Print() | correction factors are set to one";
  }
};

class MomentumCorrection : public DataCorrection {
public:
  MomentumCorrection(std::shared_ptr<ComPWA::PartList> list,
                     std::vector<ComPWA::Data::CorrectionTable> inCorr,
                     std::string t = "");
  virtual ~MomentumCorrection() {}

  virtual double correction(Event &ev);

  //virtual void print() const;

  std::string title() { return Title; }

  void setTitle(std::string t) { Title = t; }

private:
  std::shared_ptr<PartList> List;
  std::vector<ComPWA::Data::CorrectionTable> Corrections;
  std::string Title;
};
} // namespace Data
} // namespace ComPWA
#endif
