// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef DATAREADER_DATACORRECTION_HPP_
#define DATAREADER_DATACORRECTION_HPP_

#include <stdexcept>
#include <cfloat>

#include "Core/Event.hpp"
#include "Core/Logging.hpp"
#include "DataReader/CorrectionTable.hpp"

namespace ComPWA {

class DataCorrection {
public:
  virtual double correction(Event &ev) = 0;
  
  virtual void print() const = 0;
};

class UnitCorrection : public DataCorrection {
public:
  virtual double correction(Event &ev) { return 1; }
  
  virtual void print() const {
    LOG(info) << "UnitCorrection::Print() | correction factors are set to one";
  }
};

class MomentumCorrection : public DataCorrection {
public:
  MomentumCorrection(std::vector<ComPWA::DataReader::CorrectionTable> inCorr,
                     std::string t = "");

  virtual double correction(Event &ev);
  
  virtual void print() const;
  
  std::string title() { return Title; }
  
  void setTitle(std::string t) { Title = t; }

protected:
  std::vector<ComPWA::DataReader::CorrectionTable> Corrections;
  std::string Title;
};
} // ns::ComPWA
#endif
