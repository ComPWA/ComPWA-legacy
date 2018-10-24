// Copyright (c) 2013 Florian Feldbauer.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Ascii reader implementation.
///

#ifndef COMPWA_DATA_ASCIIREADER_HPP_
#define COMPWA_DATA_ASCIIREADER_HPP_

// ANSI C headers
#include <string>
#include <vector>

#include "Data/Data.hpp"
#include "Data/DataIOInterface.hpp"
// local headers
#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {
namespace Data {

///
/// \class AsciiReader
/// Reader for data in ASCII-Format. This class reads event-based data from
/// ascii-files in the same syntax. as Pawian's epemEvtReader. It implements the
/// interface of Data.hpp.
///
class AsciiReader : public DataIOInterface {
  unsigned int NumberOfParticles;

public:
  virtual ~AsciiReader();

  AsciiReader(){};

  AsciiReader(unsigned int NumberOfParticles_);

  // virtual AsciiReader *clone() const;

  // virtual AsciiReader *emptyClone() const;

  virtual void writeData(std::shared_ptr<ComPWA::Data::Data> Data,
                         const std::string &OutputFilePath) const;

  virtual std::shared_ptr<ComPWA::Data::Data>
  readData(const std::string &InputFilePath) const;
};

} // namespace Data
} // namespace ComPWA

#endif
