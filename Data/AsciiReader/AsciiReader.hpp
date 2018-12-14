// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ASCIIREADER_HPP_
#define COMPWA_DATA_ASCIIREADER_HPP_

#include <memory>
#include <vector>

#include "Core/Event.hpp"

namespace ComPWA {
namespace Data {

///
/// \class AsciiReader
/// Reader for data in ASCII-Format. This class reads event-based data from
/// ascii-files in the same syntax.
///
class AsciiReader {
  unsigned int NumberOfParticles;

public:
  virtual ~AsciiReader();

  AsciiReader(unsigned int NumberOfParticles_);

  std::shared_ptr<std::vector<ComPWA::Event>>
  readData(const std::string &InputFilePath) const;
};

} // namespace Data
} // namespace ComPWA

#endif
