// Copyright (c) 2013 Florian Feldbauer.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Ascii reader implementation.
///

#ifndef _ASCIIREADER_H_
#define _ASCIIREADER_H_

// ANSI C headers
#include <vector>
#include <string>

// local headers
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/DataPoint.hpp"


namespace ComPWA {
namespace DataReader {
namespace AsciiReader {

///
/// \class AsciiReader
/// Reader for data in ASCII-Format. This class reads event-based data from
/// ascii-files in the same syntax. as Pawian's epemEvtReader. It implements the
/// interface of Data.hpp.
///
class AsciiReader : public Data {

public:
  virtual ~AsciiReader();

  AsciiReader(){};

  AsciiReader(const std::string inConfigFile, const int particles);

  virtual AsciiReader *clone() const;

  virtual AsciiReader *emptyClone() const;

  virtual void writeData(std::string file = "", std::string trName = ""){};
};

} // ns::AsciiReader
} // ns::DataReader
} // ns::ComPWA

#endif
