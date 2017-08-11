// Copyright (c) 2013 Florian Feldbauer.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Reader for data in ASCII-Format like Pawian's epemEvtReader
/*! \class AsciiReader
 * @file AsciiReader.hpp
 * This class reads event-based data from ascii-files in the same syntax
 * as Pawian's epemEvtReader. It implements the
 * interface of Data.hpp.
 */

#ifndef _ASCII_READER_H_
#define _ASCII_READER_H_

//_____ I N C L U D E S _______________________________________________________

// ANSI C headers
#include <vector>
#include <string>

// local headers
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/DataPoint.hpp"

//_____ D E F I N I T I O N S __________________________________________________

namespace ComPWA {
namespace DataReader {
namespace AsciiReader {

class AsciiReader : public Data {

public:
  /** Destructor */
  virtual ~AsciiReader();

  AsciiReader(){};

  AsciiReader(const std::string inConfigFile, const int particles);

  virtual AsciiReader *Clone() const;

  virtual AsciiReader *EmptyClone() const;

  virtual void WriteData(std::string file = "", std::string trName = ""){};
};

} /* namespace AsciiReader */
} /* namespace DataReader */
} /* namespace ComPWA */

#endif /* _ASCII_READER_H_ */
