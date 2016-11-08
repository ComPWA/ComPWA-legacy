//-------------------------------------------------------------------------------
// Copyright (c) 2013 Florian Feldbauer.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Florian Feldbauer - initial API and implementation
//-------------------------------------------------------------------------------
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

	AsciiReader( const std::string inConfigFile, const int particles );

	virtual AsciiReader* Clone() const;

	virtual AsciiReader* EmptyClone() const;

	virtual void writeData(std::string file="", std::string trName="") { };

};

} /* namespace AsciiReader */
} /* namespace DataReader */
} /* namespace ComPWA */

#endif /* _ASCII_READER_H_ */
