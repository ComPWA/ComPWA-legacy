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
 * @file AsciiReader.h
 * This class reads event-based data from ascii-files in the same syntax as Pawian's epemEvtReader. It implements the
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

//_____ D E F I N I T I O N S __________________________________________________

class AsciiReader : public Data {

public:
  /// Default Constructor (0x0)
  AsciiReader( const std::string inConfigFile, const int particles );

  virtual const Event& getEvent( const int );
  virtual const int getBin( const int, double&, double& );

  virtual const unsigned int getNEvents() const { return EvtList_.size(); };
  virtual const unsigned int getNBins() const {return fmaxBins_;};

  /** Destructor */
  virtual ~AsciiReader();

  virtual std::shared_ptr<Data> rndSubSet(unsigned int size, std::shared_ptr<Generator> gen) {
	  return std::shared_ptr<Data>();}
protected:

private:
  std::vector<Event> EvtList_;
  unsigned int fmaxBins_;

};

#endif /* _ASCII_READER_H_ */
