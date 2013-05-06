//! Reader for data in ASCII-Format like Pawian's epemEvtReader
/*! \class PawianEpemReader
 * @file PawianEpemReader.hpp
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
  AsciiReader( const std::string inConfigFile );

  virtual const Event* getEvent( const int );
  virtual const int getBin( const int, double&, double& );

  virtual const unsigned int getNEvents() const { return EvtList.size(); };
  virtual const unsigned int getNBins() const {return fmaxBins;};

  /** Destructor */
  virtual ~AsciiReader();

protected:

private:
  std::vector<Event> EvtList_;

};

#endif /* _RootReader_HPP */
