//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------

//! Dictionary to manage information provided by modules
/*! \class Dictionary
 * @file Dictionary.hpp
 * This class gathers information from modules, namely what they need as
 * input and what they provide as output with the given input. It will only
 * be used during the setup phase, when performing a fit direct links are
 * used.
*/

#ifndef _DICTIONARY_HPP_
#define _DICTIONARY_HPP_

#include <vector>
#include <map>
#include <memory>
#include <string>

#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "DataReader/Data.hpp"
#include "Core/Amplitude.hpp"

namespace ComPWA {

struct dataInfo{
  std::string name;
  std::shared_ptr<DataReader::Data> data;
  std::vector<std::string> parProvided;

  dataInfo(std::string inN, std::shared_ptr<DataReader::Data> inD, std::vector<std::string> inP):
    name(inN), data(inD), parProvided(inP){};
};

struct ampInfo{
  std::string name;
  std::shared_ptr<Amplitude> amplitude;
  std::vector<std::string> inputNeeded;

  ampInfo(std::string inN, std::shared_ptr<Amplitude> inD, std::vector<std::string> inP):
    name(inN), amplitude(inD), inputNeeded(inP){};
};

class Dictionary
{
public:
  //! Standard constructor
   /*!
    * Standard constructor with no information provided, the dictionary
    * can't provide any information at this moment.
   */
  Dictionary();

  //! Destructor
  virtual ~Dictionary();

  //! Introduction of a DataReader
   /*!
    * Use this function to introduce a DataReader implementation. The
    * dictionary gathers information from it associates it with a name.
    * \param inData pointer to DataReader
    * \param inName unique name of the data set
    * \return associated name for this DataReader
    * \sa std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");
   */
  virtual std::string introduce(std::shared_ptr<DataReader::Data> inData, std::string inName="");

  //! Introduction of an Amplitude
   /*!
    * Use this function to introduce an Amplitude implementation. The
    * dictionary gathers information from it associates it with a name.
    * \param inAmp pointer to Amplitude
    * \param inName unique name of the data set
    * \return associated name for this Amplitude
    * \sa std::string introduce(std::shared_ptr<Data> inData, std::string inName="");
   */
  virtual std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");

  //! Check if name already used for a Data Module
   /*!
    * This function is used to check if a name is available for usage
    * \param inName name to be checked
    * \return boolean, true when name already used
    * \sa std::string introduce(std::shared_ptr<Data> inData, std::string inName="");
    * \sa std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");
   */
  virtual bool dataNameUsed(std::string inName);

  //! Check if name already used for an Amplitude
   /*!
    * This function is used to check if a name is available for usage
    * \param inName name to be checked
    * \return boolean, true when name already used
    * \sa std::string introduce(std::shared_ptr<Data> inData, std::string inName="");
    * \sa std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");
   */
  virtual bool amplitudeNameUsed(std::string inName);


protected:
  std::vector<dataInfo > mData_; /*!< List of DataReader Info */
  std::vector<ampInfo > mAmps_; /*!< List of Amplitude Info */

};

} /* namespace ComPWA */

#endif /* _DICTIONARY_HPP_ */
