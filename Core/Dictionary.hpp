//! Dictionary to manage information provided by modules
/*! \class Dictionary
 * @file Dictionary.hpp
 * This class gathers information from modules, namely what they need as
 * input and what they provide as output with the given input.
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
#include "Physics/Amplitude.hpp"

struct dataInfo{
  std::string name;
  std::shared_ptr<Data> data;
  std::vector<std::string> parProvided;

  dataInfo(std::string inN, std::shared_ptr<Data> inD, std::vector<std::string> inP):
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
    * \return associated name for this DataReader
    * \sa std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");
   */
  virtual std::string introduce(std::shared_ptr<Data> inData, std::string inName="");

  //! Introduction of an Amplitude
   /*!
    * Use this function to introduce an Amplitude implementation. The
    * dictionary gathers information from it associates it with a name.
    * \param inAmp pointer to Amplitude
    * \return associated name for this Amplitude
    * \sa std::string introduce(std::shared_ptr<Data> inData, std::string inName="");
   */
  virtual std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");

  //! Check if name already used for a Data Module
   /*!
    * This function is used to check if a name is available for usage
    * \param inAmp pointer to Amplitude
    * \return boolean, true when name already used
    * \sa std::string introduce(std::shared_ptr<Data> inData, std::string inName="");
    * \sa std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");
   */
  virtual bool dataNameUsed(std::string inName);

  //! Check if name already used for an Amplitude
   /*!
    * This function is used to check if a name is available for usage
    * \param inAmp pointer to Amplitude
    * \return boolean, true when name already used
    * \sa std::string introduce(std::shared_ptr<Data> inData, std::string inName="");
    * \sa std::string introduce(std::shared_ptr<Amplitude> inAmp, std::string inName="");
   */
  virtual bool amplitudeNameUsed(std::string inName);


protected:
  std::vector<dataInfo > mData_; /*!< List of DataReader Info */
  std::vector<ampInfo > mAmps_; /*!< List of Amplitude Info */

};

#endif /* _DICTIONARY_HPP_ */
