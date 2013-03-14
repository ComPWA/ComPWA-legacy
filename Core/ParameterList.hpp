//! Internal container representing a parameter list.
/*! \class ParameterList
 * @file ParameterList.hpp
 * This class provides a list of fit parameters which can have different types.
 * It consists of a vectors of parameters of type double, int and bool.
*/

#ifndef _PARAMETERLIST_HPP_
#define _PARAMETERLIST_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "Core/Parameter.hpp"
#include "Core/PWAExceptions.hpp"

class ParameterList
{

public:

  //! Standard constructor with empty parameter vector
  /*!
   * Standard constructor without input. The vectors of parameters are empty.
  */
  ParameterList();

  //! Standard constructor with a vector of DoubleParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vector gets
   * copied to the internal vector. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList. Non-double parameter
   * are empty
   * \param inVec input vector of double parameters
   * \sa addParameter(PWAParameter<double>&)
  */
  ParameterList(const std::vector<DoubleParameter>& inVec);

  //! Standard constructor with a vector of IntegerParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vector gets
   * copied to the internal vector. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList. Non-integer parameter
   * are empty
   * \param inVec input vector of integer parameters
   * \sa addParameter(PWAParameter<integer>&)
  */
  ParameterList(const std::vector<IntegerParameter>& inVec);

  //! Standard constructor with a vector of BoolParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vector gets
   * copied to the internal vector. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList. Non-boolean parameter
   * are empty
   * \param inVec input vector of boolean parameters
   * \sa addParameter(PWAParameter<bool>&)
  */
  ParameterList(const std::vector<BoolParameter>& inVec);

  //! Standard constructor with a vector of bool, int and double PWAParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vectors get
   * copied to the internal vectors. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList.
   * \param inD input vector of floating point parameters
   * \param inI input vector of integer parameters
   * \param inB input vector of boolean parameters
   * \sa addParameter(PWAParameter<double>&, PWAParameter<int>&, PWAParameter<bool>&)
  */
  ParameterList(const std::vector<DoubleParameter>& inD,
      const std::vector<IntegerParameter>& inI, const std::vector<BoolParameter>& inB);

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameterList which variables will be copied
  */
  ParameterList(const ParameterList& in);

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~ParameterList();

  //! Getter for number of parameter
  virtual const inline unsigned int GetNParameter() const {return (vDoublePar_.size()+vIntPar_.size()+vBoolPar_.size());}
  //! Getter for number of double parameter
  virtual const inline unsigned int GetNDouble() const {return vDoublePar_.size();}
  //! Getter for number of integer parameter
  virtual const inline unsigned int GetNInteger() const {return vIntPar_.size();}
  //! Getter for number of boolean parameter
  virtual const inline unsigned int GetNBool() const {return vBoolPar_.size();}

  //! Getter for floating point parameter
  /*!
   * Getter for floating point parameter
   * \param i input number of parameter to load
   * \return par output container for loaded parameter
  */
  virtual DoubleParameter& GetDoubleParameter(const unsigned int i) ;

  //! Getter for integer parameter
  /*!
   * Getter for integer parameter
   * \param i input number of parameter to load
   * \return par output container for loaded parameter
  */
  virtual IntegerParameter& GetIntegerParameter(const unsigned int i) ;

  //! Getter for boolean parameter
  /*!
   * Getter for boolean parameter
   * \param i input number of parameter to load
   * \return par output container for loaded parameter
  */
  virtual BoolParameter& GetBoolParameter(const unsigned int i) ;

  //! Getter for parameter value
  /*!
   * Getter for parameter value
   * \param i input number of parameter to load
   * \return par output container for loaded parameter
  */
  virtual const double GetParameterValue(const unsigned int i) const ;

  //! Setter for parameter value
  /*!
   * Setter for parameter value
   * \param i input number of parameter to load
   * \param inVal input floating value for parameter
  */
  virtual void SetParameterValue(const unsigned int i, const double inVal) ;

  //! Setter for parameter value
  /*!
   * Setter for parameter value
   * \param i input number of parameter to load
   * \param inVal input integer value for parameter
  */
  virtual void SetParameterValue(const unsigned int i, const int inVal) ;

  //! Setter for parameter value
  /*!
   * Setter for parameter value
   * \param i input number of parameter to load
   * \param inVal input boolean value for parameter
  */
  virtual void SetParameterValue(const unsigned int i, const bool inVal) ;


  //! Add floating point parameter
  /*!
   * Adds a floating point parameter to the list
   * \param par input parameter
  */
  virtual void AddParameter(DoubleParameter par);

  //! Add integer parameter
  /*!
   * Adds an integer parameter to the list
   * \param par input parameter
  */
  virtual void AddParameter(IntegerParameter par);

  //! Add boolean parameter
  /*!
   * Adds an boolean parameter to the list
   * \param par input parameter
  */
  virtual void AddParameter(BoolParameter par);

  //! A public function returning a string with parameter information
  /*!
   * This function simply returns the member string out_, which contains
   * all parameter information. The string gets created using the outstream
   * of the PWAParameter class.
   * \return string with parameter information
   * \sa operator<<
  */
  std::string const& to_str() const;

protected:
  std::vector<DoubleParameter> vDoublePar_; /*!< Vector of floating point parameters */
  std::vector<IntegerParameter> vIntPar_; /*!< Vector of integer parameters */
  std::vector<BoolParameter> vBoolPar_; /*!< Vector of boolean parameters */
  std::string out_; /*!< Output string to print information */

  //! A protected function which creates an output string for printing
  /*!
   * This function uses all available information about the parameterlist
   * to create a string which will be streamed via the stream operator <<.
   * \sa operator<<, to_str()
  */
  void make_str();

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type. The definition
   * of this class has to be outside the namespace of the class.
   * \sa make_str(), to_str()
  */
  friend std::ostream & operator<<(std::ostream &os, const ParameterList &p);

};


#endif
