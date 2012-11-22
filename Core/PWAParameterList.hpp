//! Internal container representing a parameter.
/*! \class PWAParameterList
 * @file PWAParameterList.hpp
 * This class provides a list of fit parameters which can have different types.
 * It consists of a vectors of parameters of type double, int and bool.
*/

#ifndef _PWAPARAMETERLIST_HPP_
#define _PWAPARAMETERLIST_HPP_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "PWAGenericPar.hpp"

class PWAParameterList
{

public:

  //! Standard constructor with empty parameter vector
  /*!
   * Standard constructor without input. The vectors of parameters are empty.
  */
  PWAParameterList();

  //! Standard constructor with a vector of double PWAParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vector gets
   * copied to the internal vector. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList. Non-double parameter
   * are empty
   * \param inVec input vector of double parameters
   * \sa addParameter(PWAParameter<double>&)
  */
  PWAParameterList(const std::vector<PWAGenericPar<double> >& inVec);

  //! Standard constructor with a vector of integer PWAParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vector gets
   * copied to the internal vector. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList. Non-integer parameter
   * are empty
   * \param inVec input vector of integer parameters
   * \sa addParameter(PWAParameter<integer>&)
  */
  PWAParameterList(const std::vector<PWAGenericPar<int> >& inVec);

  //! Standard constructor with a vector of boolean PWAParameter
  /*!
   * Standard constructor with list of PWAParameter provided. The vector gets
   * copied to the internal vector. To avoid copying, use the addParameter()
   * functions and empty constructor PWAParameterList. Non-boolean parameter
   * are empty
   * \param inVec input vector of boolean parameters
   * \sa addParameter(PWAParameter<bool>&)
  */
  PWAParameterList(const std::vector<PWAGenericPar<bool> >& inVec);

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
  PWAParameterList(const std::vector<PWAGenericPar<double> >& inD,
      const std::vector<PWAGenericPar<int> >& inI, const std::vector<PWAGenericPar<bool> >& inB);

  //! Copy constructor using = operator
  /*!
   * Simple copy constructor using the = operator. As this operator is not
   * overloaded in this class, c++ will copy every member variable. As this
   * is a container class, this should be fine.
   * \param in input PWAParameterList which variables will be copied
  */
  PWAParameterList(const PWAParameterList& in);

  //! Empty Destructor
  /*!
   * There is nothing to destroy :(
  */
  virtual ~PWAParameterList();


  //! Getter for floating point parameter
  /*!
   * Getter for floating point parameter
   * \param i input number of parameter to load
   * \param par output container for loaded parameter
  */
  virtual const int GetParameter(const unsigned int i, PWAGenericPar<double>& par) const ;

  //! Getter for integer parameter
  /*!
   * Getter for integer parameter
   * \param i input number of parameter to load
   * \param par output container for loaded parameter
  */
  virtual const int GetParameter(const unsigned int i, PWAGenericPar<int>& par) const ;

  //! Getter for boolean parameter
  /*!
   * Getter for boolean parameter
   * \param i input number of parameter to load
   * \param par output container for loaded parameter
  */
  virtual const int GetParameter(const unsigned int i, PWAGenericPar<bool>& par) const ;

  //! Add floating point parameter
  /*!
   * Adds a floating point parameter to the list
   * \param par input parameter
  */
  virtual void AddParameter(PWAGenericPar<double>& par);

  //! Add integer parameter
  /*!
   * Adds an integer parameter to the list
   * \param par input parameter
  */
  virtual void AddParameter(PWAGenericPar<int>& par);

  //! Add boolean parameter
  /*!
   * Adds an boolean parameter to the list
   * \param par input parameter
  */
  virtual void AddParameter(PWAGenericPar<bool>& par);

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
  std::vector<PWAGenericPar<double> > vDoublePar_; /*!< Vector of floating point parameters */
  std::vector<PWAGenericPar<int> > vIntPar_; /*!< Vector of integer parameters */
  std::vector<PWAGenericPar<bool> > vBoolPar_; /*!< Vector of boolean parameters */
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
  friend std::ostream & operator<<(std::ostream &os, const PWAParameterList &p);

};


#endif
