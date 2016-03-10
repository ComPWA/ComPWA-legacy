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
#include <memory>
#include <sstream>
#include <vector>
#include <map>

#include "boost/serialization/vector.hpp"

#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"


class ParameterList
{

public:

	//! Standard constructor with empty parameter vector
	/*!
	 * Standard constructor without input. The vectors of parameters are empty.
	 */
	ParameterList();

	//! Standard constructor with a vector of ComplexParameter
	/*!
	 * Standard constructor with list of PWAParameter provided. The vector gets
	 * copied to the internal vector. To avoid copying, use the addParameter()
	 * functions and empty constructor PWAParameterList. Non-double parameter
	 * are empty
	 * \param inVec input vector of complex parameters
	 * \sa addParameter(PWAParameter<double>&)
	 */
	ParameterList(const std::vector<std::shared_ptr<ComplexParameter> >& inVec);

	//! Standard constructor with a vector of DoubleParameter
	/*!
	 * Standard constructor with list of PWAParameter provided. The vector gets
	 * copied to the internal vector. To avoid copying, use the addParameter()
	 * functions and empty constructor PWAParameterList. Non-double parameter
	 * are empty
	 * \param inVec input vector of double parameters
	 * \sa addParameter(PWAParameter<double>&)
	 */
	ParameterList(const std::vector<std::shared_ptr<DoubleParameter> >& inVec);

	//! Standard constructor with a vector of IntegerParameter
	/*!
	 * Standard constructor with list of PWAParameter provided. The vector gets
	 * copied to the internal vector. To avoid copying, use the addParameter()
	 * functions and empty constructor PWAParameterList. Non-integer parameter
	 * are empty
	 * \param inVec input vector of integer parameters
	 * \sa addParameter(PWAParameter<integer>&)
	 */
	ParameterList(const std::vector<std::shared_ptr<IntegerParameter> >& inVec);

	//! Standard constructor with a vector of BoolParameter
	/*!
	 * Standard constructor with list of PWAParameter provided. The vector gets
	 * copied to the internal vector. To avoid copying, use the addParameter()
	 * functions and empty constructor PWAParameterList. Non-boolean parameter
	 * are empty
	 * \param inVec input vector of boolean parameters
	 * \sa addParameter(PWAParameter<bool>&)
	 */
	ParameterList(const std::vector<std::shared_ptr<BoolParameter> >& inVec);

	//! Standard constructor with a vector of bool, int, double and complex PWAParameter
	/*!
	 * Standard constructor with list of PWAParameter provided. The vectors get
	 * copied to the internal vectors. To avoid copying, use the addParameter()
	 * functions and empty constructor PWAParameterList.
	 * \param inC input vector of complex parameters
	 * \param inD input vector of floating point parameters
	 * \param inI input vector of integer parameters
	 * \param inB input vector of boolean parameters
	 * \sa addParameter(PWAParameter<double>&, PWAParameter<int>&, PWAParameter<bool>&)
	 */
	ParameterList(const std::vector<std::shared_ptr<ComplexParameter> >& inC,
			const std::vector<std::shared_ptr<DoubleParameter> >& inD,
			const std::vector<std::shared_ptr<IntegerParameter> >& inI,
			const std::vector<std::shared_ptr<BoolParameter> >& inB);

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

	//! Getter for abstract parameter
	/*!
	 * Getter for abstract parameter pointer
	 * \param i input number of parameter to load
	 * \return shared pointer to parameter
	 */
	virtual std::shared_ptr<AbsParameter> GetParameter(const unsigned int i) ;

	//! Getter for abstract parameter
	/*!
	 * Getter for abstract parameter pointer
	 * \param parname name of parameter to load
	 * \return shared pointer to parameter
	 */
	virtual std::shared_ptr<AbsParameter> GetParameter(const std::string parname) ;


	//! Getter for number of parameter
	virtual const inline unsigned int GetNParameter() const {
		return (
				vDoublePar_.size()+vIntPar_.size()+vBoolPar_.size()
				+vMultiDouble_.size()+vMultiComplex_.size()
				);
	}

	//! Getter for number of multi complex parameter
	virtual const inline unsigned int GetNMultiComplex() const {
		return vMultiComplex_.size();
	}

	//! Getter for number of multi double parameter
	virtual const inline unsigned int GetNMultiDouble() const {
		return vMultiDouble_.size();
	}

	//! Getter for number of complex parameter
	virtual const inline unsigned int GetNComplex() const {
		return vComplexPar_.size();
	}

	//! Getter for number of double parameter
	virtual const inline unsigned int GetNDouble() const {return vDoublePar_.size();
	}

	//! Getter for number of integer parameter
	virtual const inline unsigned int GetNInteger() const {
		return vIntPar_.size();
	}

	//! Getter for number of boolean parameter
	virtual const inline unsigned int GetNBool() const {
		return vBoolPar_.size();
	}

	//! Getter for complex parameter
	/*!
	 * Getter for complex parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<ComplexParameter> GetComplexParameter(const unsigned int i) const ;

	//! Getter for floating point parameter
	/*!
	 * Getter for floating point parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<DoubleParameter> GetDoubleParameter(const unsigned int i) const ;

	//! Getter for double list parameter
	/*!
	 * Getter for double list parameter MultiDouble
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiDouble> GetMultiDouble(const unsigned int i) const ;

	//! Getter for complex list parameter
	/*!
	 * Getter for complex list parameter MultiComplex
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiComplex> GetMultiComplex(const unsigned int i) const ;

	//! Getter for integer parameter
	/*!
	 * Getter for integer parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<IntegerParameter> GetIntegerParameter(const unsigned int i) const ;

	//! Getter for boolean parameter
	/*!
	 * Getter for boolean parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<BoolParameter> GetBoolParameter(const unsigned int i) const ;

	//! Getter for parameter value
	/*!
	 * Getter for parameter value
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual const double GetParameterValue(const unsigned int i) const ;

	//! Getter for complex list parameter
	/*!
	 * Getter for double list parameter MultiComplex
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiComplex> GetMultiComplex(const std::string parname) const ;

	//! Getter for double list parameter
	/*!
	 * Getter for double list parameter MultiDouble
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiDouble> GetMultiDouble(const std::string parname) const ;

	//! Getter for complex parameter
	/*!
	 * Getter for complex parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<ComplexParameter> GetComplexParameter(const std::string parname) const ;

	//! Getter for floating point parameter
	/*!
	 * Getter for floating point parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<DoubleParameter> GetDoubleParameter(const std::string parname) const ;

	//! Getter for integer parameter
	/*!
	 * Getter for integer parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<IntegerParameter> GetIntegerParameter(const std::string parname) const ;

	//! Getter for boolean parameter
	/*!
	 * Getter for boolean parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<BoolParameter> GetBoolParameter(const std::string parname) const ;


	//! Getter for complex parameter
	/*!
	 * Getter for complex parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<ComplexParameter> > GetComplexParameters() const { return vComplexPar_; }

	//! Getter for floating point parameter
	/*!
	 * Getter for floating point parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<DoubleParameter> > GetDoubleParameters() const { return vDoublePar_; }

	//! Getter for double list parameter
	/*!
	 * Getter for double list parameter MultiDouble
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<MultiDouble> > GetMultiDoubles() const { return vMultiDouble_; }

	//! Getter for complex list parameter
	/*!
	 * Getter for complex list parameter MultiComplex
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<MultiComplex> > GetMultiComplexs() const { return vMultiComplex_; }

	//! Getter for integer parameter
	/*!
	 * Getter for integer parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<IntegerParameter> > GetIntegerParameters() const { return vIntPar_; }

	//! Getter for boolean parameter
	/*!
	 * Getter for boolean parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<BoolParameter> > GetBoolParameters() const { return vBoolPar_; }

	//! Getter for parameter value
	/*!
	 * Getter for parameter value
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual const double GetParameterValue(const std::string parname) const ;

	//! Setter for parameter value
	/*!
	 * Setter for parameter value
	 * \param i input number of parameter to load
	 * \param inVal input complex value for parameter
	 */
	virtual void SetParameterValue(const unsigned int i, const std::complex<double> inVal) ;

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

	//! Add parameter via abstract pointer
	/*!
	 * Adds a parameter with to be defined type to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<AbsParameter> par);

	//! Add complex list parameter
	/*!
	 * Adds a complex parameter MultiComplex to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<MultiComplex> par);

	//! Add double list parameter
	/*!
	 * Adds a double parameter MultiDouble to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<MultiDouble> par);

	//! Add complex parameter
	/*!
	 * Adds a complex parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<ComplexParameter> par);

	//! Add floating point parameter
	/*!
	 * Adds a floating point parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<DoubleParameter> par);

	//! Add integer parameter
	/*!
	 * Adds an integer parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<IntegerParameter> par);

	//! Add boolean parameter
	/*!
	 * Adds an boolean parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<BoolParameter> par);

	//! Remove complex parameter
	/*!
	 * Remove a complex parameter from the list
	 * \param id of parameter to remove
	 */
	virtual void RemoveComplex(const unsigned int id);

	//! Remove floating point parameter
	/*!
	 * Remove a floating point parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveDouble(const unsigned int id);

	//! Remove integer parameter
	/*!
	 * Remove an integer parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveInteger(const unsigned int id);

	//! Remove boolean parameter
	/*!
	 * Remove an boolean parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveBool(const unsigned int id);

	//! Remove complex parameter
	/*!
	 * Remove a complex parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveComplex(const std::string parName);

	//! Remove floating point parameter
	/*!
	 * Remove a floating point parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveDouble(const std::string parName);

	//! Remove integer parameter
	/*!
	 * Remove an integer parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveInteger(const std::string parName);

	//! Remove boolean parameter
	/*!
	 * Remove an boolean parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveBool(const std::string parName);

	//! A public function returning a string with parameter information
	/*!
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	std::string const& to_str() ;

	//! Append ParameterList to (*this). Shared_ptr are not(!) deep copied
	virtual void Append(ParameterList& addList);

	//! Create new index
	virtual void Indexing();

protected:
	std::map<std::string,unsigned int> mMultiComplexID_; /*!< Map of complex list parameter ids */
	std::map<std::string,unsigned int> mMultiDoubleID_; /*!< Map of double list parameter ids */
	std::map<std::string,unsigned int> mComplexParID_; /*!< Map of complex parameter ids */
	std::map<std::string,unsigned int> mDoubleParID_; /*!< Map of floating point parameter ids */
	std::map<std::string,unsigned int> mIntParID_; /*!< Map of integer parameter ids */
	std::map<std::string,unsigned int> mBoolParID_; /*!< Map of boolean parameter ids */
	std::vector<std::shared_ptr<MultiComplex> > vMultiComplex_; /*!< Vector of complex parameter lists */
	std::vector<std::shared_ptr<MultiDouble> > vMultiDouble_; /*!< Vector of floating point parameter lists */
	std::vector<std::shared_ptr<ComplexParameter> > vComplexPar_; /*!< Vector of complex parameters */
	std::vector<std::shared_ptr<DoubleParameter> > vDoublePar_; /*!< Vector of floating point parameters */
	std::vector<std::shared_ptr<IntegerParameter> > vIntPar_; /*!< Vector of integer parameters */
	std::vector<std::shared_ptr<BoolParameter> > vBoolPar_; /*!< Vector of boolean parameters */
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
	 * of this function has to be outside the namespace of the class.
	 * \sa make_str(), to_str()
	 */
	friend std::ostream & operator<<(std::ostream &os, ParameterList &p);

private:
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		using namespace boost::serialization;
		ar & make_nvp("DoubleParameters",vDoublePar_); //currently only DoubleParameters can be serialized
//		ar & make_nvp("DoubleParametersID",mDoubleParID_);
		ar & make_nvp("OutString",out_);
		Indexing();
	}
};
BOOST_SERIALIZATION_SHARED_PTR(ParameterList)
BOOST_CLASS_IMPLEMENTATION( ParameterList, boost::serialization::object_serializable )
BOOST_CLASS_TRACKING( ParameterList, boost::serialization::track_never )

#include <boost/serialization/split_free.hpp>
#include <boost/unordered_map.hpp>
#include <typeinfo>

//---/ Wrapper for std::shared_ptr<> /------------------------------------------

namespace boost {
namespace serialization {

template<class Archive, class Type>
void save(Archive & archive, const std::shared_ptr<Type> & value, const unsigned int /*version*/)
{
	Type *data = value.get();
	archive << make_nvp("shared_ptr",data);
}

template<class Archive, class Type>
void load(Archive & archive, std::shared_ptr<Type> & value, const unsigned int /*version*/)
{
	Type *data;
	archive >> make_nvp("shared_ptr",data);
//	archive >>data;

	typedef std::weak_ptr<Type> WeakPtr;
	static boost::unordered_map<void*, WeakPtr> hash;

	if (hash[data].expired())
	{
		value = std::shared_ptr<Type>(data);
		hash[data] = value;
	}
	else value = hash[data].lock();
}

template<class Archive, class Type>
inline void serialize(Archive & archive, std::shared_ptr<Type> & value, const unsigned int version)
{
	split_free(archive, value, version);
}

}//ns:serialization
}//ns:boost


#endif
