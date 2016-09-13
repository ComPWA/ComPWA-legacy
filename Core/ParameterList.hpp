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

	/**! Copy constructor using = operator
	 *
	 * Simple copy constructor using the = operator. As this operator is not
	 * overloaded in this class, c++ will copy every member variable. As this
	 * is a container class, this should be fine.
	 * \param in input PWAParameterList which variables will be copied
	 */
	ParameterList(const ParameterList& in) = default;

	void DeepCopy(const ParameterList& in);

	//! Empty Destructor
	/*!
	 * There is nothing to destroy :(
	 */
	virtual ~ParameterList();

	//! Get number of parameters
	virtual const inline unsigned int GetNParameter() const {
		return (
				vDouble_.size()+vInt_.size()+vBool_.size()
				+vMultiDouble_.size()+vMultiComplex_.size()
		);
	}

	std::shared_ptr<AbsParameter> GetParameter(const unsigned int i);

	std::shared_ptr<AbsParameter> GetParameter(const std::string parname);

	/**! Remove duplicate entries
	 * If a parameter name is found multiple times only the fist occurance is kept.
	 * The parameter values of parameters with the same name are not compared.
	 */
	virtual void RemoveDuplicates();

	//! Add parameter via abstract pointer
	/*!
	 * Adds a parameter with to be defined type to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<AbsParameter> par);

	//! Get string
	std::string const& to_str() ;

	//! Append ParameterList to (*this). Shared_ptr are not(!) deep copied
	virtual void Append(const ParameterList& addList);



	//**************************************************************************
	//************* Functions to access individual parameter types *************
	//**************************************************************************

	//================= Bool Parameter ==================
	//! Getter for number of boolean parameter
	virtual const inline unsigned int GetNBool() const {
		return vBool_.size();
	}
	/**! A public function returning a string with parameter information
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	virtual std::vector<std::shared_ptr<BoolParameter> >::const_iterator
	FindBoolParameter( const std::string name ) const;

	//! Get ID of Bool parameter
	virtual unsigned int FindBoolId( const std::string name );

	/**! Add boolean parameter
	 * Adds an boolean parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<BoolParameter> par);

	//! Getter for boolean parameter
	/*!
	 * Getter for boolean parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<BoolParameter>
	GetBoolParameter(const std::string parname) const ;

	//! Getter for boolean parameter
	/*!
	 * Getter for boolean parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<BoolParameter>
	GetBoolParameter(const unsigned int i) const ;

	//! Getter for boolean parameter
	/*!
	 * Getter for boolean parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<BoolParameter> >
	GetBoolParameters() const {
		return vBool_;
	}

	//! Get value of Bool parameter by name
	virtual const bool GetBoolParameterValue(const std::string parname) const;

	//! Get value of Bool parameter by id
	virtual const bool GetBoolParameterValue(const unsigned int i) const;

	//! Set value of Bool parameter by name
	virtual void SetParameterValue(const std::string name, const bool inVal);

	/*! Setter for parameter value
	 * Setter for parameter value
	 * \param i input number of parameter to load
	 * \param inVal input boolean value for parameter
	 */
	virtual void SetParameterValue(const unsigned int i, const bool inVal);

	//! Remove boolean parameter
	/*!
	 * Remove an boolean parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveBool(const std::string name);

	//! Remove boolean parameter
	/*!
	 * Remove an boolean parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveBool(const unsigned int id);


	//================= Integer Parameter ==================
	//! Getter for number of integer parameter
	virtual const inline unsigned int GetNInteger() const {
		return vInt_.size();
	}
	/**! A public function returning a string with parameter information
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	virtual std::vector<std::shared_ptr<IntegerParameter> >::const_iterator
	FindIntegerParameter( const std::string name ) const;

	//! Get ID of Integer parameter
	virtual unsigned int FindIntegerId( const std::string name );

	/**! Add integer parameter
	 * Adds an integer parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<IntegerParameter> par);

	//! Getter for integer parameter
	/*!
	 * Getter for integer parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<IntegerParameter>
	GetIntegerParameter(const std::string parname) const ;

	//! Getter for integer parameter
	/*!
	 * Getter for integer parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<IntegerParameter>
	GetIntegerParameter(const unsigned int i) const ;

	//! Getter for integer parameter
	/*!
	 * Getter for integer parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<IntegerParameter> >
	GetIntegerParameters() const {
		return vInt_;
	}

	//! Get value of Integer parameter by name
	virtual const int GetIntegerParameterValue(const std::string parname) const;

	//! Get value of Integer parameter by id
	virtual const int GetIntegerParameterValue(const unsigned int i) const;

	//! Set value of Integer parameter by name
	virtual void SetParameterValue(const std::string name, const int inVal);

	/*! Setter for parameter value
	 * Setter for parameter value
	 * \param i input number of parameter to load
	 * \param inVal input integer value for parameter
	 */
	virtual void SetParameterValue(const unsigned int i, const int inVal);

	//! Remove integer parameter
	/*!
	 * Remove an integer parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveInteger(const std::string name);

	//! Remove integer parameter
	/*!
	 * Remove an integer parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveInteger(const unsigned int id);

	//================= Double Parameter ==================
	//! Getter for number of double parameter
	virtual const inline unsigned int GetNDouble() const {
		return vDouble_.size();
	}
	/**! A public function returning a string with parameter information
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	virtual std::vector<std::shared_ptr<DoubleParameter> >::const_iterator
	FindDoubleParameter( const std::string name ) const;

	//! Get ID of Double parameter
	virtual unsigned int FindDoubleId( const std::string name );

	/**! Add double parameter
	 * Adds an double parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<DoubleParameter> par);

	//! Getter for double parameter
	/*!
	 * Getter for double parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<DoubleParameter>
	GetDoubleParameter(const std::string parname) const ;

	//! Getter for double parameter
	/*!
	 * Getter for double parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<DoubleParameter>
	GetDoubleParameter(const unsigned int i) const ;

	//! Getter for double parameter
	/*!
	 * Getter for double parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<DoubleParameter> >
	GetDoubleParameters() const {
		return vDouble_;
	}

	//! Get value of Double parameter by name
	virtual const double
	GetDoubleParameterValue(const std::string parname) const;

	//! Get value of Double parameter by id
	virtual const double GetDoubleParameterValue(const unsigned int i) const;

	//! Set value of Double parameter by name
	virtual void SetParameterValue(const std::string name, const double inVal);

	/*! Setter for parameter value
	 * Setter for parameter value
	 * \param i input number of parameter to load
	 * \param inVal input double value for parameter
	 */
	virtual void SetParameterValue(const unsigned int i, const double inVal);

	//! Remove double parameter
	/*!
	 * Remove an double parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveDouble(const std::string name);

	//! Remove double parameter
	/*!
	 * Remove an double parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveDouble(const unsigned int id);


	//================= Complex Parameter ==================
	//! Getter for number of complex parameter
	virtual const inline unsigned int GetNComplex() const {
		return vComplex_.size();
	}
	/**! A public function returning a string with parameter information
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	virtual std::vector<std::shared_ptr<ComplexParameter> >::const_iterator
	FindComplexParameter( const std::string name ) const;

	//! Get ID of Complex parameter
	virtual unsigned int FindComplexId( const std::string name );

	/**! Add complex parameter
	 * Adds an complex parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<ComplexParameter> par);

	//! Getter for complex parameter
	/*!
	 * Getter for complex parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<ComplexParameter>
	GetComplexParameter(const std::string parname) const ;

	//! Getter for complex parameter
	/*!
	 * Getter for complex parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<ComplexParameter>
	GetComplexParameter(const unsigned int i) const ;

	//! Getter for complex parameter
	/*!
	 * Getter for complex parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<ComplexParameter> >
	GetComplexParameters() const { return vComplex_; }

	//! Get value of Complex parameter by name
	virtual const std::complex<double>
	GetComplexParameterValue(const std::string parname) const;

	//! Get value of Complex parameter by id
	virtual const std::complex<double>
	GetComplexParameterValue(const unsigned int i) const;

	//! Set value of Complex parameter by name
	virtual void
	SetParameterValue(const std::string name, const std::complex<double> inVal);

	/*! Setter for parameter value
	 * Setter for parameter value
	 * \param i input number of parameter to load
	 * \param inVal input complex value for parameter
	 */
	virtual void
	SetParameterValue(const unsigned int i, const std::complex<double> inVal);

	//! Remove complex parameter
	/*!
	 * Remove an complex parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveComplex(const std::string name);

	//! Remove complex parameter
	/*!
	 * Remove an complex parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveComplex(const unsigned int id);

	//================= MultiDouble Parameter ==================
	//! Getter for number of multi-double parameter
	virtual const inline unsigned int GetNMultiDouble() const {
		return vMultiDouble_.size();
	}
	/**! A public function returning a string with parameter information
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	virtual std::vector<std::shared_ptr<MultiDouble> >::const_iterator
	FindMultiDouble( const std::string name ) const;

	//! Get ID of MultiDouble parameter
	virtual unsigned int FindMultiDoubleId( const std::string name );

	/**! Add multi-double parameter
	 * Adds an multi-double parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<MultiDouble> par);

	//! Getter for multi-double parameter
	/*!
	 * Getter for multi-double parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiDouble>
	GetMultiDouble(const std::string parname) const ;

	//! Getter for multi-double parameter
	/*!
	 * Getter for multi-double parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiDouble>
	GetMultiDouble(const unsigned int i) const ;

	//! Getter for multi-double parameter
	/*!
	 * Getter for multi-double parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<MultiDouble> >
	GetMultiDoubles() const {
		return vMultiDouble_;
	}

	//! Remove multi-double parameter
	/*!
	 * Remove an multi-double parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveMultiDouble(const std::string name);

	//! Remove multi-double parameter
	/*!
	 * Remove an multi-double parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveMultiDouble(const unsigned int id);

	//================= MultiComplex Parameter ==================
	//! Getter for number of multi-complex parameter
	virtual const inline unsigned int GetNMultiComplex() const {
		return vMultiComplex_.size();
	}
	/**! A public function returning a string with parameter information
	 * This function simply returns the member string out_, which contains
	 * all parameter information. The string gets created using the outstream
	 * of the PWAParameter class.
	 * \return string with parameter information
	 * \sa operator<<
	 */
	virtual std::vector<std::shared_ptr<MultiComplex> >::const_iterator
	FindMultiComplex( const std::string name ) const;

	//! Get ID of MultiComplex parameter
	virtual unsigned int FindMultiComplexId( const std::string name );

	/**! Add multi-complex parameter
	 * Adds an multi-complex parameter to the list
	 * \param par input parameter
	 */
	virtual void AddParameter(std::shared_ptr<MultiComplex> par);

	//! Getter for multi-complex parameter
	/*!
	 * Getter for multi-complex parameter
	 * \param parname input name of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiComplex>
	GetMultiComplex(const std::string parname) const ;

	//! Getter for multi-complex parameter
	/*!
	 * Getter for multi-complex parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::shared_ptr<MultiComplex>
	GetMultiComplex(const unsigned int i) const ;

	//! Getter for multi-complex parameter
	/*!
	 * Getter for multi-complex parameter
	 * \param i input number of parameter to load
	 * \return par output container for loaded parameter
	 */
	virtual std::vector<std::shared_ptr<MultiComplex> >
	GetMultiComplexs() const {
		return vMultiComplex_;
	}

	//! Remove multi-complex parameter
	/*!
	 * Remove an multi-complex parameter from the list
	 * \param parName parameter name
	 */
	virtual void RemoveMultiComplex(const std::string name);

	//! Remove multi-complex parameter
	/*!
	 * Remove an multi-complex parameter from the list
	 * \param par input parameter
	 */
	virtual void RemoveMultiComplex(const unsigned int id);


protected:
	/*!< Vector of boolean parameters */
	std::vector<std::shared_ptr<BoolParameter> > vBool_;
	/*!< Vector of integer parameters */
	std::vector<std::shared_ptr<IntegerParameter> > vInt_;
	/*!< Vector of floating point parameters */
	std::vector<std::shared_ptr<DoubleParameter> > vDouble_;
	/*!< Vector of complex parameters */
	std::vector<std::shared_ptr<ComplexParameter> > vComplex_;
	/*!< Vector of floating point parameter lists */
	std::vector<std::shared_ptr<MultiDouble> > vMultiDouble_;
	/*!< Vector of complex parameter lists */
	std::vector<std::shared_ptr<MultiComplex> > vMultiComplex_;

	/*!< Output string to print information */
	std::string out_;

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
		ar & make_nvp("DoubleParameters",vDouble_); //currently only DoubleParameters can be serialized
		ar & make_nvp("OutString",out_);
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
