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
//! Base class for internal parameter.
/*! \class AbsParameter
 * @file AbsParameter.hpp
 * This class defines the internal container of a parameter.
 * For the use in the function tree, the observer pattern is used and
 * this class takes over the role of the Subject. Therefore the actual
 * implementations of AbsParameter are the ConcreteSubjects of the
 * observer pattern and the TreeNodes take the role of the observers.
 */

#ifndef _ABSPARAMETER_HPP_
#define _ABSPARAMETER_HPP_


#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include "Core/ParObserver.hpp"
#include <fstream>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/shared_ptr.hpp>

enum ParType { COMPLEX = 1, DOUBLE = 2, INTEGER = 3, BOOL = 4, UNDEFINED = 0};
enum ErrorType { SYM = 1, ASYM = 2, LHSCAN = 3, NOTDEF = 0};

template <class T> class ParError
{
public:
	ParError(ErrorType t=ErrorType::NOTDEF) : type(t){};
	virtual ~ParError() {};
	virtual ErrorType GetType() { return type; };
	virtual T GetError() =0;
	virtual T GetErrorLow() =0;
	virtual T GetErrorHigh() =0;
	operator T() { return GetError(); }
private:
	ErrorType type;
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(type);
	}

};
BOOST_CLASS_TRACKING(ParError<double>, boost::serialization::track_always)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(ParError);

template <class T> class SymError : public ParError<T>
{
public:
	SymError() : ParError<T>(ErrorType::SYM), error(0){};
	SymError(T val) : ParError<T>(ErrorType::SYM), error(val){};
	virtual T GetError() {return error;}
	friend std::ostream& operator<<( std::ostream& out, const SymError& b ){
		return out << "+-" << GetError();
	}

private:
	virtual T GetErrorLow() { return error; };
	virtual T GetErrorHigh() { return error; };
	T error;
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ParError<T>);
//		ar & boost::serialization::base_object<ParError>(*this);
		ar & BOOST_SERIALIZATION_NVP(error);
	}

};
BOOST_CLASS_IMPLEMENTATION(SymError<double>,boost::serialization::object_serializable);

template <class T> class AsymError : public ParError<T>
{
public:
	AsymError() : ParError<T>(ErrorType::ASYM), error(std::pair<T,T>(0,0)){};
	AsymError(std::pair<T,T> val) : ParError<T>(ErrorType::ASYM), error(val){};
	virtual T GetError() {return (error.first+error.second)/2;}
	virtual T GetErrorLow() {
		if(error.first<0) return (-1)*error.first;
		return error.first;
	}
	virtual T GetErrorHigh() {return error.second;}
	virtual std::pair<T,T> GetLowHigh() {return error;}
	friend std::ostream& operator<<( std::ostream& out, const AsymError& b ){
		return out << "+" << GetErrorHigh() << "-"<<GetErrorLow();
	}

private:
	std::pair<T,T> error;
	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ParError<T>);
//		ar & boost::serialization::base_object<ParError>(*this);
		ar & BOOST_SERIALIZATION_NVP(error);
	}
};
BOOST_CLASS_IMPLEMENTATION(AsymError<double> , boost::serialization::object_serializable);

class AbsParameter //: public std::enable_shared_from_this<AbsParameter>
{
public:
  //! Constructor with name of parameter and optional type
  AbsParameter(std::string name, ParType type=ParType::UNDEFINED):name_(name), type_(type){

  }

  //! Destructor
  virtual ~AbsParameter(){
	  //std::cout << "GoodBye " << name_ <<std::endl;
  }

  //! Getter for name of object
  virtual const std::string& GetName(){
    return name_;
  }

  //! Getter for type of object
  virtual const ParType type(){
    return type_;
  }

  //! Getter for typename of object, to be defined by the actual implementation
  virtual const std::string TypeName()=0;

  //Observer Pattern Functions

  //! Attaches a new TreeNode as Observer
  void Attach(std::shared_ptr<ParObserver> newObserver){
    oberservingNodes.push_back(newObserver);;
  }

  //! Removes TreeNodes not needed as Observer anymore
  void Detach(std::shared_ptr<ParObserver> obsoleteObserver){
    oberservingNodes.erase(std::remove(oberservingNodes.begin(), oberservingNodes.end(), obsoleteObserver), oberservingNodes.end());
  }

  //! Notify all observing TreeNodes that parameter changed
  void Notify(){
    for(std::vector<std::shared_ptr<ParObserver> >::const_iterator iter = oberservingNodes.begin(); iter != oberservingNodes.end(); ++iter)
    {
        if(*iter != std::shared_ptr<ParObserver>())//Ist das richtig????
        {
            (*iter)->Update();
        }
    }
  }

  //! Return shared_pointer pointing to this Parameter
  //std::shared_ptr<AbsParameter> getptr() {
  //    return shared_from_this();
  //}

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type.
   * \sa make_str(), to_str()
  */
  friend std::ostream& operator<<( std::ostream& out, std::shared_ptr<AbsParameter> b ){
    return out << b->to_str();
  }

  //! friend function to stream parameter information to output
  /*!
   * Declaring the stream-operator << as friend allows to stream parameter
   * information to the output as easily as a generic type.
   * \sa make_str(), to_str()
  */
  friend std::ostream& operator<<( std::ostream& out, AbsParameter& b ){
    return out << b.to_str();
  }

  //! A public function returning a string with parameter information
  /*!
   * This function simply returns the member string out_, which contains
   * all parameter information. The string gets rebuild with every change
   * of the parameter.
   * \return string with parameter information
   * \sa operator<<, make_str()
  */
  virtual std::string const& to_str() {
	make_str();
    return out_;
  }

  //! A public function returning a string with parameter value
  /*!
   * This function simply returns the member string outVal_, which contains
   * the parameter value. The string gets rebuild with every change
   * of the parameter.
   * \return string with parameter information
   * \sa make_str()
  */
  virtual std::string const& val_to_str() {
	make_str();
    return outVal_;
  }

protected:
	std::string out_; /*!< Output string to print information */
	std::string outVal_; /*!< Output string to print only value */
	std::string name_; /*!< internal name of the parameter */
	ParType type_; /*!< ParType enum for type of parameter */
	//  ParError error_;
	//  bool hasError_; /*!< Is an error defined for this parameter? */

	std::vector<std::shared_ptr<ParObserver> > oberservingNodes; /*!< list of observers, e.g. TreeNodes */

	//! Interface to fill output string, to be implemented by parameter implementations
	virtual void make_str() =0;
private:

	friend class boost::serialization::access;
	template<class archive>
	void serialize(archive& ar, const unsigned int version)
	{
		ar & BOOST_SERIALIZATION_NVP(type_);
		ar & BOOST_SERIALIZATION_NVP(name_);
		ar & BOOST_SERIALIZATION_NVP(outVal_);
		ar & BOOST_SERIALIZATION_NVP(out_);
	}

};
BOOST_SERIALIZATION_ASSUME_ABSTRACT(AbsParameter);

#endif
