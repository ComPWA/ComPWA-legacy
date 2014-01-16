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

enum ParType { COMPLEX = 1, DOUBLE = 2, INTEGER = 3, BOOL = 4, UNDEFINED = 0};

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
  friend std::ostream& operator<<( std::ostream& out, const AbsParameter& b ){
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
  virtual std::string const& to_str() const{
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
  virtual std::string const& val_to_str() const{
    return outVal_;
  }

protected:
  std::string out_; /*!< Output string to print information */
  std::string outVal_; /*!< Output string to print only value */
  std::string name_; /*!< internal name of the parameter */
  ParType type_; /*!< ParType enum for type of parameter */

  std::vector<std::shared_ptr<ParObserver> > oberservingNodes; /*!< list of observers, e.g. TreeNodes */

  //! Interface to fill output string, to be implemented by parameter implementations
  virtual void make_str() =0;
};

#endif
