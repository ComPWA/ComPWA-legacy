// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Parameter base class.
///

#ifndef _Parameter_HPP_
#define _Parameter_HPP_

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <fstream>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/tracking.hpp>

#include "Core/ParObserver.hpp"

namespace ComPWA {

/// Enums for the type of the parameter, should be extended if an new parameter
/// type is added
enum ParType {
  COMPLEX = 1,
  DOUBLE = 2,
  INTEGER = 3,
  BOOL = 4,
  MDOUBLE = 5,
  MCOMPLEX = 6,
  MUNSIGNEDINTEGER = 7,
  UNDEFINED = 0
};

/// Names of the parameter types, should be extended if an new parameter type is
/// added
static const char *ParNames[8] = {
    "UNDEFINED", "COMPLEX", "DOUBLE",   "INTEGER",
    "BOOL",      "MDOUBLE", "MCOMPLEX", "MUNSIGNEDINTEGER"};



///
/// \class Parameter
/// Base class for internal parameter.
/// This class defines the internal container of a parameter. For the use in
/// the function tree, the observer pattern is used and this class takes over
/// the role of the Subject. Therefore the actual implementations of
/// Parameter are the ConcreteSubjects of the observer pattern and the
/// TreeNodes take the role of the observers.
///
class Parameter {
public:
  /// Constructor with name of parameter and optional type
  Parameter(std::string name, ParType type = ParType::UNDEFINED)
      : Name(name), Type(type) {}

  virtual ~Parameter() {}

  /// Getter for name of object
  virtual std::string name() const { return Name; }
  
  /// Getter for name of object
  virtual void setName( std::string n ) { Name = n; }

  /// Getter for type of object
  virtual ParType type() const { return Type; }

  /// Getter for typename of object, to be defined by the actual implementation
  virtual std::string className() const = 0;

  // Observer Pattern Functions

  /// Attaches a new TreeNode as Observer
  void Attach(std::shared_ptr<ParObserver> newObserver) {
    oberservingNodes.push_back(newObserver);
    ;
  }

  /// Removes TreeNodes not needed as Observer anymore
  void Detach(std::shared_ptr<ParObserver> obsoleteObserver) {
    oberservingNodes.erase(std::remove(oberservingNodes.begin(),
                                       oberservingNodes.end(),
                                       obsoleteObserver),
                           oberservingNodes.end());
  }

  /// Notify all observing TreeNodes that parameter changed
  void Notify() {
    for (std::vector<std::shared_ptr<ParObserver>>::const_iterator iter =
             oberservingNodes.begin();
         iter != oberservingNodes.end(); ++iter) {
      if (*iter != std::shared_ptr<ParObserver>()) // Ist das richtig????
      {
        (*iter)->update();
      }
    }
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  std::shared_ptr<Parameter> b) {
    return out << b->to_str();
  }

  friend std::ostream &operator<<(std::ostream &out, Parameter &b) {
    return out << b.to_str();
  }

  /// A public function returning a string with parameter information
  virtual std::string to_str() const = 0;

  /// A public function returning a string with parameter value
  virtual std::string val_to_str() const = 0;
  
protected:
  /// Name of parameter
  std::string Name;
  
  /// Type of parameter (e.g. Double, Integer, ...)
  ParType Type;

  /// List of observers, e.g. TreeNodes
  std::vector<std::shared_ptr<ParObserver>> oberservingNodes;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(Name);
    ar &BOOST_SERIALIZATION_NVP(Type);
  }
};
} // ns::ComPWA

BOOST_SERIALIZATION_SHARED_PTR(Parameter);

BOOST_CLASS_IMPLEMENTATION(
    ComPWA::Parameter, boost::serialization::level_type::object_serializable)

#endif
