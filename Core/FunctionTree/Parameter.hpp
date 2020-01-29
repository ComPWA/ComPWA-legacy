// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Parameter base class.
///

#ifndef _Parameter_HPP_
#define _Parameter_HPP_

#include <algorithm>
#include <complex>
#include <fstream>

#include <memory>
#include <string>
#include <vector>

#include "Core/FunctionTree/ParObserver.hpp"

#include <boost/serialization/level.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/tracking.hpp>

namespace ComPWA {
namespace FunctionTree {

/// Enums for the type of the parameter, should be extended if an new parameter
/// type is added
enum ParType {
  UNDEFINED = 0,
  COMPLEX = 1,
  DOUBLE = 2,
  INTEGER = 3,
  MCOMPLEX = 4,
  MDOUBLE = 5,
  MINTEGER = 6
};

/// Names of the parameter types, should be extended if an new parameter type is
/// added
static const char *const ParNames[7] = {"UNDEFINED", "COMPLEX",  "DOUBLE",
                                        "INTEGER",   "MCOMPLEX", "MDOUBLE",
                                        "MINTEGER"};

/// Template functions which return above specified parameter types
template <typename T> inline ParType typeName(void) {
  return ParType::UNDEFINED;
}

template <> inline ParType typeName<std::vector<std::complex<double>>>(void) {
  return ParType::MCOMPLEX;
}
template <> inline ParType typeName<std::vector<double>>(void) {
  return ParType::MDOUBLE;
}
template <> inline ParType typeName<std::vector<int>>(void) {
  return ParType::MINTEGER;
}
template <> inline ParType typeName<std::complex<double>>(void) {
  return ParType::COMPLEX;
}
template <> inline ParType typeName<double>(void) { return ParType::DOUBLE; }
template <> inline ParType typeName<int>(void) { return ParType::INTEGER; }

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

  virtual ~Parameter() = default;

  /// Getter for name of object
  virtual std::string name() const { return Name; }

  /// Getter for name of object
  virtual void setName(std::string n) { Name = n; }

  /// Getter for type of object
  virtual ParType type() const { return Type; }

  /// Getter for typename of object, to be defined by the actual implementation
  virtual std::string className() const = 0;

  virtual bool isParameter() const { return false; }

  // Observer Pattern Functions

  /// Attaches a new TreeNode as Observer
  void attach(std::weak_ptr<ParObserver> newObserver) {
    ObservingNodes.push_back(newObserver);
  }

  /// Removes TreeNodes not needed as Observer anymore
  void detachExpired() {
    ObservingNodes.erase(std::remove_if(ObservingNodes.begin(),
                                        ObservingNodes.end(),
                                        [](auto x) { return x.expired(); }),
                         ObservingNodes.end());
  }

  /// Notify all observing TreeNodes that parameter changed
  void notify() {
    for (auto x : ObservingNodes) {
      x.lock()->update();
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

  //  template <typename T>
  //  std::shared_ptr<T> GetComponent() {
  //    return std::dynamic_pointer_cast<T>(shared_from_this());
  //  }

protected:
  /// Name of parameter
  std::string Name;

  /// Type of parameter (e.g. Double, Integer, ...)
  ParType Type;

  /// List of observers, e.g. TreeNodes
  std::vector<std::weak_ptr<ParObserver>> ObservingNodes;

private:
  friend class boost::serialization::access;
  template <class archive>
  void serialize(archive &ar, const unsigned int version) {
    ar &BOOST_SERIALIZATION_NVP(Name);
    ar &BOOST_SERIALIZATION_NVP(Type);
  }
};
} // namespace FunctionTree
} // namespace ComPWA

BOOST_SERIALIZATION_SHARED_PTR(Parameter);

BOOST_CLASS_IMPLEMENTATION(
    ComPWA::FunctionTree::Parameter,
    boost::serialization::level_type::object_serializable)

#endif
