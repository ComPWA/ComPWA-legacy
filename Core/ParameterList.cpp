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
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <map>

#include "Core/AbsParameter.hpp"
#include "Core/Exceptions.hpp"

#include "Core/ParameterList.hpp"

namespace ComPWA {

ParameterList::ParameterList() {}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<ComplexParameter>> &inVec)
    : vComplex_(inVec) {}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<DoubleParameter>> &inVec)
    : vDouble_(inVec) {}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<IntegerParameter>> &inVec)
    : vInt_(inVec) {}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<BoolParameter>> &inVec)
    : vBool_(inVec) {}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<ComplexParameter>> &inC,
    const std::vector<std::shared_ptr<DoubleParameter>> &inD,
    const std::vector<std::shared_ptr<IntegerParameter>> &inI,
    const std::vector<std::shared_ptr<BoolParameter>> &inB)
    : vBool_(inB), vInt_(inI), vDouble_(inD), vComplex_(inC) {
  // TODO check names!
}

void ParameterList::DeepCopy(const ParameterList &in) {
  out_ = in.out_;
  vMultiComplex_.clear();
  vMultiDouble_.clear();
  vDouble_.clear();
  vInt_.clear();
  vBool_.clear();
  for (unsigned int i = 0; i < in.GetNMultiComplex(); i++)
    vMultiComplex_.push_back(std::shared_ptr<MultiComplex>(
        new MultiComplex(*(in.vMultiComplex_[i]))));
  for (unsigned int i = 0; i < in.GetNMultiDouble(); i++)
    vMultiDouble_.push_back(
        std::shared_ptr<MultiDouble>(new MultiDouble(*(in.vMultiDouble_[i]))));
  for (unsigned int i = 0; i < in.GetNDouble(); i++)
    vDouble_.push_back(std::shared_ptr<DoubleParameter>(
        new DoubleParameter(*(in.vDouble_[i]))));
  for (unsigned int i = 0; i < in.GetNMultiUnsignedInteger(); i++)
    vMultiUnsignedInteger_.push_back(std::shared_ptr<MultiUnsignedInteger>(
        new MultiUnsignedInteger(*(in.vMultiUnsignedInteger_[i]))));
  for (unsigned int i = 0; i < in.GetNInteger(); i++)
    vInt_.push_back(std::shared_ptr<IntegerParameter>(
        new IntegerParameter(*(in.vInt_[i]))));
  for (unsigned int i = 0; i < in.GetNBool(); i++)
    vBool_.push_back(
        std::shared_ptr<BoolParameter>(new BoolParameter(*(in.vBool_[i]))));
}

ParameterList::~ParameterList() { /* nothing */ }

std::shared_ptr<AbsParameter>
ParameterList::GetParameter(const unsigned int i) const {
  if (i >= GetNParameter())
    throw BadParameter("ParameterList::GetParameter() | Parameter ID=" +
                       std::to_string(i) + " not in list");

  unsigned int pos = 0;
  if (i < vBool_.size())
    return vBool_.at(i - pos);
  pos += vBool_.size();
  if (i < (pos + vInt_.size()))
    return vInt_.at(i - pos);
  pos += vInt_.size();
  if (i < (pos + vDouble_.size()))
    return vDouble_.at(i - pos);
  pos += vDouble_.size();
  if (i < (pos + vComplex_.size()))
    return vComplex_.at(i - pos);
  pos += vComplex_.size();
  if (i < (pos + vMultiDouble_.size()))
    return vMultiDouble_.at(i - pos);
  pos += vMultiDouble_.size();
  if (i < (pos + vMultiComplex_.size()))
    return vMultiComplex_.at(i - pos);
  pos += vMultiComplex_.size();
  if (i < (pos + vMultiUnsignedInteger_.size()))
    return vMultiUnsignedInteger_.at(i - pos);

  throw BadParameter("ParameterList::GetParameter() | Parameter ID=" +
                     std::to_string(i) + " not in list");
}

bool ParameterList::ParameterExists(const std::string parname) const {
  try {
    GetParameter(parname);
    return true;
  } catch (std::exception &ex) {
    return false;
  }
  return false;
}

std::shared_ptr<AbsParameter>
ParameterList::GetParameter(const std::string parname) const {
  return (*FindBoolParameter(parname));
  return (*FindIntegerParameter(parname));
  return (*FindDoubleParameter(parname));
  return (*FindComplexParameter(parname));
  return (*FindMultiDouble(parname));
  return (*FindMultiComplex(parname));
  return (*FindMultiUnsignedInteger(parname));
}

void ParameterList::Append(const ParameterList &addList) {
  for (int i = 0; i < addList.GetNBool(); i++)
    AddParameter(addList.GetBoolParameter(i));
  for (int i = 0; i < addList.GetNInteger(); i++)
    AddParameter(addList.GetIntegerParameter(i));
  for (int i = 0; i < addList.GetNDouble(); i++)
    AddParameter(addList.GetDoubleParameter(i));
  for (int i = 0; i < addList.GetNComplex(); i++)
    AddParameter(addList.GetComplexParameter(i));
  for (int i = 0; i < addList.GetNMultiDouble(); i++)
    AddParameter(addList.GetMultiDouble(i));
  for (int i = 0; i < addList.GetNMultiComplex(); i++)
    AddParameter(addList.GetMultiComplex(i));
  for (int i = 0; i < addList.GetNMultiUnsignedInteger(); i++)
    AddParameter(addList.GetMultiUnsignedInteger(i));
}

void ParameterList::RemoveDuplicates() {
  // We loop over the list of parameters and delete some. Therefore the size
  // changes during to loop. We enclose GetParameter in try...catch
  for (int i = 0; i < GetNParameter(); i++) {
    std::string name;
    try {
      name = GetParameter(i)->GetName();
    } catch (...) {
      continue;
    }
    for (int j = i + 1; j < GetNParameter(); j++) {
      std::string name2;
      std::shared_ptr<AbsParameter> delPar;
      try {
        delPar = GetParameter(j);
        name2 = GetParameter(j)->GetName();
      } catch (...) {
        continue;
      }

      if (name == name2) {
        auto type = delPar->type();
        switch (type) {
        case ParType::BOOL:
          RemoveBool(j);
          break;
        case ParType::INTEGER:
          RemoveInteger(j - vBool_.size());
          break;
        case ParType::DOUBLE:
          RemoveDouble(j - vBool_.size() - vInt_.size());
          break;
        case ParType::COMPLEX:
          RemoveComplex(j - vBool_.size() - vInt_.size() - vComplex_.size());
          break;
        case ParType::MDOUBLE:
          RemoveMultiDouble(j - vBool_.size() - vInt_.size() -
                            vComplex_.size() - vMultiDouble_.size());
          break;
        case ParType::MCOMPLEX:
          RemoveMultiComplex(j - vBool_.size() - vInt_.size() -
                             vComplex_.size() - vMultiDouble_.size() -
                             vMultiComplex_.size());
          break;
        case ParType::MUNSIGNEDINTEGER:
          RemoveMultiUnsignedInteger(j - vBool_.size() - vInt_.size() -
                                     vComplex_.size() - vMultiDouble_.size() -
                                     vMultiComplex_.size() -
                                     vMultiUnsignedInteger_.size());
          break;
        case ParType::UNDEFINED:
          break;
        }
        j--; // decrement if a parameter is removed
      }
    }
  }
}

void ParameterList::AddParameter(std::shared_ptr<AbsParameter> par) {
  switch (par->type()) {
  case ParType::BOOL: {
    AddParameter(std::dynamic_pointer_cast<BoolParameter>(par));
    break;
  }
  case ParType::INTEGER: {
    AddParameter(std::dynamic_pointer_cast<IntegerParameter>(par));
    break;
  }
  case ParType::DOUBLE: {
    AddParameter(std::dynamic_pointer_cast<DoubleParameter>(par));
    break;
  }
  case ParType::COMPLEX: {
    AddParameter(std::dynamic_pointer_cast<ComplexParameter>(par));
    break;
  }
  case ParType::MDOUBLE: {
    AddParameter(std::dynamic_pointer_cast<MultiDouble>(par));
    break;
  }
  case ParType::MCOMPLEX: {
    AddParameter(std::dynamic_pointer_cast<MultiComplex>(par));
    break;
  }
  case ParType::MUNSIGNEDINTEGER: {
    AddParameter(std::dynamic_pointer_cast<MultiUnsignedInteger>(par));
    break;
  }
  default: { break; }
  }
}

std::string const &ParameterList::to_str() {
  make_str();
  return out_;
}

std::ostream &operator<<(std::ostream &os, ParameterList &p) {
  return os << p.to_str();
}

void ParameterList::make_str() {
  std::stringstream oss;

  if (GetNParameter()) {
    oss << "Parameter List:" << std::endl;
    // print list of complex, float, int and bool parameter
    if (vComplex_.size()) {
      oss << "Complex parameters: " << vComplex_.size() << std::endl;
      auto it = vComplex_.begin();
      for (; it != vComplex_.end(); ++it)
        oss << (*it)->GetName() << ": " << (*it)->GetValue() << std::endl;
    }
    if (vDouble_.size()) {
      oss << "Double parameters: " << vDouble_.size() << std::endl;
      auto it = vDouble_.begin();
      for (; it != vDouble_.end(); ++it)
        oss << (*it)->GetName() << ": " << (*it)->GetValue()
            << " fixed=" << (*it)->IsFixed() << " hasError? "
            << (*it)->HasError() << std::endl;
    }
    if (vInt_.size()) {
      oss << "Integer parameters: " << vInt_.size() << std::endl;
      auto it = vInt_.begin();
      for (; it != vInt_.end(); ++it)
        oss << (*it)->GetName() << ": " << (*it)->GetValue() << std::endl;
    }
    if (vBool_.size()) {
      oss << "Boolean parameters: " << vBool_.size() << std::endl;
      auto it = vBool_.begin();
      for (; it != vBool_.end(); ++it)
        oss << (*it)->GetName() << ": " << (*it)->GetValue() << std::endl;
    }
    if (vMultiComplex_.size()) {
      oss << "Multi complex parameters: " << vMultiComplex_.size() << std::endl;
      auto it = vMultiComplex_.begin();
      for (; it != vMultiComplex_.end(); ++it)
        oss << (*it)->GetName() << " size=" << (*it)->GetNValues() << std::endl;
    }
    if (vMultiUnsignedInteger_.size()) {
      oss << "Multi unsigned int parameters: " << vMultiUnsignedInteger_.size()
          << std::endl;
      auto it = vMultiUnsignedInteger_.begin();
      for (; it != vMultiUnsignedInteger_.end(); ++it)
        oss << (*it)->GetName() << " size=" << (*it)->GetNValues() << std::endl;
    }
    if (vMultiDouble_.size()) {
      oss << "Multi double parameters: " << vMultiDouble_.size() << std::endl;
      auto it = vMultiDouble_.begin();
      for (; it != vMultiDouble_.end(); ++it)
        oss << (*it)->GetName() << " size=" << (*it)->GetNValues() << std::endl;
    }
  } else
    oss << "Parameter List empty!";

  out_ = oss.str();
}

//******************************************************************************
//************* Functions to access individual parameter types *****************
//******************************************************************************

//---------- BOOL PARAMETER ----------
std::vector<std::shared_ptr<BoolParameter>>::const_iterator
ParameterList::FindBoolParameter(const std::string name) const {
  auto it = vBool_.begin();
  for (; it != vBool_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter("ParameterList::FindBoolParameter() | Boolean parameter " +
                     name + " can not be found in list!");
}

unsigned int ParameterList::FindBoolId(const std::string name) {
  return (FindBoolParameter(name) - vBool_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<BoolParameter> par) {
  vBool_.push_back(par);
}

std::shared_ptr<BoolParameter>
ParameterList::GetBoolParameter(const std::string parname) const {
  return (*FindBoolParameter(parname));
}

std::shared_ptr<BoolParameter>
ParameterList::GetBoolParameter(const unsigned int i) const {
  if (!(i < vBool_.size())) {
    throw BadParameter(
        "ParameterList::GetBoolParameter() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vBool_.at(i);
}

const bool
ParameterList::GetBoolParameterValue(const std::string parname) const {
  return (*FindBoolParameter(parname))->GetValue();
}

const bool ParameterList::GetBoolParameterValue(const unsigned int i) const {
  return vBool_.at(i)->GetValue();
}

void ParameterList::SetParameterValue(const std::string name,
                                      const bool inVal) {
  (*FindBoolParameter(name))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const bool inVal) {
  if (!(i < vBool_.size())) {
    throw BadParameter("Parameter not in bool list");
    return;
  }
  (vBool_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::RemoveBool(const std::string name) {
  try {
    unsigned int id = FindBoolParameter(name) - vBool_.begin();
    RemoveBool(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveBool() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveBool(const unsigned int id) {
  try {
    vBool_.erase(vBool_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveBool() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}

//******************************************************************************
//---------- INTEGER PARAMETER ----------
std::vector<std::shared_ptr<IntegerParameter>>::const_iterator
ParameterList::FindIntegerParameter(const std::string name) const {
  auto it = vInt_.begin();
  for (; it != vInt_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter(
      "ParameterList::FindIntegerParameter() | Integerean parameter " + name +
      " can not be found in list!");
}

unsigned int ParameterList::FindIntegerId(const std::string name) {
  return (FindIntegerParameter(name) - vInt_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<IntegerParameter> par) {
  vInt_.push_back(par);
}

std::shared_ptr<IntegerParameter>
ParameterList::GetIntegerParameter(const std::string parname) const {
  return (*FindIntegerParameter(parname));
}

std::shared_ptr<IntegerParameter>
ParameterList::GetIntegerParameter(const unsigned int i) const {
  if (!(i < vInt_.size())) {
    throw BadParameter(
        "ParameterList::GetIntegerParameter() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vInt_.at(i);
}

const int
ParameterList::GetIntegerParameterValue(const std::string parname) const {
  return (*FindIntegerParameter(parname))->GetValue();
}

const int ParameterList::GetIntegerParameterValue(const unsigned int i) const {
  return vInt_.at(i)->GetValue();
}

void ParameterList::SetParameterValue(const std::string name, const int inVal) {
  (*FindIntegerParameter(name))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const int inVal) {
  if (!(i < vInt_.size())) {
    throw BadParameter("Parameter not in bool list");
    return;
  }
  (vInt_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::RemoveInteger(const std::string name) {
  try {
    unsigned int id = FindIntegerParameter(name) - vInt_.begin();
    RemoveInteger(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveInteger() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveInteger(const unsigned int id) {
  try {
    vInt_.erase(vInt_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveInteger() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}

//******************************************************************************
//---------- DOUBLE PARAMETER ----------
std::vector<std::shared_ptr<DoubleParameter>>::const_iterator
ParameterList::FindDoubleParameter(const std::string name) const {
  auto it = vDouble_.begin();
  for (; it != vDouble_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter(
      "ParameterList::FindDoubleParameter() | Double parameter " + name +
      " can not be found in list!");
}

unsigned int ParameterList::FindDoubleId(const std::string name) {
  return (FindDoubleParameter(name) - vDouble_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<DoubleParameter> par) {
  vDouble_.push_back(par);
}

std::shared_ptr<DoubleParameter>
ParameterList::GetDoubleParameter(const std::string parname) const {
  return (*FindDoubleParameter(parname));
}

std::shared_ptr<DoubleParameter>
ParameterList::GetDoubleParameter(const unsigned int i) const {
  if (!(i < vDouble_.size())) {
    throw BadParameter(
        "ParameterList::GetDoubleParameter() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vDouble_.at(i);
}

const double
ParameterList::GetDoubleParameterValue(const std::string parname) const {
  return (*FindDoubleParameter(parname))->GetValue();
}

const double
ParameterList::GetDoubleParameterValue(const unsigned int i) const {
  return vDouble_.at(i)->GetValue();
}

void ParameterList::SetParameterValue(const std::string name,
                                      const double inVal) {
  (*FindDoubleParameter(name))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i,
                                      const double inVal) {
  if (!(i < vDouble_.size())) {
    throw BadParameter("Parameter not in bool list");
    return;
  }
  (vDouble_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::RemoveDouble(const std::string name) {
  try {
    unsigned int id = FindDoubleParameter(name) - vDouble_.begin();
    RemoveDouble(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveDouble() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveDouble(const unsigned int id) {
  try {
    vDouble_.erase(vDouble_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveDouble() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}

//******************************************************************************
//---------- COMPLEX PARAMETER ----------
std::vector<std::shared_ptr<ComplexParameter>>::const_iterator
ParameterList::FindComplexParameter(const std::string name) const {
  auto it = vComplex_.begin();
  for (; it != vComplex_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter(
      "ParameterList::FindComplexParameter() | Complexean parameter " + name +
      " can not be found in list!");
}

unsigned int ParameterList::FindComplexId(const std::string name) {
  return (FindComplexParameter(name) - vComplex_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<ComplexParameter> par) {
  vComplex_.push_back(par);
}

std::shared_ptr<ComplexParameter>
ParameterList::GetComplexParameter(const std::string parname) const {
  return (*FindComplexParameter(parname));
}

std::shared_ptr<ComplexParameter>
ParameterList::GetComplexParameter(const unsigned int i) const {
  if (!(i < vComplex_.size())) {
    throw BadParameter(
        "ParameterList::GetComplexParameter() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vComplex_.at(i);
}

const std::complex<double>
ParameterList::GetComplexParameterValue(const std::string parname) const {
  return (*FindComplexParameter(parname))->GetValue();
}

const std::complex<double>
ParameterList::GetComplexParameterValue(const unsigned int i) const {
  return vComplex_.at(i)->GetValue();
}

void ParameterList::SetParameterValue(const std::string name,
                                      const std::complex<double> inVal) {
  (*FindComplexParameter(name))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i,
                                      const std::complex<double> inVal) {
  if (!(i < vComplex_.size())) {
    throw BadParameter("Parameter not in bool list");
    return;
  }
  (vComplex_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::RemoveComplex(const std::string name) {
  try {
    unsigned int id = FindComplexParameter(name) - vComplex_.begin();
    RemoveComplex(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveComplex() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveComplex(const unsigned int id) {
  try {
    vComplex_.erase(vComplex_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveComplex() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}

//******************************************************************************
//---------- MULTIDOUBLE PARAMETER ----------
std::vector<std::shared_ptr<MultiDouble>>::const_iterator
ParameterList::FindMultiDouble(const std::string name) const {
  auto it = vMultiDouble_.begin();
  for (; it != vMultiDouble_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter("ParameterList::FindMultiDouble() | Doubleean parameter " +
                     name + " can not be found in list!");
}

unsigned int ParameterList::FindMultiDoubleId(const std::string name) {
  return (FindMultiDouble(name) - vMultiDouble_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<MultiDouble> par) {
  vMultiDouble_.push_back(par);
}

std::shared_ptr<MultiDouble>
ParameterList::GetMultiDouble(const std::string parname) const {
  return (*FindMultiDouble(parname));
}

std::shared_ptr<MultiDouble>
ParameterList::GetMultiDouble(const unsigned int i) const {
  if (!(i < vMultiDouble_.size())) {
    throw BadParameter(
        "ParameterList::GetMultiDouble() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vMultiDouble_.at(i);
}

void ParameterList::RemoveMultiDouble(const std::string name) {
  try {
    unsigned int id = FindMultiDouble(name) - vMultiDouble_.begin();
    RemoveMultiDouble(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveDouble() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveMultiDouble(const unsigned int id) {
  try {
    vMultiDouble_.erase(vMultiDouble_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveDouble() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}

//******************************************************************************
//---------- MULTICOMPLEX PARAMETER ----------
std::vector<std::shared_ptr<MultiComplex>>::const_iterator
ParameterList::FindMultiComplex(const std::string name) const {
  auto it = vMultiComplex_.begin();
  for (; it != vMultiComplex_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter(
      "ParameterList::FindMultiComplex() | Complexean parameter " + name +
      " can not be found in list!");
}

unsigned int ParameterList::FindMultiComplexId(const std::string name) {
  return (FindMultiComplex(name) - vMultiComplex_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<MultiComplex> par) {
  vMultiComplex_.push_back(par);
}

std::shared_ptr<MultiComplex>
ParameterList::GetMultiComplex(const std::string parname) const {
  return (*FindMultiComplex(parname));
}

std::shared_ptr<MultiComplex>
ParameterList::GetMultiComplex(const unsigned int i) const {
  if (!(i < vMultiComplex_.size())) {
    throw BadParameter(
        "ParameterList::GetMultiComplex() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vMultiComplex_.at(i);
}

void ParameterList::RemoveMultiComplex(const std::string name) {
  try {
    unsigned int id = FindMultiComplex(name) - vMultiComplex_.begin();
    RemoveMultiComplex(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveComplex() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveMultiComplex(const unsigned int id) {
  try {
    vMultiComplex_.erase(vMultiComplex_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveComplex() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}
//******************************************************************************
//---------- MULTIUNSIGNED INT PARAMETER ----------
std::vector<std::shared_ptr<MultiUnsignedInteger>>::const_iterator
ParameterList::FindMultiUnsignedInteger(const std::string name) const {
  auto it = vMultiUnsignedInteger_.begin();
  for (; it != vMultiUnsignedInteger_.end(); ++it) {
    if ((*it)->GetName() == name)
      return it;
  }
  throw BadParameter(
      "ParameterList::FindMultiUnsignedInteger() | UnsignedInteger parameter " +
      name + " can not be found in list!");
}

unsigned int ParameterList::FindMultiUnsignedIntegerId(const std::string name) {
  return (FindMultiUnsignedInteger(name) - vMultiUnsignedInteger_.begin());
}

void ParameterList::AddParameter(std::shared_ptr<MultiUnsignedInteger> par) {
  vMultiUnsignedInteger_.push_back(par);
}

std::shared_ptr<MultiUnsignedInteger>
ParameterList::GetMultiUnsignedInteger(const std::string parname) const {
  return (*FindMultiUnsignedInteger(parname));
}

std::shared_ptr<MultiUnsignedInteger>
ParameterList::GetMultiUnsignedInteger(const unsigned int i) const {
  if (!(i < vMultiUnsignedInteger_.size())) {
    throw BadParameter(
        "ParameterList::GetMultiUnsignedInteger() | Parameter not found: " +
        std::to_string((double long)i));
  }
  return vMultiUnsignedInteger_.at(i);
}

void ParameterList::RemoveMultiUnsignedInteger(const std::string name) {
  try {
    unsigned int id =
        FindMultiUnsignedInteger(name) - vMultiUnsignedInteger_.begin();
    RemoveMultiUnsignedInteger(id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveUnsignedInteger() | Can not remove"
                  "parameter "
               << name << ": " << ex.what();
    throw;
  }
}

void ParameterList::RemoveMultiUnsignedInteger(const unsigned int id) {
  try {
    vMultiUnsignedInteger_.erase(vMultiUnsignedInteger_.begin() + id);
  } catch (BadParameter &ex) {
    LOG(error) << " ParameterList::RemoveUnsignedInteger() | Can not remove"
                  "parameter with ID="
               << id << ": " << ex.what();
    throw;
  }
}
}
