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

ParameterList::ParameterList() {
  //make_str();
}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<ComplexParameter> >& inVec) :
    vComplexPar_(inVec) {
  //make_str();
}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<DoubleParameter> >& inVec) :
    vDoublePar_(inVec) {
  //make_str();
}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<IntegerParameter> >& inVec) :
    vIntPar_(inVec) {
  //make_str();
}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<BoolParameter> >& inVec) :
    vBoolPar_(inVec) {
  //make_str();
}

ParameterList::ParameterList(
    const std::vector<std::shared_ptr<ComplexParameter> >& inC,
    const std::vector<std::shared_ptr<DoubleParameter> >& inD,
    const std::vector<std::shared_ptr<IntegerParameter> >& inI,
    const std::vector<std::shared_ptr<BoolParameter> >& inB) :
    vComplexPar_(inC), vDoublePar_(inD), vIntPar_(inI), vBoolPar_(inB) {
  //TODO check names!
  //make_str();
}

ParameterList::ParameterList(const ParameterList& in) {
  mMultiComplexID_ = in.mMultiComplexID_;
  mMultiDoubleID_ = in.mMultiDoubleID_;
  mDoubleParID_ = in.mDoubleParID_;
  mIntParID_ = in.mIntParID_;
  mBoolParID_ = in.mBoolParID_;
  out_ = in.out_;
  vMultiComplex_.clear();
  vMultiDouble_.clear();
  vDoublePar_.clear();
  vIntPar_.clear();
  vBoolPar_.clear();
  for (unsigned int i = 0; i < in.GetNMultiComplex(); i++)
    vMultiComplex_.push_back(
        std::shared_ptr<MultiComplex>(
            new MultiComplex(*(in.vMultiComplex_[i]))));
  for (unsigned int i = 0; i < in.GetNMultiDouble(); i++)
    vMultiDouble_.push_back(
        std::shared_ptr<MultiDouble>(new MultiDouble(*(in.vMultiDouble_[i]))));
  for (unsigned int i = 0; i < in.GetNDouble(); i++)
    vDoublePar_.push_back(
        std::shared_ptr<DoubleParameter>(
            new DoubleParameter(*(in.vDoublePar_[i]))));
  for (unsigned int i = 0; i < in.GetNInteger(); i++)
    vIntPar_.push_back(
        std::shared_ptr<IntegerParameter>(
            new IntegerParameter(*(in.vIntPar_[i]))));
  for (unsigned int i = 0; i < in.GetNBool(); i++)
    vBoolPar_.push_back(
        std::shared_ptr<BoolParameter>(new BoolParameter(*(in.vBoolPar_[i]))));
}

ParameterList::~ParameterList() {
  /* nothing */
}

std::shared_ptr<AbsParameter> ParameterList::GetParameter(
    const unsigned int i) const {
  if (!(i
      < (vBoolPar_.size() + vIntPar_.size() + vDoublePar_.size()
          + vComplexPar_.size() + vMultiDouble_.size() + vMultiComplex_.size()))) {
    throw BadParameter("Parameter not in list");
    return std::shared_ptr<AbsParameter>();
  }
  if (i < vComplexPar_.size())    // is in complex list
    return vComplexPar_.at(i);
  else if (i < (vComplexPar_.size() + vDoublePar_.size()))    // is in double list
    return vDoublePar_.at(i - vComplexPar_.size());
  else if (i < (vComplexPar_.size() + vDoublePar_.size() + vIntPar_.size()))    // is in integer list
    return vIntPar_.at(i - vComplexPar_.size() - vDoublePar_.size());
  else if (i
      < (vComplexPar_.size() + vDoublePar_.size() + vIntPar_.size()
          + vBoolPar_.size()))    // is in boolean list
    return vBoolPar_.at(
        i - vComplexPar_.size() - vDoublePar_.size() - vIntPar_.size());
  else if (i
      < (vComplexPar_.size() + vDoublePar_.size() + vIntPar_.size()
          + vBoolPar_.size() + vMultiDouble_.size()))    // is in multidouble list
    return vMultiDouble_.at(
        i - vComplexPar_.size() - vDoublePar_.size() - vIntPar_.size()
            - vBoolPar_.size());
  else
    // is in multicomplex list
    return vMultiComplex_.at(
        i - vComplexPar_.size() - vDoublePar_.size() - vIntPar_.size()
            - vBoolPar_.size() - vMultiDouble_.size());

  throw BadParameter("Parameter not in list");
  return std::shared_ptr<AbsParameter>();
}

std::shared_ptr<AbsParameter> ParameterList::GetParameter(
    const std::string parname) const {
  int i = -1;

  try {
    i = mMultiComplexID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vMultiComplex_.at(i);

  try {
    i = mMultiDoubleID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vMultiDouble_.at(i);

  try {
    i = mBoolParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vBoolPar_.at(i);

  try {
    i = mIntParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vIntPar_.at(i);

  try {
    i = mDoubleParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vDoublePar_.at(i);

  try {
    i = mComplexParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vComplexPar_.at(i);

  throw BadParameter("Parameter not found by name: " + parname);
  return std::shared_ptr<AbsParameter>();
}

std::shared_ptr<MultiComplex> ParameterList::GetMultiComplex(
    const unsigned int i) const {
  if (!(i < vMultiComplex_.size())) {
    throw BadParameter(
        "Double Parameter not found: " + std::to_string((double long) i));
    //return 0;
  }
  return vMultiComplex_.at(i);
}

std::shared_ptr<MultiDouble> ParameterList::GetMultiDouble(
    const unsigned int i) const {
  if (!(i < vMultiDouble_.size())) {
    throw BadParameter(
        "Double Parameter not found: " + std::to_string((double long) i));
    //return 0;
  }
  return vMultiDouble_.at(i);
}

std::shared_ptr<ComplexParameter> ParameterList::GetComplexParameter(
    const unsigned int i) const {
  if (!(i < vComplexPar_.size())) {
    throw BadParameter(
        "Complex Parameter not found: " + std::to_string((double long) i));
    //return 0;
  }
  return vComplexPar_.at(i);
}

std::shared_ptr<DoubleParameter> ParameterList::GetDoubleParameter(
    const unsigned int i) const {
  if (!(i < vDoublePar_.size())) {
    throw BadParameter(
        "Double Parameter not found: " + std::to_string((double long) i));
    //return 0;
  }
  return vDoublePar_.at(i);
}

std::shared_ptr<IntegerParameter> ParameterList::GetIntegerParameter(
    const unsigned int i) const {
  if (!(i < vIntPar_.size())) {
    throw BadParameter(
        "Integer Parameter not found: " + std::to_string((double long) i));
    //return 0;
  }
  return vIntPar_.at(i);
}

std::shared_ptr<BoolParameter> ParameterList::GetBoolParameter(
    const unsigned int i) const {
  if (!(i < vBoolPar_.size())) {
    throw BadParameter(
        "Bool Parameter not found: " + std::to_string((double long) i));
    //return 0;
  }
  return vBoolPar_.at(i);
}

const double ParameterList::GetParameterValue(const unsigned int i) const {
  if (!(i
      < (vBoolPar_.size() + vIntPar_.size() + vDoublePar_.size()
          + vComplexPar_.size() + vMultiDouble_.size() + vMultiComplex_.size()))) {
    throw BadParameter("Parameter not in list");
    return 0;
  }
  if (i < vComplexPar_.size())    // is in complex list
    return vComplexPar_.at(i)->GetValue().real();
  else if (i < (vComplexPar_.size() + vDoublePar_.size()))    // is in double list
    return (double) vDoublePar_.at(i - vComplexPar_.size())->GetValue();
  else if (i < (vComplexPar_.size() + vDoublePar_.size() + vIntPar_.size()))    // is in integer list
    return (double) vIntPar_.at(i - vComplexPar_.size() - vDoublePar_.size())->GetValue();
  else if (i
      < (vComplexPar_.size() + vDoublePar_.size() + vIntPar_.size()
          + vBoolPar_.size()))    // is in boolean list
    return (double) vBoolPar_.at(
        i - vComplexPar_.size() - vDoublePar_.size() - vIntPar_.size())->GetValue();
  else if (i
      < (vComplexPar_.size() + vDoublePar_.size() + vIntPar_.size()
          + vBoolPar_.size()) + vBoolPar_.size())    // is in multidouble list
    return (double) vMultiDouble_.at(
        i - vComplexPar_.size() - vDoublePar_.size() - vIntPar_.size()
            - vBoolPar_.size())->GetValue();
  else
    // is in multicomplex list
    return (double) vMultiComplex_.at(
        i - vComplexPar_.size() - vDoublePar_.size() - vIntPar_.size()
            - vBoolPar_.size() - vMultiDouble_.size())->GetValue().real();

  throw BadParameter("Parameter not in list");
  return 0;
}

std::shared_ptr<MultiComplex> ParameterList::GetMultiComplex(
    const std::string parname) const {
  unsigned int i = 0;
  try {
    i = mMultiComplexID_.at(parname);
  }
  catch (...) {
    throw BadParameter("MultiComplex not found: " + parname);
  };
  return vMultiComplex_.at(i);
}

std::shared_ptr<MultiDouble> ParameterList::GetMultiDouble(
    const std::string parname) const {
  unsigned int i = 0;
  try {
    i = mMultiDoubleID_.at(parname);
  }
  catch (...) {
    throw BadParameter("MultiDouble not found: " + parname);
  };
  return vMultiDouble_.at(i);
}

std::shared_ptr<ComplexParameter> ParameterList::GetComplexParameter(
    const std::string parname) const {
  unsigned int i = 0;
  try {
    i = mComplexParID_.at(parname);
  }
  catch (...) {
    throw BadParameter("Complex Parameter not found: " + parname);
  };
  return vComplexPar_.at(i);
}

std::shared_ptr<DoubleParameter> ParameterList::GetDoubleParameter(
    const std::string parname) const {
  unsigned int i = 0;
  try {
    i = mDoubleParID_.at(parname);
  }
  catch (...) {
    throw BadParameter("Double Parameter not found: " + parname);
  };
  return vDoublePar_.at(i);
}

std::shared_ptr<IntegerParameter> ParameterList::GetIntegerParameter(
    const std::string parname) const {
  unsigned int i = 0;
  try {
    i = mIntParID_.at(parname);
  }
  catch (...) {
    throw BadParameter("Integer Parameter not found: " + parname);
  };
  return vIntPar_.at(i);
}

std::shared_ptr<BoolParameter> ParameterList::GetBoolParameter(
    const std::string parname) const {
  unsigned int i = 0;
  try {
    i = mBoolParID_.at(parname);
  }
  catch (...) {
    throw BadParameter("Bool Parameter not found: " + parname);
  };
  return vBoolPar_.at(i);
}

const double ParameterList::GetParameterValue(const std::string parname) const {
  int i = -1;

  try {
    i = mMultiComplexID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return (vMultiComplex_.at(i)->GetValue()).real();

  try {
    i = mMultiDoubleID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vMultiDouble_.at(i)->GetValue();

  try {
    i = mBoolParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vBoolPar_.at(i)->GetValue();

  try {
    i = mIntParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vIntPar_.at(i)->GetValue();

  try {
    i = mDoubleParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vDoublePar_.at(i)->GetValue();

  try {
    i = mComplexParID_.at(parname);
  }
  catch (...) {
    i = -1;
  };
  if (i > -1)
    return vComplexPar_.at(i)->GetValue().real();

  throw BadParameter("Parameter not found by name: " + parname);
  return 0;
}

void ParameterList::SetParameterValue(const unsigned int i,
    const std::complex<double> inVal) {
  if (!(i < vComplexPar_.size())) {
    throw BadParameter("Parameter not in complex list");
    return;
  }
  (vComplexPar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i,
    const double inVal) {
  if (!(i < vDoublePar_.size())) {
    throw BadParameter("Parameter not in double list");
    return;
  }
  (vDoublePar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const int inVal) {
  if (!(i < vIntPar_.size())) {
    throw BadParameter("Parameter not in integer list");
    return;
  }
  (vIntPar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const bool inVal) {
  if (!(i < vBoolPar_.size())) {
    throw BadParameter("Parameter not in bool list");
    return;
  }
  (vBoolPar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::AddParameter(std::shared_ptr<AbsParameter> par) {
  //TODO check names!
  switch (par->type()) {
  case ParType::COMPLEX: {
    std::shared_ptr<ComplexParameter> tmp = std::dynamic_pointer_cast<
        ComplexParameter>(par);
    AddParameter(tmp);
    break;
  }
  case ParType::MCOMPLEX: {
    std::shared_ptr<MultiComplex> tmp = std::dynamic_pointer_cast<MultiComplex>(
        par);
    AddParameter(tmp);
    //cout << "Easy\n";
    break;
  }
  case ParType::MDOUBLE: {
    std::shared_ptr<MultiDouble> tmp = std::dynamic_pointer_cast<MultiDouble>(
        par);
    AddParameter(tmp);
    //cout << "Easy\n";
    break;
  }
  case ParType::DOUBLE: {
    std::shared_ptr<DoubleParameter> tmp = std::dynamic_pointer_cast<
        DoubleParameter>(par);
    AddParameter(tmp);
    break;
  }
  case ParType::INTEGER: {
    std::shared_ptr<IntegerParameter> tmp = std::dynamic_pointer_cast<
        IntegerParameter>(par);
    AddParameter(tmp);
    break;
  }
  case ParType::BOOL: {
    std::shared_ptr<BoolParameter> tmp =
        std::dynamic_pointer_cast<BoolParameter>(par);
    AddParameter(tmp);
    break;
  }
  default: {
    //TODO exception
    //cout << "Invalid Selection\n";
    break;
  }
  }
  // vDoublePar_.push_back(par);
  //make_str();
}

void ParameterList::AddParameter(std::shared_ptr<MultiComplex> par) {
  //TODO check names!
  vMultiComplex_.push_back(par);
  mMultiComplexID_.insert(
      std::pair<std::string, unsigned int>(par->GetName(),
          vMultiDouble_.size() - 1));
  //make_str();
}

void ParameterList::AddParameter(std::shared_ptr<MultiDouble> par) {
  //TODO check names!
  vMultiDouble_.push_back(par);
  mMultiDoubleID_.insert(
      std::pair<std::string, unsigned int>(par->GetName(),
          vMultiDouble_.size() - 1));
  //make_str();
}

void ParameterList::AddParameter(std::shared_ptr<ComplexParameter> par) {
  //TODO check names!
  vComplexPar_.push_back(par);
  mComplexParID_.insert(
      std::pair<std::string, unsigned int>(par->GetName(),
          vComplexPar_.size() - 1));
  //make_str();
}

void ParameterList::AddParameter(std::shared_ptr<DoubleParameter> par) {
  //TODO check names!

  auto result = mDoubleParID_.insert(
      std::pair<std::string, unsigned int>(par->GetName(), vDoublePar_.size()));
  if (result.second)
    vDoublePar_.push_back(par);
  //make_str();
}

void ParameterList::AddParameter(std::shared_ptr<IntegerParameter> par) {
  //TODO check names!
  vIntPar_.push_back(par);
  mIntParID_.insert(
      std::pair<std::string, unsigned int>(par->GetName(),
          vIntPar_.size() - 1));
  //make_str();
}

void ParameterList::AddParameter(std::shared_ptr<BoolParameter> par) {
  //TODO check names!
  vBoolPar_.push_back(par);
  mBoolParID_.insert(
      std::pair<std::string, unsigned int>(par->GetName(),
          vBoolPar_.size() - 1));
  //make_str();
}

void ParameterList::RemoveComplex(const unsigned int id) {
  //TODO: try catch
  vComplexPar_.erase(vComplexPar_.begin() + id);
  for (std::map<std::string, unsigned int>::iterator it =
      mComplexParID_.begin(); it != mComplexParID_.end(); ++it)
    if (it->second == id)
      mComplexParID_.erase(it);
  //make_str();
}

void ParameterList::RemoveDouble(const unsigned int id) {
  //TODO: try catch
  vDoublePar_.erase(vDoublePar_.begin() + id);
  for (std::map<std::string, unsigned int>::iterator it = mDoubleParID_.begin();
      it != mDoubleParID_.end(); ++it)
    if (it->second == id)
      mDoubleParID_.erase(it);
  //make_str();
}

void ParameterList::RemoveInteger(const unsigned int id) {
  //TODO: try catch
  vIntPar_.erase(vIntPar_.begin() + id);
  for (std::map<std::string, unsigned int>::iterator it = mIntParID_.begin();
      it != mIntParID_.end(); ++it)
    if (it->second == id)
      mIntParID_.erase(it);
  //make_str();
}

void ParameterList::RemoveBool(const unsigned int id) {
  //TODO: try catch
  vBoolPar_.erase(vBoolPar_.begin() + id);
  for (std::map<std::string, unsigned int>::iterator it = mBoolParID_.begin();
      it != mBoolParID_.end(); ++it)
    if (it->second == id)
      mBoolParID_.erase(it);
  //make_str();
}

void ParameterList::RemoveComplex(const std::string parName) {
  //TODO: try catch
  unsigned int id = mComplexParID_.find(parName)->second;
  vComplexPar_.erase(vComplexPar_.begin() + id);
  mComplexParID_.erase(parName);
  //make_str();
}

void ParameterList::RemoveDouble(const std::string parName) {
  //TODO: try catch
  unsigned int id = mDoubleParID_.find(parName)->second;
  vDoublePar_.erase(vDoublePar_.begin() + id);
  mDoubleParID_.erase(parName);
  //make_str();
}

void ParameterList::RemoveInteger(const std::string parName) {
  //TODO: try catch
  unsigned int id = mIntParID_.find(parName)->second;
  vIntPar_.erase(vIntPar_.begin() + id);
  mIntParID_.erase(parName);
  //make_str();
}

void ParameterList::RemoveBool(const std::string parName) {
  //TODO: try catch
  unsigned int id = mBoolParID_.find(parName)->second;
  vBoolPar_.erase(vBoolPar_.begin() + id);
  mBoolParID_.erase(parName);
  //make_str();
}

void ParameterList::make_str() {
  std::stringstream oss;

  oss << std::endl << "Parameter List";
  if (!vDoublePar_.size() && !vIntPar_.size() && !vBoolPar_.size())
    oss << " empty" << std::endl;
  else
    oss << std::endl;
  //print list of complex, float, int and bool parameter
  if (vComplexPar_.size())
    oss << "  " << vComplexPar_.size() << " complex point parameters: "
        << std::endl;
  for (unsigned int d = 0; d < vComplexPar_.size(); d++)
    oss << vComplexPar_[d]->GetName() << ": " << vComplexPar_[d]->GetValue()
        << std::endl;
  if (vDoublePar_.size())
    oss << "  " << vDoublePar_.size() << " floating point parameters: "
        << std::endl;
  for (unsigned int d = 0; d < vDoublePar_.size(); d++)
    oss << vDoublePar_[d]->GetName() << ": " << vDoublePar_[d]->GetValue()
        << " fixed=" << vDoublePar_[d]->IsFixed() << " hasError? "
        << vDoublePar_[d]->HasError() << std::endl;
  if (vIntPar_.size())
    oss << "  " << vIntPar_.size() << " integer parameters: " << std::endl;
  for (unsigned int i = 0; i < vIntPar_.size(); i++)
    oss << vIntPar_[i]->GetName() << ": " << vIntPar_[i]->GetValue()
        << std::endl;
  if (vBoolPar_.size())
    oss << "  " << vBoolPar_.size() << " boolean parameters: " << std::endl;
  for (unsigned int b = 0; b < vBoolPar_.size(); b++)
    oss << vBoolPar_[b]->GetName() << ": " << vBoolPar_[b]->GetValue()
        << std::endl;

  out_ = oss.str();
}
std::string const& ParameterList::to_str() {
  make_str();
  return out_;
}

std::ostream & operator<<(std::ostream &os, ParameterList &p) {
  return os << p.to_str();
}
void ParameterList::Append(const ParameterList& addList) {
  for (int i = 0; i < addList.GetNParameter(); i++)
    AddParameter(addList.GetParameter(i));
}

} /* namespace ComPWA */
