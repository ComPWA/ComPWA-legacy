// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Core/Exceptions.hpp"
#include "Data/DataSet.hpp"
#include "FitParameter.hpp"

#include "ParameterList.hpp"

namespace ComPWA {
namespace FunctionTree {

ParameterList::ParameterList(const ComPWA::Data::DataSet &DataSample) {
  // Add data vector to ParameterList
  for (auto x : DataSample.Data)
    addValue(MDouble(x.first, x.second));
  // Adding weight at the end
  addValue(MDouble("Weight", DataSample.Weights));
}

void ParameterList::DeepCopy(const ParameterList &in) {
  IntValues.clear();
  DoubleValues.clear();
  ComplexValues.clear();
  MultiIntValues.clear();
  MultiDoubleValues.clear();
  MultiComplexValues.clear();
  FitParameters.clear();

  for (auto p : in.IntValues)
    IntValues.push_back(std::make_shared<Value<int>>(*p));
  for (auto p : in.DoubleValues)
    DoubleValues.push_back(std::make_shared<Value<double>>(*p));
  for (auto p : in.ComplexValues)
    ComplexValues.push_back(std::make_shared<Value<std::complex<double>>>(*p));

  for (auto p : in.MultiIntValues)
    MultiIntValues.push_back(std::make_shared<Value<std::vector<int>>>(*p));
  for (auto p : in.MultiDoubleValues)
    MultiDoubleValues.push_back(
        std::make_shared<Value<std::vector<double>>>(*p));
  for (auto p : in.MultiComplexValues)
    MultiComplexValues.push_back(
        std::make_shared<Value<std::vector<std::complex<double>>>>(*p));

  for (auto p : in.FitParameters)
    FitParameters.push_back(
        std::make_shared<ComPWA::FunctionTree::FitParameter>(*p));
}

std::size_t ParameterList::numParameters() const {
  return FitParameters.size();
}

void ParameterList::addParameters(
    std::vector<std::shared_ptr<Parameter>> pars) {
  for (auto i : pars)
    addParameter(i);
}

std::shared_ptr<FitParameter>
ParameterList::addUniqueParameter(std::shared_ptr<FitParameter> par) {
  std::shared_ptr<FitParameter> tmp;
  try {
    tmp = FindParameter(par->name(), FitParameters);
  } catch (std::exception &ex) {
    tmp = par;
    FitParameters.push_back(par);
  }

  if (*tmp != *par)
    throw BadParameter("ParameterList::addUniqueParameter() |  FitParameter " +
                       par->name() +
                       " found in list but match is not identical!");

  return tmp;
}

void ParameterList::addParameter(std::shared_ptr<FitParameter> par) {
  FitParameters.push_back(
      std::dynamic_pointer_cast<FunctionTree::FitParameter>(par));
}

void ParameterList::addParameter(std::shared_ptr<Parameter> par) {
  switch (par->type()) {
  case ParType::DOUBLE: {
    addParameter(std::dynamic_pointer_cast<FunctionTree::FitParameter>(par));
    break;
  }
  default: { break; }
  }
}

std::size_t ParameterList::numValues() const {
  return IntValues.size() + DoubleValues.size() + ComplexValues.size() +
         MultiIntValues.size() + MultiDoubleValues.size() +
         MultiComplexValues.size();
}
void ParameterList::addValues(std::vector<std::shared_ptr<Parameter>> values) {
  for (auto i : values)
    addValue(i);
}

void ParameterList::addValue(std::shared_ptr<Parameter> par) {
  switch (par->type()) {
  case ParType::INTEGER: {
    IntValues.push_back(std::dynamic_pointer_cast<Value<int>>(par));
    break;
  }
  case ParType::DOUBLE: {
    DoubleValues.push_back(std::dynamic_pointer_cast<Value<double>>(par));
    break;
  }
  case ParType::COMPLEX: {
    ComplexValues.push_back(
        std::dynamic_pointer_cast<Value<std::complex<double>>>(par));
    break;
  }
  case ParType::MINTEGER: {
    MultiIntValues.push_back(
        std::dynamic_pointer_cast<Value<std::vector<int>>>(par));
    break;
  }
  case ParType::MDOUBLE: {
    MultiDoubleValues.push_back(
        std::dynamic_pointer_cast<Value<std::vector<double>>>(par));
    break;
  }
  case ParType::MCOMPLEX: {
    MultiComplexValues.push_back(
        std::dynamic_pointer_cast<Value<std::vector<std::complex<double>>>>(
            par));
    break;
  }
  default: { break; }
  }
}

std::string ParameterList::to_str() const {
  std::stringstream s;
  if (IntValues.size()) {
    s << "Integer values [" << IntValues.size() << "]:" << std::endl;
    for (auto p : IntValues)
      s << p->to_str() << std::endl;
  }
  if (DoubleValues.size()) {
    s << "Double values [" << DoubleValues.size() << "]:" << std::endl;
    for (auto p : DoubleValues)
      s << p->to_str() << std::endl;
  }
  if (ComplexValues.size()) {
    s << "Complex values [" << ComplexValues.size() << "]:" << std::endl;
    for (auto p : ComplexValues)
      s << p->to_str() << std::endl;
  }
  if (MultiIntValues.size()) {
    s << "Multi integer values [" << MultiIntValues.size() << "]:" << std::endl;
    for (auto p : MultiIntValues)
      s << p->to_str() << std::endl;
  }
  if (MultiDoubleValues.size()) {
    s << "Multi double values [" << MultiDoubleValues.size()
      << "]:" << std::endl;
    for (auto p : MultiDoubleValues)
      s << p->to_str() << std::endl;
  }
  if (MultiComplexValues.size()) {
    s << "Multi complex values [" << MultiComplexValues.size()
      << "]:" << std::endl;
    for (auto p : MultiComplexValues)
      s << p->to_str() << std::endl;
  }
  if (FitParameters.size()) {
    s << "Fit parameters [" << FitParameters.size() << "]:" << std::endl;
    for (auto p : FitParameters)
      s << p->to_str() << std::endl;
  }
  return s.str();
};

} // namespace FunctionTree
} // namespace ComPWA
