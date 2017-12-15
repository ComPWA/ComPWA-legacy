// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <map>

#include "Core/FitParameter.hpp"
#include "Core/Exceptions.hpp"

#include "Core/ParameterList.hpp"

using namespace ComPWA;

std::size_t ParameterList::numParameters() const {
  return FitParameters.size();
}

void ParameterList::addParameters(
    std::vector<std::shared_ptr<Parameter>> pars) {
  for (auto i : pars)
    addParameter(i);
}

void ParameterList::addParameter(std::shared_ptr<Parameter> par) {
  switch (par->type()) {
  case ParType::DOUBLE: {
    FitParameters.push_back(std::dynamic_pointer_cast<FitParameter>(par));
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
    IntValues.push_back(std::dynamic_pointer_cast<ComPWA::Value<int>>(par));
    break;
  }
  case ParType::DOUBLE: {
    DoubleValues.push_back(
        std::dynamic_pointer_cast<ComPWA::Value<double>>(par));
    break;
  }
  case ParType::COMPLEX: {
    ComplexValues.push_back(
        std::dynamic_pointer_cast<ComPWA::Value<std::complex<double>>>(par));
    break;
  }
  case ParType::MINTEGER: {
    MultiIntValues.push_back(
        std::dynamic_pointer_cast<ComPWA::Value<std::vector<int>>>(par));
    break;
  }
  case ParType::MDOUBLE: {
    MultiDoubleValues.push_back(
        std::dynamic_pointer_cast<ComPWA::Value<std::vector<double>>>(par));
    break;
  }
  case ParType::MCOMPLEX: {
    MultiComplexValues.push_back(
        std::dynamic_pointer_cast<
            ComPWA::Value<std::vector<std::complex<double>>>>(par));
    break;
  }
  default: { break; }
  }
}
