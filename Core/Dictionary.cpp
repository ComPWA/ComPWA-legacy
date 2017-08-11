// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <vector>
#include <map>
#include <memory>
#include <string>

#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "DataReader/Data.hpp"

#include "Core/Dictionary.hpp"

namespace ComPWA {

Dictionary::Dictionary() {}

Dictionary::~Dictionary() { /* nothing */ }

std::string Dictionary::introduce(std::shared_ptr<DataReader::Data> inData,
                                  std::string inName) {
  if (inName == "" || dataNameUsed(inName)) {
    inName = "test"; // TODO: Generate new name, check available
  }

  // dataInfo tmp(inName, inData, inData->getVariableNames());//TODO!!

  // mData_.push_back( tmp );

  return inName;
}

std::string Dictionary::introduce(std::shared_ptr<AmpIntensity> inAmp,
                                  std::string inName) {
  if (inName == "" || amplitudeNameUsed(inName)) {
    inName = "test"; // TODO: Generate new name, check available
  }

  // mData_.insert( std::pair<std::string, std::shared_ptr<Amplitude>
  // >(inName,inData) );

  return inName;
}

bool Dictionary::dataNameUsed(std::string inName) {
  bool used = false;

  for (std::vector<dataInfo>::iterator it = mData_.begin(); it != mData_.end();
       ++it)
    if (it->name == inName)
      used = true;

  return used;
}

bool Dictionary::amplitudeNameUsed(std::string inName) {
  bool used = false;

  for (std::vector<ampInfo>::iterator it = mAmps_.begin(); it != mAmps_.end();
       ++it)
    if (it->name == inName)
      used = true;

  return used;
}

} /* namespace ComPWA */
