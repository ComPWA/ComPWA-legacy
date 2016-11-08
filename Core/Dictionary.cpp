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

#include <vector>
#include <map>
#include <memory>
#include <string>

#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "DataReader/Data.hpp"
#include "Core/Amplitude.hpp"

#include "Core/Dictionary.hpp"

namespace ComPWA {

Dictionary::Dictionary() {

}

Dictionary::~Dictionary(){
  /* nothing */
}

std::string Dictionary::introduce(std::shared_ptr<DataReader::Data> inData, std::string inName){
  if(inName=="" || dataNameUsed(inName)){
    inName="test";//TODO: Generate new name, check available
  }


  //dataInfo tmp(inName, inData, inData->getVariableNames());//TODO!!

  //mData_.push_back( tmp );

  return inName;
}

std::string Dictionary::introduce(std::shared_ptr<Amplitude> inAmp, std::string inName){
  if(inName=="" || amplitudeNameUsed(inName)){
    inName="test";//TODO: Generate new name, check available
  }

  //mData_.insert( std::pair<std::string, std::shared_ptr<Amplitude> >(inName,inData) );

  return inName;
}

bool Dictionary::dataNameUsed(std::string inName){
  bool used = false;

  for (std::vector<dataInfo >::iterator it=mData_.begin(); it!=mData_.end(); ++it)
    if( it->name == inName) used = true;

  return used;
}

bool Dictionary::amplitudeNameUsed(std::string inName){
  bool used = false;

  for (std::vector<ampInfo >::iterator it=mAmps_.begin(); it!=mAmps_.end(); ++it)
    if( it->name == inName) used = true;

  return used;
}

} /* namespace ComPWA */
