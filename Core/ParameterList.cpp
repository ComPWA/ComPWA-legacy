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

ParameterList::ParameterList(){
  make_str();
}

ParameterList::ParameterList(const std::vector<std::shared_ptr<DoubleParameter> >& inVec)
:vDoublePar_(inVec){
  make_str();
}

ParameterList::ParameterList(const std::vector<std::shared_ptr<IntegerParameter> >& inVec)
:vIntPar_(inVec){
  make_str();
}

ParameterList::ParameterList(const std::vector<std::shared_ptr<BoolParameter> >& inVec)
:vBoolPar_(inVec){
  make_str();
}

ParameterList::ParameterList(const std::vector<std::shared_ptr<DoubleParameter> >& inD,
    const std::vector<std::shared_ptr<IntegerParameter> >& inI,
    const std::vector<std::shared_ptr<BoolParameter> >& inB)
:vDoublePar_(inD), vIntPar_(inI), vBoolPar_(inB){
  //TODO check names!
  make_str();
}

ParameterList::ParameterList(const ParameterList& in){
  *this = in;
}

ParameterList::~ParameterList() {
  /* nothing */
}

std::shared_ptr<AbsParameter> ParameterList::GetParameter(const unsigned int i) {
  if( !(i < (vBoolPar_.size()+vIntPar_.size()+vDoublePar_.size()) ) ){
      throw BadParameter("Parameter not in list");
      return 0;
  }
  if( i < vDoublePar_.size() ) // is in double list
    return vDoublePar_.at(i);
  else if( i < (vDoublePar_.size()+vIntPar_.size()) ) // is in integer list
    return vIntPar_.at(i-vDoublePar_.size());
  else // is in boolean list
    return vBoolPar_.at(i-vDoublePar_.size()-vIntPar_.size());

  throw BadParameter("Parameter not in list");
  return 0;
}

std::shared_ptr<AbsParameter> ParameterList::GetParameter(const std::string parname) {
  int i=-1;

  try{
      i=mBoolParID_.at(parname);
  }
  catch(...)
  {
      i=-1;
  };
  if(i>-1)
    return vBoolPar_.at(i);

  try{
      i=mIntParID_.at(parname);
  }
  catch(...)
  {
      i=-1;
  };
  if(i>-1)
    return vIntPar_.at(i);

  try{
      i=mDoubleParID_.at(parname);
  }
  catch(...)
  {
      i=-1;
  };
  if(i>-1)
    return vDoublePar_.at(i);


  throw BadParameter("Parameter not found by name: "+parname);
  return 0;
}

std::shared_ptr<DoubleParameter> ParameterList::GetDoubleParameter(const unsigned int i) {
  if( !(i < vDoublePar_.size()) ){
      throw BadParameter("Double Parameter not found");
      //return 0;
  }
  return vDoublePar_.at(i);
}

std::shared_ptr<IntegerParameter> ParameterList::GetIntegerParameter(const unsigned int i) {
  if( !(i < vIntPar_.size()) ){
      throw BadParameter("Integer Parameter not found");
      //return 0;
  }
  return vIntPar_.at(i);
}

std::shared_ptr<BoolParameter> ParameterList::GetBoolParameter(const unsigned int i) {
  if( !(i < vBoolPar_.size()) ){
      throw BadParameter("Bool Parameter not found");
      //return 0;
  }
  return vBoolPar_.at(i);
}

const double ParameterList::GetParameterValue(const unsigned int i) const {
  if( !(i < (vBoolPar_.size()+vIntPar_.size()+vDoublePar_.size()) ) ){
      throw BadParameter("Parameter not in list");
      return 0;
  }
  if( i < vDoublePar_.size() ) // is in double list
    return vDoublePar_.at(i)->GetValue();
  else if( i < (vDoublePar_.size()+vIntPar_.size()) ) // is in integer list
    return (double) vIntPar_.at(i-vDoublePar_.size())->GetValue();
  else // is in boolean list
    return (double) vBoolPar_.at(i-vDoublePar_.size()-vIntPar_.size())->GetValue();

  throw BadParameter("Parameter not in list");
  return 0;
}

std::shared_ptr<DoubleParameter> ParameterList::GetDoubleParameter(const std::string parname) {
  unsigned int i=0;
  try{
      i=mDoubleParID_.at(parname);
  }
  catch(...)
  {
      throw BadParameter("Double Parameter not found");
  };
  return vDoublePar_.at(i);
}

std::shared_ptr<IntegerParameter> ParameterList::GetIntegerParameter(const std::string parname) {
  unsigned int i=0;
  try{
      i=mIntParID_.at(parname);
  }
  catch(...)
  {
      throw BadParameter("Integer Parameter not found");
  };
  return vIntPar_.at(i);
}

std::shared_ptr<BoolParameter> ParameterList::GetBoolParameter(const std::string parname) {
  unsigned int i=0;
  try{
      i=mBoolParID_.at(parname);
  }
  catch(...)
  {
      throw BadParameter("Bool Parameter not found");
  };
  return vBoolPar_.at(i);
}

const double ParameterList::GetParameterValue(const std::string parname) const {
  int i=-1;

  try{
      i=mBoolParID_.at(parname);
  }
  catch(...)
  {
      i=-1;
  };
  if(i>-1)
    return vBoolPar_.at(i)->GetValue();

  try{
      i=mIntParID_.at(parname);
  }
  catch(...)
  {
      i=-1;
  };
  if(i>-1)
    return vIntPar_.at(i)->GetValue();

  try{
      i=mDoubleParID_.at(parname);
  }
  catch(...)
  {
      i=-1;
  };
  if(i>-1)
    return vDoublePar_.at(i)->GetValue();


  throw BadParameter("Parameter not found by name: "+parname);
  return 0;
}

void ParameterList::SetParameterValue(const unsigned int i, const double inVal) {
  if( !(i < vDoublePar_.size() ) ){
      throw BadParameter("Parameter not in double list");
      return ;
  }
  (vDoublePar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const int inVal) {
  if( !(i < vIntPar_.size() ) ){
      throw BadParameter("Parameter not in integer list");
      return ;
  }
  (vIntPar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const bool inVal) {
  if( !(i < vBoolPar_.size() ) ){
      throw BadParameter("Parameter not in bool list");
      return ;
  }
  (vBoolPar_.at(i))->SetValue(inVal);
  return;
}

void ParameterList::AddParameter(std::shared_ptr<AbsParameter> par) {
  //TODO check names!
  switch(par->type())
  {
  case ParType::COMPLEX: //TODO
    //cout << "Easy\n";
    break;
  case ParType::DOUBLE:{
    std::shared_ptr<DoubleParameter> tmp = std::dynamic_pointer_cast<DoubleParameter>(par);
    vDoublePar_.push_back(tmp);
    mDoubleParID_.insert(std::pair<std::string,unsigned int>(par->GetName(),vDoublePar_.size()-1));
    break;}
  case ParType::INTEGER:{
    std::shared_ptr<IntegerParameter> tmp = std::dynamic_pointer_cast<IntegerParameter>(par);
    vIntPar_.push_back(tmp);
    mIntParID_.insert(std::pair<std::string,unsigned int>(par->GetName(),vIntPar_.size()-1));
    break;}
  case ParType::BOOL:{
    std::shared_ptr<BoolParameter> tmp = std::dynamic_pointer_cast<BoolParameter>(par);
    vBoolPar_.push_back(tmp);
    mBoolParID_.insert(std::pair<std::string,unsigned int>(par->GetName(),vBoolPar_.size()-1));
    break;}
  default:{
    //TODO exception
    //cout << "Invalid Selection\n";
    break;}
  }
 // vDoublePar_.push_back(par);
  make_str();
}

void ParameterList::AddParameter(std::shared_ptr<DoubleParameter> par) {
  //TODO check names!
  vDoublePar_.push_back(par);
  mDoubleParID_.insert(std::pair<std::string,unsigned int>(par->GetName(),vDoublePar_.size()-1));
  make_str();
}

void ParameterList::AddParameter(std::shared_ptr<IntegerParameter> par) {
  //TODO check names!
  vIntPar_.push_back(par);
  mIntParID_.insert(std::pair<std::string,unsigned int>(par->GetName(),vIntPar_.size()-1));
  make_str();
}

void ParameterList::AddParameter(std::shared_ptr<BoolParameter> par) {
  //TODO check names!
  vBoolPar_.push_back(par);
  mBoolParID_.insert(std::pair<std::string,unsigned int>(par->GetName(),vBoolPar_.size()-1));
  make_str();
}

void ParameterList::RemoveDouble(const unsigned int id){
  //TODO: try catch
  vDoublePar_.erase(vDoublePar_.begin()+id);
  for (std::map<std::string,unsigned int>::iterator it = mDoubleParID_.begin(); it != mDoubleParID_.end(); ++it )
    if (it->second == id)
      mDoubleParID_.erase(it);
  make_str();
}

void ParameterList::RemoveInteger(const unsigned int id){
  //TODO: try catch
  vIntPar_.erase(vIntPar_.begin()+id);
  for (std::map<std::string,unsigned int>::iterator it = mIntParID_.begin(); it != mIntParID_.end(); ++it )
    if (it->second == id)
      mIntParID_.erase(it);
  make_str();
}

void ParameterList::RemoveBool(const unsigned int id){
  //TODO: try catch
  vBoolPar_.erase(vBoolPar_.begin()+id);
  for (std::map<std::string,unsigned int>::iterator it = mBoolParID_.begin(); it != mBoolParID_.end(); ++it )
    if (it->second == id)
      mBoolParID_.erase(it);
  make_str();
}

void ParameterList::RemoveDouble(const std::string parName){
  //TODO: try catch
  unsigned int id = mDoubleParID_.find(parName)->second;
  vDoublePar_.erase(vDoublePar_.begin()+id);
  mDoubleParID_.erase(parName);
  make_str();
}

void ParameterList::RemoveInteger(const std::string parName){
  //TODO: try catch
  unsigned int id = mIntParID_.find(parName)->second;
  vIntPar_.erase(vIntPar_.begin()+id);
  mIntParID_.erase(parName);
  make_str();
}

void ParameterList::RemoveBool(const std::string parName){
  //TODO: try catch
  unsigned int id = mBoolParID_.find(parName)->second;
  vBoolPar_.erase(vBoolPar_.begin()+id);
  mBoolParID_.erase(parName);
  make_str();
}

void ParameterList::make_str() {
  std::stringstream oss;

  oss << std::endl << "Parameter List";
  if( !vDoublePar_.size() && !vIntPar_.size() && !vBoolPar_.size() )
    oss << " empty" << std::endl;
  else
    oss << std::endl;
  //print list of double, int and bool parameter
  if(vDoublePar_.size())
    oss << "  " << vDoublePar_.size() << " floating point parameters: " << std::endl;
  for(unsigned int d=0; d< vDoublePar_.size(); d++)
    oss << vDoublePar_[d] << std::endl;
  if(vIntPar_.size())
    oss << "  " << vIntPar_.size() << " integer parameters: " << std::endl;
  for(unsigned int i=0; i< vIntPar_.size(); i++)
    oss << vIntPar_[i] << std::endl;
  if(vBoolPar_.size())
    oss << "  " << vBoolPar_.size() << " boolean parameters: " << std::endl;
  for(unsigned int b=0; b< vBoolPar_.size(); b++)
    oss << vBoolPar_[b] << std::endl;

  out_ = oss.str();
}

std::string const& ParameterList::to_str() {
  make_str();
  return out_;
}

std::ostream & operator<<(std::ostream &os, ParameterList &p){
  return os << p.to_str();
}

