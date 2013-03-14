#include <vector>

#include "Core/ParameterList.hpp"

ParameterList::ParameterList(){
  make_str();
}

ParameterList::ParameterList(const std::vector<DoubleParameter>& inVec)
:vDoublePar_(inVec){
  make_str();
}

ParameterList::ParameterList(const std::vector<IntegerParameter>& inVec)
:vIntPar_(inVec){
  make_str();
}

ParameterList::ParameterList(const std::vector<BoolParameter>& inVec)
:vBoolPar_(inVec){
  make_str();
}

ParameterList::ParameterList(const std::vector<DoubleParameter>& inD,
    const std::vector<IntegerParameter>& inI, const std::vector<BoolParameter>& inB)
:vDoublePar_(inD), vIntPar_(inI), vBoolPar_(inB){
  make_str();
}

ParameterList::ParameterList(const ParameterList& in){
  *this = in;
}

ParameterList::~ParameterList() {
  /* nothing */
}

const DoubleParameter ParameterList::GetDoubleParameter(const unsigned int i) const {
  if( !(i < vDoublePar_.size()) ){
      throw BadParameter("Double Parameter not found");
      return 0;
  }
  return vDoublePar_.at(i);
}

const IntegerParameter ParameterList::GetIntegerParameter(const unsigned int i) const {
  if( !(i < vIntPar_.size()) ){
      throw BadParameter("Integer Parameter not found");
      return 0;
  }
  return vIntPar_.at(i);
}

const BoolParameter ParameterList::GetBoolParameter(const unsigned int i) const {
  if( !(i < vBoolPar_.size()) ){
      throw BadParameter("Bool Parameter not found");
      return 0;
  }
  return vBoolPar_.at(i);
}

const double ParameterList::GetParameterValue(const unsigned int i) const {
  if( !(i < (vBoolPar_.size()+vIntPar_.size()+vDoublePar_.size()) ) ){
      throw BadParameter("Parameter not in list");
      return 0;
  }
  if( i < vDoublePar_.size() ) // is in double list
    return vDoublePar_.at(i).GetValue();
  else if( i < (vDoublePar_.size()+vIntPar_.size()) ) // is in integer list
    return (double) vIntPar_.at(i-vDoublePar_.size()).GetValue();
  else // is in boolean list
    return (double) vBoolPar_.at(i-vDoublePar_.size()-vIntPar_.size()).GetValue();

  throw BadParameter("Parameter not in list");
  return 0;
}

void ParameterList::SetParameterValue(const unsigned int i, const double inVal) {
  if( !(i < vDoublePar_.size() ) ){
      throw BadParameter("Parameter not in double list");
      return ;
  }
  (vDoublePar_.at(i)).SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const int inVal) {
  if( !(i < vIntPar_.size() ) ){
      throw BadParameter("Parameter not in integer list");
      return ;
  }
  (vIntPar_.at(i)).SetValue(inVal);
  return;
}

void ParameterList::SetParameterValue(const unsigned int i, const bool inVal) {
  if( !(i < vBoolPar_.size() ) ){
      throw BadParameter("Parameter not in bool list");
      return ;
  }
  (vBoolPar_.at(i)).SetValue(inVal);
  return;
}

void ParameterList::AddParameter(DoubleParameter& par) {
  vDoublePar_.push_back(par);
  make_str();
}

void ParameterList::AddParameter(IntegerParameter& par) {
  vIntPar_.push_back(par);
  make_str();
}

void ParameterList::AddParameter(BoolParameter& par) {
  vBoolPar_.push_back(par);
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

std::string const& ParameterList::to_str() const {
  return out_;
}

std::ostream & operator<<(std::ostream &os, const ParameterList &p){
  return os << p.to_str();
}

