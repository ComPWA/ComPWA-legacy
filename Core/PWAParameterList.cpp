#include <vector>

#include "Core/PWAParameterList.hpp"

PWAParameterList::PWAParameterList(){
  make_str();
}

PWAParameterList::PWAParameterList(const std::vector<PWAGenericPar<double> >& inVec)
:vDoublePar_(inVec){
  make_str();
}

PWAParameterList::PWAParameterList(const std::vector<PWAGenericPar<int> >& inVec)
:vIntPar_(inVec){
  make_str();
}

PWAParameterList::PWAParameterList(const std::vector<PWAGenericPar<bool> >& inVec)
:vBoolPar_(inVec){
  make_str();
}

PWAParameterList::PWAParameterList(const std::vector<PWAGenericPar<double> >& inD,
    const std::vector<PWAGenericPar<int> >& inI, const std::vector<PWAGenericPar<bool> >& inB)
:vDoublePar_(inD), vIntPar_(inI), vBoolPar_(inB){
  make_str();
}

PWAParameterList::PWAParameterList(const PWAParameterList& in){
  *this = in;
}

PWAParameterList::~PWAParameterList() {
  /* nothing */
}

const int PWAParameterList::GetParameter(const unsigned int i, PWAGenericPar<double>& par) const {
  if( i < vDoublePar_.size() ){
      //TODO Exception?
      return 0;
  }
  par = vDoublePar_.at(i);
  return 1;
}

const int PWAParameterList::GetParameter(const unsigned int i, PWAGenericPar<int>& par) const {
  if( i < vIntPar_.size() ){
      //TODO Exception?
      return 0;
  }
  par = vIntPar_.at(i);
  return 1;
}

const int PWAParameterList::GetParameter(const unsigned int i, PWAGenericPar<bool>& par) const {
  if( i < vBoolPar_.size() ){
      //TODO Exception?
      return 0;
  }
  par = vBoolPar_.at(i);
  return 1;
}

void PWAParameterList::AddParameter(PWAGenericPar<double>& par) {
  vDoublePar_.push_back(par);
  make_str();
}

void PWAParameterList::AddParameter(PWAGenericPar<int>& par) {
  vIntPar_.push_back(par);
  make_str();
}

void PWAParameterList::AddParameter(PWAGenericPar<bool>& par) {
  vBoolPar_.push_back(par);
  make_str();
}

void PWAParameterList::make_str() {
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

std::string const& PWAParameterList::to_str() const {
  return out_;
}

std::ostream & operator<<(std::ostream &os, const PWAParameterList &p){
  return os << p.to_str();
}

