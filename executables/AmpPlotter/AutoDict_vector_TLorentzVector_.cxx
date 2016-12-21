#include "vector"
#include "TLorentzVector.h"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<TLorentzVector>::iterator;
#pragma link C++ operators vector<TLorentzVector>::const_iterator;
#pragma link C++ operators vector<TLorentzVector>::reverse_iterator;
#endif
#endif
