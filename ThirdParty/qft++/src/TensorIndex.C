// Author: Mike Williams 
/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "qft++/TensorIndex.h"
//_____________________________________________________________________________
/** @file TensorIndex.C
 *  @brief TensorIndex class source file.
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

TensorIndex& TensorIndex::Permute() {

  if(_rank < 1) return *this;
  int level;
  bool valid = false;
  while(valid == false){
    level = _rank - 1;
    _index += (1 << (level << 1));
    while(((*this)[level] > (int)_rank - 1) && (level >= 0)){
      this->SetIndex(level,0);
      if(level > 0){
	this->SetIndex(level-1,(*this)[level-1] + 1);
	level--; 
      }
      else{
	_index = -1;
	level--;
      }
    }
    valid = true;
    for(unsigned int i = 0; i < _rank; ++i){
      for(unsigned int j = 0; j < i; ++j){			  
	if((j < i) && (*this)[i] == (*this)[j]) valid = false;
      }
    }
  }
  return (*this); 
}
//_____________________________________________________________________________

void TensorIndex::Print(std::ostream &__os){
  __os << "(" ;
  for(unsigned int i = 0; i < _rank; i++){
    __os << (*this)[i];
    if(i < (_rank - 1)) std::cout << "," ;
  }
  __os << ")" << std::endl;
}
//_____________________________________________________________________________
    }
}
