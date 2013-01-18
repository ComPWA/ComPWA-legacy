// Template utility functions header file -*- C++ -*-
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
// Author: Mike Williams (9/15/2004)
#ifndef _TemplateUtilFuncs_H
#define _TemplateUtilFuncs_H
//_____________________________________________________________________________
// Headers:
#include <complex>
//_____________________________________________________________________________
/** @file TemplateUtilFuncs.h
 *  @brief Definitions of template utility functions.
 *  @author Mike Williams
 *
 * Some utility template functions needed by the template classes in this    
 * package. A specialization for a class should be added to that class'      
 * header file.                                                              
 *                                                                           
 * <b> Example Specialization: </b> 
 * <!--                                                        
 *  template <> inline SomeType conj(const SomeType &__x){...return ...} 
 * 
 * Adding whatever code is appropriate to return the complex conjugate for a 
 * SomeType object.                                                         
 * -->
 * \include TemplateUtilFuncs.ex
 */
//_____________________________________________________________________________

using namespace std;
//_____________________________________________________________________________

// Some utility functions templates:

/// Return @a zero for the type (defaults to numeric 0)
template <typename _Tp> 
inline _Tp zero(const _Tp &__var){return 0;}

/// Return the conjugate for the type (defualts to the variable itself)
template <typename _Tp> inline _Tp conj(const _Tp &__var){return __var;}

/// Return the imaginary part of the type (defaults to 0)
template <typename _Tp> inline _Tp imag(const _Tp &__var){return 0;}

/// Return @a unity for the type (defaults to numeric 1)
template <typename _Tp> 
inline _Tp unity(const _Tp &__var){return 1;}
//_____________________________________________________________________________

#endif /* _TemplateUtilFuncs_H */
