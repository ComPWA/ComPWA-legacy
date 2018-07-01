// Type class definition file -*- C++ -*-
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
#ifndef _Type_H
#define _Type_H
//_____________________________________________________________________________
/** @file Type.h
 *  @brief Type class definition file.
 */
namespace ComPWA {
    namespace QFT {

//_____________________________________________________________________________
/** @class Type
 *  @author Mike Williams
 *
 *  @brief Internal template class used to determine type properties.
 *
 * Type is a template class used to determine properties of the type used to
 * instantiate it. Each instantiation contains the folloing: <br><br>
 *
 * The typedef Type<_Tp>::ParamType is how _Tp should be passed as an argument.
 * For basic types it's the type (Type<float>::ParamType is float), 
 * for pointers it's the pointer, for objects it's a constant reference to the 
 * object. So Type<SomeObject>::ParamType would be <em> const SomeObject& 
 * </em>. 
 *
 * Any type without a specialization will have Type<_Tp>::ParamType as 
 * const _Tp&. 
 *                                                                           
 * <b> Exmaple Use: </b> 
 * <!--
 * template<typename _Tp> void foo(typename Type<_Tp>::ParamType __val){
 *    ...
 *   };                                                            
 * -->                                                        
 * \include Type_ParamType.ex
 *
 * This class also defines a boolean data member, @a IsPointer, which is 
 * @a true when the class is instantized with a pointer and @a false otherwise.
 *                                                                           
 * <b> Specifications Provided: </b><br>
 *                         float,double,long double,int,long int,short int,  
 *                         signed char,unsigned char,unsigned short int,     
 *                         unsigned int,unsigned long int,bool               
 *                                     
 * Any non-pointer/non-basic type will instantiate the default template. Thus,
 * it'll have @a IsPointer set as @a false and @a ParamType will be a constant
 * reference to the type. If this isn't the desired behavior, then a 
 * specification should be defined for the type.
 *
 */                                                                           
//_____________________________________________________________________________

//_____________________________________________________________________________

// default template: 
template <typename _Tp> class Type {
  
 public:

  typedef const _Tp& ParamType; ///< How _Tp should be passed as an argument

  const static bool IsPointer = false; ///< Is _Tp a pointer?

};
//_____________________________________________________________________________
#ifndef DOXYGEN_SKIP_THIS

// primitive specifications:

template<> class Type<float> {
  
 public:
  
  typedef float ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<double> {
  
 public:
  
  typedef double ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<long double> {
  
 public:
  
  typedef long double ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<int> {
  
 public:
  
  typedef int ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<long int> {
  
 public:
  
  typedef long int ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<short int> {
  
 public:
  
  typedef short int ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<signed char> {
  
 public:
  
  typedef signed char ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<unsigned char> {
  
 public:
  
  typedef unsigned char ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<unsigned short int> {
  
 public:
  
  typedef unsigned short int ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<unsigned int> {
  
 public:
  
  typedef unsigned int ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<unsigned long int> {
  
 public:
  
  typedef unsigned long int ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

template<> class Type<bool> {

 public:

  typedef bool ParamType;

  const static bool IsPointer = false;

};
//_____________________________________________________________________________

// partial specialization for pointers:
template<typename _Tp>
class Type<_Tp*> {
 public:
 
  typedef _Tp* ParamType;
  
  const static bool IsPointer = true;

};
//_____________________________________________________________________________
#endif /* DOXYGEN_SKIP_THIS */

    }}
#endif /* _Type_H */
