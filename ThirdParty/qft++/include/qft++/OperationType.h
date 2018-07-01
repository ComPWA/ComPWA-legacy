// OperationType classes definition file. -*- C++ -*-
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
#ifndef _OperationType_H
#define _OperationType_H
//_____________________________________________________________________________
/** @file OperationType.h
 *  @brief OperationType classes definition file (MultType,AddType,etc...)
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

//_____________________________________________________________________________

/** @class MultType
 *  @author Mike Williams
 *
 *  @brief Internal template class used to determine type of V * U.
 *
 * <b> Exmaple Use </b> 
 * \include OperationType.ex                              
 */                                                                           
//_____________________________________________________________________________

template <typename _V,typename _U> class MultType {  

protected:
  static _V _v_type;
  static _U _u_type;

public:
  typedef decltype(_v_type * _u_type) Type; ///< Return type of V * U
};
//_____________________________________________________________________________

/** @class DivType
 *  @author Mike Williams
 *
 *  @brief Internal template class used to determine type of V / U.
 *
 * <b> Exmaple Use </b> 
 * \include OperationType.ex                              
 */                                                                           
//_____________________________________________________________________________

template <typename _V,typename _U> class DivType {

protected:
  static _V _v_type;
  static _U _u_type;
 
public:
  typedef decltype(_v_type / _u_type) Type; ///< Return type of V / U
};
//_____________________________________________________________________________

/** @class AddType
 *  @author Mike Williams
 *
 *  @brief Internal template class used to determine type of V + U.
 *
 * <b> Exmaple Use </b> 
 * \include OperationType.ex                              
 */                                                                           
//_____________________________________________________________________________

template <typename _V,typename _U> class AddType {

protected:
  static _V _v_type;
  static _U _u_type;
 
public:
  typedef decltype(_v_type + _u_type) Type; ///< Return type of V + U
};
//_____________________________________________________________________________

/** @class SubType
 *  @author Mike Williams
 *
 *  @brief Internal template class used to determine type of V - U.
 *
 * <b> Exmaple Use </b> 
 * \include OperationType.ex                              
 */                                                                           
//_____________________________________________________________________________

template <typename _V,typename _U> class SubType {

protected:
  static _V _v_type;
  static _U _u_type;
 
public:
  typedef decltype(_v_type - _u_type) Type; ///< Return type of V - U
};
//_____________________________________________________________________________

    }
}
#endif /* _OperationType_H */
