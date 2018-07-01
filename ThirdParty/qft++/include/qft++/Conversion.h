// Conversion class definition file -*- C++ -*-
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
// Author: Mike Williams 
#ifndef _Conversion_H
#define _Conversion_H
//_____________________________________________________________________________
/** @file Conversion.h
 *  @brief Conversion class definition file (also defines IsDerived macro).
 */
//_____________________________________________________________________________
/** @class Conversion
 *  @author Mike Williams
 *
 *  @brief Internal template class used to determine type convertibility.
 *
 * Note: This class follows very closey the work of Andrei Alexandrescu from
 *       <em>Modern C++ Design</em>, Ch. 2, Sec. 7.
 *
 * Conversion<T,U> determines whether or not a converion from T to U exists.
 * If such a conversion does exist, then Conversion<T,U>::exists is @a true,
 * otherwise it's @a false. If T and U are the same type, then 
 * Conversion<T,U>::same_type is @a true, otherwise it's @a false.
 *
 * <b> Example Usage </b>
 *
 * \include Conversion.ex
 */
//_____________________________________________________________________________

//_____________________________________________________________________________

namespace ComPWA {
    namespace QFT {
// default template:
template <typename _T,typename _U> class Conversion {

private:

#ifndef DOXYGEN_SKIP_THIS
  class _big {char dummy[2];};
#endif  
  typedef char _small_type; ///< a type with size 1
  typedef _big _big_type;  ///< a type with size > 1

  static _small_type _Test(_U); ///< return _small_type
  static _big_type _Test(...); ///< return _big_type
  
  static _T _MakeT(); ///< make an object of type _T

public:

  /// Does a converion of @a T to @a U exist?
  static const bool exists = sizeof(_Test(_MakeT())) == sizeof(_small_type);

  /// Are @a T and @a U the same type?
  static const bool same_type = false;
};
//_____________________________________________________________________________
#ifndef DOXYGEN_SKIP_THIS

// specialization for _T = _U:

template <typename _T> class Conversion<_T,_T> {

public:

  static const bool exists = true;
  static const bool same_type = true;
};
#endif /* DOXYGEN_SKIP_THIS */
//_____________________________________________________________________________
  
/// Does @a Derived inherit from @a Base? (note IsDerived(T,T) is @a true)
#define IsDerived(Derived,Base) (Conversion<const Derived*,const Base*>::exists && !Conversion<const Derived*,const void*>::same_type)

#include "qft++/Type.h"
/// Should @a T be passed by value?
#define PassByValue(__T) Conversion<__T,typename Type<__T>::ParamType>::same_type
//_____________________________________________________________________________

}    
}
#endif /* _Conversion_H */
