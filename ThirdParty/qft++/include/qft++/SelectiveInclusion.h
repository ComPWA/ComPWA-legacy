// Templates for selectively including overloaded template functions. -*-C++-*-
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
#ifndef _SelectiveInclusion_H
#define _SelectiveInclusion_H
//_____________________________________________________________________________
/** @file SelectiveInclusion.h
 *  @brief Utilities for selectively including overloaded template functions.
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

//_____________________________________________________________________________
/** @class DisableIf
 *  @author Mike Williams
 *
 *  @brief Used to disable overloaded template function given a condition.
 *
 * DisableIf<condition,T> defines the type Type only if condition is 
 * @a false. Thus, this template can be used to disable an overloaded template
 * function when the condition is @a true.
 *
 * <b> Example </b>
 *
 * \include DisableIf.ex
 */
//_____________________________________________________________________________

// default template: (this will be used for condition = false)

template <bool _condition, typename _Tp = void> class DisableIf {

public:
  /// Only defined if condition is false
  typedef _Tp Type; 
};
//_____________________________________________________________________________
#ifndef DOXYGEN_SKIP_THIS

// specification for true:

template <typename _Tp> class DisableIf<true,_Tp>{};

#endif /* DOXYGEN_SKIP_THIS */
//_____________________________________________________________________________
/** @class EnableIf
 *  @author Mike Williams
 *
 *  @brief Used to enable overloaded template function given a condition.
 *
 * EnableIf<condition,T> defines the type Type only if condition is 
 * @a true. Thus, this template can be used to enable an overloaded template
 * function when the condition is @a true.
 *
 * <b> Example </b>
 *
 * \include EnableIf.ex
 */
//_____________________________________________________________________________

// default template: (this will be used for condition = true)

template <bool _condition, typename _Tp = void> class EnableIf {

public:
  /// Only defined if condition is true
  typedef _Tp Type; 
};
//_____________________________________________________________________________
#ifndef DOXYGEN_SKIP_THIS

// specification for false:

template <typename _Tp> class EnableIf<false,_Tp>{};

#endif /* DOXYGEN_SKIP_THIS */
//_____________________________________________________________________________
}
}
#endif /* _SelectiveInclusion_H */
