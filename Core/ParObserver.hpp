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
//! Base class parameter observer.
/*! \class ParObserver
 * @file ParObserver.hpp
 * For the use in the function tree, the observer pattern is used. 
 * This class takes the role of the Observer. It's implemented by the 
 * TreeNode class, which then are able to observe a parameter and note 
 * changes.
*/

#ifndef _PAROBSERVER_HPP_
#define _PAROBSERVER_HPP_

namespace ComPWA {

class ParObserver
{
public:
    //! This function gets called by the parameter to inform observing TreeNodes
    virtual void Update() = 0;
};

} /* namespace ComPWA */

#endif
