// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// ParObserver class
///

#ifndef _PAROBSERVER_HPP_
#define _PAROBSERVER_HPP_

namespace ComPWA {
namespace FunctionTree {

///
/// ParObserver Base class parameter observer.
/// For the use in the function tree, the observer pattern is used.
/// This class takes the role of the Observer. It's implemented by the
/// TreeNode class, which then are able to observe a parameter and note
/// changes.
///
class ParObserver {
public:
  /// Call this function to mark the observing node as modified.
  virtual void update() = 0;
};

} // namespace FunctionTree
} // namespace ComPWA

#endif
