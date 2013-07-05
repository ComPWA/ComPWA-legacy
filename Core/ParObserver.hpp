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

class ParObserver
{
public:
    //! This function gets called by the parameter to inform observing TreeNodes
    virtual void Update() = 0;
};

#endif
