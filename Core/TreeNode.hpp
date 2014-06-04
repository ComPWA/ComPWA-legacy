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
//! TreeNode is the interface for elements of the FunctionTree
/*! \class TreeNode
 * @file TreeNode.hpp
 * This class acts as a container for a parameter in a function tree. It has a
 * Strategy to calculate its value and a unique name.
*/

#ifndef _TREENODE_HPP_
#define _TREENODE_HPP_

#include <string>
#include <complex>
#include <memory>

#include "Core/Functions.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParObserver.hpp"

class TreeNode : public std::enable_shared_from_this<TreeNode>, public ParObserver {
public:
  //! Standard constructor
   /*!
    * Standard constructor for one dimensional tree using a string the identifying name
    * /param name unique name of this node
    * /param intResult start value of node
    * /param strat strategy how this node is calculated
    * /param parent pointer to connected upper level node
   */
  TreeNode(std::string name, std::shared_ptr<AbsParameter> intResult, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent);

  //! Standard constructor
   /*!
    * Standard constructor for multi-dimensional tree using a string the identifying name
    * /param name unique name of this node
    * /param intResult start value of node
    * /param strat strategy how this node is calculated
    * /param parent pointer to connected upper level node
   */
  TreeNode(std::string name, std::vector<std::shared_ptr<AbsParameter>>& intResult, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent);

  //! Destructor
  ~TreeNode();

  //! Add this node to parents children-list
  void linkParents(){
    for(unsigned int i=0; i<_parents.size(); i++)
      _parents[i]->_children.push_back(shared_from_this());
  }

  //! Check if recalculation is needed
  inline bool needsCalculation(){
    return _changed;
  };

  //inline void changeVal(const double& newVal){
  //  _value=newVal;
  //  for(unsigned int i=0; i<_parents.size(); i++)
  //    _parents[i]->Update();
    //changed=true;
  //};

  //! Trigger flagged to be recalculated
  void Update();

  //! Trigger recalculation
  void recalculate();

  //! Get pointer to node value
  //const std::shared_ptr<AbsParameter> getValue(){
  //  return _value[0];
  //};

  //! Get pointer to node value
  const std::shared_ptr<AbsParameter> getValue(unsigned int ele=0){
    return _value[ele];
  };

  //! Get dimension
  const unsigned int getDim(){
    return _value.size();
  };


  //! Get node name
  const std::string& getName(){
    return _name;
  };

  //! Add link to children list
  void addChild(std::shared_ptr<TreeNode> newChild){
    _children.push_back(newChild);
  };

  //! Add link to parents list
  void addParent(std::shared_ptr<TreeNode> newParent){
    _parents.push_back(newParent);
    newParent->_children.push_back(shared_from_this());
  };

  //! return parents names
  void getParentNames(std::vector<std::string>& names){
    for(unsigned int i=0; i<_parents.size(); i++)
      names.push_back(_parents[i]->getName());
  };

  //! return children names
  void getChildrenNames(std::vector<std::string>& names){
    for(unsigned int i=0; i<_children.size(); i++)
      names.push_back(_children[i]->getName());
  };

  //! return children pointer
  const std::vector<std::shared_ptr<TreeNode> > getChildren(){
    return _children;
  };

  //! delete children & parent pointer
  const void deleteLinks(){
    _children.clear();
    _parents.clear();
    for(unsigned int i=0; i<_value.size(); i++){
    	_value[i]->Detach(shared_from_this());
    }
  };

  //! Get value of this node
   /*!
    * pure virtual function inheriting classes must implement to provide a
    * value calculated or set in this node
    * /return complex number for this node-value
   */
  //virtual const std::complex<double> getNodeValue() =0; //TODO: complex? Template?

  //! String used to display tree
  std::string to_str(std::string beginning);

  //! Stream-Operator used to display tree
  friend std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p);

protected:
  std::vector<std::shared_ptr<TreeNode> > _parents; /*!< Link to parents */
  std::vector<std::shared_ptr<TreeNode> > _children; /*!< Link to children */

  std::vector<std::shared_ptr<AbsParameter>> _value; /*!< Value of this node */
  std::string _name; /*!< Unique name of this node */
  bool _changed; /*!< flag if node needs recalculation */
  //std::string childOP;

  std::shared_ptr<Strategy> _strat; /*!< Strategy how node calculates its value */
};

//std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p){
//  return os << p->to_str();
//}

#endif /* _TREENODE_HPP_ */
