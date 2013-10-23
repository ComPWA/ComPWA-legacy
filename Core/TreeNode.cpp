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
#include <string>
#include <complex>
#include <memory>

#include "Core/TreeNode.hpp"
#include "Core/Functions.hpp"

  TreeNode::TreeNode(std::string name, std::shared_ptr<AbsParameter> intResult, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent)
    :_value(intResult),_name(name),_changed(true),_strat(strat){
    if(parent){
        _parents.push_back(parent);
        //parent->children.push_back(shared_from_this());
    }
  };

  //! Destructor
  TreeNode::~TreeNode(){
    //;
  }

  void TreeNode::Update(){ //sollte nur von kindern oder observed objects aufgerufen werden!
    //std::vector<double> newVals;
    //for(unsigned int i=0; i<children.size(); i++){
    //    newVals.push_back(children[i]->value);
    //}  //end children-loop
    //changeVal(myStrat->execute(newVals));
    for(unsigned int i=0; i<_parents.size(); i++)
      _parents[i]->Update();
    _changed=true;
  }; //end update()

  void TreeNode::recalculate(){
    if(_children.size()<1){
      _changed=false;
      return;
    }
    ParameterList newVals;
    for(unsigned int i=0; i<_children.size(); i++){
        if(_children[i]->needsCalculation())
          _children[i]->recalculate();
        //std::shared_ptr<AbsParameter> para = _children[i]->getValue();
        //para->type();
        newVals.AddParameter(_children[i]->getValue());
    }  //end children-loop
    _value = _strat->execute(newVals);
    _changed=false;
  }; //end update()

  std::string TreeNode::to_str(std::string beginning = ""){
    std::stringstream oss;
    if(_changed && _children.size())
      oss << beginning << _name << " = ?";
    else if(_value)
      oss << beginning << _name << " = " << _value->val_to_str();
    if(_children.size())
      oss << " with " << _children.size() << " children" << std::endl;
    else
      oss << std::endl;

    for(unsigned int i=0; i<_children.size(); i++){
      //oss << " -> ";
      oss << _children[i]->to_str(beginning+" -> ");
    }
    return oss.str();
  };


std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p){
  return os << p->to_str();
}
