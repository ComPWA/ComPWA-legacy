#include <string>
#include <complex>
#include <memory>

#include "Core/TreeNode.hpp"
#include "Core/Functions.hpp"

  TreeNode::TreeNode(std::string name, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent)
    :_value(0.),_name(name),_changed(true),_strat(strat){
    if(parent){
        _parents.push_back(parent);
        //parent->children.push_back(shared_from_this());
    }
  };

  //! Destructor
  TreeNode::~TreeNode(){
    //;
  }

  void TreeNode::update(){ //darf nur von kindern aufgerufen werden!
    //std::vector<double> newVals;
    //for(unsigned int i=0; i<children.size(); i++){
    //    newVals.push_back(children[i]->value);
    //}  //end children-loop
    //changeVal(myStrat->execute(newVals));
    for(unsigned int i=0; i<_parents.size(); i++)
      _parents[i]->update();
    _changed=true;
  }; //end update()

  void TreeNode::recalculate(){
    if(_children.size()<1){
      _changed=false;
      return;
    }
    std::vector<double> newVals;
    for(unsigned int i=0; i<_children.size(); i++){
        if(_children[i]->needsCalculation())
          _children[i]->recalculate();
        newVals.push_back(_children[i]->getValue());
    }  //end children-loop
    _value = _strat->execute(newVals);
    _changed=false;
  }; //end update()

  std::string TreeNode::to_str(std::string beginning = ""){
    std::stringstream oss;
    if(_changed && _children.size())
      oss << beginning << _name << " = ?";
    else
      oss << beginning << _name << " = " << _value;
    if(_children.size())
      oss << " with " << _children.size() << " children" << std::endl;
    else
      oss << std::endl;

    for(unsigned int i=0; i<_children.size(); i++){
      //oss << " -> ";
      oss << beginning << _children[i]->to_str(" -> ");
    }
    return oss.str();
  };


std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p){
  return os << p->to_str();
}
