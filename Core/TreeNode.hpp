//! TreeNode is the interface for elements of the FunctionTree
/*! \class TreeNode
 * @file TreeNode.hpp
 * This class provides the interface for all kind of elements which later can
 * be used as nodes in the FunctionTree. There needs to be a unique naming
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
    * Standard constructor using a string the identifying name
    * /param inName unique name of this node
    * /param strat strategy how this node is calculated
    * /param parent pointer to connected upper level node
   */
  TreeNode(std::string name, std::shared_ptr<AbsParameter> intResult, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent);

  //! Destructor
  ~TreeNode();

  void linkParents(){
    for(unsigned int i=0; i<_parents.size(); i++)
      _parents[i]->_children.push_back(shared_from_this());
  }

  //!Check if recalculation is needed
  inline bool needsCalculation(){
    return _changed;
  };

  //inline void changeVal(const double& newVal){
  //  _value=newVal;
  //  for(unsigned int i=0; i<_parents.size(); i++)
  //    _parents[i]->Update();
    //changed=true;
  //};

  void Update();

  void recalculate();

  const std::shared_ptr<AbsParameter> getValue(){
    return _value;
  };

  const std::string& getName(){
    return _name;
  };

  void addChild(std::shared_ptr<TreeNode> newChild){
    _children.push_back(newChild);
  };

  //! Get value of this node
   /*!
    * pure virtual function inheriting classes must implement to provide a
    * value calculated or set in this node
    * /return complex number for this node-value
   */
  //virtual const std::complex<double> getNodeValue() =0; //TODO: complex? Template?

  std::string to_str(std::string beginning);

  friend std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p);

protected:
  std::vector<std::shared_ptr<TreeNode> > _parents;
  std::vector<std::shared_ptr<TreeNode> > _children;

  std::shared_ptr<AbsParameter> _value;
  std::string _name; /*!< Unique name of this node */
  bool _changed;
  //std::string childOP;

  std::shared_ptr<Strategy> _strat;
};

//std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p){
//  return os << p->to_str();
//}

#endif /* _TREENODE_HPP_ */
