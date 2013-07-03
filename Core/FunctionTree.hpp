//! FunctionTree for the actual Optimization
/*! \class FunctionTree
 * @file FunctionTree.hpp
 * This class can be used to store a function in a tree-like structure. This
 * reduces the amount of recalculation needed when performing an optimization,
 * as most often only parts need to be recalculated in one iteration.
*/

#ifndef _FUNCTIONTREE_HPP_
#define _FUNCTIONTREE_HPP_

#include <vector>
#include <memory>

#include "Core/TreeNode.hpp"

class FunctionTree
{
public:
  //! Standard constructor
   /*!
    * Standard constructor with the top node provided
   */
  FunctionTree(std::shared_ptr<TreeNode> head):head_(head){
    nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(head->getName(),head));
  }

  //! Destructor
  virtual ~FunctionTree(){
    //;
  }

  //! Add node to fcnTree
   /*!
    * Add a Node to the function tree
    * Adds Top-Down linking
    * \param TreeNode to be added
   */
  virtual void addNode(std::shared_ptr<TreeNode> newNode){
    //TODO: check existence, throw exception
    nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(newNode->getName(),newNode));
    newNode->linkParents();
  }

  //! return the head of the tree
   /*!
    * Access to the head element of the tree
    * \return FuntionTreeNode at head of tree
   */
  virtual std::shared_ptr<TreeNode> head(){
    //fcnTree_.
    return head_;
    //TODO: return double? ;
  }

  void recalculate(){
    head_->recalculate();
  }


protected:
  std::shared_ptr<TreeNode> head_;
  std::map<std::string, std::shared_ptr<TreeNode> > nodes_; /*!< boost directed graph used for our function */

};

#endif /* _FUNCTIONTREE_HPP_ */
