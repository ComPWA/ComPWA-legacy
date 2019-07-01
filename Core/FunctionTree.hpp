// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// FunctionTree class
///

#ifndef _FUNCTIONTREE_HPP_
#define _FUNCTIONTREE_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Core/FitParameter.hpp"
#include "Core/Functions.hpp"
#include "Core/TreeNode.hpp"

namespace ComPWA {

///
/// \class FunctionTree for automatically cashing of mathematical expressions.
/// This class can be used to store a function in a tree-like structure. Parts
/// of the tree that were calculated before and have not been changed are
/// cashed. This reduces the amount of recalculation at evaluation time.
/// This class is mainly a container for TreeNode's. Most funtionality is
/// implemented in TreeNode.
///
class FunctionTree {
public:
  //  FunctionTree(){};

  /// Create FunctionTree with head node.
  FunctionTree(std::string name, std::shared_ptr<ComPWA::Parameter> parameter,
               std::shared_ptr<ComPWA::Strategy> strategy);

  /// Create FunctionTree with head leaf.
  FunctionTree(std::string name, std::shared_ptr<ComPWA::Parameter> parameter);

  /// Create FunctionTree with head leaf.
  FunctionTree(std::string name, double value);

  /// Create FunctionTree with head leaf.
  FunctionTree(std::string name, std::complex<double> value);

  /// Constructor with the top node provided.
  /// \p head is used as FunctionTree's head node
  FunctionTree(std::shared_ptr<ComPWA::TreeNode> head);

  virtual ~FunctionTree();

  /// Add an existing node to FunctionTree. Can be used to link tree's to each
  /// other: simply insert head node of tree A to tree B
  virtual void insertNode(std::shared_ptr<ComPWA::TreeNode> node,
                          std::string parent);

  /// Insert an existing FunctionTree as TreeNode
  virtual void insertTree(std::shared_ptr<ComPWA::FunctionTree> tree,
                          std::string parent);

  /// Create and add a node to the function tree. Adds Top-Down-Linking to the
  /// node. For this node the caching is disabled.
  virtual void createNode(std::string name,
                          std::shared_ptr<ComPWA::Strategy> strategy,
                          std::string parent);

  /// Create and add a node to the function tree. Adds Top-Down-Linking to the
  /// node.
  virtual void createNode(std::string name,
                          std::shared_ptr<ComPWA::Parameter> parameter,
                          std::shared_ptr<ComPWA::Strategy> strategy,
                          std::string parent);

  /// Create a leaf in FunctionTree. If leaf does not exist it is created and
  /// the corresponding linking is added. If it already exists a link to the
  /// existing element is added to \p parent.
  virtual void createLeaf(std::string name,
                          std::shared_ptr<ComPWA::Parameter> parameter,
                          std::string parent);

  /// Create a leaf in FunctionTree. A DoubleParamter is created and added.
  virtual void createLeaf(std::string name, double value, std::string parent);

  /// Create a leaf in FunctionTree. A ComplexParameter is created and added.
  virtual void createLeaf(std::string name, std::complex<double> value,
                          std::string parent);

  /// Get the head of FunctionTree
  virtual std::shared_ptr<ComPWA::TreeNode> head() const { return Head; }

  /// Recalculate those parts of the tree that have been changed.
  virtual std::shared_ptr<ComPWA::Parameter> parameter() {
    return Head->parameter();
  }

  /// Check if FunctionTree is properly linked and some further checks.
  virtual bool sanityCheck();

  /// Fill ParameterList with parameters. The function is intended to be filled
  /// with fit parameters, so we add only FitParameters.
  virtual void fillParameters(ComPWA::ParameterList &list);

  /// Streaming operator
  friend std::ostream &operator<<(std::ostream &out,
                                  const ComPWA::FunctionTree &b) {
    return out << b.head();
  }

  /// Streaming operator
  friend std::ostream &operator<<(std::ostream &out,
                                  std::shared_ptr<ComPWA::FunctionTree> b) {
    return out << b->head();
  }

  /// Print structure of tree and its values.
  /// Print down to \p level, level = -1 print the whole tree
  virtual std::string print(unsigned int level = -1) {
    return Head->print(level);
  };
  /// Helper function to set all nodes to status changed
  virtual void UpdateAll(std::shared_ptr<ComPWA::TreeNode> startNode);

protected:
  // List of child tree's
  // We need to store the childTreee that were added via insertTree() here.
  // Because otherwise the
  // destructor of these tree's would delete the linking of the tree nodes.
  std::vector<std::shared_ptr<ComPWA::FunctionTree>> ChildTrees;

  // Head node storing the absolute result
  std::shared_ptr<ComPWA::TreeNode> Head;

  // Store the TreeNodes is std::map. TreeNodes of childTrees in \p _childTrees
  // are not included here
  std::map<std::string, std::shared_ptr<ComPWA::TreeNode>> Nodes;

  /// Recursive function to get all used NodeNames
  void GetNamesDownward(std::shared_ptr<ComPWA::TreeNode> start,
                        std::vector<std::string> &childNames,
                        std::vector<std::string> &parentNames);

  /// Helper function to recursively add child nodes of a new tree
  virtual void AddChildNodes(std::shared_ptr<ComPWA::TreeNode> startNode);


};

class FunctionTreeInterface {
public:
  virtual ~FunctionTreeInterface(){};
  virtual std::shared_ptr<FunctionTree>
  createFunctionTree(const ParameterList &DataSample,
                     const std::string &suffix) const = 0;
};

} // namespace ComPWA
#endif
