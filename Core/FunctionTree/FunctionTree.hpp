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

#include "Core/FunctionTree/Functions.hpp"
#include "Core/FunctionTree/TreeNode.hpp"

namespace ComPWA {
namespace FunctionTree {

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
  /// Create FunctionTree with head node.
  FunctionTree(std::string name, std::shared_ptr<Parameter> parameter,
               std::shared_ptr<Strategy> strategy);

  /// Create FunctionTree with head leaf.
  FunctionTree(std::string name, std::shared_ptr<Parameter> parameter);

  /// Create FunctionTree with head leaf.
  FunctionTree(std::string name, double value);

  /// Create FunctionTree with head leaf.
  FunctionTree(std::string name, std::complex<double> value);

  /// Constructor with the top node provided.
  /// \p head is used as FunctionTree's head node
  FunctionTree(std::shared_ptr<TreeNode> head);

  virtual ~FunctionTree();

  /// Add an existing node to FunctionTree. Can be used to link tree's to each
  /// other: simply insert head node of tree A to tree B
  virtual void insertNode(std::shared_ptr<TreeNode> node, std::string parent);

  /// Insert an existing FunctionTree as TreeNode
  virtual void
  insertTree(std::shared_ptr<ComPWA::FunctionTree::FunctionTree> tree,
             std::string parent);

  /// Create and add a node to the function tree. Adds Top-Down-Linking to the
  /// node. For this node the caching is disabled.
  virtual void createNode(std::string name, std::shared_ptr<Strategy> strategy,
                          std::string parent);

  /// Create and add a node to the function tree. Adds Top-Down-Linking to the
  /// node.
  virtual void createNode(std::string name,
                          std::shared_ptr<Parameter> parameter,
                          std::shared_ptr<Strategy> strategy,
                          std::string parent);

  /// Create a leaf in FunctionTree. If leaf does not exist it is created and
  /// the corresponding linking is added. If it already exists a link to the
  /// existing element is added to \p parent.
  virtual void createLeaf(std::string name,
                          std::shared_ptr<Parameter> parameter,
                          std::string parent);

  /// Create a leaf in FunctionTree. A DoubleParamter is created and added.
  virtual void createLeaf(std::string name, double value, std::string parent);

  /// Create a leaf in FunctionTree. A ComplexParameter is created and added.
  virtual void createLeaf(std::string name, std::complex<double> value,
                          std::string parent);

  /// Recalculate those parts of the tree that have been changed.
  virtual std::shared_ptr<Parameter> parameter() { return Head->parameter(); }

  std::shared_ptr<TreeNode> Head;

protected:
  /// DummyNode is inserted artifically increase as parent of the head node. This
  /// is a workaround to ensure that a node (and its childs) are only unlinked of
  /// no other FunctionTree points to it.
  std::shared_ptr<TreeNode> DummyNode;

  /// Recursive function to get all used NodeNames
  void GetNamesDownward(std::shared_ptr<TreeNode> start,
                        std::vector<std::string> &childNames,
                        std::vector<std::string> &parentNames);
};

} // namespace FunctionTree
} // namespace ComPWA
#endif
