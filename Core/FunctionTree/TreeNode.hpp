// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// TreeNode class
///

#ifndef COMPWA_FUNCTIONTREE_TREENODE_HPP_
#define COMPWA_FUNCTIONTREE_TREENODE_HPP_

#include <complex>
#include <memory>
#include <string>

#include "Core/FunctionTree/Functions.hpp"
#include "Core/FunctionTree/ParObserver.hpp"
#include "Core/FunctionTree/ParameterList.hpp"

namespace ComPWA {
namespace FunctionTree {

/// TreeNode is the basic building block of the FunctionTree. A FunctionTree
/// is merely a collection of TreeNodes that are connected to each other in a
/// specific way via addNode(). The FunctionTree represents an arbitrary
/// function. structure. Parts of the tree that were calculated before and have
/// not been changed are cached. This reduces the amount of recalculation at
/// evaluation time.
///
/// There are normal TreeNodes which perform a calculation.
/// They have a Strategy to calculate its value.
/// The other type of TreeNodes are leaves. They simply need an output parameter
/// and acts as a data source. Leaves can be created via createLeaf().
class TreeNode : public std::enable_shared_from_this<TreeNode>,
                 public ParObserver {
public:
  /// Constructor for tree using a \p strategy. This will not cache the output
  /// value
  TreeNode(std::shared_ptr<Strategy> strategy);
  /// Constructor for tree using a \p parameter and a \p strategy. The output
  /// value will be cached.
  TreeNode(std::shared_ptr<Parameter> parameter,
           std::shared_ptr<Strategy> strategy = nullptr);

  virtual ~TreeNode();

  void removeExpiredParents();

  void addNode(std::shared_ptr<TreeNode> node);
  void addNodes(std::vector<std::shared_ptr<TreeNode>> nodes);
  void addParent(std::weak_ptr<TreeNode> node);

  /// Obtain parameter of node. In case child nodes have changed, child nodes
  /// are recalculated and Parameter is updated
  std::shared_ptr<Parameter> parameter();

  /// Fill ParameterList with parameters. The function is intended to be filled
  /// with fit parameters, so we add only FitParameters.
  void fillParameters(ParameterList &list);

  /// Flags the node as modified. Should only be called from its child nodes.
  void update();

  friend std::ostream &operator<<(std::ostream &os,
                                  std::shared_ptr<TreeNode> p) {
    return os << p->print(-1);
  }

  /// Print node and its child nodes to std::string. The recursion goes down
  /// until \p level is reached. A \p prefix can be added inorder to create
  /// a tree like output.
  std::string print(int level = -1, std::string prefix = "");

private:
  /// List of parent nodes
  std::vector<std::weak_ptr<TreeNode>> Parents;

  //// List of child nodes
  std::vector<std::shared_ptr<TreeNode>> ChildNodes;

  /// (cached) node value
  std::shared_ptr<ComPWA::FunctionTree::Parameter> OutputParameter;

  /// Node has changed and needs to call recalculate()
  bool HasChanged;

  /// Node strategy. Strategy defines how the node value calculated given its
  /// child nodes and child leafs.
  std::shared_ptr<Strategy> Strat;

  /// Obtain parameter of node. In case child nodes have changed, child nodes
  /// are recalculated and Parameter is updated
  std::shared_ptr<ComPWA::FunctionTree::Parameter> recalculate();
};

std::shared_ptr<TreeNode> createLeaf(std::shared_ptr<Parameter> parameter);

/// helper function to create TreeNode leaves which are constants
template <typename T,
          typename = std::enable_if_t<
              (std::is_same<int, T>::value || std::is_same<double, T>::value ||
               std::is_same<std::complex<double>, T>::value ||
               std::is_same<std::vector<int>, T>::value ||
               std::is_same<std::vector<double>, T>::value ||
               std::is_same<std::vector<std::complex<double>>, T>::value)>>
std::shared_ptr<TreeNode> createLeaf(const T &value, std::string name = "") {
  auto leaf =
      std::make_shared<TreeNode>(std::make_shared<Value<T>>(name, value));
  return leaf;
}

} // namespace FunctionTree
} // namespace ComPWA

#endif
