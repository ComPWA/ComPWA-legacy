// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// TreeNode class
///

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

namespace ComPWA {

// Forward decalaration since we want both classes to be friends
class FunctionTree;

///
/// \class TreeNode is the interface for elements of the FunctionTree
/// This class acts as a container for a parameter in a function tree. It has a
/// Strategy to calculate its value and a name.
///
class TreeNode : public std::enable_shared_from_this<ComPWA::TreeNode>,
                 public ComPWA::ParObserver {
  /// We add FunctionTree as a friend so we can declare some functionts as
  /// protected which deal with the linking between parents and children.
  friend class ComPWA::FunctionTree;

public:
  /// Constructor for tree using a \p name, a \p parameter, a \p strategy and
  /// an identifier for the parent node.
  TreeNode(std::string name, std::shared_ptr<ComPWA::AbsParameter> parameter,
           std::shared_ptr<ComPWA::Strategy> strategy,
           std::shared_ptr<ComPWA::TreeNode> parent);

  /// Constructor for multidimensional TreeNode using a \p name, a vector of
  /// \p parameters, a \p strategy and an identifier for the parent node.
  TreeNode(std::string name,
           std::vector<std::shared_ptr<ComPWA::AbsParameter>> &parameters,
           std::shared_ptr<ComPWA::Strategy> strategy,
           std::shared_ptr<ComPWA::TreeNode> parent);

  virtual ~TreeNode();

  /// Check if recalculation is needed (obsolete)
  virtual bool Modified() const { return _changed; };

  /// Flags the node as modified. Should only be called from its child nodes.
  virtual void Update();

  /// Recalculate node and all of its child nodes
  virtual void Recalculate();

  /// Get child parameter at \p position
  virtual std::shared_ptr<ComPWA::AbsParameter> Parameter(unsigned int position = 0);

  /// Get list of child parameters
  virtual std::vector<std::shared_ptr<ComPWA::AbsParameter>> &Parameters();
  
  /// Fill ParameterList with parameters. The function is intended to be filled
  /// with fit parameters, so we add only DoubleParameters.
  virtual void FillParameters(ComPWA::ParameterList &list);

  /// Dimension of node (= number of parameters)
  virtual std::size_t Dimension() const { return _parameters.size(); };

  virtual std::string Name() const { return _name; };

  /// Add link to children list. This function is intended to be used in
  /// debugging and testing.
  virtual void AddChild(std::shared_ptr<ComPWA::TreeNode> childNode);

  /// Add link to parents list. This function is intended to be used in
  /// debugging and testing.
  virtual void AddParent(std::shared_ptr<ComPWA::TreeNode> parentNode);

  /// Get list of child nodes
  virtual std::vector<std::shared_ptr<ComPWA::TreeNode>> &GetChildNodes();

  /// Find node with \p name within all downstream nodes. The first match is
  /// returned. In case no node exisits a NULL pointer is returned.
  virtual std::shared_ptr<ComPWA::TreeNode> FindChildNode(std::string name) const;

  //  /// Find node \p name are return its parameter. The first match is
  //  /// returned. In case no node exisits a NULL pointer is returned.
  //  virtual std::shared_ptr<ComPWA::AbsParameter> ChildValue(std::string name) const;
  //
  //  /** Return value of certain child node
  //   * We go recursively through out tree to find the specified node and
  //   return
  //   * its value. In case
  //   * of a node with multiple values we return the first one. Currently we
  //   assume
  //   * that the variable
  //   * is a std::complex<double> or can be converted to it.
  //   *
  //   * @param name node specifier
  //   * @return current value of node
  //   */
  //  virtual std::complex<double> getChildSingleValue(std::string name) const;

  /// Streaming operator
  friend std::ostream &operator<<(std::ostream &out,
                                  const ComPWA::TreeNode &p) {
    return out << p.Print(-1);
  }

  /// Stream operator
  friend std::ostream &operator<<(std::ostream &os,
                                  std::shared_ptr<ComPWA::TreeNode> p) {
    return os << p->Print(-1);
  }

  /// Print node and its child nodes to std::string. The recursion goes down
  /// until \p level is reached. A \p prefix can be added inorder to create
  /// a tree like output.
  virtual std::string Print(int level = -1, std::string prefix = "") const;

protected:
  /// List of parent nodes
  std::vector<std::shared_ptr<ComPWA::TreeNode>> _parents;

  //// List of child nodes
  std::vector<std::shared_ptr<ComPWA::TreeNode>> _children;

  /// List of child leafes
  std::vector<std::shared_ptr<ComPWA::AbsParameter>> _parameters;

  std::string _name;

  /// Node has changed and needs to call Recalculate()
  bool _changed;

  /// Node strategy. Strategy defines how the node value calculated given its
  /// child nodes and child leafs.
  std::shared_ptr<ComPWA::Strategy> _strat;
  
  /// Add this node to parents children-list
  virtual void LinkParents();

  /// Delete links to child and parent nodes
  virtual void DeleteLinks();
  
  /// Fill list with names of parent nodes
  virtual void FillParentNames(std::vector<std::string> &names) const;

  /// Fill list with names of children nodes
  virtual void FillChildNames(std::vector<std::string> &names) const;
  
};

} // namespace ComPWA

#endif
