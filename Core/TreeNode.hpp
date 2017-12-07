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
  TreeNode(std::string name, std::shared_ptr<ComPWA::Parameter> parameter,
           std::shared_ptr<ComPWA::Strategy> strategy,
           std::shared_ptr<ComPWA::TreeNode> parent);

  /// Constructor for multidimensional TreeNode using a \p name, a vector of
  /// \p parameters, a \p strategy and an identifier for the parent node.
  TreeNode(std::string name,
           std::vector<std::shared_ptr<ComPWA::Parameter>> &parameters,
           std::shared_ptr<ComPWA::Strategy> strategy,
           std::shared_ptr<ComPWA::TreeNode> parent);

  virtual ~TreeNode();

  /// Check if recalculation is needed (obsolete)
  virtual bool isModified() const { return _changed; };

  /// Flags the node as modified. Should only be called from its child nodes.
  virtual void update();

  /// Recalculate node and all of its child nodes
  virtual void recalculate();

  /// Get child parameter at \p position
  virtual std::shared_ptr<ComPWA::Parameter> parameter(unsigned int position);

  /// Get list of child parameters
  virtual std::vector<std::shared_ptr<ComPWA::Parameter>> &parameters();

  /// Fill ParameterList with parameters. The function is intended to be filled
  /// with fit parameters, so we add only DoubleParameters.
  virtual void fillParameters(ComPWA::ParameterList &list);

  /// Dimension of node (= number of parameters)
  virtual std::size_t dimension() const { return _parameters.size(); };

  virtual std::string name() const { return _name; };

  virtual void setName(std::string name) { _name = name; };

  /// Add link to children list. This function is intended to be used in
  /// debugging and testing.
  virtual void addChild(std::shared_ptr<ComPWA::TreeNode> childNode);

  /// Add link to parents list. This function is intended to be used in
  /// debugging and testing.
  virtual void addParent(std::shared_ptr<ComPWA::TreeNode> parentNode);

  /// Get list of child nodes
  virtual std::vector<std::shared_ptr<ComPWA::TreeNode>> &childNodes();

  /// Find node with \p name within all downstream nodes. The first match is
  /// returned. In case no node exisits a NULL pointer is returned.
  virtual std::shared_ptr<ComPWA::TreeNode>
  findChildNode(std::string name) const;

  friend std::ostream &operator<<(std::ostream &out,
                                  const ComPWA::TreeNode &p) {
    return out << p.print(-1);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  std::shared_ptr<ComPWA::TreeNode> p) {
    return os << p->print(-1);
  }

  /// Print node and its child nodes to std::string. The recursion goes down
  /// until \p level is reached. A \p prefix can be added inorder to create
  /// a tree like output.
  virtual std::string print(int level = -1, std::string prefix = "") const;

protected:
  /// List of parent nodes
  std::vector<std::shared_ptr<ComPWA::TreeNode>> _parents;

  //// List of child nodes
  std::vector<std::shared_ptr<ComPWA::TreeNode>> _children;

  /// List of child leafes
  std::vector<std::shared_ptr<ComPWA::Parameter>> _parameters;

  std::string _name;

  /// Node has changed and needs to call recalculate()
  bool _changed;

  /// Node strategy. Strategy defines how the node value calculated given its
  /// child nodes and child leafs.
  std::shared_ptr<ComPWA::Strategy> _strat;

  /// Add this node to parents children-list
  virtual void linkParents();

  /// Delete links to child and parent nodes
  virtual void deleteLinks();

  /// Fill list with names of parent nodes
  virtual void fillParentNames(std::vector<std::string> &names) const;

  /// Fill list with names of children nodes
  virtual void fillChildNames(std::vector<std::string> &names) const;
};

} // ns::ComPWA

#endif
