// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <complex>
#include <memory>
#include <string>

#include "Functions.hpp"
#include "TreeNode.hpp"

namespace ComPWA {
namespace FunctionTree {

TreeNode::TreeNode(std::string name, std::shared_ptr<Parameter> parameter,
                   std::shared_ptr<Strategy> strategy,
                   std::shared_ptr<TreeNode> parent)
    : Name(name), OutputParameter(parameter), HasChanged(true),
      Strat(strategy) {
  if (!parameter && !strategy)
    throw std::runtime_error(
        "TreeNode::TreeNode() | Neither strategy nor parameter given!");

  if (parent) {
    Parents.push_back(parent);
  }
}

TreeNode::TreeNode(std::string name, std::shared_ptr<Strategy> strategy,
                   std::shared_ptr<TreeNode> parent)
    : Name(name), OutputParameter(std::shared_ptr<Parameter>()),
      HasChanged(true), Strat(strategy) {

  if (!strategy)
    throw std::runtime_error(
        "TreeNode::TreeNode() | Neither strategy nor parameter given!");

  if (parent) {
    Parents.push_back(parent);
  }
}

TreeNode::~TreeNode() {}

void TreeNode::update() {
  for (unsigned int i = 0; i < Parents.size(); i++)
    Parents.at(i)->update();
  HasChanged = true;
};

std::shared_ptr<Parameter> TreeNode::parameter() {
  if (!OutputParameter && !ChildNodes.size())
    throw std::runtime_error("TreeNode::parameter() | Caching is disabled but "
                             "Node is a lead node!");

  // has been changed or is lead node -> return Parameter
  if (OutputParameter && (!HasChanged || !ChildNodes.size()))
    return OutputParameter;

  auto result = recalculate();

  if (OutputParameter) {
    OutputParameter = result;
    HasChanged = false;
  }

  return result;
}

std::shared_ptr<Parameter> TreeNode::recalculate() const {
  // has been changed or is lead node -> return Parameter
  if (OutputParameter && (!HasChanged || !ChildNodes.size()))
    return OutputParameter;

  std::shared_ptr<Parameter> result;
  if (OutputParameter)
    result = OutputParameter;

  ParameterList newVals;
  for (auto ch : ChildNodes) {
    auto p = ch->parameter();
    if (p->isParameter())
      newVals.addParameter(p);
    else
      newVals.addValue(p);
  }
  try {
    Strat->execute(newVals, result);
  } catch (std::exception &ex) {
    LOG(INFO) << "TreeNode::Recalculate() | Strategy " << Strat
              << " failed on node " << name() << ": " << ex.what();
    throw;
  }

  return result;
}

void TreeNode::fillParameters(ParameterList &list) {
  for (auto ch : ChildNodes) {
    ch->fillParameters(list);
  }
  list.addParameter(parameter());
}

std::shared_ptr<TreeNode> TreeNode::findNode(std::string name) {
  if (Name == name)
    return shared_from_this();

  for (unsigned int i = 0; i < ChildNodes.size(); i++) {
    auto n = ChildNodes.at(i)->findNode(name);
    if (n)
      return n;
  }
  return std::shared_ptr<TreeNode>();
}

std::string TreeNode::print(int level, std::string prefix) const {
  std::stringstream oss;
  oss << prefix << Name;

  auto p = recalculate();
  if (!ChildNodes.size()) { // Print leaf nodes
    if (p->name() != "")
      oss << " [" << p->name() << "]";
    oss << " = " << p->val_to_str() << std::endl;
  } else { // Print non-leaf nodes
    oss << " [";
    if (!OutputParameter)
      oss << "-, ";
    oss << ChildNodes.size() << "]";
    if (OutputParameter && HasChanged)
      oss << " = ?";
    else
      oss << " = " << p->val_to_str() << std::endl;
  }

  // Abort recursion
  if (level == 0)
    return oss.str();

  for (auto ch : ChildNodes) {
    oss << ch->print(level - 1, prefix + ". ");
  }
  return oss.str();
}

void TreeNode::fillParentNames(std::vector<std::string> &names) const {
  for (auto i : Parents) {
    names.push_back(i->name());
  }
}

void TreeNode::deleteParentLinks(std::shared_ptr<TreeNode> parent) {
  // Remove parent node from this node's Parents
  auto r = std::find(Parents.begin(), Parents.end(), parent);
  if (r != Parents.end())
    Parents.erase(r);
  
  // If this node does not have any remaining parents we need to delete links
  // down to tree to ensure that shared_ptr's reference count goes to zero.
  // We do this recursively until we arrive at a node which still has other
  // parents.
  if (!Parents.size()) {
    for (auto ch : childNodes()) {
      ch->deleteParentLinks(shared_from_this());
    }
    if (OutputParameter)
      this->parameter()->Detach(shared_from_this());
  }
}

std::vector<std::shared_ptr<TreeNode>> &TreeNode::childNodes() {
  return ChildNodes;
}

void TreeNode::fillChildNames(std::vector<std::string> &names) const {
  for (auto ch : ChildNodes)
    names.push_back(ch->name());
}

} // namespace FunctionTree
} // namespace ComPWA
