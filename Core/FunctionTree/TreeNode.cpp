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

TreeNode::TreeNode(std::shared_ptr<Strategy> strategy)
    : OutputParameter(nullptr), HasChanged(true), Strat(strategy) {
  if (!strategy)
    throw std::runtime_error("TreeNode::TreeNode() | No strategy given!");
}

TreeNode::TreeNode(std::shared_ptr<Parameter> parameter,
                   std::shared_ptr<Strategy> strategy)
    : OutputParameter(parameter), HasChanged(false), Strat(strategy) {
  if (!parameter)
    throw std::runtime_error(
        "TreeNode::TreeNode() | No output parameter given!");
  if (strategy) {
    // if a strategy is given its no leaf
    HasChanged = true;
  }
}

TreeNode::~TreeNode() {
  for (auto x : ChildNodes) {
    x->removeExpiredParents();
  }
  if (0 == ChildNodes.size()) { // if leaf
    OutputParameter->detachExpired();
  }
}

void TreeNode::removeExpiredParents() {
  Parents.erase(std::remove_if(Parents.begin(), Parents.end(),
                               [](auto x) { return x.expired(); }),
                Parents.end());
}

void TreeNode::addNode(std::shared_ptr<TreeNode> node) {
  node->addParent(shared_from_this());
  ChildNodes.push_back(node);
}

void TreeNode::addNodes(std::vector<std::shared_ptr<TreeNode>> nodes) {
  auto ThisNode = shared_from_this();
  for (auto node : nodes) {
    node->addParent(ThisNode);
    ChildNodes.push_back(node);
  }
}

void TreeNode::addParent(std::weak_ptr<TreeNode> node) {
  Parents.push_back(node);
}

void TreeNode::update() {
  for (auto x : Parents)
    x.lock()->update();
  HasChanged = true;
};

std::shared_ptr<Parameter> TreeNode::parameter() {
  if (!OutputParameter && 0 == ChildNodes.size())
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

std::shared_ptr<Parameter> TreeNode::recalculate() {
  // has been changed or is leaf node -> return Parameter
  if (OutputParameter && (!HasChanged || 0 == ChildNodes.size()))
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
              << " failed with input ParameterList:\n"
              << newVals << ": " << ex.what();
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

std::string TreeNode::print(int level, std::string prefix) {
  std::stringstream oss;
  oss << prefix;

  auto p = recalculate();
  if (!ChildNodes.size()) { // Print leaf nodes
    if (p->name() != "")
      oss << "[" << p->name() << "]";
    oss << " = " << p->val_to_str() << std::endl;
  } else { // Print non-leaf nodes
    oss << Strat->str() << " [";
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

std::shared_ptr<TreeNode> createLeaf(std::shared_ptr<Parameter> parameter) {
  auto leaf = std::make_shared<TreeNode>(parameter);
  parameter->attach(leaf);
  return leaf;
}

} // namespace FunctionTree
} // namespace ComPWA
