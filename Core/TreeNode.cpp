// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <string>
#include <complex>
#include <memory>

#include "Core/TreeNode.hpp"
#include "Core/Functions.hpp"

using namespace ComPWA;

TreeNode::TreeNode(std::string name,
                   std::shared_ptr<ComPWA::Parameter> parameter,
                   std::shared_ptr<Strategy> strategy,
                   std::shared_ptr<TreeNode> parent)
    : Name(name), Parameter(parameter), HasChanged(true), UseCache(true),
      Strat(strategy) {
  if (!parameter)
    UseCache = false;
  if (!parameter && !strategy)
    throw std::runtime_error(
        "TreeNode::TreeNode() | Neither strategy nor parameter given!");

  if (parent) {
    Parents.push_back(parent);
  }
}

TreeNode::TreeNode(std::string name, std::shared_ptr<Strategy> strategy,
                   std::shared_ptr<TreeNode> parent)
    : Name(name), Parameter(std::shared_ptr<ComPWA::Parameter>()),
      HasChanged(true), UseCache(false), Strat(strategy) {

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

std::shared_ptr<ComPWA::Parameter> TreeNode::parameter() {
  if (UseCache && !Parameter)
    throw std::runtime_error("TreeNode::parameter() | Caching is requested but "
                             "Parameter is not initialized!");

  if (!UseCache && !ChildNodes.size())
    throw std::runtime_error("TreeNode::parameter() | Caching is disabled but "
                             "Node is a lead node!");

  // has been changed or is lead node -> return Parameter
  if (Parameter && (!HasChanged || !ChildNodes.size()))
    return Parameter;

  auto result = recalculate();
  
  if(UseCache){
    Parameter = result;
    HasChanged = false;
  }

  return result;
}

std::shared_ptr<ComPWA::Parameter> TreeNode::recalculate() const {
  // has been changed or is lead node -> return Parameter
  if (Parameter && (!HasChanged || !ChildNodes.size()))
    return Parameter;

  std::shared_ptr<ComPWA::Parameter> result;
  if (Parameter)
    result = Parameter;

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

void TreeNode::fillParameters(ComPWA::ParameterList &list) {
  for (auto ch : ChildNodes) {
    ch->fillParameters(list);
  }
  list.addParameter(parameter());
}

std::shared_ptr<TreeNode> TreeNode::findChildNode(std::string name) const {
  std::shared_ptr<TreeNode> node;
  if (!ChildNodes.size())
    node = std::shared_ptr<TreeNode>();
  for (unsigned int i = 0; i < ChildNodes.size(); i++) {
    if (ChildNodes.at(i)->name() == name) {
      return ChildNodes.at(i);
    } else
      node = ChildNodes.at(i)->findChildNode(name);
    if (node)
      return node;
  }
  return node;
}

std::string TreeNode::print(int level, std::string prefix) const {
  std::stringstream oss;
  oss << prefix << Name;

  auto p = recalculate();
  if (!ChildNodes.size()) { // Print leaf nodes
    if ( p->name() != "" )
      oss << " [" << p->name() << "]";
    oss << " = " << p->val_to_str() << std::endl;
  } else { // Print non-leaf nodes
          oss << " [";
    if (!UseCache)
      oss << "-, ";
    oss << ChildNodes.size() << "]";
    if (UseCache && HasChanged)
      oss << " = ?";
    else
      oss << " = " << p->val_to_str() << std::endl;
  }

  // Abort recursion
  if (level == 0)
    return oss.str();
  
  for( auto ch : ChildNodes ) {
    oss << ch->print(level - 1, prefix + ". ");
  }
  return oss.str();
}

void TreeNode::addChild(std::shared_ptr<TreeNode> childNode) {
  ChildNodes.push_back(childNode);
}

void TreeNode::addParent(std::shared_ptr<TreeNode> parentNode) {
  Parents.push_back(parentNode);
  parentNode->ChildNodes.push_back(shared_from_this());
}

void TreeNode::fillParentNames(std::vector<std::string> &names) const {
  for (auto i : Parents) {
    names.push_back(i->name());
  }
}

void TreeNode::linkParents() {
  for (auto p : Parents)
    p->ChildNodes.push_back(shared_from_this());
}

void TreeNode::deleteLinks() {
  ChildNodes.clear();
  Parents.clear();
  if (Parameter)
    this->parameter()->Detach(shared_from_this());
}

std::vector<std::shared_ptr<TreeNode>> &TreeNode::childNodes() {
  return ChildNodes;
}

void TreeNode::fillChildNames(std::vector<std::string> &names) const {
  for (auto ch : ChildNodes)
    names.push_back(ch->name());
}
