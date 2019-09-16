// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "FunctionTree.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace FunctionTree {


///
/// \class Dummy
/// Strategy which does not do anything. Used to create a DummyNode in
/// FunctionTree.
///
class Dummy : public Strategy {
public:
  Dummy(ParType in) : Strategy(in, "Dummy"){};
  virtual ~Dummy(){};
  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out) {};
};

FunctionTree::FunctionTree(std::string name,
                           std::shared_ptr<Parameter> parameter,
                           std::shared_ptr<Strategy> strategy) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  createNode(name, parameter, strategy, "");
}

FunctionTree::FunctionTree(std::string name,
                           std::shared_ptr<Parameter> parameter) {
                             DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  createLeaf(name, parameter, "");
}

FunctionTree::FunctionTree(std::string name, double value) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  createLeaf(name, value, "");
}

FunctionTree::FunctionTree(std::string name, std::complex<double> value) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  createLeaf(name, value, "");
}

FunctionTree::FunctionTree(std::shared_ptr<TreeNode> head) : Head(head) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  head->addParent(DummyNode);
  Nodes.insert(
      std::pair<std::string, std::shared_ptr<TreeNode>>(head->name(), head));
}

FunctionTree::~FunctionTree() {
  Head->deleteParentLinks(DummyNode);
}

void FunctionTree::insertNode(std::shared_ptr<TreeNode> node,
                              std::string parent) {
  std::string n = node->name();
  bool ret = Nodes
                 .insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
                     node->name(), node))
                 .second;
  if (!ret)
    throw TreeBuildError("FunctionTree::insertNode | Cannot insert node " +
                         node->name() + ", already exists!");

  // Assign new parent to head of new tree, create links
  std::shared_ptr<TreeNode> parentNode;
  try {
    parentNode = Nodes.at(parent);
  } catch (std::out_of_range &ex) {
    LOG(ERROR) << "FunctionTree::insertNode() | Parent node " << parent
               << " not found in FunctionTree!";
    throw;
  }

  // In case of an existing node, it is possible that this node
  // already have parents. Do need to consider this here?
  node->addParent(parentNode);
  
  parentNode->update();
  // Subtree already linked, but need to be added to list of nodes
  AddChildNodes(node);
}

void FunctionTree::insertTree(std::shared_ptr<FunctionTree> tree,
                              std::string parent) {
  insertNode(tree->head(), parent);
  return;
}

void FunctionTree::createNode(std::string name,
                              std::shared_ptr<Parameter> parameter,
                              std::shared_ptr<Strategy> strategy,
                              std::string parent) {

  if (parent == "" && Head)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> newNode, parentNode;
  if (parent == "") // is this a head node?
    parentNode = DummyNode;
  else {
    try {
      parentNode = Nodes.at(parent);
    } catch (std::out_of_range &ex) {
      LOG(ERROR) << "FunctionTree::createNode() | Parent node " << parent
                 << " not found in FunctionTree for node " << name;
      throw;
    }
  }
  newNode = std::shared_ptr<TreeNode>(
      new TreeNode(name, parameter, strategy, parentNode));

  Nodes.insert(
      std::pair<std::string, std::shared_ptr<TreeNode>>(name, newNode));
  newNode->linkParents();
  if (parentNode)
    parentNode->update();
  if (parent == "")
    Head = newNode; // if we created a head redirect pointer
}

void FunctionTree::createNode(std::string name,
                              std::shared_ptr<Strategy> strategy,
                              std::string parent) {
  createNode(name, std::shared_ptr<Parameter>(), strategy, parent);
}

void FunctionTree::createLeaf(std::string name,
                              std::shared_ptr<Parameter> parameter,
                              std::string parent) {
  if (parent == "" && Head)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else {
    try {
      parentNode = Nodes.at(parent);
    } catch (std::out_of_range &ex) {
      LOG(ERROR) << "FunctionTree::createLeaf() | Parent node " << parent
                 << " not found in FunctionTree!";
      throw;
    }
  }

  std::shared_ptr<TreeNode> leaf;

  // check if Leaf already exists
  bool exists = true;
  if (Nodes.find(name) == Nodes.end()) {
    exists = false;
  } else {
    leaf = Nodes.at(name);
  }

  // setup connections
  if (exists) {
    leaf->addParent(parentNode);
    if (parentNode)
      parentNode->update();
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, parameter, std::shared_ptr<Strategy>(), parentNode));
    Nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    parameter->Attach(leaf);
    if (parentNode)
      parentNode->update();
    if (parent == "")
      Head = leaf; // if we created a head, redirect pointer
  }
}

void FunctionTree::createLeaf(std::string name, double value,
                              std::string parent) {
  createLeaf(name, std::make_shared<Value<double>>("", value), parent);
}

void FunctionTree::createLeaf(std::string name, std::complex<double> value,
                              std::string parent) {
  createLeaf(name, std::make_shared<Value<std::complex<double>>>("", value),
             parent);
}

bool FunctionTree::sanityCheck() {
  if (!Head)
    throw std::runtime_error("FunctionTree::sanityCheck() | "
                             "This tree has no head!");

  // collect all children and parent names
  std::vector<std::string> childNames, parentNames;
  GetNamesDownward(Head, childNames, parentNames);

  // check if matches with available nodeNames
  std::vector<std::string> missedChild, missedParent;
  for (auto const &childName : childNames) {
    try {
      Nodes.at(childName);
    } catch (std::out_of_range &ex) {
      missedChild.push_back(childName);
    }
  }
  for (auto const &parentName : parentNames) {
    try {
      Nodes.at(parentName);
    } catch (std::out_of_range &ex) {
      if (parentName == "DummyNode")
        continue;
      missedParent.push_back(parentName);
    }
  }
  if (missedChild.size()) {
    for (auto const &child : missedChild)
      LOG(DEBUG) << "This tree misses a child: " << child;
    return false;
  }
  if (missedParent.size()) {
    for (auto const &parent : missedParent)
      LOG(DEBUG) << "This tree misses a parent: " << parent;
    return false;
  }

  return true;
}

void FunctionTree::fillParameters(ParameterList &list) {
  Head->fillParameters(list);
}

void FunctionTree::GetNamesDownward(std::shared_ptr<TreeNode> start,
                                    std::vector<std::string> &childNames,
                                    std::vector<std::string> &parentNames) {

  start->fillParentNames(parentNames);
  start->fillChildNames(childNames);

  std::vector<std::shared_ptr<TreeNode>> childs = start->childNodes();
  for (size_t i = 0; i < childs.size(); i++)
    GetNamesDownward(childs.at(i), childNames, parentNames);
  return;
}

void FunctionTree::AddChildNodes(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren = startNode->childNodes();

  for (size_t i = 0; i < newChildren.size(); i++) {
    Nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
        newChildren.at(i)->name(), newChildren.at(i)));
    //        if( !ret )
    //            throw TreeBuildError("FunctionTree::addChildNodes | "
    //                                 "Can not insert node
    //                                 "+inNode->getName()+"!");

    AddChildNodes(newChildren.at(i));
  }
  return;
}

void FunctionTree::UpdateAll(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren = startNode->childNodes();
  for (size_t i = 0; i < newChildren.size(); i++) {
    newChildren.at(i)->update();
    UpdateAll(newChildren[i]);
  }
  return;
}

} // namespace FunctionTree
} // namespace ComPWA
