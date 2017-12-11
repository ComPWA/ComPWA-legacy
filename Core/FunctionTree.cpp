// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/FunctionTree.hpp"
#include "Core/Logging.hpp"

using namespace ComPWA;

FunctionTree::FunctionTree(std::string name,
                           std::shared_ptr<ComPWA::Parameter> parameter,
                           std::shared_ptr<ComPWA::Strategy> strategy) {
  createNode(name, parameter, strategy, "");
}

FunctionTree::FunctionTree(std::string name,
                           std::shared_ptr<ComPWA::Parameter> parameter) {
  createLeaf(name, parameter, "");
}

FunctionTree::FunctionTree(std::string name, double value) {
  createLeaf(name, value, "");
}

FunctionTree::FunctionTree(std::string name, std::complex<double> value) {
  createLeaf(name, value, "");
}

FunctionTree::FunctionTree(std::shared_ptr<ComPWA::TreeNode> head)
    : Head(head) {
  Nodes.insert(std::pair<std::string, std::shared_ptr<ComPWA::TreeNode>>(
      head->name(), head));
}

FunctionTree::~FunctionTree() {
  // We have to delete the links. Otherwise, we have a circular reference
  // on shared pointers and those can not be deleted.
  auto iter = Nodes.begin();
  for (; iter != Nodes.end(); ++iter) {
    iter->second->deleteLinks();
  }
}

//void FunctionTree::createHead(std::string name, double value) {
//  if (_head) // if head exists throw exception
//    throw std::runtime_error(
//        "FunctionTree::createNode() | head node already exists!");
//  createLeaf(name, value, "");
//}
//
//void FunctionTree::createHead(std::string name,
//                              std::shared_ptr<ComPWA::Parameter> parameter) {
//  if (_head) // if head exists throw exception
//    throw std::runtime_error(
//        "FunctionTree::createNode() | head node already exists!");
//  createLeaf(name, parameter, "");
//}
//
//void FunctionTree::createHead(std::string name,
//                              std::shared_ptr<ComPWA::Parameter> parameter,
//                              std::shared_ptr<Strategy> strategy) {
//  if (_head) // if head exists throw exception
//    throw std::runtime_error("FunctionTree::createNode() | "
//                             "Head node already exists!");
//
//  createNode(name, parameter, strategy, "");
//}

void FunctionTree::insertNode(std::shared_ptr<TreeNode> node,
                              std::string parent) {
  std::string n = node->name();
  bool ret = Nodes
                 .insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
                     node->name(), node))
                 .second;
  if (!ret)
    throw TreeBuildError("FunctionTree::insertNode | Can not insert node " +
                         node->name() + "!");

  // Assign new parent to head of new tree, create links
  std::shared_ptr<TreeNode> parentNode;
  try {
    parentNode = Nodes.at(parent);
  } catch (std::out_of_range &ex) {
    LOG(error) << "FunctionTree::insertNode() | Parent node " << parent
               << " not found in FunctionTree!";
    throw;
  }

  // In case of an existing node, it is possible that this node
  // already have parents. Do need to consider this here?
  node->addParent(parentNode);
  //  inNode->linkParents();

  // Subtree already linked, but need to be added to list of nodes
  AddChildNodes(node);
}

void FunctionTree::insertTree(std::shared_ptr<FunctionTree> tree,
                              std::string parent) {
  ChildTrees.push_back(tree);
  insertNode(tree->head(), parent);
  return;
}

void FunctionTree::addNode(std::shared_ptr<ComPWA::TreeNode> newNode) {
  // TODO: check existence, throw exception
  Nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
      newNode->name(), newNode));
  newNode->linkParents();
}

void FunctionTree::createNode(std::string name,
                              std::shared_ptr<ComPWA::Parameter> parameter,
                              std::shared_ptr<ComPWA::Strategy> strategy,
                              std::string parent) {

  if (parent == "" && Head)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> newNode, parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = Nodes.at(parent);

  newNode = std::shared_ptr<TreeNode>(
      new TreeNode(name, parameter, strategy, parentNode));

  Nodes.insert(
      std::pair<std::string, std::shared_ptr<TreeNode>>(name, newNode));
  newNode->linkParents();
  if (parent == "")
    Head = newNode; // if we created a head redirect pointer
}

void FunctionTree::createNode(std::string name,
                              std::shared_ptr<ComPWA::Strategy> strategy,
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
  else
    parentNode = Nodes.at(parent);
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
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, parameter, std::shared_ptr<Strategy>(), parentNode));
    Nodes.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    parameter->Attach(leaf);
    if (parent == "")
      Head = leaf; // if we created a head, redirect pointer
  }
}

void FunctionTree::createLeaf(std::string name, double value,
                              std::string parent) {
  createLeaf(name, std::make_shared<Value<double>>("ParOfNode_" + name, value),
             parent);
}

void FunctionTree::createLeaf(std::string name,
                              std::complex<double> value,
                              std::string parent) {
  createLeaf(name, std::make_shared<Value<std::complex<double>>>(
                       "ParOfNode" + name, value),
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
  for (int i = 0; i < childNames.size(); i++) {
    try {
      Nodes.at(childNames[i]);
    } catch (std::out_of_range &ex) {
      missedChild.push_back(childNames[i]);
    }
  }
  for (int i = 0; i < parentNames.size(); i++) {
    try {
      Nodes.at(parentNames[i]);
    } catch (std::out_of_range &ex) {
      missedParent.push_back(parentNames[i]);
    }
  }
  if (missedChild.size()) {
    for (int i = 0; i < missedChild.size(); i++)
      LOG(debug) << "This tree misses a child: " << missedChild[i];
    return false;
  }
  if (missedParent.size()) {
    for (int i = 0; i < missedParent.size(); i++)
      LOG(debug) << "This tree misses a parent: " << missedParent[i];
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
