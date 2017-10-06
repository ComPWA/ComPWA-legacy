// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/FunctionTree.hpp"
#include "Core/Logging.hpp"

using namespace ComPWA;

FunctionTree::FunctionTree(std::shared_ptr<ComPWA::TreeNode> head)
    : _head(head) {
  _nodes.insert(std::pair<std::string, std::shared_ptr<ComPWA::TreeNode>>(
      head->Name(), head));
}

FunctionTree::~FunctionTree() {
  // We have to delete the links. Otherwise, we have a circular reference
  // on shared pointers and those can not be deleted.
  auto iter = _nodes.begin();
  for (; iter != _nodes.end(); ++iter) {
    iter->second->DeleteLinks();
  }
}

void FunctionTree::CreateHead(const std::string &name, const double value) {
  if (_head) // if head exists throw exception
    throw std::runtime_error(
        "FunctionTree::CreateNode() | head node already exists!");
  CreateLeaf(name, value, "");
}

void FunctionTree::CreateHead(const std::string &name,
                              std::shared_ptr<ComPWA::AbsParameter> parameter) {
  if (_head) // if head exists throw exception
    throw std::runtime_error(
        "FunctionTree::CreateNode() | head node already exists!");
  CreateLeaf(name, parameter, "");
}

void FunctionTree::CreateHead(const std::string &name,
                              std::shared_ptr<Strategy> strategy,
                              unsigned int dimension, bool useVec) {
  if (_head) // if head exists throw exception
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  CreateNode(name, strategy, "", dimension, useVec);
}

void FunctionTree::InsertNode(std::shared_ptr<TreeNode> node,
                              std::string parent) {
  std::string n = node->Name();
  bool ret = _nodes
                 .insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
                     node->Name(), node))
                 .second;
  if (!ret)
    throw TreeBuildError("FunctionTree::insertNode | "
                         "Can not insert node " +
                         node->Name() + "!");

  // Assign new parent to head of new tree, create links
  std::shared_ptr<TreeNode> parentNode;
  try {
    parentNode = _nodes.at(parent);
  } catch (std::out_of_range &ex) {
    LOG(error) << "FunctionTree::insertNode() | Parent node " << parent
               << " not found in FunctionTree!";
    throw;
  }

  // In case of an existing node, it is possible that this node
  // already have parents. Do need to consider this here?
  node->AddParent(parentNode);
  //  inNode->linkParents();

  // Subtree already linked, but need to be added to list of nodes
  AddChildNodes(node);
}

void FunctionTree::InsertTree(std::shared_ptr<FunctionTree> tree,
                              std::string parent) {
  _childTrees.push_back(tree);
  InsertNode(tree->Head(), parent);
  return;
}

void FunctionTree::AddNode(std::shared_ptr<ComPWA::TreeNode> newNode) {
  // TODO: check existence, throw exception
  _nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
      newNode->Name(), newNode));
  newNode->LinkParents();
}

void FunctionTree::CreateNode(const std::string &name,
                              std::shared_ptr<Strategy> strategy,
                              std::string parent, unsigned int dimension,
                              bool useVec) {
  if (dimension == 0)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Zero dimension! ");
  if (parent == "" && _head)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::vector<std::shared_ptr<AbsParameter>> inter;

  switch (strategy->OutType()) {
  case ParType::MCOMPLEX: {
    inter.push_back(std::shared_ptr<AbsParameter>(new MultiComplex(
        "par" + name, std::vector<std::complex<double>>(
                          dimension, std::complex<double>(0., 0.)))));
    break;
  }
  case ParType::MDOUBLE: {
    inter.push_back(std::shared_ptr<AbsParameter>(
        new MultiDouble("par" + name, std::vector<double>(dimension, 0.))));
    break;
  }
  case ParType::COMPLEX: {
    std::complex<double> start(0., 0.);
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<AbsParameter>(
            new ComplexParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<AbsParameter>(
          new ComplexParameter("par" + name, start)));
    }
    break;
  }
  case ParType::DOUBLE: {
    double start(0.);
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<AbsParameter>(
            new DoubleParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<AbsParameter>(
          new DoubleParameter("par" + name, start)));
    }
    break;
  }
  case ParType::INTEGER: {
    int start(0);
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<AbsParameter>(
            new IntegerParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<AbsParameter>(
          new IntegerParameter("par" + name, start)));
    }
    break;
  }
  case ParType::BOOL: {
    bool start = false;
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<AbsParameter>(
            new BoolParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<AbsParameter>(
          new BoolParameter("par" + name, start)));
    }
    break;
  }
  default: {
    throw BadParameter("FunctionTree::createNode() | Bad parameter type " +
                       std::to_string(strategy->OutType()));
  }
  } // end switch

  std::shared_ptr<TreeNode> newNode, parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = _nodes.at(parent);

  if (dimension == 1 && !useVec)
    newNode = std::shared_ptr<TreeNode>(
        new TreeNode(name, inter[0], strategy, parentNode));
  else
    newNode = std::shared_ptr<TreeNode>(
        new TreeNode(name, inter, strategy, parentNode));

  _nodes.insert(
      std::pair<std::string, std::shared_ptr<TreeNode>>(name, newNode));
  newNode->LinkParents();
  if (parent == "")
    _head = newNode; // if we created a head redirect pointer
}

void FunctionTree::CreateLeaf(const std::string name,
                              std::shared_ptr<AbsParameter> parameter,
                              std::string parent) {
  if (parent == "" && _head)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = _nodes.at(parent);
  std::shared_ptr<TreeNode> leaf;

  // check if Leaf already exists
  bool exists = true;
  if (_nodes.find(name) == _nodes.end()) {
    exists = false;
  } else {
    leaf = _nodes.at(name);
  }

  // setup connections
  if (exists) {
    leaf->AddParent(parentNode);
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, parameter, std::shared_ptr<Strategy>(), parentNode));
    _nodes.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->LinkParents();
    parameter->Attach(leaf);
    if (parent == "")
      _head = leaf; // if we created a head, redirect pointer
  }
}

void FunctionTree::CreateLeaf(const std::string name, const double value,
                              std::string parent) {

  std::shared_ptr<DoubleParameter> staticVal(
      new DoubleParameter("ParOfNode_" + name, value));
  CreateLeaf(name, staticVal, parent);
}

void FunctionTree::CreateLeaf(const std::string name,
                              const std::complex<double> value,
                              std::string parent) {

  std::shared_ptr<ComplexParameter> staticVal(
      new ComplexParameter("ParOfNode_" + name, value));
  CreateLeaf(name, staticVal, parent);
}

void FunctionTree::CreateLeaf(
    const std::string name,
    std::vector<std::shared_ptr<AbsParameter>> &parameters,
    std::string parent) {

  for (auto i : parameters) {
    CreateLeaf(name, i, parent);
  }
}

bool FunctionTree::SanityCheck() {
  bool isSane = true;
  if (!_head)
    throw std::runtime_error("FunctionTree::sanityCheck() | "
                             "This tree has no head!");

  // collect all children and parent names
  std::vector<std::string> childNames, parentNames;
  GetNamesDownward(_head, childNames, parentNames);

  // check if matches with available nodeNames
  std::vector<std::string> missedChild, missedParent;
  for (unsigned int i = 0; i < childNames.size(); i++) {
    try {
      _nodes.at(childNames[i]);
    } catch (std::out_of_range &ex) {
      missedChild.push_back(childNames[i]);
    }
  }
  for (unsigned int i = 0; i < parentNames.size(); i++) {
    try {
      _nodes.at(parentNames[i]);
    } catch (std::out_of_range &ex) {
      missedParent.push_back(parentNames[i]);
    }
  }
  if (missedChild.size()) {
    for (unsigned int i = 0; i < missedChild.size(); i++)
      LOG(debug) << "This tree misses a child: " << missedChild[i];
    return false;
  }
  if (missedParent.size()) {
    for (unsigned int i = 0; i < missedParent.size(); i++)
      LOG(debug) << "This tree misses a parent: " << missedParent[i];
    return false;
  }

  return isSane;
}

void FunctionTree::FillParameters(ParameterList &list) {
  _head->FillParameters(list);
}

void FunctionTree::GetNamesDownward(std::shared_ptr<TreeNode> start,
                                    std::vector<std::string> &childNames,
                                    std::vector<std::string> &parentNames) {

  start->FillParentNames(parentNames);
  start->FillChildNames(childNames);

  const std::vector<std::shared_ptr<TreeNode>> childs = start->GetChildNodes();
  for (unsigned int i = 0; i < childs.size(); i++)
    GetNamesDownward(childs[i], childNames, parentNames);
  return;
}

void FunctionTree::AddChildNodes(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren =
      startNode->GetChildNodes();

  for (unsigned int i = 0; i < newChildren.size(); i++) {
    _nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
        newChildren.at(i)->Name(), newChildren.at(i)));
    //        if( !ret )
    //            throw TreeBuildError("FunctionTree::addChildNodes | "
    //                                 "Can not insert node
    //                                 "+inNode->getName()+"!");

    AddChildNodes(newChildren.at(i));
  }
  return;
}

void FunctionTree::UpdateAll(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren =
      startNode->GetChildNodes();
  for (unsigned int i = 0; i < newChildren.size(); i++) {
    newChildren[i]->Update();
    UpdateAll(newChildren[i]);
  }
  return;
}
