// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/FunctionTree.hpp"
#include "Core/Logging.hpp"

using namespace ComPWA;
/// Create FunctionTree with head node.
FunctionTree::FunctionTree(const std::string &name,
                           std::shared_ptr<ComPWA::Strategy> strategy,
                           unsigned int dimension, bool useVec) {
  createHead(name, strategy, dimension, useVec);
}

FunctionTree::FunctionTree(std::shared_ptr<ComPWA::TreeNode> head)
    : _head(head) {
  _nodes.insert(std::pair<std::string, std::shared_ptr<ComPWA::TreeNode>>(
      head->name(), head));
}

FunctionTree::~FunctionTree() {
  // We have to delete the links. Otherwise, we have a circular reference
  // on shared pointers and those can not be deleted.
  auto iter = _nodes.begin();
  for (; iter != _nodes.end(); ++iter) {
    iter->second->deleteLinks();
  }
}

void FunctionTree::createHead(const std::string &name, const double value) {
  if (_head) // if head exists throw exception
    throw std::runtime_error(
        "FunctionTree::createNode() | head node already exists!");
  createLeaf(name, value, "");
}

void FunctionTree::createHead(const std::string &name,
                              std::shared_ptr<ComPWA::Parameter> parameter) {
  if (_head) // if head exists throw exception
    throw std::runtime_error(
        "FunctionTree::createNode() | head node already exists!");
  createLeaf(name, parameter, "");
}

void FunctionTree::createHead(const std::string &name,
                              std::shared_ptr<Strategy> strategy,
                              unsigned int dimension, bool useVec) {
  if (_head) // if head exists throw exception
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  createNode(name, strategy, "", dimension, useVec);
}

void FunctionTree::insertNode(std::shared_ptr<TreeNode> node,
                              std::string parent) {
  std::string n = node->name();
  bool ret = _nodes
                 .insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
                     node->name(), node))
                 .second;
  if (!ret)
    throw TreeBuildError("FunctionTree::insertNode | Can not insert node " +
                         node->name() + "!");

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
  node->addParent(parentNode);
  //  inNode->linkParents();

  // Subtree already linked, but need to be added to list of nodes
  AddChildNodes(node);
}

void FunctionTree::insertTree(std::shared_ptr<FunctionTree> tree,
                              std::string parent) {
  _childTrees.push_back(tree);
  insertNode(tree->head(), parent);
  return;
}

void FunctionTree::addNode(std::shared_ptr<ComPWA::TreeNode> newNode) {
  // TODO: check existence, throw exception
  _nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
      newNode->name(), newNode));
  newNode->linkParents();
}

void FunctionTree::createNode(const std::string &name,
                              std::shared_ptr<Strategy> strategy,
                              std::string parent, unsigned int dimension,
                              bool useVec) {
  if (dimension == 0)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Zero dimension! ");
  if (parent == "" && _head)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::vector<std::shared_ptr<Parameter>> inter;

  switch (strategy->OutType()) {
  case ParType::MCOMPLEX: {
    inter.push_back(std::shared_ptr<Parameter>(new MultiComplex(
        "par" + name, std::vector<std::complex<double>>(
                          dimension, std::complex<double>(0., 0.)))));
    break;
  }
  case ParType::MDOUBLE: {
    inter.push_back(std::shared_ptr<Parameter>(
        new MultiDouble("par" + name, std::vector<double>(dimension, 0.))));
    break;
  }
  case ParType::COMPLEX: {
    std::complex<double> start(0., 0.);
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<Parameter>(
            new ComplexParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<Parameter>(
          new ComplexParameter("par" + name, start)));
    }
    break;
  }
  case ParType::DOUBLE: {
    double start(0.);
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<Parameter>(
            new DoubleParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<Parameter>(
          new DoubleParameter("par" + name, start)));
    }
    break;
  }
  case ParType::INTEGER: {
    int start(0);
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<Parameter>(
            new IntegerParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<Parameter>(
          new IntegerParameter("par" + name, start)));
    }
    break;
  }
  case ParType::BOOL: {
    bool start = false;
    if (useVec) {
      for (unsigned int i = 0; i < dimension; i++)
        inter.push_back(std::shared_ptr<Parameter>(
            new BoolParameter("par" + name, start)));
    } else {
      inter.push_back(std::shared_ptr<Parameter>(
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
  newNode->linkParents();
  if (parent == "")
    _head = newNode; // if we created a head redirect pointer
}

void FunctionTree::createLeaf(const std::string name,
                              std::shared_ptr<Parameter> parameter,
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
    leaf->addParent(parentNode);
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, parameter, std::shared_ptr<Strategy>(), parentNode));
    _nodes.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    parameter->Attach(leaf);
    if (parent == "")
      _head = leaf; // if we created a head, redirect pointer
  }
}

void FunctionTree::createLeaf(const std::string name, const double value,
                              std::string parent) {

  std::shared_ptr<DoubleParameter> staticVal(
      new DoubleParameter("ParOfNode_" + name, value));
  createLeaf(name, staticVal, parent);
}

void FunctionTree::createLeaf(const std::string name,
                              const std::complex<double> value,
                              std::string parent) {

  std::shared_ptr<ComplexParameter> staticVal(
      new ComplexParameter("ParOfNode_" + name, value));
  createLeaf(name, staticVal, parent);
}

void FunctionTree::createLeaf(
    const std::string name,
    std::vector<std::shared_ptr<Parameter>> &parameters,
    std::string parent) {

  for (auto i : parameters) {
    createLeaf(name, i, parent);
  }
}

bool FunctionTree::sanityCheck() {
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

  return true;
}

void FunctionTree::fillParameters(ParameterList &list) {
  _head->fillParameters(list);
}

void FunctionTree::GetNamesDownward(std::shared_ptr<TreeNode> start,
                                    std::vector<std::string> &childNames,
                                    std::vector<std::string> &parentNames) {

  start->fillParentNames(parentNames);
  start->fillChildNames(childNames);

  const std::vector<std::shared_ptr<TreeNode>> childs = start->childNodes();
  for (unsigned int i = 0; i < childs.size(); i++)
    GetNamesDownward(childs.at(i), childNames, parentNames);
  return;
}

void FunctionTree::AddChildNodes(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren =
      startNode->childNodes();

  for (unsigned int i = 0; i < newChildren.size(); i++) {
    _nodes.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
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
  std::vector<std::shared_ptr<TreeNode>> newChildren =
      startNode->childNodes();
  for (unsigned int i = 0; i < newChildren.size(); i++) {
    newChildren.at(i)->update();
    UpdateAll(newChildren[i]);
  }
  return;
}
