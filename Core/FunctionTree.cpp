/*
 * FunctionTree.cpp
 *
 *  Created on: Jul 20, 2015
 *      Author: weidenka
 */

#include "Core/FunctionTree.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
void FunctionTree::insertTree(std::shared_ptr<FunctionTree> inTree,
                              std::string parent) {
  _childTrees.push_back(inTree);
  insertNode(inTree->head(), parent);
  return;
}

void FunctionTree::insertNode(std::shared_ptr<TreeNode> inNode,
                              std::string parent) {
  std::string n = inNode->getName();
  bool ret = nodes_
                 .insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
                     inNode->getName(), inNode))
                 .second;
  if (!ret)
    throw TreeBuildError("FunctionTree::insertNode | "
                         "Can not insert node " +
                         inNode->getName() + "!");

  // Assign new parent to head of new tree, create links
  std::shared_ptr<TreeNode> parentNode;
  try {
    parentNode = nodes_.at(parent);
  } catch (std::out_of_range &ex) {
    LOG(error) << "FunctionTree::insertNode() | Parent node " << parent
               << " not found in FunctionTree!";
    throw;
  }

  /* In case of an existing node, it is possible that this node
   * already have parents. Do need to consider this here? */
  inNode->addParent(parentNode);
  //	inNode->linkParents();

  // Subtree already linked, but need to be added to list of nodes
  addChildNodes(inNode);
}

void FunctionTree::addChildNodes(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren = startNode->getChildren();

  for (unsigned int i = 0; i < newChildren.size(); i++) {
    nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode>>(
        newChildren.at(i)->getName(), newChildren.at(i)));
    //        if( !ret )
    //            throw TreeBuildError("FunctionTree::addChildNodes | "
    //                                 "Can not insert node
    //                                 "+inNode->getName()+"!");

    addChildNodes(newChildren.at(i));
  }
  return;
}

void FunctionTree::UpdateAll(std::shared_ptr<TreeNode> startNode) {
  std::vector<std::shared_ptr<TreeNode>> newChildren = startNode->getChildren();
  for (unsigned int i = 0; i < newChildren.size(); i++) {
    newChildren[i]->Update();
    UpdateAll(newChildren[i]);
  }
  return;
}

void FunctionTree::getNamesDownward(std::shared_ptr<TreeNode> start,
                                    std::vector<std::string> &childNames,
                                    std::vector<std::string> &parentNames) {

  start->getParentNames(parentNames);
  start->getChildrenNames(childNames);

  const std::vector<std::shared_ptr<TreeNode>> childs = start->getChildren();
  for (unsigned int i = 0; i < childs.size(); i++)
    getNamesDownward(childs[i], childNames, parentNames);
  return;
}

void FunctionTree::createHead(const std::string &name,
                              std::shared_ptr<Strategy> strat, unsigned int dim,
                              bool useVec) {
  if (head_) // if head exists throw exception
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  createNode(name, strat, "", dim, useVec);
}

void FunctionTree::createNode(const std::string &name,
                              std::shared_ptr<Strategy> strat,
                              std::string parent, unsigned int dim,
                              bool useVec) {
  if (dim == 0)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Zero dimension! ");
  if (parent == "" && head_)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::vector<std::shared_ptr<AbsParameter>> inter;

  switch (strat->OutType()) {
  case ParType::MCOMPLEX: {
    inter.push_back(std::shared_ptr<AbsParameter>(new MultiComplex(
        "par" + name,
        std::vector<std::complex<double>>(dim, std::complex<double>(0., 0.)))));
    break;
  }
  case ParType::MDOUBLE: {
    inter.push_back(std::shared_ptr<AbsParameter>(
        new MultiDouble("par" + name, std::vector<double>(dim, 0.))));
    break;
  }
  case ParType::COMPLEX: {
    std::complex<double> start(0., 0.);
    if (useVec) {
      for (unsigned int i = 0; i < dim; i++)
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
      for (unsigned int i = 0; i < dim; i++)
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
      for (unsigned int i = 0; i < dim; i++)
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
      for (unsigned int i = 0; i < dim; i++)
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
                       std::to_string(strat->OutType()));
  }
  } // end switch

  std::shared_ptr<TreeNode> newNode, parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = nodes_.at(parent);

  if (dim == 1 && !useVec)
    newNode = std::shared_ptr<TreeNode>(
        new TreeNode(name, inter[0], strat, parentNode));
  else
    newNode =
        std::shared_ptr<TreeNode>(new TreeNode(name, inter, strat, parentNode));

  nodes_.insert(
      std::pair<std::string, std::shared_ptr<TreeNode>>(name, newNode));
  newNode->linkParents();
  if (parent == "")
    head_ = newNode; // if we created a head redirect pointer
}

void FunctionTree::createLeaf(const std::string name, const double extPar,
                              std::string parent) {
  if (parent == "" && head_)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = nodes_.at(parent);
  //		std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
  std::shared_ptr<TreeNode> leaf;
  std::shared_ptr<DoubleParameter> staticVal(
      new DoubleParameter("ParOfNode_" + name, extPar));

  // check if Leaf already exists
  bool exists = true;
  if (nodes_.find(name) == nodes_.end()) {
    exists = false;
  } else {
    leaf = nodes_.at(name);
  }

  // setup connections
  if (exists) {
    leaf->addParent(parentNode);
    // TODO: check if also static?
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, staticVal, std::shared_ptr<Strategy>(), parentNode));
    nodes_.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    if (parent == "")
      head_ = leaf; // if we created a head redirect pointer
  }
}

void FunctionTree::createLeaf(const std::string name,
                              const std::complex<double> extPar,
                              std::string parent) {
  if (parent == "" && head_)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = nodes_.at(parent);
  //		std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
  std::shared_ptr<TreeNode> leaf;
  std::shared_ptr<ComplexParameter> staticVal(
      new ComplexParameter("ParOfNode_" + name, extPar));

  // check if Leaf already exists
  bool exists = true;
  if (nodes_.find(name) == nodes_.end()) {
    exists = false;
  } else {
    leaf = nodes_.at(name);
  }

  // setup connections
  if (exists) {
    leaf->addParent(parentNode);
    // TODO: check if also static?
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, staticVal, std::shared_ptr<Strategy>(), parentNode));
    nodes_.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    if (parent == "")
      head_ = leaf; // if we created a head redirect pointer
  }
}

void FunctionTree::createLeaf(const std::string name,
                              std::shared_ptr<AbsParameter> extPar,
                              std::string parent) {
  if (parent == "" && head_)
    throw std::runtime_error("FunctionTree::createNode() | "
                             "Head node already exists!");

  std::shared_ptr<TreeNode> parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = nodes_.at(parent);
  std::shared_ptr<TreeNode> leaf;

  // check if Leaf already exists
  bool exists = true;
  if (nodes_.find(name) == nodes_.end()) {
    exists = false;
  } else {
    leaf = nodes_.at(name);
  }

  // setup connections
  if (exists) {
    leaf->addParent(parentNode);
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, extPar, std::shared_ptr<Strategy>(), parentNode));
    nodes_.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    extPar->Attach(leaf);
    if (parent == "")
      head_ = leaf; // if we created a head, redirect pointer
  }
}

void FunctionTree::createLeaf(
    const std::string name, std::vector<std::shared_ptr<AbsParameter>> &extPar,
    std::string parent) {
  std::shared_ptr<TreeNode> parentNode;
  if (parent == "") // is this a head node?
    parentNode = std::shared_ptr<TreeNode>();
  else
    parentNode = nodes_.at(parent);
  //		std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
  std::shared_ptr<TreeNode> leaf;

  // check if Leaf already exists
  bool exists = true;
  if (nodes_.find(name) == nodes_.end()) {
    exists = false;
  } else {
    leaf = nodes_.at(name);
  }

  // setup connections
  if (exists) {
    leaf->addParent(parentNode);
  } else {
    leaf = std::shared_ptr<TreeNode>(
        new TreeNode(name, extPar, std::shared_ptr<Strategy>(), parentNode));
    nodes_.insert(
        std::pair<std::string, std::shared_ptr<TreeNode>>(name, leaf));
    leaf->linkParents();
    for (unsigned int i = 0; i < extPar.size(); i++)
      extPar[i]->Attach(leaf);
  }
}

bool FunctionTree::sanityCheck() {
  bool isSane = true;
  if (!head_)
    throw std::runtime_error("FunctionTree::sanityCheck() | "
                             "This tree has no head!");

  // collect all children and parent names
  std::vector<std::string> childNames, parentNames;
  getNamesDownward(head_, childNames, parentNames);

  // check if matches with available nodeNames
  std::vector<std::string> missedChild, missedParent;
  for (unsigned int i = 0; i < childNames.size(); i++) {
    try {
      nodes_.at(childNames[i]);
    } catch (std::out_of_range &ex) {
      missedChild.push_back(childNames[i]);
    }
  }
  for (unsigned int i = 0; i < parentNames.size(); i++) {
    try {
      nodes_.at(parentNames[i]);
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
} /* End ComPWA namespace */
