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
  virtual void execute(ParameterList &paras, std::shared_ptr<Parameter> &out){};
};

FunctionTree::FunctionTree(std::string name,
                           std::shared_ptr<Parameter> parameter,
                           std::shared_ptr<Strategy> strategy) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  Head = std::shared_ptr<TreeNode>(
      new TreeNode(name, parameter, strategy, DummyNode));
}

FunctionTree::FunctionTree(std::string name,
                           std::shared_ptr<Parameter> parameter) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  Head = std::shared_ptr<TreeNode>(
      new TreeNode(name, parameter, std::shared_ptr<Strategy>(), DummyNode));
  parameter->Attach(Head);
}

FunctionTree::FunctionTree(std::string name, double value) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  Head =
      std::make_shared<TreeNode>(name, std::make_shared<Value<double>>(value),
                                 std::shared_ptr<Strategy>(), DummyNode);
}

FunctionTree::FunctionTree(std::string name, std::complex<double> value) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  //  createLeaf(name, value, "");
  Head = std::make_shared<TreeNode>(
      name, std::make_shared<Value<std::complex<double>>>(value),
      std::shared_ptr<Strategy>(), DummyNode);
}

FunctionTree::FunctionTree(std::shared_ptr<TreeNode> head) : Head(head) {
  DummyNode = std::make_shared<TreeNode>(
      "DummyNode", std::make_shared<Dummy>(ParType::UNDEFINED),
      std::shared_ptr<TreeNode>());
  head->Parents.push_back(DummyNode);
  Head = head;
}

FunctionTree::~FunctionTree() { Head->deleteParentLinks(DummyNode); }

void FunctionTree::insertNode(std::shared_ptr<TreeNode> node,
                              std::string parent) {

  auto parentNode = Head->findNode(parent);
  if (!parentNode)
    throw TreeBuildError("FunctionTree::insertNode | Cannot insert node " +
                         node->name() + ", already exists!");

  // Reuse existing node
  auto oldNode = Head->findNode(node->name());
  if (oldNode) {
    oldNode->Parents.push_back(parentNode);
    parentNode->update();
    return;
  }

  node->Parents.push_back(parentNode);
  parentNode->ChildNodes.push_back(node);
  parentNode->update();
}

void FunctionTree::insertTree(std::shared_ptr<FunctionTree> tree,
                              std::string parent) {
  insertNode(tree->Head, parent);
  return;
}

void FunctionTree::createNode(std::string name,
                              std::shared_ptr<Parameter> parameter,
                              std::shared_ptr<Strategy> strategy,
                              std::string parent) {
  auto leaf = std::make_shared<TreeNode>(name, parameter, strategy,
                                         std::shared_ptr<TreeNode>());
  if (parameter)
    parameter->Attach(leaf); // Can not be set in TreeNode constructor
  insertNode(leaf, parent);
}

void FunctionTree::createNode(std::string name,
                              std::shared_ptr<Strategy> strategy,
                              std::string parent) {
  createNode(name, std::shared_ptr<Parameter>(), strategy, parent);
}

void FunctionTree::createLeaf(std::string name,
                              std::shared_ptr<Parameter> parameter,
                              std::string parent) {
  createNode(name, parameter, std::shared_ptr<Strategy>(), parent);
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

} // namespace FunctionTree
} // namespace ComPWA
