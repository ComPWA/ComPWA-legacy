// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <string>
#include <complex>
#include <memory>

#include "Core/TreeNode.hpp"
#include "Core/Functions.hpp"

using namespace ComPWA;

TreeNode::TreeNode(std::string name, std::shared_ptr<AbsParameter> parameter,
                   std::shared_ptr<Strategy> strategy,
                   std::shared_ptr<TreeNode> parent)
    : _name(name), _changed(true), _strat(strategy) {
  _parameters.push_back(parameter);
  if (parent) {
    _parents.push_back(parent);
  }
}

TreeNode::TreeNode(std::string name,
                   std::vector<std::shared_ptr<AbsParameter>> &parameters,
                   std::shared_ptr<Strategy> strategy,
                   std::shared_ptr<TreeNode> parent)
    : _name(name), _changed(true), _strat(strategy) {
  for (unsigned int i = 0; i < parameters.size(); i++) {
    _parameters.push_back(parameters.at(i));
  }
  if (parent) {
    _parents.push_back(parent);
  }
}

TreeNode::~TreeNode() {}

void TreeNode::LinkParents() {
  for (unsigned int i = 0; i < _parents.size(); i++)
    _parents.at(i)->_children.push_back(shared_from_this());
}

void TreeNode::Update() {
  for (unsigned int i = 0; i < _parents.size(); i++)
    _parents.at(i)->Update();
  _changed = true;
};

void TreeNode::Recalculate() {
  if (_changed == false)
    return;

  if (_children.size() < 1) {
    _changed = false;
    return;
  }

  if (_parameters.size() == 1) { // Single dimension
    ParameterList newVals;
    for (auto ch : _children) {
      ch->Recalculate();
      for (auto p : ch->Parameters()) {
        newVals.AddParameter(p);
      }
    }

    try {
      _strat->execute(newVals, _parameters.at(0));
    } catch (std::exception &ex) {
      LOG(error) << "TreeNode::Recalculate() | Strategy " << _strat
                 << " failed on node " << Name() << ": " << ex.what();
      throw;
    }
  } else { // Multi dimensions
    for (unsigned int ele = 0; ele < _parameters.size(); ele++) {
      ParameterList newVals;

      for (auto ch : _children) {
        ch->Recalculate();
        if (ch->Dimension() == 1)
          newVals.AddParameter(ch->Parameter(0));
        else if (ch->Dimension() != _parameters.size())
          newVals.AddParameter(ch->Parameter(ele));
        else
          throw std::runtime_error("TreeNode::Recalculate() | Dimension of "
                                   "child nodes does not match");
      }

      try {
        _strat->execute(newVals, _parameters.at(ele));
      } catch (std::exception &ex) {
        LOG(error) << "TreeNode::recalculate() | Strategy " << _strat
                   << " failed on node " << Name() << ": " << ex.what();
        throw;
      }
    }
  }
  _changed = false;
}

std::shared_ptr<AbsParameter> TreeNode::Parameter(unsigned int position) {
  return _parameters.at(position);
}

std::vector<std::shared_ptr<AbsParameter>> &TreeNode::Parameters() {
  return _parameters;
}

void TreeNode::FillParameters(ComPWA::ParameterList &list) {
  for (auto ch : _children) {
    ch->FillParameters(list);
  }
  for (auto i : _parameters) {
    if (i->type() == ComPWA::ParType::DOUBLE)
      list.AddParameter(i);
  }
}

std::shared_ptr<TreeNode> TreeNode::FindChildNode(std::string name) const {
  std::shared_ptr<TreeNode> node;
  if (!_children.size())
    node = std::shared_ptr<TreeNode>();
  for (unsigned int i = 0; i < _children.size(); i++) {
    if (_children.at(i)->Name() == name) {
      return _children.at(i);
    } else
      node = _children.at(i)->FindChildNode(name);
    if (node)
      return node;
  }
  return node;
}

// std::shared_ptr<AbsParameter> TreeNode::ChildValue(std::string name) const {
//  std::shared_ptr<TreeNode> node = FindChildNode(name);
//  if (node)
//    return node->Parameter();
//
//  return std::shared_ptr<AbsParameter>();
//}
//
// std::complex<double> TreeNode::getChildSingleValue(std::string name) const {
//  std::shared_ptr<TreeNode> node = std::shared_ptr<TreeNode>();
//  node = FindChildNode(name);
//  if (node) {
//    std::shared_ptr<AbsParameter> val = node->Parameter();
//    if (val->type() == ParType::DOUBLE)
//      return std::complex<double>(
//          (std::dynamic_pointer_cast<DoubleParameter>(val))->GetValue(), 0);
//    if (val->type() == ParType::COMPLEX)
//      return std::complex<double>(
//          (std::dynamic_pointer_cast<ComplexParameter>(val))->GetValue());
//    if (val->type() == ParType::INTEGER)
//      return std::complex<double>(
//          (std::dynamic_pointer_cast<IntegerParameter>(val))->GetValue(), 0);
//    if (val->type() == ParType::BOOL)
//      return std::complex<double>(
//          (std::dynamic_pointer_cast<BoolParameter>(val))->GetValue(), 0);
//    if (val->type() == ParType::MDOUBLE)
//      return std::complex<double>(
//          (std::dynamic_pointer_cast<MultiDouble>(val))->GetValue(0), 0);
//    if (val->type() == ParType::MCOMPLEX)
//      return std::complex<double>(
//          (std::dynamic_pointer_cast<MultiComplex>(val))->GetValue(0));
//  }
//  return std::complex<double>(-999, 0);
//}

std::string TreeNode::Print(int level, std::string prefix) const {
  std::stringstream oss;
  if (_changed && _children.size()) {
    oss << prefix << _name << " = ?";
  } else {
    oss << prefix << _name;
    auto it = _parameters.begin();
    for (; it != _parameters.end(); ++it) {
      if (!_children.size()) // print parameter name for leafs
        oss << " [" << (*it)->GetName() << "]";
      oss << " = " << (*it)->val_to_str();
      if (it != _parameters.end())
        oss << ", ";
    }
  }

  if (_children.size())
    oss << " (" << _children.size() << " children/" << _parameters.size()
        << " values)" << std::endl;
  else
    oss << std::endl;

  if (level == 0)
    return oss.str();
  for (unsigned int i = 0; i < _children.size(); i++) {
    oss << _children.at(i)->Print(level - 1, prefix + ". ");
  }
  return oss.str();
}

void TreeNode::AddChild(std::shared_ptr<TreeNode> childNode) {
  _children.push_back(childNode);
}

void TreeNode::AddParent(std::shared_ptr<TreeNode> parentNode) {
  _parents.push_back(parentNode);
  parentNode->_children.push_back(shared_from_this());
}

void TreeNode::FillParentNames(std::vector<std::string> &names) const {
  for (auto i : _parents) {
    names.push_back(i->Name());
  }
}

void TreeNode::FillChildNames(std::vector<std::string> &names) const {
  for (unsigned int i = 0; i < _children.size(); i++)
    names.push_back(_children.at(i)->Name());
}

std::vector<std::shared_ptr<TreeNode>> &TreeNode::GetChildNodes() {
  return _children;
}

void TreeNode::DeleteLinks() {
  _children.clear();
  _parents.clear();
  for (unsigned int i = 0; i < _parameters.size(); i++) {
    _parameters.at(i)->Detach(shared_from_this());
  }
}
