//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <string>
#include <complex>
#include <memory>

#include "Core/TreeNode.hpp"
#include "Core/Functions.hpp"

TreeNode::TreeNode(std::string name, std::shared_ptr<AbsParameter> intResult, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent)
:_name(name),_changed(true),_strat(strat){
	_value.push_back(intResult);
	if(parent){
		_parents.push_back(parent);
		//parent->children.push_back(shared_from_this());
	}
};

TreeNode::TreeNode(std::string name, std::vector<std::shared_ptr<AbsParameter>>& intResult, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent)
:_name(name),_changed(true),_strat(strat){
	for(unsigned int i=0; i<intResult.size(); i++){
		_value.push_back(intResult[i]);
	}
	if(parent){
		_parents.push_back(parent);
		//parent->children.push_back(shared_from_this());
	}
};

//! Destructor
TreeNode::~TreeNode(){
	//;
}

void TreeNode::Update(){ //sollte nur von kindern oder observed objects aufgerufen werden!
	//std::vector<double> newVals;
	//for(unsigned int i=0; i<children.size(); i++){
	//    newVals.push_back(children[i]->value);
	//}  //end children-loop
	//changeVal(myStrat->execute(newVals));
	for(unsigned int i=0; i<_parents.size(); i++)
		_parents[i]->Update();
	_changed=true;
}; //end update()

void TreeNode::recalculate(){
	if(_children.size()<1){
		_changed=false;
		return;
	}

	if(_value.size()==1){ //i have just one dim, merge everything
		ParameterList newVals;
		for(unsigned int i=0; i<_children.size(); i++){ //all children
			for(unsigned int ele=0; ele<_children[i]->getDim(); ele++){ //all dims
				if(_children[i]->needsCalculation())
					_children[i]->recalculate();
				std::shared_ptr<AbsParameter> para = _children[i]->getValue(ele);
				if(!para) std::cout << this->getName() << " child failed: " << i << std::endl;
				//para->type();
				newVals.AddParameter(para);
			}  //end children-loop
		}
		_changed=false;
		//std::cout << "Values: " << newVals.GetNDouble() << " " << _strat << std::endl;
		if(!_strat->execute(newVals, _value[0])) std::cout << this->getName() << " strat failed: " << _strat << std::endl;
		//std::cout << "New Val: " << _value[0] << std::endl;

	}else{ //i have a certain dim, children must fill it

		for(unsigned int ele=0; ele<_value.size(); ele++){
			ParameterList newVals;
			for(unsigned int i=0; i<_children.size(); i++){
				if(!(_children[i]->getDim() == _value.size() || _children[i]->getDim() == 1))
					continue; //TODO: exception;
				std::shared_ptr<AbsParameter> para;

				if(_children[i]->needsCalculation())
					_children[i]->recalculate();

				if(_children[i]->getDim()==1){
					para = _children[i]->getValue();
				}else{
					para = _children[i]->getValue(ele);
				}

				if(!para) std::cout << this->getName() << " " << i << std::endl;
				//para->type();
				newVals.AddParameter(para);
			}  //end children-loop
			_strat->execute(newVals, _value[ele]);
			_changed=false;
		}//end loop dims

	}//end which dim of this node
}; //end update()

std::string TreeNode::to_str(int lv, std::string beginning){
	std::stringstream oss;
	if(_changed && _children.size())
		oss << beginning << _name << " = ?";
	else{
		oss << beginning << _name << " = ";
		for(unsigned int i=0; i<_value.size()-1; i++)
			if(_value[i])
				oss << _value[i]->val_to_str() << ", ";
		if(_value[_value.size()-1])
			oss << _value[_value.size()-1]->val_to_str();
	}
	if(_children.size())
		oss << " (" << _children.size() << " children/"<<_value.size()<<" values)" << std::endl;
	else
		oss << std::endl;

	if(lv==0) return oss.str();
	for(unsigned int i=0; i<_children.size(); i++){
		//oss << " -> ";
		oss << _children.at(i)->to_str(lv-1,beginning+". ");
	}
	return oss.str();
};

std::string TreeNode::print(unsigned int lv){
	return to_str(lv);
}

const void TreeNode::deleteLinks(){
	_children.clear();
	_parents.clear();
	for(unsigned int i=0; i<_value.size(); i++){
		_value[i]->Detach(shared_from_this());
	}
};

std::shared_ptr<TreeNode> TreeNode::getChildNode(std::string name) const{
	std::shared_ptr<TreeNode> node;
	if(!_children.size()) node = std::shared_ptr<TreeNode>();
	for(unsigned int i=0; i<_children.size(); i++){
		if(_children.at(i)->getName()==name){
			return _children.at(i);
		} else
			node = _children.at(i)->getChildNode(name);
		if(node) return node;
	}
	return node;
}

std::shared_ptr<AbsParameter> TreeNode::getChildValue(std::string name) const{
	std::shared_ptr<TreeNode> node = std::shared_ptr<TreeNode>();
	node = getChildNode(name);
	if(node) return node->getValue();
	return std::shared_ptr<AbsParameter>();

	//	std::shared_ptr<AbsParameter> ret = std::shared_ptr<AbsParameter>();
	//	for(unsigned int i=0; i<_children.size(); i++){
	//		if(_children[i]->getName()==name)
	//			return _children[i]->getValue(0);
	//		else {
	//			ret = _children[i]->getChildValue(name);
	//			if(ret) return ret;
	//			else continue;
	//		}
	//	}
	//	return ret;
}

std::complex<double> TreeNode::getChildSingleValue(std::string name) const{
	std::shared_ptr<TreeNode> node = std::shared_ptr<TreeNode>();
	node = getChildNode(name);
	if(node) {
		std::shared_ptr<AbsParameter> val = node->getValue();
		if(val->type() == ParType::DOUBLE)
			return std::complex<double>(
					(std::dynamic_pointer_cast<DoubleParameter>(val))->GetValue(),0);
		if(val->type() == ParType::COMPLEX)
			return std::complex<double>(
					(std::dynamic_pointer_cast<ComplexParameter>(val))->GetValue());
		if(val->type() == ParType::INTEGER)
			return std::complex<double>(
					(std::dynamic_pointer_cast<IntegerParameter>(val))->GetValue(),0);
		if(val->type() == ParType::BOOL)
			return std::complex<double>(
					(std::dynamic_pointer_cast<BoolParameter>(val))->GetValue(),0);
		if(val->type() == ParType::MDOUBLE)
			return std::complex<double>(
					(std::dynamic_pointer_cast<MultiDouble>(val))->GetValue(0),0);
		if(val->type() == ParType::MCOMPLEX)
			return std::complex<double>(
					(std::dynamic_pointer_cast<MultiComplex>(val))->GetValue(0));
	}
	return std::complex<double>(-999,0);
}

std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p){
	return os << p->to_str(-1);
}
