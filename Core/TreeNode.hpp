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
//! TreeNode is the interface for elements of the FunctionTree
/*! \class TreeNode
 * @file TreeNode.hpp
 * This class acts as a container for a parameter in a function tree. It has a
 * Strategy to calculate its value and a unique name.
 */

#ifndef _TREENODE_HPP_
#define _TREENODE_HPP_

#include <string>
#include <complex>
#include <memory>

#include "Core/Functions.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParObserver.hpp"

class TreeNode : public std::enable_shared_from_this<TreeNode>, public ParObserver {
public:
	//! Standard constructor
	/*!
	 * Standard constructor for one dimensional tree using a string the identifying name
	 * /param name unique name of this node
	 * /param intResult start value of node
	 * /param strat strategy how this node is calculated
	 * /param parent pointer to connected upper level node
	 */
	TreeNode(std::string name, std::shared_ptr<AbsParameter> intResult,
			std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent);

	//! Standard constructor
	/*!
	 * Standard constructor for multi-dimensional tree using a string the identifying name
	 * /param name unique name of this node
	 * /param intResult start value of node
	 * /param strat strategy how this node is calculated
	 * /param parent pointer to connected upper level node
	 */
	TreeNode(std::string name,
			std::vector<std::shared_ptr<AbsParameter>>& intResult,
			std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent);

	//! Destructor
	~TreeNode();

	//! Add this node to parents children-list
	void linkParents(){
		for(unsigned int i=0; i<_parents.size(); i++)
			_parents[i]->_children.push_back(shared_from_this());
	}

	//! Check if recalculation is needed
	inline bool needsCalculation(){
		return _changed;
	};

	//! Trigger flagged to be recalculated
	void Update();

	//! Trigger recalculation
	void recalculate();

	//! Get pointer to node value
	const std::shared_ptr<AbsParameter> getValue(unsigned int ele=0){
		return _value[ele];
	};

	//! Get dimension
	const unsigned int getDim(){
		return _value.size();
	};

	//! Get node name
	const std::string& getName(){
		return _name;
	};

	//! Add link to children list
	void addChild(std::shared_ptr<TreeNode> newChild){
		_children.push_back(newChild);
	};

	//! Add link to parents list
	void addParent(std::shared_ptr<TreeNode> newParent){
		_parents.push_back(newParent);
		newParent->_children.push_back(shared_from_this());
	};

	//! return parents names
	void getParentNames(std::vector<std::string>& names){
		for(unsigned int i=0; i<_parents.size(); i++)
			names.push_back(_parents[i]->getName());
	};

	//! return children names
	void getChildrenNames(std::vector<std::string>& names){
		for(unsigned int i=0; i<_children.size(); i++)
			names.push_back(_children[i]->getName());
	};

	//! return children pointer
	const std::vector<std::shared_ptr<TreeNode> > getChildren(){
		return _children;
	};

	//! delete children & parent pointer
	const void deleteLinks();

	/*! String used to display tree
	 *\param lv print down to level lv, if lv=-1 print to whole tree, if lv=0 print current node only
	 */
	std::string to_str(int lv=-1, std::string beginning="");

	/** Return Node of tree
	 * We go recursively through the tree to find the specified node and return a shared_ptr of it
	 *
	 * @param name node specifier
	 * @return pointer to node
	 */
	virtual std::shared_ptr<TreeNode> getChildNode(std::string name) const;

	/** Return value of certain child node
	 * We go recursively through out tree to find the specified node and return its value. In case
	 * of a node with multiple values we return the first one. Currently we assume that the variable
	 * is a std::complex<double> or can be converted to it.
	 *
	 * @param name node specifier
	 * @return current value of node
	 */
	virtual std::complex<double> getChildSingleValue(std::string name) const;

	/** Return vector of values of certain child node
	 * We go recursively through out tree to find the specified node and the vector of its contents.
	 * In cast the nodes doesn't have multiple values with return a vector of the size 1.
	 *
	 * @param name node specifier
	 * @return current vector of values of node
	 */
	virtual std::shared_ptr<AbsParameter> getChildValue(std::string name) const;

	//! Stream-Operator used to display tree
	friend std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p);

	/*! Print structure of tree and its values
	 * \param lv Print down to level lv, lv=-1 print the whole tree
	 */
	std::string print(unsigned int lv=-1);

protected:
	std::vector<std::shared_ptr<TreeNode> > _parents; /*!< Link to parents */
	std::vector<std::shared_ptr<TreeNode> > _children; /*!< Link to children */

	std::vector<std::shared_ptr<AbsParameter>> _value; /*!< Value of this node */
	std::string _name; /*!< Unique name of this node */
	bool _changed; /*!< flag if node needs recalculation */

	std::shared_ptr<Strategy> _strat; /*!< Strategy how node calculates its value */
};

#endif /* _TREENODE_HPP_ */
