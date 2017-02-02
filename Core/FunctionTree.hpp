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
//! FunctionTree for the actual Optimization
/*! \class FunctionTree
 * @file FunctionTree.hpp
 * This class can be used to store a function in a tree-like structure. This
 * reduces the amount of recalculation needed when performing an optimization,
 * as most often only parts need to be recalculated in one iteration.
 */

#ifndef _FUNCTIONTREE_HPP_
#define _FUNCTIONTREE_HPP_

#include <vector>
#include <memory>
#include <string>
#include <map>

#include <boost/log/trivial.hpp>

#include "Core/Logging.hpp"
#include "Core/Functions.hpp"
#include "Core/TreeNode.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"

//#include "Optimizer/ControlParameter.hpp"

using namespace boost::log;

namespace ComPWA {

class FunctionTree    //: public ControlParametr
{
public:
	//! Standard constructor
	/*!
	 * Standard constructor for empty tree
	 */
	FunctionTree(){

	}

	//! Standard constructor
	/*!
	 * Standard constructor with the top node provided
	 * /param head first node to be used as head
	 */
	FunctionTree(std::shared_ptr<TreeNode> head):head_(head){
		nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(head->getName(),head));
	}

	//! Destructor
	virtual ~FunctionTree(){
		std::map<std::string, std::shared_ptr<TreeNode> >::iterator iter = nodes_.begin();
		for( ; iter != nodes_.end(); ++iter){
			iter->second->deleteLinks();
		}
	}
  
  virtual void recursiveAddDaughters(std::shared_ptr<TreeNode> newNode) {
    addNode(newNode);
    auto children = newNode->getChildren();
    for (auto child : children) {
      recursiveAddDaughters(child);
    }
  }

	//! Add node to FcnTree
	/*!
	 * Add a node to the function tree
	 * Adds Top-Down-Linking to the node
	 * \param newNode Node to be added
	 */
	virtual void addNode(std::shared_ptr<TreeNode> newNode){
		//TODO: check existence, throw exception
		nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(newNode->getName(),newNode));
		newNode->linkParents();
	}

	//! Create head node of FcnTree
	/*!
	 * Add top node to the function tree
	 * \param name identifier of node
	 * \param strat Strategy with which the node calculates its value
	 * \param dim dimension of the head node
	 * \param useVec true if value is a vector
	 * \sa addNode(), createNode(), createLeaf()
	 */
	virtual void createHead(const std::string& name, std::shared_ptr<Strategy> strat,
			unsigned int dim=1, bool useVec=false);

	//! Create head node of FcnTree
	/*!
	 * Create head node as a leaf. Constant parameter! We use this to generate an empty tree.
	 * \param name identifier of node
	 * \param extPar value of the head
	 * \sa addNode(), createNode(), createLeaf()
	 */
	virtual void createHead(const std::string& name, const double extPar){
		if( head_ )//if head exists throw exception
			throw std::runtime_error("FunctionTree::createNode() | head node already exists!");
		createLeaf(name,extPar,"");
	}

	//! Create head node of FcnTree
	/*!
	 * Create head node as a leaf.
	 * \param name identifier of node
	 * \param extPar parameter which holds the value of the head
	 * \sa addNode(), createNode(), createLeaf()
	 */
	virtual void createHead(const std::string& name, std::shared_ptr<AbsParameter> extPar){
		if( head_ )//if head exists throw exception
			throw std::runtime_error("FunctionTree::createNode() | head node already exists!");
		createLeaf(name,extPar,"");
	}

	//! Add an existing node to FunctionTree
	/*!
	 * Add an existing node to FunctionTree. Can be used to link tree's to each other: simply insert
	 * head node of tree A to tree B
	 *
	 * \param inNode Node to insert in FunctionTree
	 * \param parent the parent of this node (for linking)
	 * \sa addNode(), createHead(), createLeaf(), insertTree
	 */
	virtual void insertNode(std::shared_ptr<TreeNode> inNode, std::string parent);

	//! Insert an existing FunctionTree as TreeNode
	/*!
	 * Insert head of a existing FunctionTree to current tree.
	 *
	 * \param inTree Tree to insert in FunctionTree
	 * \param parent the parent of this node (for linking)
	 * \sa addNode(), createHead(), createLeaf(), insertNode()
	 */
	virtual void insertTree(std::shared_ptr<FunctionTree> inTree, std::string parent);

	//! Create a node for the FcnTree
	/*!
	 * Create and add a node to the function tree
	 * Adds Top-Down-Linking to the node
	 * \param name identifier of node
	 * \param strat Strategy with which the node calculates its value
	 * \param parent the parent of this node (for linking)
	 * \param dim dimension of the head node
     * \param useVec true if value is a vector
	 * \sa addNode(), createHead(), createLeaf()
	 */
	virtual void createNode(const std::string& name, std::shared_ptr<Strategy> strat,
			std::string parent, unsigned int dim=1, bool useVec=false);

	//! Create a leaf for the FcnTree
	/*!
	 * Create and add a static node to the function tree if not existing yet
	 * Adds Top-Down-Linking to the node
	 * \param name identifier of node
	 * \param extPar the parameter this node represents
	 * \param parent the parent of this node (for linking)
	 * \sa addNode(), createHead(), createNode()
	 */
	virtual void createLeaf(const std::string name, const double extPar, std::string parent);

	//! Create a leaf for the FcnTree
	/*!
	 * Create and add a node to the function tree if not existing yet
	 * Adds Top-Down-Linking to the node
	 * Attaches the Node as Observer to the external parameter
	 * \param name identifier of node
	 * \param extPar the parameter this node represents
	 * \param parent the parent of this node (for linking)
	 * \sa addNode(), createHead(), createNode()
	 */
	virtual void createLeaf(const std::string name, std::shared_ptr<AbsParameter> extPar,
			std::string parent);

	//! Create a leaf for the FcnTree with higher dimension
	/*!
	 * Create and add a node to the function tree if not existing yet
	 * Adds Top-Down-Linking to the node
	 * Attaches the Node as Observer to the external parameter
	 * \param name identifier of node
	 * \param extPar list of parameters this node represents
	 * \param parent the parent of this node (for linking)
	 * \sa addNode(), createHead(), createNode()
	 */
	virtual void createLeaf(const std::string name,
			std::vector<std::shared_ptr<AbsParameter> >& extPar, std::string parent);

	//! return the head of the tree
	/*!
	 * Access to the head element of the tree
	 * \return FuntionTreeNode at head of tree
	 */
	virtual const std::shared_ptr<TreeNode> head() const {
		return head_;
		//TODO: return double? ;
	}

	//! trigger calculation
	void recalculate(){
		head_->recalculate();
	}

	//! check if Tree functions, create some debug messages if not
	bool sanityCheck();

	//! friend function to stream parameter information to output
	/*!
	 * Declaring the stream-operator << as friend allows to stream parameter
	 * information to the output as easily as a generic type.
	 * \sa make_str(), to_str()
	 */
	friend std::ostream& operator<<( std::ostream& out, const FunctionTree& b ){
		return out << b.head();
	}

	//! friend function to stream parameter information to output
	/*!
	 * Declaring the stream-operator << as friend allows to stream parameter
	 * information to the output as easily as a generic type.
	 * \param out the ostream the object is piped to
	 * \param b the piped object
	 * \sa make_str(), to_str()
	 */
	friend std::ostream& operator<<( std::ostream& out, std::shared_ptr<FunctionTree> b ){
		return out << b->head();
	}
	/*! Print structure of tree and its values
	 * \param lv Print down to level lv, lv=-1 print the whole tree
	 */
	std::string print(unsigned int lv=-1) { return head_->print(lv); };

	//! Get number of nodes
	int GetNumberOfNodes() { return nodes_.size(); }

protected:
	/** List of child tree's
	 * We need to store the childTreee that were added via insertTree() here. Because otherwise the
	 * destructor of these tree's would delete the linking of the tree nodes.
	 * */
	std::vector<std::shared_ptr<FunctionTree> > _childTrees;

	std::shared_ptr<TreeNode> head_; /*!< the head node storing the absolute result */
	std::map<std::string, std::shared_ptr<TreeNode> > nodes_; /*!< map to store the nodes */
	//! recursive function to get all used NodeNames
	void getNamesDownward(std::shared_ptr<TreeNode> start, std::vector<std::string>& childNames,
			std::vector<std::string>& parentNames);
	//! Helper function to recursively add child nodes of a new tree
	virtual void addChildNodes(std::shared_ptr<TreeNode> startNode);
	//! Helper function to set all nodes to status changed
	virtual void UpdateAll(std::shared_ptr<TreeNode> startNode);
};

} /* End ComPWA namespace */
#endif /* _FUNCTIONTREE_HPP_ */
