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

#include "Core/Functions.hpp"
#include "Core/TreeNode.hpp"
#include "Core/AbsParameter.hpp"
#include "Core/Parameter.hpp"

#include "Optimizer/ControlParameter.hpp"

using namespace boost::log;

class FunctionTree //: public ControlParametr
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
	  for(std::map<std::string, std::shared_ptr<TreeNode> >::iterator iter = nodes_.begin(); iter != nodes_.end(); ++iter){
		  iter->second->deleteLinks();
	  }

    //;
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
    * \sa addNode(), createNode(), createLeaf()
   */
  virtual void createHead(const std::string& name, std::shared_ptr<Strategy> strat){
    //TODO: (type of) parameter from strategy!
    //std::shared_ptr<AbsParameter> inter = strat->GetResultContainer();

    std::shared_ptr<DoubleParameter> inter(new DoubleParameter("par"+name,0.));
    std::shared_ptr<TreeNode> newNode(new TreeNode(name, inter, strat, std::shared_ptr<TreeNode>()));
    head_ = newNode;
    nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(name,newNode));
  }

  //! Create a node for the FcnTree
   /*!
    * Create and add a node to the function tree
    * Adds Top-Down-Linking to the node
    * \param name identifier of node
    * \param strat Strategy with which the node calculates its value
    * \param parent the parent of this node (for linking)
    * \sa addNode(), createHead(), createLeaf()
   */
  virtual void createNode(const std::string& name, std::shared_ptr<Strategy> strat, std::string parent, unsigned int dim=1, bool useVec=false){
    //TODO: (type of) parameter from strategy!
    //std::shared_ptr<AbsParameter> inter = strat->GetResultContainer();
	if(dim==0) return; //TODO: exception

	std::vector<std::shared_ptr<AbsParameter>> inter;
	switch(strat->OutType()){
      case ParType::MCOMPLEX:{
        //TODO: error if dim>1
        std::vector<std::complex<double> > start(dim, std::complex<double>(0.,0.));
        inter.push_back(std::shared_ptr<AbsParameter>(new MultiComplex("par"+name,start)));
        break;
      }//end multi complex

      case ParType::MDOUBLE:{
        //TODO: error if dim>1
        std::vector<double> start(dim, 0.);
        inter.push_back(std::shared_ptr<AbsParameter>(new MultiDouble("par"+name,start)));
        break;
      }//end multi double

      case ParType::COMPLEX:{
        std::complex<double> start(0.,0.);
        if(useVec){
          for(unsigned int i; i<dim; i++)
            inter.push_back(std::shared_ptr<AbsParameter>(new ComplexParameter("par"+name,start)));
        }else{
          inter.push_back(std::shared_ptr<AbsParameter>(new ComplexParameter("par"+name,start)));
        }
        break;
      }//end complex

      case ParType::DOUBLE:{
        double start(0.);
        if(useVec){
          for(unsigned int i; i<dim; i++)
            inter.push_back(std::shared_ptr<AbsParameter>(new DoubleParameter("par"+name,start)));
        }else{
          inter.push_back(std::shared_ptr<AbsParameter>(new DoubleParameter("par"+name,start)));
        }
        break;
      }//end double

      case ParType::INTEGER:{
        int start(0);
        if(useVec){
          for(unsigned int i; i<dim; i++)
            inter.push_back(std::shared_ptr<AbsParameter>(new IntegerParameter("par"+name,start)));
        }else{
          inter.push_back(std::shared_ptr<AbsParameter>(new IntegerParameter("par"+name,start)));
        }
        break;
      }//end int

      case ParType::BOOL:{
        bool start=false;
        if(useVec){
          for(unsigned int i; i<dim; i++)
            inter.push_back(std::shared_ptr<AbsParameter>(new BoolParameter("par"+name,start)));
        }else{
          inter.push_back(std::shared_ptr<AbsParameter>(new BoolParameter("par"+name,start)));
        }
        break;
      }//end bool

      default:{
        //TODO: exception output partype wrong
        return;
      }
	}//end switch

    std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
	if(dim==1 && !useVec){
      std::shared_ptr<TreeNode> newNode(new TreeNode(name, inter[0], strat, parentNode));
      nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(name,newNode));
      newNode->linkParents();
	}else{
      std::shared_ptr<TreeNode> newNode(new TreeNode(name, inter, strat, parentNode));
      nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(name,newNode));
      newNode->linkParents();
	}

  }

  //! Create a leaf for the FcnTree
   /*!
    * Create and add a static node to the function tree if not existing yet
    * Adds Top-Down-Linking to the node
    * \param name identifier of node
    * \param extPar the parameter this node represents
    * \param parent the parent of this node (for linking)
    * \sa addNode(), createHead(), createNode()
   */
  virtual void createLeaf(const std::string name, const double extPar, std::string parent){
    std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
    std::shared_ptr<TreeNode> leaf;
    std::shared_ptr<DoubleParameter> staticVal(new DoubleParameter("ParOfNode_"+name,extPar));

    //check if Leaf already exists
    bool exists = true;
    if ( nodes_.find( name ) == nodes_.end() ) {
      exists=false;
    } else {
      leaf = nodes_.at(name);
    }

    //setup connections
    if(exists){
      leaf->addParent(parentNode);
      //TODO: check if also static?
    }else{
      leaf = std::shared_ptr<TreeNode>(new TreeNode(name, staticVal, std::shared_ptr<Strategy>(), parentNode));
      nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(name,leaf));
      leaf->linkParents();
    }
  }

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
  virtual void createLeaf(const std::string name, std::shared_ptr<AbsParameter> extPar, std::string parent){
    std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
    std::shared_ptr<TreeNode> leaf;

    //check if Leaf already exists
    bool exists = true;
    if ( nodes_.find( name ) == nodes_.end() ) {
      exists=false;
    } else {
      leaf = nodes_.at(name);
    }

    //setup connections
    if(exists){
      leaf->addParent(parentNode);
    }else{
      leaf = std::shared_ptr<TreeNode>(new TreeNode(name, extPar, std::shared_ptr<Strategy>() , parentNode));
      nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(name,leaf));
      leaf->linkParents();
      extPar->Attach(leaf);
    }
  }

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
  virtual void createLeaf(const std::string name, std::vector<std::shared_ptr<AbsParameter>>& extPar, std::string parent){
    std::shared_ptr<TreeNode> parentNode = nodes_.at(parent);
    std::shared_ptr<TreeNode> leaf;

    //check if Leaf already exists
    bool exists = true;
    if ( nodes_.find( name ) == nodes_.end() ) {
      exists=false;
    } else {
      leaf = nodes_.at(name);
    }

    //setup connections
    if(exists){
      leaf->addParent(parentNode);
    }else{
      leaf = std::shared_ptr<TreeNode>(new TreeNode(name, extPar, std::shared_ptr<Strategy>() , parentNode));
      nodes_.insert(std::pair<std::string, std::shared_ptr<TreeNode> >(name,leaf));
      leaf->linkParents();
      for(unsigned int i=0; i<extPar.size(); i++)
        extPar[i]->Attach(leaf);
    }
  }

  //! return the head of the tree
   /*!
    * Access to the head element of the tree
    * \return FuntionTreeNode at head of tree
   */
  virtual const std::shared_ptr<TreeNode> head() const {
    return head_;
    //TODO: return double? ;
  }

 /* virtual double controlParameter(ParameterList& minPar){
    ParType resultType = head_->getValue()->type();
    switch(par->type())
      {
      case ParType::COMPLEX: //TODO
        //cout << "Easy\n";
        break;
      case ParType::DOUBLE:{
        DoubleParameter* tmp = (DoubleParameter*) par.get(); //TODO: saver
        return tmp->GetValue();
        break;}
      case ParType::INTEGER:{
        IntegerParameter* tmp = (IntegerParameter*) par.get();
        return tmp->GetValue();
        break;}
      case ParType::BOOL:{
        BoolParameter* tmp = (BoolParameter*) par.get();
        return tmp->GetValue();
        break;}
      default:{
        //TODO exception
        //cout << "Invalid Selection\n";
        break;}
      }
    return 0;
  }*/

  //! trigger calculation
  void recalculate(){
    head_->recalculate();
  }

  //! check if Tree functions, create some debug messages if not
  bool sanityCheck(){
    bool isSane=true;

    //first fo all: is there a tree?
    if(!head_){
      BOOST_LOG_TRIVIAL(error)<<"This tree has no head!";
      return false;
    }

    //some basic infos next
    BOOST_LOG_TRIVIAL(debug)<<"# of Nodes: "<<nodes_.size();

    //collect all children and parent names
    std::vector<std::string> childNames, parentNames;
    getNamesDownward(head_, childNames, parentNames);

    //check if matches with available nodeNames
    std::vector<std::string> missedChild, missedParent;
    for(unsigned int i=0; i<childNames.size();i++){
      try{
        nodes_.at(childNames[i]);
      }catch(std::out_of_range()){
        missedChild.push_back(childNames[i]);
      }
    }
    for(unsigned int i=0; i<parentNames.size();i++){
      try{
        nodes_.at(parentNames[i]);
      }catch(std::out_of_range()){
        missedParent.push_back(parentNames[i]);
      }
    }
    if(missedChild.size()){
      for(unsigned int i=0; i<missedChild.size();i++)
        BOOST_LOG_TRIVIAL(debug)<<"This tree misses a child: "<<missedChild[i];
      return false;
    }
    if(missedParent.size()){
      for(unsigned int i=0; i<missedParent.size();i++)
        BOOST_LOG_TRIVIAL(debug)<<"This tree misses a parent: "<<missedParent[i];
      return false;
    }


    return isSane;
  }

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


protected:
  std::shared_ptr<TreeNode> head_; /*!< the head node storing the absolute result */
  std::map<std::string, std::shared_ptr<TreeNode> > nodes_; /*!< map to store the nodes */

  //! recursive function to get all used NodeNames
  void getNamesDownward(std::shared_ptr<TreeNode> start, std::vector<std::string>& childNames, std::vector<std::string>& parentNames){

    start->getParentNames(parentNames);
    start->getChildrenNames(childNames);

    const std::vector<std::shared_ptr<TreeNode> > childs = start->getChildren();
    for(unsigned int i=0; i<childs.size(); i++)
      getNamesDownward(childs[i], childNames, parentNames);

  }

};

#endif /* _FUNCTIONTREE_HPP_ */
