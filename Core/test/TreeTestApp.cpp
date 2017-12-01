//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Test-Application to show usage of the FunctionTree.
/*!
 * @file TreeTestApp.cpp
 * This application test uses the ComPWA FunctionTree to calculate a simple
 * equation with cached intermediate results. The Tree can be set up in two
 * different way's, one using the functionality of the FunctionTree and the
 * other setting up the nodes and links manually. The second method is shown
 * mainly for a better understanding on how the FunctionTree internally works.
*/

#include <iostream>
#include <memory>
#include <string>
#include <regex>
#include <algorithm>
#include <vector>
#include <map>

#include "Core/Functions.hpp"
#include "Core/TreeNode.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Parameter.hpp"

using namespace ComPWA;

// This function will be called from a thread

int main(int argc, char **argv) {
  std::cout << "  ComPWA Copyright (C) 2013  Mathias Michel " << std::endl;
  std::cout << "  This program comes with ABSOLUTELY NO WARRANTY; for details "
               "see license.txt"
            << std::endl;
  std::cout << std::endl;

  bool autonodes = true;

  //------------SetUp some operations for R = a * ( b + c * d)-----------
  std::shared_ptr<Strategy> add =
      std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE));
  std::shared_ptr<Strategy> mult =
      std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE));
  std::shared_ptr<Strategy> madd =
      std::shared_ptr<Strategy>(new AddAll(ParType::MDOUBLE));
  std::shared_ptr<Strategy> mmult =
      std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE));

  //------------SetUp the parameters for R = a * ( b + c * d)-----------
  std::shared_ptr<DoubleParameter> parA(new DoubleParameter("parA", 5.));
  std::shared_ptr<DoubleParameter> parB(new DoubleParameter("parB", 2.));
  std::shared_ptr<DoubleParameter> parC(new DoubleParameter("parC", 3.));
  std::shared_ptr<DoubleParameter> parD(new DoubleParameter("parD", 1.));

  // std::cout << parA << std::endl;
  // std::cout << " ValToStr: " << parA->val_to_str() << std::endl;

  std::shared_ptr<FunctionTree> myTree, subTree;

  if (autonodes) { // Let FunctionTree manage the creation of the tree

    myTree = std::shared_ptr<FunctionTree>(new FunctionTree());
    myTree->CreateHead("R", mult);
    myTree->CreateLeaf("a", parA, "R");
    myTree->CreateNode("bcd", add, "R");
    myTree->CreateLeaf("b", parB, "bcd");

    subTree = std::shared_ptr<FunctionTree>(new FunctionTree());
    subTree->CreateHead("cd", mult);
    subTree->CreateLeaf("c", parC, "cd");
    subTree->CreateLeaf("d", parD, "cd");

    myTree->InsertTree(subTree, "bcd");

  } else { // create Tree manually

    // parameter container for intermediate results and final result
    std::shared_ptr<DoubleParameter> finalA(new DoubleParameter("finalA", 0.));
    std::shared_ptr<DoubleParameter> interBCD(
        new DoubleParameter("interBCD", 0.));
    std::shared_ptr<DoubleParameter> interCD(
        new DoubleParameter("interCD", 0.));

    //------------SetUp some nodes for R = a * ( b + c * d)----------------
    std::shared_ptr<TreeNode> R = std::shared_ptr<TreeNode>(
        new TreeNode("R", finalA, mult, std::shared_ptr<TreeNode>()));
    std::shared_ptr<TreeNode> A = std::shared_ptr<TreeNode>(
        new TreeNode("a", parA, std::shared_ptr<Strategy>(), R));
    parA->Attach(A);
    std::shared_ptr<TreeNode> BCD =
        std::shared_ptr<TreeNode>(new TreeNode("bcd", interBCD, add, R));
    std::shared_ptr<TreeNode> B = std::shared_ptr<TreeNode>(
        new TreeNode("b", parB, std::shared_ptr<Strategy>(), BCD));
    parB->Attach(B);
    std::shared_ptr<TreeNode> CD =
        std::shared_ptr<TreeNode>(new TreeNode("cd", interCD, add, BCD));
    std::shared_ptr<TreeNode> C = std::shared_ptr<TreeNode>(
        new TreeNode("c", parC, std::shared_ptr<Strategy>(), CD));
    parC->Attach(C);
    std::shared_ptr<TreeNode> D = std::shared_ptr<TreeNode>(
        new TreeNode("d", parD, std::shared_ptr<Strategy>(), CD));
    parD->Attach(D);
    myTree = std::shared_ptr<FunctionTree>(new FunctionTree(R));
    myTree->AddNode(A);
    myTree->AddNode(BCD);
    myTree->AddNode(B);
    myTree->AddNode(CD);
    myTree->AddNode(C);
    myTree->AddNode(D);
  }

  //------------Finished SetUp, now check Tree----------------

  std::cout << "R = a * ( b + c * d) Tree set up, not calculated" << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  //------------Trigger Calculation----------------
  myTree->Recalculate();

  std::cout << "Tree calculated" << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  parD->setValue(2.);

  std::cout << "Changed d from 1 to 2 " << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  myTree->Recalculate();

  std::cout << "Changed d from 1 to 2 and recalculated " << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  std::cout << std::endl << std::endl;

  //------------new Tree with more dimensions----------------

  unsigned int nElements = 10;
  std::shared_ptr<FunctionTree> myTreeMult;

  //------------SetUp the parameters for R = Sum of (a * b)-----------
  std::vector<std::shared_ptr<Parameter>> nVecParA;
  std::shared_ptr<Parameter> nParB(new DoubleParameter("parB", 2));
  for (unsigned int i = 0; i < nElements; i++) {
    std::shared_ptr<DoubleParameter> tmpA(
        new DoubleParameter("parA_" + std::to_string(i), i + 1));
    nVecParA.push_back(tmpA);
  }

  //   if(autonodes){ //Let FunctionTree manage the creation of the tree
  myTreeMult = std::shared_ptr<FunctionTree>(new FunctionTree());
  myTreeMult->CreateHead("R", add);
  myTreeMult->CreateNode("ab", mult, "R", nElements, true);
  myTreeMult->CreateLeaf("a", nVecParA, "ab");
  myTreeMult->CreateLeaf("b", nParB, "ab");

  // }

  //------------Trigger Calculation----------------
  myTreeMult->Recalculate();

  std::cout
      << std::endl
      << "Multitree Setup and calculated R = Sum[a*b] with one leaf containing "
      << nElements << " elements" << std::endl;
  std::cout << std::endl << myTreeMult << std::endl << std::endl;

  //------------new Tree with multiDouble Par----------------
  nElements = 5;
  std::shared_ptr<FunctionTree> myTreeMultD;

  //------------SetUp the parameters for R = Sum of (a * b)-----------
  std::vector<double> nMasses, nPhsp;
  std::shared_ptr<Parameter> mParB(new DoubleParameter("parB", 2));
  std::shared_ptr<Parameter> mParD(new DoubleParameter("parD", 3));
  for (unsigned int i = 0; i < nElements; i++) {
    // std::shared_ptr<DoubleParameter> tmpA(new
    // DoubleParameter("parA_"+i,i+1));
    nMasses.push_back(i + 1);
    nPhsp.push_back(2 * i + 1);
    nPhsp.push_back(2 * i + 2);
  }
  std::shared_ptr<MultiDouble> mParA(new MultiDouble("parA", nMasses));
  std::shared_ptr<MultiDouble> mParC(new MultiDouble("parC", nPhsp));

  myTreeMultD = std::shared_ptr<FunctionTree>(new FunctionTree());
  myTreeMultD->CreateHead("R", add);
  myTreeMultD->CreateNode("Rmass", mmult, "R", nElements);
  myTreeMultD->CreateNode("ab", mmult, "Rmass", nElements, false);
  myTreeMultD->CreateLeaf("a", mParA, "ab");
  myTreeMultD->CreateLeaf("b", mParB, "ab");
  myTreeMultD->CreateNode("Rphsp", add, "Rmass");
  myTreeMultD->CreateNode("cd", mmult, "Rphsp", nElements * 2, false);
  myTreeMultD->CreateLeaf("c", mParC, "cd");
  myTreeMultD->CreateLeaf("d", mParD, "cd");

  //------------Trigger Calculation----------------
  myTreeMultD->Recalculate();

  std::cout << std::endl
            << "MultiDouble Setup and calculated R = Sum[a*b] with one leaf "
               "containing "
            << nElements << " elements" << std::endl;
  std::cout << std::endl << myTreeMultD << std::endl << std::endl;

  return 0;
}

/*struct TreeNode{
  TreeNode(double inValue, std::string inName, std::shared_ptr<Strategy> strat,
std::shared_ptr<TreeNode> parent)
    :value(inValue),name(inName),changed(true),myStrat(strat){
    if(parent){
        parents.push_back(parent);
        //parent->children.push_back(shared_from_this());
    }
  };

  inline bool needsCalculation(){
    return changed;
  };

  inline void changeVal(double newVal){
    value=newVal;
    for(unsigned int i=0; i<parents.size(); i++)
      parents[i]->update();
    //changed=true;
  };

  void update(){ //darf nur von kindern aufgerufen werden!
    //std::vector<double> newVals;
    //for(unsigned int i=0; i<children.size(); i++){
    //    newVals.push_back(children[i]->value);
    //}  //end children-loop
    //changeVal(myStrat->execute(newVals));
    for(unsigned int i=0; i<parents.size(); i++)
      parents[i]->update();
    changed=true;
  }; //end update()

  void recalculate(){
    if(children.size()<1){
      changed=false;
      return;
    }
    std::vector<double> newVals;
    for(unsigned int i=0; i<children.size(); i++){
        if(children[i]->needsCalculation())
          children[i]->recalculate();
        newVals.push_back(children[i]->value);
    }  //end children-loop
    value = myStrat->execute(newVals);
    changed=false;
  }; //end update()

  std::string to_str(std::string beginning = ""){
    std::stringstream oss;
    if(changed && children.size())
      oss << beginning << name << " = ?";
    else
      oss << beginning << name << " = " << value;
    if(children.size())
      oss << " with " << children.size() << " children" << std::endl;
    else
      oss << std::endl;

    for(unsigned int i=0; i<children.size(); i++){
      //oss << " -> ";
      oss << beginning << children[i]->to_str(" -> ");
    }
    return oss.str();
  };

  friend std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode>
p);

  std::vector<std::shared_ptr<TreeNode> > parents;
  std::vector<std::shared_ptr<TreeNode> > children;

  double value;
  std::string name;
  bool changed;
  //std::string childOP;

  std::shared_ptr<Strategy> myStrat;
};

std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p){
  return os << p->to_str();
}*/
