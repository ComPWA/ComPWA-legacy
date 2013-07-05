//! Test-Application to check c++11 support.
/*!
 * @file CX0TestApp.cpp
 * This tiny application test some features of the new c++11 standard. One can
 * compile and run this to check the c++11 support of the used compiler/system.
 * The features tested are: shared_pointer, raw strings, auto type's, lamdas,
 * threading and regular expressions.
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

/*struct TreeNode{
  TreeNode(double inValue, std::string inName, std::shared_ptr<Strategy> strat, std::shared_ptr<TreeNode> parent)
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

  friend std::ostream & operator<<(std::ostream &os, std::shared_ptr<TreeNode> p);

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

/*class FcnTree{
public:
  FcnTree(){};

  void addHead(std::shared_ptr<TreeNode> inHead){
    nodeMap.add(inHead->name, inHead);
    head = inHead;
  };

  void addNode(std::shared_ptr<TreeNode> inNode){
    nodeMap.add(inNode->name, inNode);
  };

  std::shared_ptr<TreeNode> getNode(std::string nodeName){
    return nodeMap.get(nodeName);
  };

protected:
  std::shared_ptr<TreeNode> head;
  std::map<std::string, std::shared_ptr<TreeNode> > nodeMap;

};*/

//using namespace std;

//This function will be called from a thread

int main(int argc, char **argv) {
  bool autonodes=true;

  //------------SetUp some operations for R = a * ( b + c * d)-----------
  std::shared_ptr<Strategy> add = std::shared_ptr<Strategy>(new AddAll());
  std::shared_ptr<Strategy> mult = std::shared_ptr<Strategy>(new MultAll());

  //------------SetUp the parameters for R = a * ( b + c * d)-----------
  std::shared_ptr<DoubleParameter> parA(new DoubleParameter("parA",5.));
  std::shared_ptr<DoubleParameter> parB(new DoubleParameter("parB",2.));
  std::shared_ptr<DoubleParameter> parC(new DoubleParameter("parC",3.));
  std::shared_ptr<DoubleParameter> parD(new DoubleParameter("parD",1.));

  std::shared_ptr<FunctionTree> myTree;

  if(autonodes){ //Let FunctionTree manage the creation of the tree

    myTree = std::shared_ptr<FunctionTree>(new FunctionTree());
    myTree->createHead("R", mult);
    myTree->createLeaf("a", parA, "R");
    myTree->createNode("bcd", add, "R");
    myTree->createLeaf("b", parB, "bcd");
    myTree->createNode("cd", mult, "bcd");
    myTree->createLeaf("c", parC, "cd");
    myTree->createLeaf("d", parD, "cd");

  }else{ //create Tree manually

    //parameter container for intermediate results and final result
    std::shared_ptr<DoubleParameter> finalA(new DoubleParameter("finalA",0.));
    std::shared_ptr<DoubleParameter> interBCD(new DoubleParameter("interBCD",0.));
    std::shared_ptr<DoubleParameter> interCD(new DoubleParameter("interCD",0.));

    //------------SetUp some nodes for R = a * ( b + c * d)----------------
    std::shared_ptr<TreeNode> R =
        std::shared_ptr<TreeNode>(new TreeNode("R", finalA, mult, NULL));
    std::shared_ptr<TreeNode> A =
        std::shared_ptr<TreeNode>(new TreeNode("a", parA, NULL, R));
    parA->Attach(A);
    std::shared_ptr<TreeNode> BCD =
        std::shared_ptr<TreeNode>(new TreeNode("bcd", interBCD, add, R));
    std::shared_ptr<TreeNode> B =
        std::shared_ptr<TreeNode>(new TreeNode("b", parB, NULL, BCD));
    parB->Attach(B);
    std::shared_ptr<TreeNode> CD =
        std::shared_ptr<TreeNode>(new TreeNode("cd", interCD, add, BCD));
    std::shared_ptr<TreeNode> C =
        std::shared_ptr<TreeNode>(new TreeNode("c", parC, NULL, CD));
    parC->Attach(C);
    std::shared_ptr<TreeNode> D =
        std::shared_ptr<TreeNode>(new TreeNode("d", parD, NULL, CD));
    parD->Attach(D);
    myTree = std::shared_ptr<FunctionTree>(new FunctionTree(R));
    myTree->addNode(A);
    myTree->addNode(BCD);
    myTree->addNode(B);
    myTree->addNode(CD);
    myTree->addNode(C);
    myTree->addNode(D);

  }

  //------------Finished SetUp, now check Tree----------------

  std::cout << "Tree set up, not calculated" << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  //------------Trigger Calculation----------------
  myTree->recalculate();

  std::cout << "Tree calculated" << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  parB->SetValue(3.);

  std::cout << "Changed b from 2 to 3 " << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;

  myTree->recalculate();

  std::cout << "Changed b from 2 to 3 and recalculated " << std::endl;
  std::cout << std::endl << myTree << std::endl << std::endl;


  return 0;
}
