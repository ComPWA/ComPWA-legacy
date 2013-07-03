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

  //------------SetUp some operations for Amp = a * ( b + c)-----------
  std::shared_ptr<Strategy> add = std::shared_ptr<Strategy>(new AddAll());
  std::shared_ptr<Strategy> mult = std::shared_ptr<Strategy>(new MultAll());

  //------------SetUp some nodes for Amp = a * ( b + c)----------------
  std::shared_ptr<TreeNode> Amplitude =
      std::shared_ptr<TreeNode>(new TreeNode("Amplitude", mult, NULL));
  std::shared_ptr<TreeNode> A =
      std::shared_ptr<TreeNode>(new TreeNode("a", NULL, Amplitude));
  A->changeVal(5);
  std::shared_ptr<TreeNode> BC =
      std::shared_ptr<TreeNode>(new TreeNode("bc", add, Amplitude));
  //Amplitude->addChild(A);
  //Amplitude->addChild(BC);
  std::shared_ptr<TreeNode> B =
      std::shared_ptr<TreeNode>(new TreeNode("b", NULL, BC));
  B->changeVal(2);
  std::shared_ptr<TreeNode> C =
      std::shared_ptr<TreeNode>(new TreeNode("c", NULL, BC));
  C->changeVal(3);
  //BC->addChild(B);
  //BC->addChild(C);
  FunctionTree test(Amplitude);
  test.addNode(A);
  test.addNode(BC);
  test.addNode(B);
  test.addNode(C);

  std::cout << "Tree set up, not calculated" << std::endl;
  std::cout << std::endl << test.head() << std::endl << std::endl;

  //------------Trigger Calculation----------------
  test.recalculate();

  std::cout << "Tree calculated" << std::endl;
  std::cout << std::endl << test.head() << std::endl << std::endl;

  A->changeVal(1.);

  std::cout << "Changed a from 5 to 1 " << std::endl;
  std::cout << std::endl << test.head() << std::endl << std::endl;

  test.recalculate();

  std::cout << "Changed a from 5 to 1 and recalculated " << std::endl;
  std::cout << std::endl << test.head() << std::endl << std::endl;


  return 0;
}
