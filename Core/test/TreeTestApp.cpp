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

///
/// \file
/// This application test uses the ComPWA FunctionTree to calculate a simple
/// equation with cached intermediate results. The Tree can be set up in two
/// different way's, one using the functionality of the FunctionTree and the
/// other setting up the nodes and links manually. The second method is shown
/// mainly for a better understanding on how the FunctionTree internally works.
///

#define BOOST_TEST_MODULE Core

#include <iostream>
#include <memory>
#include <string>
#include <regex>
#include <algorithm>
#include <vector>
#include <map>

#include <boost/test/unit_test.hpp>

#include "Core/Functions.hpp"
#include "Core/TreeNode.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/FitParameter.hpp"
#include "Core/Value.hpp"

using namespace ComPWA;

BOOST_AUTO_TEST_SUITE(FunctionTreeTest);

BOOST_AUTO_TEST_CASE(SubTree) {
  ComPWA::Logging log("", "trace");

  std::shared_ptr<FitParameter> parA(new FitParameter("parA", 5.));
  parA->fixParameter(0);
  std::shared_ptr<FitParameter> parB(new FitParameter("parB", 2.));
  parB->fixParameter(0);
  std::shared_ptr<FitParameter> parC(new FitParameter("parC", 3.));
  parC->fixParameter(0);
  std::shared_ptr<FitParameter> parD(new FitParameter("parD", 1.));
  parD->fixParameter(0);

  // Calculate R = a * ( b + c * d)
  auto result = std::make_shared<Value<double>>();
  auto myTree = std::make_shared<FunctionTree>(
      "R", result, std::make_shared<MultAll>(ParType::DOUBLE));
  myTree->createLeaf("a", parA, "R");
  myTree->createNode("bcd", std::make_shared<AddAll>(ParType::DOUBLE),
                     "R"); // node is not cached since no parameter is passed
  myTree->createLeaf("b", parB, "bcd");

  auto subTree = std::make_shared<FunctionTree>(
      "cd", std::make_shared<Value<double>>("cd"),
      std::make_shared<MultAll>(ParType::DOUBLE));
  subTree->createLeaf("c", parC, "cd");
  subTree->createLeaf("d", parD, "cd");

  myTree->insertTree(subTree, "bcd");

  myTree->parameter(); // (re-)calculation
  LOG(INFO) << "Initial tree:" << std::endl << myTree;
  BOOST_CHECK_EQUAL(result->value(), 25);
  LOG(INFO) << "Changing parameter parD to 2!";
  parD->setValue(2.);
  myTree->parameter();
  LOG(INFO) << "Changed and recalculated tree:" << std::endl << myTree;
  BOOST_CHECK_EQUAL(result->value(), 40);
}

BOOST_AUTO_TEST_CASE(SingleParameters) {
  size_t nElements = 10;

  // Calculate R = Sum of (a * b)
  auto result = std::make_shared<Value<double>>();
  auto myTreeMult = std::make_shared<FunctionTree>(
      "R", result, std::make_shared<AddAll>(ParType::DOUBLE));
  myTreeMult->createNode("ab", std::make_shared<Value<double>>(),
                         std::make_shared<MultAll>(ParType::DOUBLE), "R");

  std::shared_ptr<Parameter> nParB(new FitParameter("parB", 2));
  myTreeMult->createLeaf("b", nParB, "ab");
  for (unsigned int i = 0; i < 10; i++) {
    std::string n = "parA_" + std::to_string(i);
    myTreeMult->createLeaf(n, std::make_shared<FitParameter>(n, i + 1),
                           "ab");
  }
  myTreeMult->parameter(); // Trigger recalculation
  BOOST_CHECK_EQUAL(result->value(), 7.2576e+06);
  
  LOG(INFO) << "Multitree Setup and calculated R = Sum[a*b] with one leaf "
               "containing "
            << nElements << " elements";
  LOG(INFO) << std::endl << myTreeMult;
}
BOOST_AUTO_TEST_CASE(MultiParameters) {

  //------------new Tree with multiDouble Par----------------
  size_t nElements = 5;
  std::shared_ptr<FunctionTree> myTreeMultD;

  //------------SetUp the parameters for R = Sum of (a * b)-----------
  std::vector<double> nMasses, nPhsp;
  std::shared_ptr<Parameter> mParB(new FitParameter("parB", 2));
  std::shared_ptr<Parameter> mParD(new FitParameter("parD", 3));
  for (unsigned int i = 0; i < nElements; i++) {
    // std::shared_ptr<FitParameter> tmpA(new
    // FitParameter("parA_"+i,i+1));
    nMasses.push_back(i + 1);
    nPhsp.push_back(2 * i + 1);
    nPhsp.push_back(2 * i + 2);
  }
  auto mParA = std::make_shared<Value<std::vector<double>>>("parA", nMasses);
  auto mParC = std::make_shared<Value<std::vector<double>>>("parC", nPhsp);
  auto result = std::make_shared<Value<double>>();
  
  myTreeMultD = std::make_shared<FunctionTree>(
      "R", result, std::make_shared<AddAll>(ParType::DOUBLE));
  myTreeMultD->createNode("Rmass", MDouble("par_Rnass", nElements),
                          std::make_shared<MultAll>(ParType::MDOUBLE), "R");
  myTreeMultD->createNode("ab", MDouble("par_ab", nElements),
                          std::make_shared<MultAll>(ParType::MDOUBLE), "Rmass");
  myTreeMultD->createLeaf("a", mParA, "ab");
  myTreeMultD->createLeaf("b", mParB, "ab");
  myTreeMultD->createNode("Rphsp", std::make_shared<AddAll>(ParType::DOUBLE),
                          "Rmass"); // this node will not be cached
  myTreeMultD->createNode("cd", MDouble("par_cd", 2 * nElements),
                          std::make_shared<MultAll>(ParType::MDOUBLE), "Rphsp");
  myTreeMultD->createLeaf("c", mParC, "cd");
  myTreeMultD->createLeaf("d", mParD, "cd");

  //------------Trigger Calculation----------------
  myTreeMultD->parameter();

  BOOST_CHECK_EQUAL(result->value(), 4950);

  LOG(INFO) << "MultiDouble Setup and calculated R = Sum[a*b] with one leaf "
               "containing "
            << nElements << " elements";
  LOG(INFO) << std::endl << myTreeMultD;
}

BOOST_AUTO_TEST_SUITE_END();
