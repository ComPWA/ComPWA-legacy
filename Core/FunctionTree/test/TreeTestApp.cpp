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

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "Core/FunctionTree/FitParameter.hpp"
#include "Core/FunctionTree/Functions.hpp"
#include "Core/FunctionTree/TreeNode.hpp"
#include "Core/FunctionTree/Value.hpp"

using namespace ComPWA::FunctionTree;

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
  auto myTree = std::make_shared<TreeNode>(
      result, std::make_shared<MultAll>(ParType::DOUBLE));

  auto CTimesD =
      std::make_shared<TreeNode>(std::make_shared<Value<double>>(),
                                 std::make_shared<MultAll>(ParType::DOUBLE));
  CTimesD->addNode(createLeaf(parC));
  CTimesD->addNode(createLeaf(parD));
  auto BPlusCD =
      std::make_shared<TreeNode>(std::make_shared<Value<double>>(),
                                 std::make_shared<AddAll>(ParType::DOUBLE));
  BPlusCD->addNode(createLeaf(parB));
  BPlusCD->addNode(CTimesD);
  myTree->addNode(createLeaf(parA));
  myTree->addNode(BPlusCD);
  myTree->parameter(); // (re-)calculation
  myTree->print();
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
  auto myTreeMult = std::make_shared<TreeNode>(
      result, std::make_shared<AddAll>(ParType::DOUBLE));

  auto ab =
      std::make_shared<TreeNode>(std::make_shared<Value<double>>(),
                                 std::make_shared<MultAll>(ParType::DOUBLE));
  myTreeMult->addNode(ab);
  auto b = createLeaf(std::make_shared<FitParameter>("parB", 2));
  ab->addNode(b);
  for (unsigned int i = 0; i < 10; i++) {
    std::string n = "parA_" + std::to_string(i);
    ab->addNode(createLeaf(std::make_shared<FitParameter>(n, i + 1)));
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
  std::shared_ptr<TreeNode> myTreeMultD;

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

  myTreeMultD = std::make_shared<TreeNode>(
      result, std::make_shared<AddAll>(ParType::DOUBLE));
  auto Rmass =
      std::make_shared<TreeNode>(MDouble("par_Rnass", nElements),
                                 std::make_shared<MultAll>(ParType::MDOUBLE));
  auto ab =
      std::make_shared<TreeNode>(MDouble("par_ab", nElements),
                                 std::make_shared<MultAll>(ParType::MDOUBLE));
  ab->addNodes({createLeaf(mParA), createLeaf(mParB)});

  auto Rphsp = std::make_shared<TreeNode>(std::make_shared<AddAll>(
      ParType::DOUBLE)); // this node will not be cached
  auto cd =
      std::make_shared<TreeNode>(MDouble("par_cd", 2 * nElements),
                                 std::make_shared<MultAll>(ParType::MDOUBLE));
  cd->addNodes({createLeaf(mParC), createLeaf(mParD)});
  Rphsp->addNode(cd);
  Rmass->addNodes({ab, Rphsp});
  myTreeMultD->addNode(Rmass);

  //------------Trigger Calculation----------------
  myTreeMultD->parameter();

  BOOST_CHECK_EQUAL(result->value(), 4950);

  LOG(INFO) << "MultiDouble Setup and calculated R = Sum[a*b] with one leaf "
               "containing "
            << nElements << " elements";
  LOG(INFO) << std::endl << myTreeMultD;
}

BOOST_AUTO_TEST_SUITE_END();
