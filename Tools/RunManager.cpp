// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "DataReader/Data.hpp"
#include "Core/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Event.hpp"
#include "Core/Generator.hpp"
#include "Core/ProgressBar.hpp"
#include "Tools/Integration.hpp"
#include "Tools/Generate.hpp"

#include "Tools/RunManager.hpp"

using namespace ComPWA::Tools;

RunManager::RunManager(std::shared_ptr<DataReader::Data> data,
                       std::shared_ptr<AmpIntensity> intens,
                       std::shared_ptr<Optimizer::Optimizer> optimizer)
    : sampleData_(data), opti_(optimizer), intens_(intens) {}

RunManager::RunManager(unsigned int size, std::shared_ptr<AmpIntensity> intens,
                       std::shared_ptr<Generator> gen)
    : gen_(gen), intens_(intens) {}

RunManager::~RunManager() {
  if (gen_)
    LOG(debug) << "~RunManager: Last seed: " << gen_->GetSeed();
}

std::shared_ptr<ComPWA::FitResult> RunManager::Fit(ParameterList &inPar) {
  LOG(info) << "RunManager::startFit() | Starting minimization.";

  // MINIMIZATION
  std::shared_ptr<FitResult> result = opti_->exec(inPar);

  LOG(info) << "RunManager::startFit() | Minimization finished!"
               " Result = "
            << result->GetResult() << ".";

  return result;
}

void RunManager::SetPhspSample(std::shared_ptr<ComPWA::DataReader::Data> phsp,
                               std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && truePhsp->GetNEvents() != phsp->GetNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");
  samplePhsp_ = phsp;
  sampleTruePhsp_ = truePhsp;
}

void RunManager::SetTruePhspSample(std::shared_ptr<DataReader::Data> truePhsp) {
  if (truePhsp && samplePhsp_ &&
      truePhsp->GetNEvents() != samplePhsp_->GetNEvents())
    throw std::runtime_error(
        "RunManager::setPhspSample() | "
        "Reconstructed sample and true sample have not the same size!");

  sampleTruePhsp_ = truePhsp;
}

bool RunManager::Generate(std::shared_ptr<Kinematics> kin, int number) {
  LOG(info) << "RunManager::generate() | "
               "Generating "
            << number << " signal events!";

  return ComPWA::Tools::Generate(number, kin, gen_, intens_, sampleData_,
                                 samplePhsp_, sampleTruePhsp_);
}

bool RunManager::GeneratePhsp(int nEvents) {
  return ComPWA::Tools::GeneratePhsp(nEvents, gen_, samplePhsp_);
}
