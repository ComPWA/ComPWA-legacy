// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
//#include <math.h>
#include "Core/Value.hpp"
#include "Physics/DecayDynamics/AmpFlatteRes.hpp"
#include <limits>

using namespace ComPWA::Physics::DecayDynamics;

AmpFlatteRes::AmpFlatteRes(std::string name,
                           std::pair<std::string, std::string> daughters,
                           std::shared_ptr<ComPWA::PartList> partL) {

  LOG(TRACE) << "AmpFlatteRes::Factory() | Construction of " << name << ".";
  setName(name);

  // All further information on the decay is stored in a ParticleProperty list
  auto partProp = partL->find(name)->second;

  SetMassParameter(std::make_shared<FitParameter>(partProp.GetMassPar()));

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "flatte")
    throw std::runtime_error(
        "AmpFlatteRes::Factory() | Decay type does not match! ");

  auto spin = partProp.GetSpinQuantumNumber("Spin");
  SetSpin(spin);
  // in default, using spin J as Orbital Angular Momentum
  // update by calling SetOrbitalAngularMomentum() before any further process
  // after RelBW is created by calling of constructor
  SetOrbitalAngularMomentum(spin);

  auto ffType = formFactorType(decayTr.get<int>("FormFactor.<xmlattr>.Type"));
  SetFormFactorType(ffType);

  std::pair<double, double> daughterMasses(
      partL->find(daughters.first)->second.GetMass(),
      partL->find(daughters.second)->second.GetMass());

  SetDecayMasses(daughterMasses);
  SetDecayNames(daughters);

  // Read parameters
  std::vector<Coupling> vC;
  for (const auto &v : decayTr.get_child("")) {
    if (v.first != "Parameter")
      continue;
    std::string type = v.second.get<std::string>("<xmlattr>.Type");
    if (type == "Coupling") {
      vC.push_back(Coupling(partL, v.second));
    } else if (type == "MesonRadius") {
      auto mesonRadius = FitParameter();
      mesonRadius.load(v.second);
      SetMesonRadiusParameter(std::make_shared<FitParameter>(mesonRadius));
    } else {
      throw std::runtime_error("AmpFlatteRes::Factory() | Parameter of type " +
                               type + " is unknown.");
    }
  }
  SetCouplings(vC);

  LOG(TRACE) << "AmpFlatteRes::Factory() | Construction of the decay "
             << partProp.name() << " -> " << daughters.first << " + "
             << daughters.second;
}

AmpFlatteRes::~AmpFlatteRes() {}

bool AmpFlatteRes::isModified() const {
  if (GetMass() != Current_mass || Couplings.at(0).value() != Current_g ||
      Couplings.at(1).value() != Current_gHidden ||
      Couplings.at(2).value() != Current_gHidden2) {
    return true;
  }
  return false;
}

void AmpFlatteRes::setModified(bool b) {
  if (b) {
    Current_mass = std::numeric_limits<double>::quiet_NaN();
    Current_g = std::numeric_limits<double>::quiet_NaN();
    Current_gHidden = std::numeric_limits<double>::quiet_NaN();
    Current_gHidden2 = std::numeric_limits<double>::quiet_NaN();
  } else {
    Current_mass = Mass->value();
    Current_g = Couplings.at(0).value();
    Current_gHidden = Couplings.at(1).value();
    Current_gHidden2 = Couplings.at(2).value();
  }
}

std::complex<double> AmpFlatteRes::evaluate(const DataPoint &point,
                                            int pos) const {

  std::complex<double> result;
  try {
    result = dynamicalFunction(
        point.value(pos), Mass->value(), Couplings.at(0).GetMassA(),
        Couplings.at(0).GetMassB(), Couplings.at(0).value(),
        Couplings.at(1).GetMassA(), Couplings.at(1).GetMassB(),
        Couplings.at(1).value(), Couplings.at(2).GetMassA(),
        Couplings.at(2).GetMassB(), Couplings.at(2).value(), (double)L,
        MesonRadius->value(), FormFactorType);
  } catch (std::exception &ex) {
    LOG(ERROR) << "AmpFlatteRes::evaluate() | "
                  "Dynamical function can not be evaluated: "
               << ex.what();
    throw;
  }
  assert(!std::isnan(result.real()) && !std::isnan(result.imag()));
  return result;
}

std::complex<double> AmpFlatteRes::dynamicalFunction(
    double mSq, double mR, double gA, std::complex<double> termA,
    std::complex<double> termB, std::complex<double> termC) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  std::complex<double> denom = std::complex<double>(mR * mR - mSq, 0) +
                               (-1.0) * i * sqrtS * (termA + termB + termC);

  std::complex<double> result = std::complex<double>(gA, 0) / denom;

#ifndef NDEBUG
  if (std::isnan(result.real()) || std::isnan(result.imag())) {
    std::cout << "AmpFlatteRes::dynamicalFunction() | " << mR << " " << mSq
              << " " << termA << " " << termB << " " << termC << std::endl;
    return 0;
  }
#endif

  return result;
}

/// Helper function to calculate the coupling terms for the Flatte formular.
inline std::complex<double> flatteCouplingTerm(double sqrtS, double mR,
                                               double coupling, double massA,
                                               double massB, unsigned int J,
                                               double mesonRadius,
                                               formFactorType ffType) {
  auto qR = qValue(mR, massA, massB);
  auto phspR = phspFactor(sqrtS, massA, massB);
  auto ffR = FormFactor(qR, J, mesonRadius, ffType);
  auto barrierA = FormFactor(sqrtS, massA, massB, J, mesonRadius, ffType) / ffR;

  // Calculate normalized vertex functions vtxA(s_R)
  std::complex<double> vtxA(1, 0); // spin==0
  if (J > 0 || ffType == formFactorType::CrystalBarrel) {
    vtxA = ffR * std::pow(qR, J);
  }
  auto width = couplingToWidth(mR, coupling, vtxA, phspR);
  // Including the factor qTermA, as suggested by PDG 2014, Chapter 47.2,
  // leads to an amplitude that doesn't converge.
  //  qTermA = qValue(sqrtS,massA1,massA2) / qValue(mR,massA1,massA2);
  //  termA = gammaA * barrierA * barrierA * std::pow(qTermA, (double)2 * J +
  //  1);

  return (width * barrierA * barrierA);
}

std::complex<double>
AmpFlatteRes::dynamicalFunction(double mSq, double mR, double massA1,
                                double massA2, double gA, double massB1,
                                double massB2, double couplingB, double massC1,
                                double massC2, double couplingC, unsigned int L,
                                double mesonRadius, formFactorType ffType) {
  double sqrtS = sqrt(mSq);

  // channel A - signal channel
  auto termA =
      flatteCouplingTerm(sqrtS, mR, gA, massA1, massA2, L, mesonRadius, ffType);
  // channel B - hidden channel
  auto termB = flatteCouplingTerm(sqrtS, mR, couplingB, massB1, massB2, L,
                                  mesonRadius, ffType);

  // channel C - hidden channel
  std::complex<double> termC;
  if (couplingC != 0.0) {
    termC = flatteCouplingTerm(sqrtS, mR, couplingC, massC1, massC2, L,
                               mesonRadius, ffType);
  }
  return dynamicalFunction(mSq, mR, gA, termA, termB, termC);
}

std::shared_ptr<ComPWA::FunctionTree>
AmpFlatteRes::tree(const ParameterList &sample, int pos, std::string suffix) {

  size_t sampleSize = sample.mDoubleValue(pos)->values().size();
  auto tr = std::make_shared<FunctionTree>(
      "Flatte" + suffix, MComplex("", sampleSize),
      std::make_shared<FlatteStrategy>(""));

  tr->createLeaf("Mass", Mass, "Flatte" + suffix);
  for (int i = 0; i < Couplings.size(); i++) {
    tr->createLeaf("g_" + std::to_string(i) + "_massA",
                   Couplings.at(i).GetMassA(), "Flatte" + suffix);
    tr->createLeaf("g_" + std::to_string(i) + "_massB",
                   Couplings.at(i).GetMassB(), "Flatte" + suffix);
    tr->createLeaf("g_" + std::to_string(i),
                   Couplings.at(i).GetValueParameter(), "Flatte" + suffix);
  }
  tr->createLeaf("OrbitalAngularMomentum", (double)L, "Flatte" + suffix);
  tr->createLeaf("MesonRadius", MesonRadius, "Flatte" + suffix);
  tr->createLeaf("FormFactorType", FormFactorType, "Flatte" + suffix);
  //_daughterMasses actually not used here. But we put it in as a cross check.
  tr->createLeaf("MassA", DaughterMasses.first, "Flatte" + suffix);
  tr->createLeaf("MassB", DaughterMasses.second, "Flatte" + suffix);
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 sample.mDoubleValue(pos), "Flatte" + suffix);

  return tr;
}

void FlatteStrategy::execute(ParameterList &paras,
                             std::shared_ptr<Parameter> &out) {
  if (out && checkType != out->type())
    throw BadParameter("FlatteStrategy::execute() | Parameter type mismatch!");

#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("FlatteStrategy::execute() | "
                       "Output parameter is of type " +
                       std::string(ParNames[out->type()]) +
                       " and conflicts with expected type " +
                       std::string(ParNames[checkType])));

  // How many parameters do we expect?
  size_t check_nInt = 0;
  size_t nInt = paras.intValues().size();
  size_t check_nDouble = 15;
  size_t nDouble = paras.doubleValues().size();
  nDouble += paras.doubleParameters().size();
  size_t check_nComplex = 0;
  size_t nComplex = paras.complexValues().size();
  size_t check_nMInteger = 0;
  size_t nMInteger = paras.mIntValues().size();
  size_t check_nMDouble = 1;
  size_t nMDouble = paras.mDoubleValues().size();
  size_t check_nMComplex = 0;
  size_t nMComplex = paras.mComplexValues().size();

  // Check size of parameter list
  if (nInt != check_nInt)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(nInt) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (nDouble != check_nDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of FitParameters does not match: " +
                       std::to_string(nDouble) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (nComplex != check_nComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(nComplex) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (nMInteger != check_nMInteger)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiInt does not match: " +
                       std::to_string(nMInteger) + " given but " +
                       std::to_string(check_nMInteger) + " expected."));
  if (nMDouble != check_nMDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(nMDouble) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (nMComplex != check_nMComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(nMComplex) + " given but " +
                       std::to_string(check_nMComplex) + " expected."));
#endif

  size_t n = paras.mDoubleValue(0)->values().size();
  if (!out)
    out = MComplex("", n);

  auto par =
      std::static_pointer_cast<Value<std::vector<std::complex<double>>>>(out);
  auto &results = par->values(); // reference

  // calc function for each point
  for (size_t ele = 0; ele < n; ele++) {
    try {
      // Generally we need to add a factor q^{2J+1} to each channel term.
      // But since Flatte resonances are usually J=0 we neglect it here.
      results.at(ele) = AmpFlatteRes::dynamicalFunction(
          paras.mDoubleValue(0)->values().at(ele),
          paras.doubleParameter(0)->value(), // mass
          paras.doubleValue(0)->value(),     // g1_massA
          paras.doubleValue(1)->value(),     // g1_massB
          paras.doubleParameter(1)->value(), // g1
          paras.doubleValue(2)->value(),     // g2_massA
          paras.doubleValue(3)->value(),     // g2_massB
          paras.doubleParameter(2)->value(), // g2
          paras.doubleValue(4)->value(),     // g3_massA
          paras.doubleValue(5)->value(),     // g3_massB
          paras.doubleParameter(3)->value(), // g3
          paras.doubleValue(6)->value(),     // OrbitalAngularMomentum
          paras.doubleParameter(4)->value(), // mesonRadius
          formFactorType(paras.doubleValue(7)->value()) // ffType
          );
    } catch (std::exception &ex) {
      LOG(ERROR) << "FlatteStrategy::execute() | " << ex.what();
      throw(std::runtime_error("FlatteStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
}

void AmpFlatteRes::SetCouplings(std::vector<Coupling> vC) {
  if (vC.size() != 2 && vC.size() != 3)
    throw std::runtime_error(
        "AmpFlatteRes::SetCouplings() | Vector with "
        "couplings has a wrong size. We expect either 2 or 3 couplings.");

  Couplings = vC;

  if (Couplings.size() == 2)
    Couplings.push_back(Coupling(0.0, 0.0, 0.0));
  // Check if one of the  coupling match the final state (_daughterMasses)
  auto mm = GetDecayMasses();
  if (mm == std::pair<double, double>(-999, -999))
    LOG(INFO)
        << "AmpFlatteRes::SetCouplings() | Masses of decay products not set. "
           " Can not determine if correct couplings were set.";

  bool ok = false;
  for (auto i : Couplings) {
    if (i.GetMassA() == mm.first && i.GetMassB() == mm.second)
      ok = true;
    if (i.GetMassB() == mm.first && i.GetMassA() == mm.second)
      ok = true;
  }
  if (!ok)
    throw std::runtime_error("AmpFlatteRes::SetCouplings() | No couplings "
                             "for the current decay particles set!");
}

void AmpFlatteRes::parameters(ParameterList &list) {
  AbstractDynamicalFunction::parameters(list);
  // We check for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
  for (auto i : Couplings) {
    if (i.value() == 0.0)
      continue;
    i.GetValueParameter() = list.addUniqueParameter(i.GetValueParameter());
  }

  MesonRadius = list.addUniqueParameter(MesonRadius);
}

void AmpFlatteRes::updateParameters(const ParameterList &list) {
  // Try to update mesonRadius
  try {
    auto rad = FindParameter(MesonRadius->name(), list);
    MesonRadius->updateParameter(rad);
  } catch (BadParameter &ex) {
  }

  // Try to update Couplings
  for (auto i : Couplings) {
    try {
      auto g = FindParameter(i.GetValueParameter()->name(), list);
      i.GetValueParameter()->updateParameter(g);
    } catch (BadParameter &ex) {
    }
  }
  return;
}
