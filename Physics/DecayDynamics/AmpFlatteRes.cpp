// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cmath>
#include <math.h>
#include "Physics/DecayDynamics/AmpFlatteRes.hpp"

using namespace ComPWA::Physics::DecayDynamics;

AmpFlatteRes::~AmpFlatteRes() {}

std::shared_ptr<AbstractDynamicalFunction>
AmpFlatteRes::Factory(std::shared_ptr<PartList> partL,
                      const boost::property_tree::ptree &pt) {
  auto obj = std::make_shared<AmpFlatteRes>();

  std::string name = pt.get<std::string>("DecayParticle.<xmlattr>.Name");
  LOG(trace) << "AmpFlatteRes::Factory() | Construction of " << name << ".";
  obj->setName(name);

  // All further information on the decay is stored in a ParticleProperty list
  auto partProp = partL->find(name)->second;

  obj->SetMassParameter(
      std::make_shared<DoubleParameter>(partProp.GetMassPar()));

  auto decayTr = partProp.GetDecayInfo();
  if (partProp.GetDecayType() != "flatte")
    throw std::runtime_error(
        "AmpFlatteRes::Factory() | Decay type does not match! ");

  auto spin = partProp.GetSpinQuantumNumber("Spin");
  obj->SetSpin(spin);

  auto ffType = formFactorType(decayTr.get<int>("FormFactor.<xmlattr>.Type"));
  obj->SetFormFactorType(ffType);

  // Get masses of decay products
  auto decayProducts = pt.get_child("SubSystem.DecayProducts");
  if (decayProducts.size() != 2)
    throw boost::property_tree::ptree_error(
        "AmpWignerD::Factory() | Expect exactly two decay products (" +
        std::to_string(decayProducts.size()) + " given)!");

  auto firstItr = decayProducts.begin();
  auto secondItr = decayProducts.begin();
  secondItr++;

  // Mass of daughter particles needs to be set before the coupling constants
  std::pair<std::string, std::string> daughterNames(
      firstItr->second.get<std::string>("<xmlattr>.Name"),
      secondItr->second.get<std::string>("<xmlattr>.Name"));
  std::pair<double, double> daughterMasses(
      partL->find(daughterNames.first)->second.GetMass(),
      partL->find(daughterNames.second)->second.GetMass());

  obj->SetDecayMasses(daughterMasses);
  obj->SetDecayNames(daughterNames);

  // Read parameters
  std::vector<Coupling> vC;
  for (const auto &v : decayTr.get_child("")) {
    if (v.first != "Parameter")
      continue;
    std::string type = v.second.get<std::string>("<xmlattr>.Type");
    if (type == "Coupling") {
      vC.push_back(Coupling(partL, v.second));
    } else if (type == "MesonRadius") {
      auto mesonRadius = DoubleParameter();
      mesonRadius.load(v.second);
      obj->SetMesonRadiusParameter(
          std::make_shared<DoubleParameter>(mesonRadius));
    } else {
      throw std::runtime_error("AmpFlatteRes::Factory() | Parameter of type " +
                               type + " is unknown.");
    }
  }
  obj->SetCouplings(vC);

  LOG(trace)
      << "RelativisticBreitWigner::Factory() | Construction of the decay "
      << partProp.name() << " -> " << daughterNames.first << " + "
      << daughterNames.second;

  return std::static_pointer_cast<AbstractDynamicalFunction>(obj);
}

bool AmpFlatteRes::CheckModified() const {
  if (AbstractDynamicalFunction::CheckModified())
    return true;
  if (_g.at(0).value() != _current_g ||
      _g.at(1).value() != _current_gHidden ||
      _g.at(2).value() != _current_gHidden2) {
    SetModified();
    const_cast<double &>(_current_g) = _g.at(0).value();
    const_cast<double &>(_current_gHidden) = _g.at(1).value();
    const_cast<double &>(_current_gHidden2) = _g.at(2).value();
    return true;
  }
  return false;
}

std::complex<double> AmpFlatteRes::evaluate(const DataPoint &point,
                                            int pos) const {

  std::complex<double> result;
  try {
    result = dynamicalFunction(
        point.value(pos), _mass->value(), _g.at(0).GetMassA(),
        _g.at(0).GetMassB(), _g.at(0).value(), _g.at(1).GetMassA(),
        _g.at(1).GetMassB(), _g.at(1).value(), _g.at(2).GetMassA(),
        _g.at(2).GetMassB(), _g.at(2).value(), (double)_spin,
        _mesonRadius->value(), _ffType);
  } catch (std::exception &ex) {
    LOG(error) << "AmpFlatteRes::EvaluateAmp() | "
                  "Dynamical function can not be evalutated: "
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

std::complex<double>
AmpFlatteRes::dynamicalFunction(double mSq, double mR, double massA1,
                                double massA2, double gA, double massB1,
                                double massB2, double couplingB, double massC1,
                                double massC2, double couplingC, unsigned int J,
                                double mesonRadius, formFactorType ffType) {
  std::complex<double> i(0, 1);
  double sqrtS = sqrt(mSq);

  // channel A - signal channel
  std::complex<double> gammaA, qTermA, termA;
  double barrierA;
  // break-up momentum
  barrierA = FormFactor(sqrtS, massA1, massA2, J, mesonRadius, ffType) /
             FormFactor(mR, massA1, massA2, J, mesonRadius, ffType);
  // convert coupling to partial width of channel A
  gammaA = couplingToWidth(mSq, mR, gA, massA1, massA2, J, mesonRadius, ffType);
  // including the factor qTermA, as suggested by PDG, leads to an amplitude
  // that doesn't converge.
  //		qTermA = Kinematics::qValue(sqrtS,massA1,massA2) /
  // Kinematics::qValue(mR,massA1,massA2);
  qTermA = std::complex<double>(1, 0);
  termA = gammaA * barrierA * barrierA * std::pow(qTermA, (double)2 * J + 1);

  // channel B - hidden channel
  std::complex<double> gammaB, qTermB, termB;
  double barrierB, gB;
  // break-up momentum
  barrierB = FormFactor(sqrtS, massB1, massB2, J, mesonRadius, ffType) /
             FormFactor(mR, massB1, massB2, J, mesonRadius, ffType);
  gB = couplingB;
  // convert coupling to partial width of channel B
  gammaB = couplingToWidth(mSq, mR, gB, massB1, massB2, J, mesonRadius, ffType);
  //		qTermB = Kinematics::qValue(sqrtS,massB1,massB2) /
  // Kinematics::qValue(mR,massB1,massB2);
  qTermB = std::complex<double>(1, 0);
  termB = gammaB * barrierB * barrierB * std::pow(qTermB, (double)2 * J + 1);

  // channel C - hidden channel
  std::complex<double> gammaC, qTermC, termC;
  double barrierC, gC;
  if (couplingC != 0.0) {
    // break-up momentum
    barrierC = FormFactor(sqrtS, massC1, massC2, J, mesonRadius, ffType) /
               FormFactor(mR, massC1, massC2, J, mesonRadius, ffType);
    gC = couplingC;
    // convert coupling to partial width of channel C
    gammaC =
        couplingToWidth(mSq, mR, gC, massC1, massC2, J, mesonRadius, ffType);
    //		qTermC = Kinematics::qValue(sqrtS,massC1,massC2) /
    // Kinematics::qValue(mR,massC1,massC2);
    qTermC = std::complex<double>(1, 0);
    termC = gammaC * barrierC * barrierC * std::pow(qTermC, (double)2 * J + 1);
  }

  return dynamicalFunction(mSq, mR, gA, termA, termB, termC);
}

std::shared_ptr<ComPWA::FunctionTree>
AmpFlatteRes::GetTree(const ParameterList &sample, int pos,
                      std::string suffix) {

  //  int sampleSize = sample.GetMultiDouble(0)->numValues();

  std::shared_ptr<FunctionTree> tr(new FunctionTree());
  //  tr->createHead("DynamicalFunction",
  //                 std::shared_ptr<Strategy>(new MultAll(ParType::MCOMPLEX)));
  tr->createHead("Flatte" + suffix,
                 std::shared_ptr<Strategy>(new FlatteStrategy("")));

  //  tr->createNode("Flatte", std::shared_ptr<Strategy>(new
  //  FlatteStrategy("")),
  //                 "DynamicalFunction", sampleSize);
  tr->createLeaf("Mass", _mass, "Flatte" + suffix);
  for (int i = 0; i < _g.size(); i++) {
    tr->createLeaf("g_" + std::to_string(i) + "_massA", _g.at(i).GetMassA(),
                   "Flatte" + suffix);
    tr->createLeaf("g_" + std::to_string(i) + "_massB", _g.at(i).GetMassB(),
                   "Flatte" + suffix);
    tr->createLeaf("g_" + std::to_string(i), _g.at(i).GetValueParameter(),
                   "Flatte" + suffix);
  }
  tr->createLeaf("Spin", (double)_spin, "Flatte" + suffix);
  tr->createLeaf("MesonRadius", _mesonRadius, "Flatte" + suffix);
  tr->createLeaf("FormFactorType", _ffType, "Flatte" + suffix);
  //_daughterMasses actually not used here. But we put it in as a cross check.
  tr->createLeaf("MassA", _daughterMasses.first, "Flatte" + suffix);
  tr->createLeaf("MassB", _daughterMasses.second, "Flatte" + suffix);
  tr->createLeaf("Data_mSq[" + std::to_string(pos) + "]",
                 sample.GetMultiDouble(pos), "Flatte" + suffix);

  return tr;
}

bool FlatteStrategy::execute(ParameterList &paras,
                             std::shared_ptr<Parameter> &out) {
#ifndef NDEBUG
  // Check parameter type
  if (checkType != out->type())
    throw(WrongParType("FlatteStrategy::execute() | "
                       "Output parameter is of type " +
                       std::string(ParNames[out->type()]) +
                       " and conflicts with expected type " +
                       std::string(ParNames[checkType])));

  // How many parameters do we expect?
  int check_nBool = 0;
  int check_nInt = 0;
  int check_nComplex = 0;
  int check_nDouble = 15;
  int check_nMDouble = 1;
  int check_nMComplex = 0;

  // Check size of parameter list
  if (paras.GetNBool() != check_nBool)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of BoolParameters does not match: " +
                       std::to_string(paras.GetNBool()) + " given but " +
                       std::to_string(check_nBool) + " expected."));
  if (paras.GetNInteger() != check_nInt)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of IntParameters does not match: " +
                       std::to_string(paras.GetNInteger()) + " given but " +
                       std::to_string(check_nInt) + " expected."));
  if (paras.GetNDouble() != check_nDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of DoubleParameters does not match: " +
                       std::to_string(paras.GetNDouble()) + " given but " +
                       std::to_string(check_nDouble) + " expected."));
  if (paras.GetNComplex() != check_nComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of ComplexParameters does not match: " +
                       std::to_string(paras.GetNComplex()) + " given but " +
                       std::to_string(check_nComplex) + " expected."));
  if (paras.GetNMultiDouble() != check_nMDouble)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiDoubles does not match: " +
                       std::to_string(paras.GetNMultiDouble()) + " given but " +
                       std::to_string(check_nMDouble) + " expected."));
  if (paras.GetNMultiComplex() != check_nMComplex)
    throw(BadParameter("FlatteStrategy::execute() | "
                       "Number of MultiComplexes does not match: " +
                       std::to_string(paras.GetNMultiComplex()) +
                       " given but " + std::to_string(check_nMComplex) +
                       " expected."));
#endif

  // invariant masses
  std::vector<double> mSq = paras.GetMultiDouble(0)->values();

  std::vector<std::complex<double>> results(mSq.size(),
                                            std::complex<double>(0., 0.));
  // TODO make sure all vectors have to same size

  // calc function for each point
  for (unsigned int ele = 0; ele < mSq.size(); ele++) {
    try {
      /* Generally we need to add a factor q^{2J+1} to each channel term.
       * But since Flatte resonances are usually J=0 we neglect it here.*/
      results.at(ele) = AmpFlatteRes::dynamicalFunction(
          mSq.at(ele),
          paras.GetDoubleParameter(0)->value(),  // mass
          paras.GetDoubleParameter(1)->value(),  // g1_massA
          paras.GetDoubleParameter(2)->value(),  // g1_massB
          paras.GetDoubleParameter(3)->value(),  // g1
          paras.GetDoubleParameter(4)->value(),  // g2_massA
          paras.GetDoubleParameter(5)->value(),  // g2_massB
          paras.GetDoubleParameter(6)->value(),  // g2
          paras.GetDoubleParameter(7)->value(),  // g3_massA
          paras.GetDoubleParameter(8)->value(),  // g3_massB
          paras.GetDoubleParameter(9)->value(),  // g3
          paras.GetDoubleParameter(10)->value(), // Spin
          paras.GetDoubleParameter(11)->value(), // mesonRadius
          formFactorType(paras.GetDoubleParameter(12)->value()) // ffType
          );
    } catch (std::exception &ex) {
      LOG(error) << "FlatteStrategy::execute() | " << ex.what();
      throw(std::runtime_error("FlatteStrategy::execute() | "
                               "Evaluation of dynamic function failed!"));
    }
  }
  out =
      std::shared_ptr<Parameter>(new MultiComplex(out->name(), results));
  return true;
}

void AmpFlatteRes::SetCouplings(std::vector<Coupling> vC) {
  if (vC.size() != 2 && vC.size() != 3)
    throw std::runtime_error(
        "AmpFlatteRes::SetCouplings() | Vector with "
        "couplings has a wrong size. We expect either 2 or 3 couplings.");

  _g = vC;

  if (_g.size() == 2)
    _g.push_back( Coupling(0.0, 0.0, 0.0) );
  // Check if one of the  coupling match the final state (_daughterMasses)
  auto mm = GetDecayMasses();
  if (mm == std::pair<double, double>(-999, -999))
    LOG(info)
        << "AmpFlatteRes::SetCouplings() | Masses of decay products not set. "
           " Can not determine if correct couplings were set.";

  bool ok = false;
  for (auto i : _g) {
    if (i.GetMassA() == mm.first && i.GetMassB() == mm.second)
      ok = true;
    if (i.GetMassB() == mm.first && i.GetMassA() == mm.second)
      ok = true;
  }
  if (!ok)
    throw std::runtime_error("AmpFlatteRes::SetCouplings() | No couplings "
                             "for the current decay particles set!");
}

void AmpFlatteRes::GetParameters(ParameterList &list) {
  AbstractDynamicalFunction::GetParameters(list);
  // We check for each parameter if a parameter of the same name exists in
  // list. If so we check if both are equal and set the local parameter to the
  // parameter from the list. In this way we connect parameters that occur on
  // different positions in the amplitude.
  std::shared_ptr<DoubleParameter> tmp;
  for (auto i : _g) {
    if (i.value() == 0.0)
      continue;
    try { // catch BadParameter
      tmp = list.GetDoubleParameter(i.GetValueParameter()->name());
      try { // catch and throw std::runtime_error due to failed parameter
            // comparisson
        if (*tmp == *i.GetValueParameter())
          i.GetValueParameter() = tmp;
      } catch (std::exception &ex) {
        throw;
      }
    } catch (BadParameter &ex) {
      list.AddParameter(i.GetValueParameter());
    }
  }

  try { // catch BadParameter
    tmp = list.GetDoubleParameter(GetMesonRadiusParameter()->name());
    try { // catch and throw std::runtime_error due to failed parameter
          // comparisson
      if (*tmp == *_mesonRadius)
        SetMesonRadiusParameter(tmp);
    } catch (std::exception &ex) {
      throw;
    }
  } catch (BadParameter &ex) {
    list.AddParameter(GetMesonRadiusParameter());
  }
}

void AmpFlatteRes::updateParameters(const ParameterList &par) {

  // Try to update mesonRadius
  std::shared_ptr<DoubleParameter> rad;
  try {
    rad = par.GetDoubleParameter(_mesonRadius->name());
  } catch (std::exception &ex) {
  }
  if (rad)
    _mesonRadius->updateParameter(rad);

  // Try to update Couplings
  for (auto i : _g) {
    std::shared_ptr<DoubleParameter> c;
    try {
      c = par.GetDoubleParameter(i.GetValueParameter()->name());
    } catch (std::exception &ex) {
    }
    if (c)
      i.GetValueParameter()->updateParameter(rad);
  }

  return;
}
