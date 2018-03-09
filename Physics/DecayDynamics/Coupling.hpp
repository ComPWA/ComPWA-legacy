// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef Coupling_h
#define Coupling_h

namespace ComPWA {
namespace Physics {
namespace DecayDynamics {

inline std::complex<double> couplingToWidth(double mSq, double mR, double g,
                                            double ma, double mb, double orbitL,
                                            double mesonRadius,
                                            formFactorType type,
                                            std::complex<double> phspFactor) {
  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // orbitL==0

  if (orbitL > 0 || type == formFactorType::CrystalBarrel) {
    std::complex<double> qV = qValue(mR, ma, mb);
    double ffR = FormFactor(mR, ma, mb, orbitL, mesonRadius, qV, type);
    std::complex<double> qR = std::pow(qV, orbitL);
    gammaA = ffR * qR;
  }

  // calculate phsp factor
  std::complex<double> res = std::norm(gammaA) * g * g * phspFactor / mR;

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::couplingToWidth() | "
                             "Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::couplingToWidth() | "
                             "Result is inf!");
#endif

  return res;
}

inline std::complex<double> couplingToWidth(double mSq, double mR, double g,
                                            double ma, double mb, double orbitL,
                                            double mesonRadius,
                                            formFactorType type) {
  double sqrtM = sqrt(mSq);
  std::complex<double> phspF = phspFactor(sqrtM, ma, mb);

  return couplingToWidth(mSq, mR, g, ma, mb, orbitL, mesonRadius, type, phspF);
}

inline std::complex<double> widthToCoupling(double mSq, double mR, double width,
                                            double ma, double mb, double orbitL,
                                            double mesonRadius,
                                            formFactorType type) {
  double sqrtS = sqrt(mSq);

  // calculate gammaA(s_R)
  std::complex<double> gammaA(1, 0); // orbitL==0
  std::complex<double> qV;
  if (orbitL > 0) {
    qV = qValue(mR, ma, mb);
    double ffR = FormFactor(mR, ma, mb, orbitL, mesonRadius, qV, type);
    std::complex<double> qR = std::pow(qV, orbitL);
    gammaA = ffR * qR;
  }

  // calculate phsp factor
  std::complex<double> rho = phspFactor(sqrtS, ma, mb);

  std::complex<double> denom = gammaA * sqrt(rho);
  std::complex<double> res = std::complex<double>(sqrt(mR * width), 0) / denom;

#ifndef NDEBUG
  // check for NaN
  if (std::isnan(res.real()) || std::isnan(res.imag()))
    throw std::runtime_error("AmpAbsDynamicalFunction::widthToCoupling() | "
                             "Result is NaN!");
  // check for inf
  if (std::isinf(res.real()) || std::isinf(res.imag()))
    throw std::runtime_error("AbstractDynamicalFunction::widthToCoupling() | "
                             "Result is inf!");
#endif

  return res;
}

class Coupling {
public:
  Coupling(double c, double massA, double massB)
      : _g(new FitParameter("", c)), _massA(massA), _massB(massB){};

  Coupling(std::shared_ptr<PartList> partL,
           const boost::property_tree::ptree tr) {
    _g = std::make_shared<FitParameter>();
    _g->load(tr.get_child(""));
    std::string nameA = tr.get<std::string>("ParticleA");
    std::string nameB = tr.get<std::string>("ParticleB");
    _massA = partL->find(nameA)->second.GetMass();
    _massB = partL->find(nameB)->second.GetMass();
  };

  void SetValueParameter(std::shared_ptr<FitParameter> g) { _g = g; }

  std::shared_ptr<FitParameter> GetValueParameter() { return _g; }

  double value() const { return _g->value(); }

  double GetMassA() const { return _massA; }

  double GetMassB() const { return _massB; }

  void SetMassA(double m) { _massA = m; }

  void SetMassB(double m) { _massB = m; }

protected:
  std::shared_ptr<FitParameter> _g;

  double _massA;

  double _massB;
};

} // ns::DecayDynamics
} // ns::Physics
} // ns::ComPWA

#endif
