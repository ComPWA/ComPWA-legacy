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
//****************************************************************************
// Abstract base class for dynamical functions.
//****************************************************************************

// --CLASS DESCRIPTION [MODEL] --
// Abstract base class for dynamical functions.

#ifndef AMP_ABS_DYNAMICAL_FUNCTION
#define AMP_ABS_DYNAMICAL_FUNCTION

#include <vector>
#include <complex>

#include <boost/property_tree/ptree.hpp>

#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/PhysConst.hpp"
#include "Core/Resonance.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

namespace ComPWA {
namespace Physics {
namespace AmplitudeSum {

class AmpAbsDynamicalFunction : public Resonance {
public:
  AmpAbsDynamicalFunction(normStyle nS = normStyle::one, int calls = 30000);

  AmpAbsDynamicalFunction(
      const char *name, unsigned int varIdA, unsigned int varIdB,
      std::shared_ptr<DoubleParameter> mag,
      std::shared_ptr<DoubleParameter> phase,
      std::shared_ptr<DoubleParameter> mass, Spin spin, Spin m, Spin n, int P,
      int C, std::string mother, std::string particleA, std::string particleB,
      std::shared_ptr<DoubleParameter> mesonR,  //  meson radius
      std::shared_ptr<DoubleParameter> motherR, //  mother radius
      formFactorType type = formFactorType::BlattWeisskopf, int nCalls = 30000,
      normStyle nS = normStyle::one);

  AmpAbsDynamicalFunction(const char *name, unsigned int varIdA,
                          unsigned int varIdB,
                          std::shared_ptr<DoubleParameter> mag,
                          std::shared_ptr<DoubleParameter> phase,
                          std::shared_ptr<DoubleParameter> mass, Spin spin,
                          Spin m, Spin n, int P, int C, std::string mother,
                          std::string particleA, std::string particleB,
                          formFactorType type = formFactorType::BlattWeisskopf,
                          int nCalls = 30000, normStyle nS = normStyle::one);

  virtual ~AmpAbsDynamicalFunction();

  //! Configure resonance from ptree
  virtual void Configure(boost::property_tree::ptree::value_type const &v,
                         ParameterList &list);
  //! Save resonance to ptree
  virtual void Save(boost::property_tree::ptree &) = 0;
  //! Implementation of interface for streaming info about the strategy
  virtual std::string to_str() const;
  //! value of resonance at \param point
  virtual std::complex<double> Evaluate(dataPoint &point);
  //! value of dynamical amplitude at \param point
  virtual std::complex<double> EvaluateAmp(dataPoint &point) = 0;
  //! value of angular distribution at \param point
  virtual double EvaluateWignerD(dataPoint &point) {
    return _wignerD.evaluate(point);
  };
  //! value of angular distribution at \param point
  virtual double EvaluateAngular(dataPoint &point) {
    return EvaluateWignerD(point);
  };

  //! Calculation integral |dynamical amplitude|^2
  //	virtual double GetIntegral() const;

  /** Calculation integral |c * dynamical amplitude * WignerD|^2
   * Used to check the correct normalization of the amplitude. Should always be
   * 1.
   * @return
   */
  virtual double GetTotalIntegral() const;

  /**! Get current normalization.
   * In case that resonance parameters has change, it is recalculated.
   */
  virtual double GetNormalization();

  //! Set normalization style
  virtual void SetNormalizationStyle(normStyle n) { _normStyle = n; };
  //! Get normalization style
  virtual normStyle GetNormalizationStyle() const { return _normStyle; };
  //! Trigger recalculation of normalization
  virtual void SetModified() { _modified = 1; }
  //! Trigger recalculation of normalization
  virtual void CheckModified();
  //! Enable/Disable resonance
  virtual void SetEnable(bool s) { _enable = s; }
  //! Get Enable/Disable status
  virtual bool GetEnable() const { return _enable; }
  //! Get resonance name
  virtual std::string GetName() const { return _name; }
  //! Set resonance name
  virtual void SetName(std::string n) { _name = n; }

  //! Set prefactor
  virtual void SetPrefactor(std::complex<double> pre) { _prefactor = pre; }

  //! Get prefactor
  virtual std::complex<double> GetPrefactor() const { return _prefactor; }

  //! Get coefficient
  virtual std::complex<double> GetCoefficient() const;
  //! Get magnitude
  virtual double GetMagnitude() const { return _mag->GetValue(); };
  //! Get magnitude
  virtual std::shared_ptr<DoubleParameter> GetMagnitudePar() { return _mag; };
  //! Get phase
  virtual double GetPhase() const { return _phase->GetValue(); };
  //! Get phase
  virtual std::shared_ptr<DoubleParameter> GetPhasePar() { return _phase; };
  //! Get resonance mass
  virtual double GetMass() const { return _mass->GetValue(); };
  //! Get resonance mass
  virtual std::shared_ptr<DoubleParameter> GetMassPar() { return _mass; };
  //! Get resonance width
  virtual double GetWidth() const = 0;
  //! Get resonance spin
  virtual double GetSpin() const { return _spin; }
  virtual double GetM() const { return _m; };
  virtual double GetN() const { return _n; };
  //! Get mass of daughter A
  virtual double GetMassA() const { return _mass1; };
  //! Get mass of daughter B
  virtual double GetMassB() const { return _mass2; };
  //! Get resonance meson radius
  virtual double GetMesonRadius() const { return _mesonRadius->GetValue(); }
  //! Set resonance mother meson radius
  virtual void SetMesonRadius(double r) { _mesonRadius->SetValue(r); }
  //! Get resonance mother meson radius
  virtual double GetMotherRadius() const { return _motherRadius->GetValue(); }
  //! Set resonance mother meson radius
  virtual void SetMotherRadius(double r) { _motherRadius->SetValue(r); }

  virtual unsigned int GetVarIdA() const { return _subSys; };

  virtual unsigned int GetVarIdB() const { return _wignerD.GetVarId(); };

  virtual void SetVarIdA(unsigned int id) { _subSys = id; };

  virtual void SetVarIdB(unsigned int id) { _wignerD.SetVarId(id); };

  virtual std::shared_ptr<FunctionTree> SetupTree(ParameterList &sample,
                                                  ParameterList &toySample,
                                                  std::string suffix) = 0;

  /** Convert width of resonance to coupling
   *
   * Implementation of Eq.47-21 of PDG2014. Only valid for narrow, isolated
   * resonances.
   * @param mSq invariant mass
   * @param mR mass of resonance
   * @param width width of resonance in channel [a,b]
   * @param ma mass of particle a
   * @param mb mass of particle b
   * @param spin spin of the resonance
   * @param mesonRadius MesonRadius
   * @param type formfactor type
   * @return
   */
  static std::complex<double>
  widthToCoupling(double mSq, double mR, double width, double ma, double mb,
                  double spin, double mesonRadius,
                  formFactorType type = formFactorType::BlattWeisskopf);

  /** Convert coupling to width
   *
   * Convert coupling to channel (ma ,mb) to partial width. Only valid for
   * narrow, isolated resonances. Implementation of inverted Eq.47-21
   * of PDG2014.
   * @param mSq invariant mass
   * @param mR mass of resonance
   * @param g coupling to channel [a,b]
   * @param ma mass of particle a
   * @param mb mass of particle b
   * @param spin Spin of resonance
   * @param mesonRadius Meson radius of resonance
   * @param type Type of barrier factor
   * @return
   */
  static std::complex<double>
  couplingToWidth(double mSq, double mR, double g, double ma, double mb,
                  double spin, double mesonRadius,
                  formFactorType type = formFactorType::BlattWeisskopf);

  /** Convert coupling to width
   *
   * Convert coupling to channel (ma ,mb) to partial width. Only valid for
   * narrow, isolated resonances. Implementation of inverted
   * Eqs.47-21 of PDG2014.
   * @param mSq invariant mass
   * @param mR mass of resonance
   * @param g coupling to channel [a,b]
   * @param ma mass of particle a
   * @param mb mass of particle b
   * @param spin Spin of resonance
   * @param mesonRadius Meson radius of resonance
   * @param type Type of barrier factor
   * @param phspFactor Phase-space factor
   * @return
   */
  static std::complex<double> couplingToWidth(double mSq, double mR, double g,
                                              double ma, double mb, double spin,
                                              double mesonRadius,
                                              formFactorType type,
                                              std::complex<double> phspFactor);

protected:
  virtual void put(boost::property_tree::ptree &pt);

  void initialize();

  //! Calculation integral |dynamical amplitude|^2
  virtual double integral() const;

  /** Calculation integral |c * dynamical amplitude * WignerD|^2
   * Used to check the correct normalization of the amplitude. Should
   * always be 1.
   * @return
   */
  virtual double totalIntegral() const;

  //! Name of resonance
  std::string _name;
  //! enable/disable resonance
  bool _enable;

  //! Resonance shape was modified (recalculate the normalization)
  bool _modified;

  //! Integral value (temporary)
  double tmp_integral;

  //! Precision of MC integration
  int _nCalls;

  //! Pre factor
  std::complex<double> _prefactor;

  //! Mass of mother particle
  double _M;
  std::string _nameMother;

  //! Masses of daughter particles
  double _mass1, _mass2;
  //! Name of daughter particles
  std::string _name1, _name2;

  //! Resonance magnitude
  std::shared_ptr<DoubleParameter> _mag;
  bool _mag_writeByName;
  //! Resonance phase
  std::shared_ptr<DoubleParameter> _phase;
  bool _phase_writeByName;

  //! Resonance mass
  std::shared_ptr<DoubleParameter> _mass;
  double tmp_mass;
  bool _mass_writeByName;

  //! Type of resonance normalization
  normStyle _normStyle;

  //! Form factor type
  formFactorType _ffType;

  //! Barrier radi for resonance and mother particle
  std::shared_ptr<DoubleParameter> _mesonRadius, _motherRadius;
  bool _mesonRadius_writeByName;
  bool _motherRadius_writeByName;

  //! Resonance sub system
  unsigned int _subSys;

  //! Resonance spin
  ComPWA::Spin _spin, _m, _n;
  //! Parity of resonance +-
  int _parity;
  //! Charge parity of resonance +-0
  int _cparity;

  //! Angular distribution
  AmpWigner2 _wignerD;
};

class couplingToWidthStrat : public Strategy {
public:
  couplingToWidthStrat() : Strategy(ParType::MCOMPLEX) {}

  virtual const std::string to_str() const { return ("coupling to width "); }

  static std::shared_ptr<FunctionTree>
  SetupTree(std::shared_ptr<MultiDouble> mSq,
            std::shared_ptr<DoubleParameter> mR,
            std::shared_ptr<DoubleParameter> g, double ma, double mb, Spin spin,
            std::shared_ptr<DoubleParameter> mesonRadius, formFactorType type,
            std::string suffix = "");

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out);
};

class phspFactorStrat : public Strategy {
public:
  phspFactorStrat() : Strategy(ParType::MCOMPLEX) {}

  virtual const std::string to_str() const { return ("phspFactor"); }

  static std::shared_ptr<FunctionTree>
  SetupTree(std::shared_ptr<MultiDouble> mSq, double ma, double mb,
            std::string suffix = "") {

    std::shared_ptr<phspFactorStrat> thisStrat(new phspFactorStrat);

    std::string stratName = "phspFactorStrat" + suffix;
    //------------Setup Tree---------------------
    std::shared_ptr<FunctionTree> newTree(new FunctionTree());

    newTree->createHead(stratName, thisStrat, mSq->GetNValues());
    newTree->createLeaf("massA", ma, stratName);
    newTree->createLeaf("massB", mb, stratName);
    newTree->createLeaf("mSq", mSq, stratName);

    return newTree;
  }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out) {

#ifndef NDEBUG
    // Check parameter type
    if (checkType != out->type())
      throw(WrongParType("phspFactorStrat::execute() | "
                         "Output parameter is of type " +
                         std::string(ParNames[out->type()]) +
                         " and conflicts with expected type " +
                         std::string(ParNames[checkType])));

    // How many parameters do we expect?
    int check_nBool = 0;
    int check_nInt = 0;
    int check_nComplex = 0;
    int check_nDouble = 2;
    int check_nMComplex = 0;
    int check_nMDouble = 1;

    // Check size of parameter list
    if (paras.GetNBool() != check_nBool)
      throw(BadParameter("phspFactorStrat::execute() | "
                         "Number of BoolParameters does not match: " +
                         std::to_string(paras.GetNBool()) + " given but " +
                         std::to_string(check_nBool) + " expected."));
    if (paras.GetNInteger() != check_nInt)
      throw(BadParameter("phspFactorStrat::execute() | "
                         "Number of IntParameters does not match: " +
                         std::to_string(paras.GetNInteger()) + " given but " +
                         std::to_string(check_nInt) + " expected."));
    if (paras.GetNDouble() != check_nDouble)
      throw(BadParameter("phspFactorStrat::execute() | "
                         "Number of DoubleParameters does not match: " +
                         std::to_string(paras.GetNDouble()) + " given but " +
                         std::to_string(check_nDouble) + " expected."));
    if (paras.GetNComplex() != check_nComplex)
      throw(BadParameter("phspFactorStrat::execute() | "
                         "Number of ComplexParameters does not match: " +
                         std::to_string(paras.GetNComplex()) + " given but " +
                         std::to_string(check_nComplex) + " expected."));
    if (paras.GetNMultiDouble() != check_nMDouble)
      throw(BadParameter("phspFactorStrat::execute() | "
                         "Number of MultiDoubles does not match: " +
                         std::to_string(paras.GetNMultiDouble()) +
                         " given but " + std::to_string(check_nMDouble) +
                         " expected."));
    if (paras.GetNMultiComplex() != check_nMComplex)
      throw(BadParameter("phspFactorStrat::execute() | "
                         "Number of MultiComplexes does not match: " +
                         std::to_string(paras.GetNMultiComplex()) +
                         " given but " + std::to_string(check_nMComplex) +
                         " expected."));
#endif

    /** Get parameters from ParameterList:
     * We use the same order of the parameters as was used during tree
     * construction */
    double ma = paras.GetDoubleParameter(0)->GetValue();
    double mb = paras.GetDoubleParameter(1)->GetValue();

    std::vector<double> mp = paras.GetMultiDouble(0)->GetValues();

    std::vector<std::complex<double>> results(mp.size(),
                                              std::complex<double>(0., 0.));
    // calc function for each point
    for (unsigned int ele = 0; ele < mp.size(); ele++) {
      try {
        results.at(ele) = Kinematics::phspFactor(std::sqrt(mp.at(ele)), ma, mb);
      } catch (std::exception &ex) {
        BOOST_LOG_TRIVIAL(error) << "phspFactorStrat::execute() | "
                                 << ex.what();
        throw(std::runtime_error("phspFactorStrat::execute() | "
                                 "Evaluation of dynamic function failed!"));
      }
    }
    out = std::shared_ptr<AbsParameter>(
        new MultiComplex(out->GetName(), results));
    return true;
  }
};

class barrierStrat : public Strategy {
public:
  barrierStrat() : Strategy(ParType::MDOUBLE) {}

  virtual const std::string to_str() const { return ("barrierFactor"); }

  static std::shared_ptr<FunctionTree>
  SetupTree(std::shared_ptr<MultiDouble> mSq,
            std::shared_ptr<DoubleParameter> mR, double ma, double mb,
            Spin spin, std::shared_ptr<DoubleParameter> mesonRadius,
            formFactorType type, std::string suffix = "") {

    std::shared_ptr<barrierStrat> thisStrat(new barrierStrat);

    std::string stratName = "barrierFactor" + suffix;
    //------------Setup Tree---------------------
    std::shared_ptr<FunctionTree> newTree(new FunctionTree());
    newTree->createHead(stratName, thisStrat, mSq->GetNValues());
    newTree->createLeaf("mass", mR, stratName);
    newTree->createLeaf("massA1", ma, stratName);
    newTree->createLeaf("massA2", mb, stratName);
    newTree->createLeaf("spin", (double)spin, stratName);
    newTree->createLeaf("mesonRadius", mesonRadius, stratName);
    newTree->createLeaf("formFactorType", type, stratName);
    newTree->createLeaf("mSq", mSq, stratName);

    return newTree;
  }

  virtual bool execute(ParameterList &paras,
                       std::shared_ptr<AbsParameter> &out) {

#ifndef NDEBUG
    // Check parameter type
    if (checkType != out->type())
      throw(WrongParType("barrierStrat::execute() | "
                         "Output parameter is of type " +
                         std::string(ParNames[out->type()]) +
                         " and conflicts with expected type " +
                         std::string(ParNames[checkType])));

    // How many parameters do we expect?
    int check_nBool = 0;
    int check_nInt = 0;
    int check_nComplex = 0;
    int check_nDouble = 6;
    int check_nMDouble = 1;
    int check_nMComplex = 0;

    // Check size of parameter list
    if (paras.GetNBool() != check_nBool)
      throw(BadParameter("barrierStrat::execute() | "
                         "Number of BoolParameters does not match: " +
                         std::to_string(paras.GetNBool()) + " given but " +
                         std::to_string(check_nBool) + " expected."));
    if (paras.GetNInteger() != check_nInt)
      throw(BadParameter("barrierStrat::execute() | "
                         "Number of IntParameters does not match: " +
                         std::to_string(paras.GetNInteger()) + " given but " +
                         std::to_string(check_nInt) + " expected."));
    if (paras.GetNDouble() != check_nDouble)
      throw(BadParameter("barrierStrat::execute() | "
                         "Number of DoubleParameters does not match: " +
                         std::to_string(paras.GetNDouble()) + " given but " +
                         std::to_string(check_nDouble) + " expected."));
    if (paras.GetNComplex() != check_nComplex)
      throw(BadParameter("barrierStrat::execute() | "
                         "Number of ComplexParameters does not match: " +
                         std::to_string(paras.GetNComplex()) + " given but " +
                         std::to_string(check_nComplex) + " expected."));
    if (paras.GetNMultiDouble() != check_nMDouble)
      throw(BadParameter("barrierStrat::execute() | "
                         "Number of MultiDoubles does not match: " +
                         std::to_string(paras.GetNMultiDouble()) +
                         " given but " + std::to_string(check_nMDouble) +
                         " expected."));
    if (paras.GetNMultiComplex() != check_nMComplex)
      throw(BadParameter("barrierStrat::execute() | "
                         "Number of MultiComplexes does not match: " +
                         std::to_string(paras.GetNMultiComplex()) +
                         " given but " + std::to_string(check_nMComplex) +
                         " expected."));

#endif

    /** Get parameters from ParameterList:
     * We use the same order of the parameters as was used during tree
     * construction
     */
    double mR = paras.GetDoubleParameter(0)->GetValue();
    double ma = paras.GetDoubleParameter(1)->GetValue();
    double mb = paras.GetDoubleParameter(2)->GetValue();
    double spin = paras.GetDoubleParameter(3)->GetValue();
    double mesonRadius = paras.GetDoubleParameter(4)->GetValue();
    formFactorType type =
        formFactorType(paras.GetDoubleParameter(5)->GetValue());

    std::vector<double> mp = paras.GetMultiDouble(0)->GetValues();

    // Initialize results with one
    std::vector<double> results(mp.size(), 1.);

    // If form factors are one anyway, we skip the loop
    if (spin == 0 && type != formFactorType::CrystalBarrel) {
      out = std::shared_ptr<AbsParameter>(
          new MultiDouble(out->GetName(), results));
      return true;
    }
    // calc function for each point
    for (unsigned int ele = 0; ele < mp.size(); ele++) {
      double s = mp.at(ele);
      double sqrtS = sqrt(s);
      std::complex<double> qValue = Kinematics::qValue(sqrtS, ma, mb);
      std::complex<double> qRValue = Kinematics::qValue(mR, ma, mb);

      try {
        double nom = Kinematics::FormFactor(sqrtS, ma, mb, spin, mesonRadius,
                                            qValue, type);
        double denom = Kinematics::FormFactor(mR, ma, mb, spin, mesonRadius,
                                              qRValue, type);

        results.at(ele) = nom / denom * nom / denom;
      } catch (std::exception &ex) {
        BOOST_LOG_TRIVIAL(error) << "barrierStrat::execute() | " << ex.what();
        throw(std::runtime_error("barrierStrat::execute() | "
                                 "Evaluation of dynamic function failed!"));
      }
    }
    out =
        std::shared_ptr<AbsParameter>(new MultiDouble(out->GetName(), results));
    return true;
  }
};

} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
