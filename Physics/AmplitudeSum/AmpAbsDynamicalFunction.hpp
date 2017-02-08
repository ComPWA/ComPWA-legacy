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


class AmpAbsDynamicalFunction : public Resonance
{
public:
	AmpAbsDynamicalFunction( normStyle nS=normStyle::one, int calls=30000);

	AmpAbsDynamicalFunction(const char *name,
			unsigned int varIdA, unsigned int varIdB,
			std::shared_ptr<DoubleParameter> mag,
			std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass,
			Spin spin, Spin m, Spin n, int P, int C,
			std::string mother, std::string particleA, std::string particleB,
			std::shared_ptr<DoubleParameter> mesonR, //  meson radius
			std::shared_ptr<DoubleParameter> motherR, //  mother radius
			formFactorType type = formFactorType::BlattWeisskopf,
			int nCalls=30000, normStyle nS=normStyle::one);

	AmpAbsDynamicalFunction(const char *name,
			unsigned int varIdA, unsigned int varIdB,
			std::shared_ptr<DoubleParameter> mag,
			std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass,
			Spin spin, Spin m, Spin n, int P, int C,
			std::string mother, std::string particleA, std::string particleB,
			formFactorType type = formFactorType::BlattWeisskopf,
			int nCalls=30000, normStyle nS=normStyle::one);

	virtual ~AmpAbsDynamicalFunction();



	//! Configure resonance from ptree
	virtual void Configure(boost::property_tree::ptree::value_type const& v,
			ParameterList& list);
	//! Save resonance to ptree
	virtual void Save(boost::property_tree::ptree&) = 0;
	//! Implementation of interface for streaming info about the strategy
	virtual std::string to_str() const;
	//! value of resonance at \param point
	virtual std::complex<double> Evaluate(dataPoint& point);
	//! value of resonance's angular distribution at \param point
	virtual double EvaluateAngular(dataPoint& point);
	//! value of dynamical amplitude at \param point
	virtual std::complex<double> EvaluateAmp(dataPoint& point) = 0;
	//! value of angular distribution at \param point
	virtual double EvaluateWignerD(dataPoint& point) {
		return _wignerD.evaluate(point);
	};

	//! Calculation integral |dynamical amplitude|^2
	//	virtual double GetIntegral() const;

	/** Calculation integral |c * dynamical amplitude * WignerD|^2
	 * Used to check the correct normalization of the amplitude. Should always be 1.
	 * @return
	 */
	virtual double GetTotalIntegral() const;

	/**! Get current normalization.
	 * In case that resonance parameters has change, it is recalculated.
	 */
	virtual double GetNormalization();

	//! Set normalization style
	virtual void SetNormalizationStyle(normStyle n){ _normStyle=n; };
	//! Get normalization style
	virtual normStyle GetNormalizationStyle() const { return _normStyle; };
	//! Trigger recalculation of normalization
	virtual void SetModified() { _modified=1; }
	//! Trigger recalculation of normalization
	virtual void CheckModified();
	//! Enable/Disable resonance
	virtual void SetEnable(bool s){ _enable=s; }
	//! Get Enable/Disable status
	virtual bool GetEnable() const { return _enable; }
	//! Get resonance name
	virtual std::string GetName() const { return _name; }
	//! Set resonance name
	virtual void SetName(std::string n){ _name = n; }

	//!Set prefactor
	virtual void SetPrefactor(std::complex<double> pre) { _prefactor = pre; }

	//!Get prefactor
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
	virtual double GetMass() const {return _mass->GetValue();};
	//! Get resonance mass
	virtual std::shared_ptr<DoubleParameter> GetMassPar() {return _mass;};
	//! Get resonance width
	virtual double GetWidth() const = 0;
	//! Get resonance spin
	virtual double GetSpin() const { return _spin.Val(); }
	virtual double GetM() const {return _m.Val();};
	virtual double GetN() const {return _n.Val();};
	//! Get mass of daughter A
	virtual double GetMassA() const {return _mass1;};
	//! Get mass of daughter B
	virtual double GetMassB() const {return _mass2;};
	//! Get resonance meson radius
	virtual double GetMesonRadius() const { return _mesonRadius->GetValue(); }
	//! Set resonance mother meson radius
	virtual void SetMesonRadius(double r) {	_mesonRadius->SetValue(r); }
	//! Get resonance mother meson radius
	virtual double GetMotherRadius() const { return _motherRadius->GetValue(); }
	//! Set resonance mother meson radius
	virtual void SetMotherRadius(double r) { _motherRadius->SetValue(r); }

	virtual unsigned int GetVarIdA() const { return _subSys; };

	virtual unsigned int GetVarIdB() const { return _wignerD.GetVarId(); };

	virtual void SetVarIdA(unsigned int id) { _subSys = id; };

	virtual void SetVarIdB(unsigned int id) { _wignerD.SetVarId(id); };

	virtual std::shared_ptr<FunctionTree> SetupTree(
			ParameterList& sample, ParameterList& toySample,
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
	static std::complex<double> widthToCoupling(double mSq, double mR,
			double width, double ma, double mb, double spin, double mesonRadius,
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
	static std::complex<double> couplingToWidth(double mSq, double mR, double g,
			double ma, double mb, double spin, double mesonRadius,
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
			double ma, double mb, double spin, double mesonRadius,
			formFactorType type, std::complex<double> phspFactor);

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

	//!Form factor type
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

class couplingToWidthStrat: public Strategy
{
public:
	couplingToWidthStrat() :
		Strategy(ParType::MCOMPLEX) { }

	virtual const std::string to_str() const {
		return ("coupling to width "); }

	static std::shared_ptr<FunctionTree> SetupTree(
			std::shared_ptr<MultiDouble> mSq,
			std::shared_ptr<DoubleParameter> mR,
			std::shared_ptr<DoubleParameter> g, double ma, double mb,
			Spin spin, std::shared_ptr<DoubleParameter> mesonRadius,
			formFactorType type);

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out);
};
} /* namespace AmplitudeSum */
} /* namespace Physics */
} /* namespace ComPWA */

#endif
