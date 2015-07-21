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

#include "Core/Parameter.hpp"
#include "Core/DataPoint.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/PhysConst.hpp"
#include "Physics/AmplitudeSum/AmpWigner2.hpp"

enum normStyle {
	none, /*!< no normaliztion between Amplitudes. */
	/*!< all amplitudes are normalized to one.
	 *  The normalization factor is \f$ 1/\sqrt(\int |A|^2)\f$ */
	one
};

class AmpAbsDynamicalFunction
{
public:
	AmpAbsDynamicalFunction(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass, int subsys, Spin spin, Spin m, Spin n,
			std::shared_ptr<DoubleParameter> mesonR, //  meson radius
			std::shared_ptr<DoubleParameter> motherR, //  mother radius
			int nCalls=30000, normStyle nS=normStyle::one);

	AmpAbsDynamicalFunction(const char *name,
			std::shared_ptr<DoubleParameter> mag, std::shared_ptr<DoubleParameter> phase,
			std::shared_ptr<DoubleParameter> mass, int subsys, Spin spin, Spin m, Spin n,
			int nCalls=30000, normStyle nS=normStyle::one);

	virtual ~AmpAbsDynamicalFunction();

	//! Implementation of interface for streaming info about the strategy
	virtual const std::string to_str() const { return (_name); }
	//! value of resonance at \param point
	virtual std::complex<double> evaluate(dataPoint& point);
	//! value of dynamical amplitude at \param point
	virtual std::complex<double> evaluateAmp(dataPoint& point) = 0;
	//! value of angular distribution at \param point
	virtual double evaluateWignerD(dataPoint& point) {
		return _wignerD.evaluate(point);
	};
	//! Calculation integral |dynamical amplitude|^2
	virtual double integral();
	/** Calculation integral |c * dynamical amplitude * WignerD|^2
	 * Used to check the correct normalization of the amplitude. Should always be 1.
	 * @return
	 */
	virtual double totalIntegral() const;

	//! Get current normalization. In case that resonance parameters has change, it is recalculated.
	virtual double GetNormalization();
	//! Set normalization manually. Setting to values <0 disables normalization
	virtual void SetNormalization(double n){ _norm=n; };
	//! Set normalization style
	virtual void SetNormalizationStyle(normStyle n){ _normStyle=n; };
	//! Get normalization style
	virtual normStyle GetNormalizationStyle(){ return _normStyle; };
	//! Trigger recalculation of normalization
	virtual void SetModified() { modified=1; }

	//! Get resonance name
	virtual std::string GetName(){ return _name; }
	//! Set resonance name
	virtual void SetName(std::string n){ _name = n; }
	//! Get magnitude
	virtual double GetMagnitude() { return _mag->GetValue(); };
	//! Get phase
	virtual double GetPhase() { return _phase->GetValue(); };
	//! Get resonance mass
	double GetMass() {return _mass->GetValue();};
	//! Get resonance spin
	virtual double getSpin() { return _spin; }
	double GetM() {return _m;};
	double GetN() {return _n;};
	//! Get resonance meson radius
	double GetMesonRadius() {return _mesonRadius->GetValue();};
	//! Set resonance mother meson radius
	void SetMesonRadius(double r) {	_mesonRadius->SetValue(r); }
	//! Get resonance mother meson radius
	double GetMotherRadius() {return _motherRadius->GetValue();};
	//! Set resonance mother meson radius
	void SetMotherRadius(double r) {	_motherRadius->SetValue(r); }

	inline virtual bool isSubSys(const unsigned int subSys) const{ return (subSys==_subSys); };
	//	virtual std::shared_ptr<FunctionTree> setupTree(
	//			allMasses& theMasses,allMasses& toyPhspSample,std::string suffix,
	//			ParameterList params) { return std::shared_ptr<FunctionTree>(); };
	virtual std::shared_ptr<FunctionTree> setupTree(
			allMasses& theMasses,allMasses& toyPhspSample,std::string suffix,
			ParameterList& params) = 0;

	/** Convert width of resonance to coupling
	 *
	 * Implementation of Eq.47-21 of PDG2014. Only valid for narrow, isolated resonances.
	 * @param mSq invariant mass
	 * @param mR mass of resonance
	 * @param width width of resonance in channel [a,b]
	 * @param ma mass of particle a
	 * @param mb mass of particle b
	 * @return
	 */
	static std::complex<double> widthToCoupling(double mSq, double mR, double width, double ma, double mb, double spin, double mesonRadius);
	/** Convert coupling to width
	 *
	 * Convert coupling to channel (@ma,@mb) to partial width. Only valid for narrow, isolated resonances.
	 * Implementation of inverted Eq.47-21 of PDG2014.
	 * @param mSq invariant mass
	 * @param mR mass of resonance
	 * Implementation of inverted Eq.47-21 of PDG2014. Only valid for narrow, isolated resonances.
	 * @param mSq invariant mass
	 * @param mR mass of resonance
	 * @param g coupling to channel [a,b]
	 * @param ma mass of particle a
	 * @param mb mass of particle b
	 * @return
	 */
	static std::complex<double> couplingToWidth(double mSq, double mR, double g, double ma, double mb, double spin, double mesonRadius);

protected:
	void initialize();
	//! Name of resonance
	std::string _name;
	//! Precision of MC integration
	int _nCalls;
	normStyle _normStyle;
	std::shared_ptr<DoubleParameter> _mag;
	std::shared_ptr<DoubleParameter> _phase;
	double _ma, _mb, _mc, _M;
	//! Resonance mass
	std::shared_ptr<DoubleParameter> _mass;
	double tmp_mass;
	//! Resonance sub system
	unsigned int _subSys;
	//! Resonance spin
	Spin _spin, _m, _n;
	//! Barrier radi for resonance and mother particle
	std::shared_ptr<DoubleParameter> _mesonRadius, _motherRadius;

	AmpWigner2 _wignerD;
	double _norm;
	bool modified;

};

#endif
