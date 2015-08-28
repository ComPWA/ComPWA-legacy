/*
 * CorrectionTabletematics.hpp
 *
 *  Created on: Feb 27, 2015
 *      Author: weidenka
 */

#ifndef MOMENTUMSYSTEMATICS_HPP_
#define MOMENTUMSYSTEMATICS_HPP_
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>

class CorrectionTable
{
public:
	CorrectionTable(std::string t="") : title(t){ }

	~CorrectionTable() { }

	/** Get Data/MC difference
	 * Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0 averages over both values.
	 * @param charge Charge of particle
	 * @param momentum Particle momentum
	 * @return
	 */
	double GetValue(int charge, double momentum) const;
	/** Get error of Data/MC difference
	 * Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0 averages over both values.
	 * @param charge Charge of particle
	 * @param momentum Particle momentum
	 * @return
	 */
	double GetError(int charge, double momentum) const;
	/** Get MC/Data difference
	 * Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0 averages over both values.
	 * @param charge Charge of particle
	 * @param momentum Particle momentum
	 * @return
	 */
	double GetInvValue(int charge, double momentum) const;
	/** Get error of MC/Data difference
	 * Choose a charge for particle(=1) or anti-particle(=-1). A charge of 0 averages over both values.
	 * @param charge Charge of particle
	 * @param momentum Particle momentum
	 * @return
	 */
	double GetInvError(int charge, double momentum) const;

	/** Add momentum bin with efficiency correction
	 * The efficiency is given as (epsilon_data/epsilon_mc - 1)
	 *
	 * @param binMin lower bin boundary
	 * @param binMax higher bin boundary
	 * @param s efficiency difference for particles
	 * @param sError uncertainty of efficiency difference for particles
	 * @param antiS efficiency difference for anti-particles
	 * @param antiSerror uncertainty of efficiency difference for anti-particles
	 */
	void Add(double binMin, double binMax, double s, double sError,
			double antiS=-999, double antiSerror=-999);

	/** Add momentum bin with efficiency correction
	 * The efficiency is given as (epsilon_mc/epsilon_data - 1)
	 *
	 * @param binMin lower bin boundary
	 * @param binMax higher bin boundary
	 * @param s efficiency difference for particles
	 * @param sError uncertainty of efficiency difference for particles
	 * @param antiS efficiency difference for anti-particles
	 * @param antiSerror uncertainty of efficiency difference for anti-particles
	 */
	void AddInv(double binMin, double binMax, double s, double sError,
			double antiS=-999, double antiSerror=-999);
	//! Get binning vector
	std::vector<std::pair<double, double> > GetBins(){ return bins; }
	//! Set binning vector
	void SetBins(std::vector<std::pair<double, double> > b){ bins = b; }
	//! Get number of bins
	int GetNbins() {return bins.size(); }
	//! Get systematics vector
	std::vector<double> GetSystematics(){ return sys; }
	//! Get systematics error vector
	std::vector<double> GetSystematicsError(){ return sysError; }
	//! Set systematics vector
	void SetSystematics(std::vector<double> b, std::vector<double> bError=std::vector<double>() );
	//! Set systematics error vector
	void SetSystematicsError(std::vector<double> bError);
	//! Get anti particle systematics vector
	std::vector<double> GetAntiSystematics(){ return antiSys; }
	//! Get anti particle systematics error vector
	std::vector<double> GetAntiSystematicsError(){ return antiSysError; }
	//! Set anti particle systematics vector
	void SetAntiSystematics(std::vector<double> b, std::vector<double> bError=std::vector<double>() );
	void SetAntiSystematicsError(std::vector<double> bError);
	/** Count total systematics internally
	 * Add systematics of track to total systematic error
	 * @param charge
	 * @param momentum
	 */
	void AddToTotalError(int charge, double momentum);

	/** Calculated total systematic uncertainty
	 * The weighted mean of the uncertainties of the single tracks is calculated
	 *
	 * @return total systematic uncertainty
	 */
	double GetTotalSystematics(bool useArithmeticMean=0);

	/** Calculated error of total systematic uncertainty
	 * The weighted error of the uncertainties of the single tracks is calculated
	 *
	 * @return error of total systematic uncertainty
	 */
	double GetTotalSystematicsError(bool useArithmeticMean=0);

	//! Print table
	void Print() const;
	//! Get title
	std::string GetTitle() { return title; }
	//! Set title
	void SetTitle( std::string t ) { title=t; }

protected:
	//! Title
	std::string title;
	//! invert e_mc/e_data-1 to e_data/e_mc-1
	static double inverse(double x);
	//! Calculate error for inversion e_mc/e_data-1 to e_data/e_mc-1
	static double inverseError(double x, double xErr);
	//! find bin position
	int findBin(double momentum) const;
	//! check if there is a bin defined that overlaps with [min,max]
	int findBin(double min, double max) const;
	//! add bin
	void addBin(double min, double max);
	//! check for consistency
	bool check() const;
	std::vector<std::pair<double,double> > bins;
	//! Data/MC difference in momentum bins for particle
	std::vector<double> sys;
	std::vector<double> sysError;
	//! Data/MC difference in momentum bins for anti-particle
	std::vector<double> antiSys;
	std::vector<double> antiSysError;

	std::vector<double> totalSys, totalSysError;
};
#endif /* MOMENTUMSYSTEMATICS_HPP_ */
