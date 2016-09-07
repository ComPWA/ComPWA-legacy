/*
 * AdvancedStrategies.hpp
 *
 *  Created on: 6 Sep 2016
 *      Author: weidenka
 */

#ifndef PHYSICS_AMPLITUDESUM_ADVANCEDSTRATEGIES_HPP_
#define PHYSICS_AMPLITUDESUM_ADVANCEDSTRATEGIES_HPP_

class phspFactorStrat : public Strategy
{
public:
	phspFactorStrat() : Strategy(ParType::MCOMPLEX) { }

	virtual const std::string to_str() const { return ("phspFactor"); }

	static std::shared_ptr<FunctionTree> SetupTree(
			std::shared_ptr<MultiDouble> mSq, double ma, double mb){

		std::shared_ptr<phspFactorStrat> thisStrat(
				new phspFactorStrat);

		std::string stratName = "phspFactorStrat";
		//------------Setup Tree---------------------
		std::shared_ptr<FunctionTree> newTree(new FunctionTree());

		newTree->createHead(stratName, thisStrat, mSq->GetNValues());
		newTree->createLeaf("massA", ma, stratName);
		newTree->createLeaf("massB", mb, stratName);
		newTree->createLeaf("mSq", mSq ,stratName);

		return newTree;

	}

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out){

		//Check parameter type
		if( checkType != out->type() )
			throw( WrongParType("phspFactorStrat::execute() | "
					"Output parameter is of type "
					+ std::string(ParNames[out->type()])
		+ " and conflicts with expected type "
		+ std::string(ParNames[checkType]) )
			);

		//How many parameters do we expect?
		int check_nBool = 0;
		int check_nInt = 0;
		int check_nComplex = 0;
		int check_nDouble = 2;
		int check_nMDouble = 1;
		int check_nMComplex = 0;

		//Check size of parameter list
		if( paras.GetNBool() != check_nBool )
			throw( BadParameter("phspFactorStrat::execute() | "
					"Number of BoolParameters does not match: "
					+std::to_string(paras.GetNBool())+" given but "
					+std::to_string(check_nBool)+ " expected.")
			);
		if( paras.GetNInteger() != check_nInt )
			throw( BadParameter("phspFactorStrat::execute() | "
					"Number of IntParameters does not match: "
					+std::to_string(paras.GetNInteger())+" given but "
					+std::to_string(check_nInt)+ " expected.")
			);
		if( paras.GetNDouble() != check_nDouble )
			throw( BadParameter("phspFactorStrat::execute() | "
					"Number of DoubleParameters does not match: "
					+std::to_string(paras.GetNDouble())+" given but "
					+std::to_string(check_nDouble)+ " expected.")
			);
		if( paras.GetNComplex() != check_nComplex )
			throw( BadParameter("phspFactorStrat::execute() | "
					"Number of ComplexParameters does not match: "
					+std::to_string(paras.GetNComplex())+" given but "
					+std::to_string(check_nComplex)+ " expected.")
			);
		if( paras.GetNMultiDouble() != check_nMDouble )
			throw( BadParameter("phspFactorStrat::execute() | "
					"Number of MultiDoubles does not match: "
					+std::to_string(paras.GetNMultiDouble())+" given but "
					+std::to_string(check_nMDouble)+ " expected.")
			);
		if( paras.GetNMultiComplex() != check_nMComplex )
			throw( BadParameter("phspFactorStrat::execute() | "
					"Number of MultiComplexes does not match: "
					+std::to_string(paras.GetNMultiComplex())+" given but "
					+std::to_string(check_nMComplex)+ " expected.")
			);

		/** Get parameters from ParameterList:
		 * We use the same order of the parameters as was used during tree
		 * construction */
		double ma = paras.GetDoubleParameter(0)->GetValue();
		double mb = paras.GetDoubleParameter(1)->GetValue();

		std::shared_ptr<MultiDouble> mp = paras.GetMultiDouble(0);
		std::vector<std::complex<double> > results(mp->GetNValues(),
				std::complex<double>(0.));
		//calc function for each point
		for(unsigned int ele=0; ele<mp->GetNValues(); ele++){
			try{
				results.at(ele) = Kinematics::phspFactor(
						mp->GetValue(ele), ma, mb );
			} catch (std::exception& ex) {
				BOOST_LOG_TRIVIAL(error) << "phspFactorStrat::execute() | "
						<<ex.what();
				std::cout<<mp->GetValue(ele)<<" "<<ma<<" "<<mb<<std::endl;
				throw( std::runtime_error("phspFactorStrat::execute() | "
						"Evaluation of dynamic function failed!")
				);
			}
		}
		out = std::shared_ptr<AbsParameter>(
				new MultiComplex(out->GetName(),results));
		return true;
	}
};

class barrierStrat : public Strategy
{
public:
	barrierStrat() : Strategy(ParType::MDOUBLE) { }

	virtual const std::string to_str() const { return ("barrierFactor"); }

	static std::shared_ptr<FunctionTree> SetupTree(
			std::shared_ptr<MultiDouble> mSq,
			std::shared_ptr<DoubleParameter> mR, double ma, double mb,
			Spin spin, std::shared_ptr<DoubleParameter> mesonRadius,
			formFactorType type){

		std::shared_ptr<barrierStrat> thisStrat( new barrierStrat );

		std::string stratName = "barrierFactor";
		//------------Setup Tree---------------------
		std::shared_ptr<FunctionTree> newTree(new FunctionTree());

		newTree->createHead(stratName, thisStrat, mSq->GetNValues());
		newTree->createLeaf("mass", mR, stratName);
		newTree->createLeaf("massA1", ma, stratName);
		newTree->createLeaf("massA2", mb, stratName);
		newTree->createLeaf("spin", spin, stratName);
		newTree->createLeaf("mesonRadius", mesonRadius, stratName);
		newTree->createLeaf("formFactorType", type, stratName);
		newTree->createLeaf("mSq", mSq ,stratName);

		return newTree;

	}

	virtual bool execute(ParameterList& paras,
			std::shared_ptr<AbsParameter>& out){

		//Check parameter type
		if( checkType != out->type() )
			throw( WrongParType("barrierStrat::execute() | "
					"Output parameter is of type "
					+ std::string(ParNames[out->type()])
		+ " and conflicts with expected type "
		+ std::string(ParNames[checkType]) )
			);

		//How many parameters do we expect?
		int check_nBool = 0;
		int check_nInt = 0;
		int check_nComplex = 0;
		int check_nDouble = 6;
		int check_nMDouble = 1;
		int check_nMComplex = 0;

		//Check size of parameter list
		if( paras.GetNBool() != check_nBool )
			throw( BadParameter("barrierStrat::execute() | "
					"Number of BoolParameters does not match: "
					+std::to_string(paras.GetNBool())+" given but "
					+std::to_string(check_nBool)+ " expected.")
			);
		if( paras.GetNInteger() != check_nInt )
			throw( BadParameter("barrierStrat::execute() | "
					"Number of IntParameters does not match: "
					+std::to_string(paras.GetNInteger())+" given but "
					+std::to_string(check_nInt)+ " expected.")
			);
		if( paras.GetNDouble() != check_nDouble )
			throw( BadParameter("barrierStrat::execute() | "
					"Number of DoubleParameters does not match: "
					+std::to_string(paras.GetNDouble())+" given but "
					+std::to_string(check_nDouble)+ " expected.")
			);
		if( paras.GetNComplex() != check_nComplex )
			throw( BadParameter("barrierStrat::execute() | "
					"Number of ComplexParameters does not match: "
					+std::to_string(paras.GetNComplex())+" given but "
					+std::to_string(check_nComplex)+ " expected.")
			);
		if( paras.GetNMultiDouble() != check_nMDouble )
			throw( BadParameter("barrierStrat::execute() | "
					"Number of MultiDoubles does not match: "
					+std::to_string(paras.GetNMultiDouble())+" given but "
					+std::to_string(check_nMDouble)+ " expected.")
			);
		if( paras.GetNMultiComplex() != check_nMComplex )
			throw( BadParameter("barrierStrat::execute() | "
					"Number of MultiComplexes does not match: "
					+std::to_string(paras.GetNMultiComplex())+" given but "
					+std::to_string(check_nMComplex)+ " expected.")
			);

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

		std::shared_ptr<MultiDouble> mp = paras.GetMultiDouble(0);
		std::vector<double> results(mp->GetNValues(),0.);

		//calc function for each point
		for(unsigned int ele=0; ele<mp->GetNValues(); ele++){
			double s = mp->GetValue(ele);
			double sqrtS = sqrt(s);
			std::complex<double> qValue = Kinematics::qValue(sqrtS, ma, mb);
			try{
				double nom = Kinematics::FormFactor(
						sqrtS, ma, mb, spin, mesonRadius,
						qValue,	type);
				double denom = Kinematics::FormFactor(
						mR, ma, mb, spin, mesonRadius,
						qValue,	type);

				results.at(ele) = nom/denom*nom/denom;
			} catch (std::exception& ex) {
				BOOST_LOG_TRIVIAL(error) << "barrierStrat::execute() | "
						<<ex.what();
				throw( std::runtime_error("barrierStrat::execute() | "
						"Evaluation of dynamic function failed!")
				);
			}
		}
		out = std::shared_ptr<AbsParameter>(
				new MultiDouble(out->GetName(),results));
		return true;
	}
};

#endif /* PHYSICS_AMPLITUDESUM_ADVANCEDSTRATEGIES_HPP_ */
