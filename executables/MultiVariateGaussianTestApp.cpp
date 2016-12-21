/*
 * MulitVariateGaussianTest.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: weidenka
 */

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>

#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Optimizer/Minuit2/MinuitResult.hpp"
#include <TH2D.h>

using namespace ComPWA;

int main(){
	size_t n = 3;
	const size_t nRnd = 5000;

	gsl_vector* v = gsl_vector_alloc(n);
	for(size_t i=0; i<n; i++)
		gsl_vector_set(v,i,0);
	gsl_matrix* m = gsl_matrix_alloc(n,n);
	/* We assume three variables x,y,z with sigma^2 = 0.01, 0.04, 0.09
	 * and correlations xy = 0.1, xz = 0.2 and yz = 0.3. The corresponding
	 * covariance matrix is then given by:
	 * 0.01	0.002	0.006
	 * 		0.04	0.018
	 * 				0.09
	 */
	double sigmaX=0.1, sigmaY=0.2, sigmaZ=0.3;
	double corrXY=0.4, corrXZ=0.5, corrYZ=0.6;
	gsl_matrix_set(m,0,0,sigmaX*sigmaX);
	gsl_matrix_set(m,0,1,sigmaX*sigmaY*corrXY);
	gsl_matrix_set(m,1,0,sigmaX*sigmaY*corrXY);
	gsl_matrix_set(m,0,2,sigmaX*sigmaZ*corrXZ);
	gsl_matrix_set(m,2,0,sigmaX*sigmaZ*corrXZ);
	gsl_matrix_set(m,1,1,sigmaY*sigmaY);
	gsl_matrix_set(m,1,2,sigmaY*sigmaZ*corrYZ);
	gsl_matrix_set(m,2,1,sigmaY*sigmaZ*corrYZ);
	gsl_matrix_set(m,2,2,sigmaZ*sigmaZ);
	std::cout<<"COVARIANCE MATRIX:"<<std::endl;
	Optimizer::Minuit2::gsl_matrix_print(m);

	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	gsl_rng* rnd = gsl_rng_alloc (T);

	TH2D* xy = new TH2D("xy","xy",100,-2,2,100,-2,2);
	TH2D* xz = new TH2D("xz","xz",100,-2,2,100,-2,2);
	TH2D* yz = new TH2D("yz","yz",100,-2,2,100,-2,2);
	double x[nRnd];
	double y[nRnd];
	double z[nRnd];
	std::ofstream myfile;
	myfile.open("sample.txt");
	myfile <<"x/F:y:z"<<std::endl;
	for(size_t i=0; i<nRnd; i++){
		gsl_vector* tmpV = gsl_vector_alloc(n);
		Optimizer::Minuit2::multivariateGaussian(rnd, n, v, m, tmpV);
		x[i] = gsl_vector_get(tmpV,0);
		y[i] = gsl_vector_get(tmpV,1);
		z[i] = gsl_vector_get(tmpV,2);
		xy->Fill(x[i],y[i]);
		xz->Fill(x[i],z[i]);
		yz->Fill(y[i],z[i]);
		myfile<<x[i]<<" "<<y[i]<<" "<<z[i]<<std::endl;
		gsl_vector_free(tmpV);
	}
	//gsl function does not give correct results
//	std::cout<<"CORRELATIONS:"<<std::endl;
//	std::cout<<" x = y "<<gsl_stats_correlation(x,1,y,1,n)<<std::endl;
//	std::cout<<" x = z "<<gsl_stats_correlation(x,1,z,1,n)<<std::endl;
//	std::cout<<" y = z "<<gsl_stats_correlation(y,1,z,1,n)<<std::endl;
	std::cout<<"CORRELATIONS(true):"<<std::endl;
	std::cout<<" x = y "<<corrXY<<std::endl;
	std::cout<<" x = z "<<corrXZ<<std::endl;
	std::cout<<" y = z "<<corrYZ<<std::endl;
	std::cout<<"CORRELATIONS:"<<std::endl;
	std::cout<<" x = y "<<xy->GetCorrelationFactor()<<std::endl;
	std::cout<<" x = z "<<xz->GetCorrelationFactor()<<std::endl;
	std::cout<<" y = z "<<yz->GetCorrelationFactor()<<std::endl;
	myfile.close();
	delete xy, xz, yz;
	return 0;
}




