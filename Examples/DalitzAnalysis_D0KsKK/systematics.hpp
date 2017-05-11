/*
 * systematics.hpp
 *
 *  Created on: Aug 28, 2015
 *      Author: weidenka
 */

#ifndef PWA_SYSTEMATICS_HPP_
#define PWA_SYSTEMATICS_HPP_

using namespace ComPWA;

MomentumCorrection* getTrackingCorrection(){

	//momentum depended efficiency difference for charged tracks
	CorrectionTable chargedTrackingSys("Tracking systematics charged tracks");
	//momentum depended efficiency difference KS candidates
	CorrectionTable ksSystematics("K_S^0 reconstruction systematics");

	//KS systematics (MA Tians study)
	//(min, max, sys, sysError, antiSys, antiSysError)
	ksSystematics.Add(0.00, 0.18, 0.0645,  0.071,  0.0715, 0.0745);
	ksSystematics.Add(0.18, 0.24, 0.0498,  0.0625, 0.0372, 0.0621);
	ksSystematics.Add(0.24, 0.32, 0.0131,  0.0439, 0.0195, 0.0444);
	ksSystematics.Add(0.32, 0.39, -0.0143, 0.0283, 0.0470, 0.0297);
	ksSystematics.Add(0.39, 0.44, -0.0224, 0.0294, -0.0111,0.0300);
	ksSystematics.Add(0.44, 0.50, 0.0303,  0.0292, 0.0216, 0.0292);
	ksSystematics.Add(0.50, 0.55, 0.0182,  0.0299, -0.0155,0.0295);
	ksSystematics.Add(0.55, 0.60, 0.0203,  0.0297, -0.0172,0.0297);
	ksSystematics.Add(0.60, 0.67, 0.0445,  0.0199, -0.0055,0.0196);
	ksSystematics.Add(0.67, 0.74, -0.0119, 0.0177, 0.0205, 0.0186);
	ksSystematics.Add(0.74, 0.81, 0.0214,  0.0181, 0.0237, 0.0185);
	ksSystematics.Add(0.81, 0.87, 0.0213,  0.0188, -0.0002,0.0190);
	ksSystematics.Add(0.87, 0.92, -0.0037, 0.0198, 0.0163, 0.0206);
	ksSystematics.Add(0.92, 0.98, 0.0127,  0.0185, -0.0068,0.0184);
	ksSystematics.Add(0.98, 1.04, 0.0262,  0.0195, 0.0175, 0.0191);
	ksSystematics.Add(1.04, 1.10, 0.0041,  0.0196, 0.001,  0.0198);
	ksSystematics.Add(1.10, 1.17, 0.0518,  0.0204, -0.0045,0.0193);
	ksSystematics.Add(1.17, 1.26, -0.0164, 0.0191, 0.0123, 0.0198);
	ksSystematics.Add(1.26, 1.40, -0.0015, 0.0198, 0.0109, 0.0204);

	//charged kaon track systematics ( K/pi l nu memo )
	//(min, max, sys, sysError)
	chargedTrackingSys.Add(0.0, 0.1, 0.3600, sqrt(0.13*0.13*0.04*0.04) );
	chargedTrackingSys.Add(0.1, 0.2,-0.0070, sqrt(0.0017*0.0017+0.0008*0.0008) );//very small error shifts the result to low values
	chargedTrackingSys.Add(0.2, 0.3, 0.0056, sqrt(0.0063*0.0063+0.0027*0.0027) );
	chargedTrackingSys.Add(0.3, 0.4, 0.0100, sqrt(0.0040*0.0040+0.0029*0.0029) );
	chargedTrackingSys.Add(0.4, 0.5, 0.0177, sqrt(0.0032*0.0032+0.0036*0.0036) );
	chargedTrackingSys.Add(0.5, 0.6, 0.0121, sqrt(0.0033*0.0033+0.0020*0.0020) );
	chargedTrackingSys.Add(0.6, 0.7, 0.0089, sqrt(0.0040*0.0040+0.0059*0.0059) );
	chargedTrackingSys.Add(0.7, 0.8, 0.0050, sqrt(0.0025*0.0025+0.0022*0.0022) );
	chargedTrackingSys.Add(0.8, 0.9, 0.0053, sqrt(0.0023*0.0023+0.0007*0.0007) );
	chargedTrackingSys.Add(0.9, 1.0, 0.0009, sqrt(0.0021*0.0021+0.0014*0.0014) );
	chargedTrackingSys.Add(1.0, 1.1,-0.0079, sqrt(0.0050*0.0050+0.0050*0.0050) );
	std::vector<CorrectionTable> vecTrkSys;
//	vecTrkSys.push_back(ksTrkSys);
//	vecTrkSys.push_back(kaonTrkSys);
//	vecTrkSys.push_back(kaonTrkSys);
	vecTrkSys.push_back(ksSystematics);
	vecTrkSys.push_back(chargedTrackingSys);
	vecTrkSys.push_back(chargedTrackingSys);
	MomentumCorrection* trkSys = new MomentumCorrection(vecTrkSys,"Tracking corrections");
	return trkSys;
}
MomentumCorrection* getPidCorrection(){
	//momentum depended efficiency difference for charged tracks
	CorrectionTable chargedPidSys("PID systematics charged tracks");
	//charged kaon PID systematics ( K/pi l nu memo )
	//(min, max, sys, sysError)
	//chargedPidSys.Add(0.0, 0.1, 0.0, 1.0 );
	chargedPidSys.Add(0.1, 0.2,-0.0263, 0.0056 );
	chargedPidSys.Add(0.2, 0.3,-0.0038, 0.0014 );
	chargedPidSys.Add(0.3, 0.4,-0.0021, 0.0006 );
	chargedPidSys.Add(0.4, 0.5,-0.0017, 0.0005 );
	chargedPidSys.Add(0.5, 0.6,-0.0020, 0.0007 );
	chargedPidSys.Add(0.6, 0.7, 0.0008, 0.0010 );
	chargedPidSys.Add(0.7, 0.8, 0.0046, 0.0010 );
	chargedPidSys.Add(0.8, 0.9, 0.0082, 0.0014 );
	chargedPidSys.Add(0.9, 1.0, 0.0028, 0.0017 );
	chargedPidSys.Add(1.0, 1.1, 0.0071, 0.0036 );
	//charged kaon PID systematics ( ZhilingPID )
	//(min, max, sys, sysError)
	//	chargedPidSys.Add(0.0, 0.1, 0.0204, 1.0 );
	//	chargedPidSys.Add(0.1, 0.2, 0.0045, 1.0 );
	//	chargedPidSys.Add(0.2, 0.3, 0.0060, 1.0 );
	//	chargedPidSys.Add(0.3, 0.4, 0.0083, 1.0 );
	//	chargedPidSys.Add(0.4, 0.5, 0.0053, 1.0 );
	//	chargedPidSys.Add(0.5, 0.6, 0.0023, 1.0 );
	//	chargedPidSys.Add(0.6, 0.7, 0.0022, 1.0 );
	//	chargedPidSys.Add(0.7, 0.8, 0.0021, 1.0 );
	//	chargedPidSys.Add(0.8, 0.9, 0.0005, 1.0 );
	//	chargedPidSys.Add(0.9, 1.0, 0.0027, 1.0 );
	std::vector<CorrectionTable> vecPidSys;
	vecPidSys.push_back(CorrectionTable("No corrections for K_S^0"));
	vecPidSys.push_back(chargedPidSys);
	vecPidSys.push_back(chargedPidSys);
	MomentumCorrection* pidSys = new MomentumCorrection(vecPidSys,"PID corrections");
	return pidSys;
}

#endif /* PWA_SYSTEMATICS_HPP_ */
