// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PWA_SYSTEMATICS_HPP_
#define PWA_SYSTEMATICS_HPP_

namespace ComPWA {
namespace Data {

/// Systematic uncertainties for the momentum depended tracking efficiency.
/// External input from a collaboration and usually not public available. We use
/// Sample values which are of the same order as the values of the BESIII
/// experiment.
MomentumCorrection getTrackingCorrection(const ParticleList &list) {

  // momentum depended efficiency difference KS candidates
  CorrectionTable ksSystematics("K_S^0 reconstruction systematics");

  //(min momentum, max momentum, sys, sysError (particle)
  // sys, sysError (anti-particle)
  ksSystematics.Add(0.00, 0.18, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.18, 0.24, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.24, 0.32, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.32, 0.39, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.39, 0.44, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.44, 0.50, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.50, 0.55, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.55, 0.60, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.60, 0.67, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.67, 0.74, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.74, 0.81, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.81, 0.87, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.87, 0.92, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.92, 0.98, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(0.98, 1.04, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(1.04, 1.10, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(1.10, 1.17, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(1.17, 1.26, 0.01, 0.001, 0.01, 0.001);
  ksSystematics.Add(1.26, 1.40, 0.01, 0.001, 0.01, 0.001);

  // Momentum depended efficiency difference for charged tracks
  CorrectionTable chargedTrackingSys("Tracking systematics charged tracks");
  //(min momentum, max momentum, sys, sysError)
  chargedTrackingSys.Add(0.0, 0.1, 0.01, 0.001);
  chargedTrackingSys.Add(0.1, 0.2, 0.01, 0.001);
  chargedTrackingSys.Add(0.2, 0.3, 0.01, 0.001);
  chargedTrackingSys.Add(0.3, 0.4, 0.01, 0.001);
  chargedTrackingSys.Add(0.4, 0.5, 0.01, 0.001);
  chargedTrackingSys.Add(0.5, 0.6, 0.01, 0.001);
  chargedTrackingSys.Add(0.6, 0.7, 0.01, 0.001);
  chargedTrackingSys.Add(0.7, 0.8, 0.01, 0.001);
  chargedTrackingSys.Add(0.8, 0.9, 0.01, 0.001);
  chargedTrackingSys.Add(0.9, 1.0, 0.01, 0.001);
  chargedTrackingSys.Add(1.0, 1.1, 0.01, 0.001);

  std::vector<CorrectionTable> vecTrkSys;
  vecTrkSys.push_back(ksSystematics);
  vecTrkSys.push_back(chargedTrackingSys);
  vecTrkSys.push_back(chargedTrackingSys);
  MomentumCorrection trkSys(list, vecTrkSys, "Tracking corrections");
  return trkSys;
}
/// Systematic uncertainties for the momentum depended pid efficiency.
/// External input from a collaboration and usually not public available. We use
/// Sample values which are of the same order as the values of the BESIII
/// experiment.
MomentumCorrection getPidCorrection(const ParticleList &list) {
  // momentum depended efficiency difference for charged tracks
  CorrectionTable chargedPidSys("PID systematics charged tracks");

  //(min, max, sys, sysError)
  // chargedPidSys.Add(0.0, 0.1, 0.0, 1.0 );
  chargedPidSys.Add(0.1, 0.2, 0.01, 0.001);
  chargedPidSys.Add(0.2, 0.3, 0.01, 0.001);
  chargedPidSys.Add(0.3, 0.4, 0.01, 0.001);
  chargedPidSys.Add(0.4, 0.5, 0.01, 0.001);
  chargedPidSys.Add(0.5, 0.6, 0.01, 0.001);
  chargedPidSys.Add(0.6, 0.7, 0.01, 0.001);
  chargedPidSys.Add(0.7, 0.8, 0.01, 0.001);
  chargedPidSys.Add(0.8, 0.9, 0.01, 0.001);
  chargedPidSys.Add(0.9, 1.0, 0.01, 0.001);
  chargedPidSys.Add(1.0, 1.1, 0.01, 0.001);

  std::vector<CorrectionTable> vecPidSys;
  vecPidSys.push_back(CorrectionTable("No corrections for K_S^0"));
  vecPidSys.push_back(chargedPidSys);
  vecPidSys.push_back(chargedPidSys);
  MomentumCorrection pidSys(list, vecPidSys, "PID corrections");
  return pidSys;
}
} // namespace Data
} // namespace ComPWA

#endif
