// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef _FOURMOMENTUM_HPP_
#define _FOURMOMENTUM_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <boost/property_tree/ptree.hpp>

namespace ComPWA {

/// ComPWA four momentum class.
class FourMomentum {

public:
  FourMomentum() : FourMomentum(0.0, 0.0, 0.0, 0.0){};
  FourMomentum(double E) : FourMomentum(0.0, 0.0, 0.0, E){};

  FourMomentum(double Px, double Py, double Pz, double E)
      : FourMomentum(std::array<double, 4>{{Px, Py, Pz, E}}) {}

  // main constructor
  FourMomentum(std::array<double, 4> P4_) : P4(P4_) {}

  FourMomentum(std::vector<double> P4_)
      : FourMomentum([&P4_]() {
          if (P4_.size() != 4)
            throw std::runtime_error(
                "FourMomentum::Fourmomentum() | Size of vector not equal 4!");
          return std::array<double, 4>{{P4_[0], P4_[1], P4_[2], P4_[3]}};
        }()) {}

  double px() const { return P4[0]; }
  double py() const { return P4[1]; }
  double pz() const { return P4[2]; }
  double e() const { return P4[3]; }

  FourMomentum operator+(const FourMomentum &pB) const {
    FourMomentum newP(*this);
    newP += pB;
    return newP;
  }

  void operator+=(const FourMomentum &pB) {
    std::transform(P4.begin(), P4.end(), pB.P4.begin(), P4.begin(),
                   std::plus<double>());
  }

  bool operator==(const FourMomentum &pB) const { return P4 == pB.P4; }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const FourMomentum &p4) {
    stream << "(" << p4.px() << "," << p4.py() << "," << p4.e() << ")";
    return stream;
  }

  double invariantMassSquared() const {
    return ((-1) *
            (P4[0] * P4[0] + P4[1] * P4[1] + P4[2] * P4[2] - P4[3] * P4[3]));
  }

  double invariantMass() const { return std::sqrt(invariantMassSquared()); }

  double threeMomentumSquared() const {
    return (P4[0] * P4[0] + P4[1] * P4[1] + P4[2] * P4[2]);
  }

private:
  std::array<double, 4> P4;
};

} // namespace ComPWA
#endif
