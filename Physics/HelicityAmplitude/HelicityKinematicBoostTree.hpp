//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//--------------------------------------------------------------------------------

#ifndef HELICITYKINEMATICBOOSTTREE_HPP_
#define HELICITYKINEMATICBOOSTTREE_HPP_

#include "HelicityDecayTree.hpp"

#include "Core/DataPoint.hpp"

namespace HelicityFormalism {

class HelicityKinematicBoostTree {
  std::vector<Vector4<double> > boosted_4vectors_;

public:
  HelicityKinematicBoostTree();
  virtual ~HelicityKinematicBoostTree();

  void constructBoostingTree(const HelicityDecayTree& decay_tree,
      const dataPoint& event_data);

  const Vector4<double>& getBoosted4Vector(unsigned int node_index) const;
};

} /* namespace HelicityFormalism */

#endif /* HELICITYKINEMATICBOOSTTREE_HPP_ */
