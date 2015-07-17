//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//----------------------------------------------------------------------------------

#ifndef HELICITYAMPLITUDETREE_HPP_
#define HELICITYAMPLITUDETREE_HPP_

#include "HelicityAmplitude.hpp"

#include <vector>

namespace HelicityFormalism {

typedef std::vector<std::shared_ptr<HelicityAmplitude> > SequentialDecayAmplitude;

class HelicityAmplitudeTree {
  friend class HelicityAmplitudeTreeFactory;

  std::vector<SequentialDecayAmplitude> sequential_decay_amplitude_list_;
 // std::vector<IndexPair> links_;

  /*struct HelicityComparator {
    std::shared_ptr<HelicityAmplitude> reference_;

    HelicityComparator(std::shared_ptr<HelicityAmplitude> reference) :
        reference_(reference) {
    }

    bool linkForAmplitude(std::shared_ptr<HelicityAmplitude> test_element) {
      return test_element.get() == reference_.get();
    }
  };*/

 // std::complex<double> evaluateIntermediateNode(unsigned int node_index,
 //     const std::vector<HelicityAngles>& helicity_angles) const;

public:
  HelicityAmplitudeTree();
  virtual ~HelicityAmplitudeTree();

  std::complex<double> evaluate(
      const std::vector<HelicityAngles>& helicity_angles) const;
};

} /* namespace HelicityFormalism */

#endif /* HELICITYAMPLITUDETREE_HPP_ */
