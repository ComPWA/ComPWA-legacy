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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "Core/Efficiency.hpp"
#include "Physics/HelicityFormalism/CoherentIntensity.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

double CoherentIntensity::Intensity(const dataPoint &point) const {
  return IntensityNoEff(point) * _eff->Evaluate(point);
}

double CoherentIntensity::IntensityNoEff(const dataPoint &point) const {
  std::complex<double> result(0., 0.);
  for (auto i : _seqDecays)
    result += i->Evaluate(point);
  return GetStrengthValue()*std::norm(result);
};

std::shared_ptr<CoherentIntensity>
CoherentIntensity::Factory(const boost::property_tree::ptree &pt) {
  LOG(trace) << " CoherentIntensity::Factory() | Construction....";
  auto obj = std::make_shared<CoherentIntensity>();
  obj->SetName(pt.get<std::string>("<xmlattr>.Name"));

//  boost::property_tree::xml_writer_settings<char> settings('\t', 1);
  write_xml(std::cout,pt);
  
  auto ptCh = pt.get_child_optional("Strength");
  if (ptCh) {
    auto strength = ComPWA::DoubleParameterFactory(ptCh.get());
    obj->SetStrength(std::make_shared<DoubleParameter>(strength));
  } else {
    obj->SetStrength(std::make_shared<ComPWA::DoubleParameter>("", 1.0));
  }

  for (const auto &v : pt.get_child("")) {
    if( v.first == "Amplitude" )
    obj->Add(
        ComPWA::Physics::HelicityFormalism::SequentialTwoBodyDecay::Factory(
            v.second));
  }
  return obj;
}

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */
