#include <boost/archive/xml_iarchive.hpp>

#include "Core/Properties.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/ParticleList.hpp"

#include "Physics/HelicityFormalism/CreateHelicityModel.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

double intensityMaximum(std::shared_ptr<IIntensity> amp, int size) {
  return 1.0;
}

std::pair<std::shared_ptr<ComPWA::IIntensity>, std::shared_ptr<IKinematics>>
createHelicityModel(const char* modelXMLFile, int seed, std::vector<int> initialState,
                std::vector<int> finalState, const char* particleXMLFile) {
  std::string modelStr;
  if (modelXMLFile) {
    modelStr = std::string(modelXMLFile);
  } else {
    //TODO: error
  }

  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, ComPWA::Physics::defaultParticleList);

  if (particleXMLFile) {
    std::string particleStr(particleXMLFile);
    boost::property_tree::ptree particles;
    boost::property_tree::xml_parser::read_xml(particleXMLFile, particles);
    ReadParticles(partL, particles);
  }

  auto kin =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, initialState, finalState);

  boost::property_tree::ptree model;
  boost::property_tree::xml_parser::read_xml(modelXMLFile, model);

  auto intens = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
      partL, kin, model.get_child("Intensity"));
  return std::make_pair(intens, kin);
}

}  // namespace HelicityFormalism
}  // namespace Physics
}  // namespace ComPWA
