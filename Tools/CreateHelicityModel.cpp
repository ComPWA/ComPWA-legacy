#include <boost/archive/xml_iarchive.hpp>

#include "Core/Logging.hpp"
#include "Core/Properties.hpp"
#include "Physics/HelicityFormalism/HelicityKinematics.hpp"
#include "Physics/IncoherentIntensity.hpp"
#include "Physics/ParticleList.hpp"

#include "Tools/CreateHelicityModel.hpp"
#include "Tools/RootGenerator.hpp"
#include "Tools/Generate.hpp"

namespace ComPWA {
namespace Tools {

double intensityMaximum(std::shared_ptr<IIntensity> amp, int size) {
  return 1.0;
}

std::pair<std::shared_ptr<ComPWA::IIntensity>, std::shared_ptr<IKinematics>>
createHelicityModel(const char* modelXMLFile, int seed, int mcPrecision, const char* logLv) {
  Logging log(std::string(logLv));

  //---------------------------------------------------
  // Create particle list
  //
  // Initially filled with a default list, optionally
  // entry from the XML model file are added.
  //---------------------------------------------------
  auto partL = std::make_shared<ComPWA::PartList>();
  ReadParticles(partL, ComPWA::Physics::defaultParticleList);

  if (!modelXMLFile) {
    throw std::runtime_error(
        "createHelicityModel() | No model file is provided.");
  }

  // Optinally add particles from the XML model file. If a
  // particle exists it gets overwritten by the ones defined in
  // the model file
  std::string modelStr = std::string(modelXMLFile);
  boost::property_tree::ptree modelTree;
  boost::property_tree::xml_parser::read_xml(modelStr, modelTree);
  ReadParticles(partL, modelTree);

  //---------------------------------------------------
  // Create Kinematics object
  //---------------------------------------------------
  auto kin =
      std::make_shared<ComPWA::Physics::HelicityFormalism::HelicityKinematics>(
          partL, modelTree.get_child("HelicityKinematics"));

  //---------------------------------------------------
  // 2) Generate a large phase space sample
  //---------------------------------------------------
  auto gen = std::make_shared<ComPWA::Tools::RootGenerator>(partL, kin, seed);
  std::shared_ptr<ComPWA::Data::Data> phspSample(
      ComPWA::Tools::generatePhsp(mcPrecision, gen));

  boost::property_tree::ptree model;
  boost::property_tree::xml_parser::read_xml(modelXMLFile, model);

  auto intens = std::make_shared<ComPWA::Physics::IncoherentIntensity>(
      partL, kin, model.get_child("Intensity"));

  // Pass phsp sample to intensity for normalization.
  // Convert to dataPoints first.
  auto phspPoints =
      std::make_shared<std::vector<DataPoint>>(phspSample->dataPoints(kin));
  intens->setPhspSample(phspPoints, phspPoints);

  return std::make_pair(intens, kin);
}

}  // namespace Tools
}  // namespace ComPWA
