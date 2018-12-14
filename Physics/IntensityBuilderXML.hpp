
#ifndef COMPWA_PHYSICS_INTENSITYBUILDER_HPP_
#define COMPWA_PHYSICS_INTENSITYBUILDER_HPP_

#include <memory>
#include <tuple>

#include "Physics/ParticleStateTransitionKinematicsInfo.hpp"

#include <boost/property_tree/ptree_fwd.hpp>

namespace ComPWA {

class Kinematics;
class Intensity;

namespace Physics {

class Amplitude;

namespace HelicityFormalism {
class HelicityKinematics;
}

namespace IntensityBuilderXML {

std::tuple<std::shared_ptr<Intensity>,
           std::shared_ptr<HelicityFormalism::HelicityKinematics>>
createIntensityAndKinematics(const boost::property_tree::ptree &pt);

std::shared_ptr<HelicityFormalism::HelicityKinematics>
createHelicityKinematics(std::shared_ptr<PartList> partL,
                         const boost::property_tree::ptree &pt);

std::shared_ptr<HelicityFormalism::HelicityKinematics>
createHelicityKinematics(std::shared_ptr<PartList> partL,
                         const boost::property_tree::ptree &pt);

ParticleStateTransitionKinematicsInfo
createKinematicsInfo(std::shared_ptr<PartList> partL,
                     const boost::property_tree::ptree &pt);

FourMomentum createFourMomentum(const boost::property_tree::ptree &pt);

std::shared_ptr<ComPWA::Intensity>
createIntensity(std::shared_ptr<PartList> partL,
                std::shared_ptr<Kinematics> kin,
                const boost::property_tree::ptree &pt);

std::shared_ptr<Intensity>
createIncoherentIntensity(std::shared_ptr<PartList> partL,
                          std::shared_ptr<Kinematics> kin,
                          const boost::property_tree::ptree &pt);
std::shared_ptr<Intensity>
createCoherentIntensity(std::shared_ptr<PartList> partL,
                        std::shared_ptr<Kinematics> kin,
                        const boost::property_tree::ptree &pt);

std::shared_ptr<Intensity>
createStrengthIntensity(std::shared_ptr<PartList> partL,
                        std::shared_ptr<Kinematics> kin,
                        const boost::property_tree::ptree &pt);

std::shared_ptr<Amplitude>
createAmplitude(std::shared_ptr<PartList> partL,
                std::shared_ptr<Kinematics> kin,
                const boost::property_tree::ptree &pt);

std::shared_ptr<Amplitude>
createCoefficientAmplitude(std::shared_ptr<PartList> partL,
                           std::shared_ptr<Kinematics> kin,
                           const boost::property_tree::ptree &pt);

std::shared_ptr<Amplitude>
createSequentialAmplitude(std::shared_ptr<PartList> partL,
                          std::shared_ptr<Kinematics> kin,
                          const boost::property_tree::ptree &pt);

std::shared_ptr<Amplitude>
createHelicityDecay(std::shared_ptr<PartList> partL,
                    std::shared_ptr<HelicityFormalism::HelicityKinematics> kin,
                    const boost::property_tree::ptree &pt);

}; // namespace IntensityBuilderXML

} // namespace Physics
} // namespace ComPWA

#endif
