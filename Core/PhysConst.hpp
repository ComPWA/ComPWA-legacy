
//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

//! PhysConst provides particle information to the framework
/*!
 * @file PhysConst.hpp
 *\class PhysConst
 *      PhysConst reads particle properties from an xml file. And provides
 *      it to the framework via a singleton structure. The structure of
 *      the xml file is defined in particleSchema.xsd. It is also intended
 *      to store other physical constants in this way.
 */

#ifndef PHYSCONST_HPP_
#define PHYSCONST_HPP_

#include <iostream>
#include <vector>

#include <boost/property_tree/ptree.hpp>

#include <Core/Parameter.hpp>
#include <Core/Spin.hpp>
#include <Core/Properties.hpp>

namespace ComPWA {

class PhysConst {
  
public:
  static PhysConst *CreateInstance(std::string file = "./particles.xml") {
    if (_inst)
      throw std::runtime_error("PhysConst::CreateInstance() | Instance already "
                               "exists. Use Instance() to access it!");
    _inst = new PhysConst(file);
    return _inst;
  }

  static PhysConst *CreateInstance(boost::property_tree::ptree pt) {
    if (_inst)
      throw std::runtime_error("PhysConst::CreateInstance() | Instance already "
                               "exists. Use Instance() to access it!");

    _inst = new PhysConst(pt);
    return _inst;
  }

  static PhysConst *Instance() {
    if (!_inst) {
      throw std::runtime_error("No instance of PhysConst created! "
                               "Create one first!");
    }
    return _inst;
  }

  ~PhysConst();

  PhysConst(PhysConst const &) = delete;

  void operator=(PhysConst const &) = delete;

  const Constant &FindConstant(const std::string &name) const;

  const ParticleProperties &FindParticle(const std::string &name) const;

  const ParticleProperties &FindParticle(pid id) const;

  bool ParticleExists(const std::string &name) const;

protected:
  PhysConst(std::string filePath = "./particles.xml");

  PhysConst(boost::property_tree::ptree pt);

  //! Singleton stuff
  static PhysConst *_inst;

  void readTree(boost::property_tree::ptree pt);

  //! Input file name
  std::string _particleFileName;
  //! Input file name
  std::string _constantFileName;

  std::vector<ParticleProperties> _partList;
  std::vector<Constant> _constList;
};

} /* namespace ComPWA */

#endif /* PHYSCONST_HPP_ */
