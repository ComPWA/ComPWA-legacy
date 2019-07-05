#include "UpdatePTreeParameter.hpp"

namespace ComPWA {
namespace Tools {

void updateParameterRangeByType(boost::property_tree::ptree &Tree,
                                const std::string &ParameterType, double Min,
                                double Max) {
  updateParameter(Tree, "Type", ParameterType,
                  // dummy args real args  not update   update
                  0, false, Min, Max, false, false, true);
}

void updateParameterRangeByName(boost::property_tree::ptree &Tree,
                                const std::string &ParameterName, double Min,
                                double Max) {
  updateParameter(Tree, "Name", ParameterName,
                  // dummy args real args  not--update  update
                  0, false, Min, Max, false, false, true);
}

void updateParameterValue(boost::property_tree::ptree &Tree,
                          const std::string ParameterName, double Value) {
  updateParameter(Tree, "Name", ParameterName,
                  // real dummy ------ args, update not---update
                  Value, false, -999, -999, true, false, false);
}

void fixParameter(boost::property_tree::ptree &Tree,
                  const std::string ParameterName, double Value) {
  if ((int)Value == -999) {
    updateParameter(Tree, "Name", ParameterName,
                    // dummy real, dummy args, notup update not-update
                    -999, true, -999, -999, false, true, false);
  } else {
    updateParameter(Tree, "Name", ParameterName,
                    // real args, dummy args, up----date not-update
                    Value, true, -999, -999, true, true, false);
  }
}

void releaseParameter(boost::property_tree::ptree &Tree,
                      const std::string ParameterName, double Value) {
  if ((int)Value == -999) {
    updateParameter(Tree, "Name", ParameterName,
                    // dummy real,  dummy args, not-up update not-up
                    -999, false, -999, -999, false, true, false);
  } else {
    updateParameter(Tree, "Name", ParameterName,
                    // real args, dummy args, not-up update not-up
                    Value, false, -999, -999, false, true, false);
  }
}

void updateParameter(boost::property_tree::ptree &Tree,
                     const std::string &KeyType, const std::string &KeyValue,
                     double Value, bool Fix, double Min, double Max,
                     bool UpdateValue, bool UpdateFix, bool UpdateRange) {

  for (auto &Node : Tree.get_child("")) {
    // If not not a parameter node, recursively update this node.
    if (Node.first != "Parameter") {
      if (Node.first != "DecayParticle" && Node.first != "DecayProducts" &&
          Node.first != "CanonicalSum") {
        updateParameter(Node.second, KeyType, KeyValue, Value, Fix, Min, Max,
                        UpdateValue, UpdateFix, UpdateRange);
      }
      continue;
    }
    // If it is a parameter node, and it's target parameter,
    // update parameter properties.
    if (KeyValue != Node.second.get<std::string>("<xmlattr>." + KeyType)) {
      continue;
    }
    if (UpdateValue)
      Node.second.put("Value", Value);
    if (UpdateFix)
      Node.second.put("Fix", Fix);
    if (UpdateRange) {
      Node.second.put("Min", Min);
      Node.second.put("Max", Max);
    }
  }
}

void updateParameter(boost::property_tree::ptree &Tree,
                     const FitParameterList &FitParameters) {
  for (const auto &FitPar : FitParameters) {
    updateParameter(Tree, "Name", FitPar.Name, FitPar.Value, FitPar.IsFixed,
                    FitPar.Bounds.first, FitPar.Bounds.second, true, true,
                    FitPar.HasBounds);
  }
}

} // end namespace Tools
} // end namespace ComPWA
