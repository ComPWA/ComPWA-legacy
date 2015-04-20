#ifndef HELICITYAMPLITUDE_HPP_
#define HELICITYAMPLITUDE_HPP_

#include <HelicityStateDefinitions.hpp>
#include <Core/Amplitude.hpp>

#include <complex>

namespace HelicityFormalism {

  class HelicityAmplitude: public Amplitude {
    SphericalWaveTwoParticleState initial_state;
    PlaneWaveTwoParticleState final_state;

  public:
    HelicityAmplitude();
    virtual ~HelicityAmplitude();

    bool validate() const;

    std::complex<double> evaluate() const;

    virtual const double integral() =0;
    virtual const double integral(ParameterList& par) =0;
    virtual const double normalization() =0;
    virtual const double normalization(ParameterList& par) =0;
    virtual double getMaxVal(ParameterList& par,
        std::shared_ptr<Generator> gen) = 0;
    virtual double getMaxVal(std::shared_ptr<Generator> gen) = 0;

    virtual const ParameterList& intensity(dataPoint& point,
        ParameterList& par) =0;
    virtual const ParameterList& intensity(dataPoint& point) =0;
    virtual const ParameterList& intensityNoEff(dataPoint& point) =0;
    virtual const ParameterList& intensity(std::vector<double> point,
        ParameterList& par) =0;

    virtual const bool fillStartParVec(ParameterList& outPar) =0;
    virtual void setParameterList(ParameterList& par) =0;

    virtual void printAmps() = 0;
    virtual void printFractions() = 0;

    virtual Amplitude* Clone() = 0;
  };

} /* namespace HelicityFormalism */
#endif /* HELICITYAMPLITUDE_HPP_ */
