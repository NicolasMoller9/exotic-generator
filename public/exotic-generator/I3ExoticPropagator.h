#ifndef EXOTIC_PROPAGATOR_I3EXOTICPROPAGATOR_H_INCLUDED
#define EXOTIC_PROPAGATOR_I3EXOTICPROPAGATOR_H_INCLUDED
/**
 * class: I3ExoticPropagator.h
 * Copyright (c) 2008 IceCube Collaboration
 * Version $Id: I3ExoticPropagator.h 124417 2014-10-10 15:16:00Z jacobi $
 *
 * Date 06 Feb 2008
 * @version $Revision: 124417 $
 * @date $Date: 2014-10-10 10:16:00 -0500 (Fr, 10. Okt 2014) $
 * @author Brian Christy <bchristy@icecube.umd.edu>
 * @author Alex Olivas <olivas@icecube.umd.edu>
 *
 * @brief A module to convert exotics into a chain of propagated particles through the ice
 *
 */

#include "icetray/I3Module.h"
#include <string>
#include "dataclasses/physics/I3MCTreeUtils.h"

class I3ExoticPropagator : public I3Module
{
 public:

  I3ExoticPropagator(const I3Context& ctx);
  ~I3ExoticPropagator();

  void Configure();
  void DAQ(I3FramePtr frame);
  void Finish();

  double CalculateMeanFreePath(int partID, double sig0, double b);

 private:
  I3ExoticPropagator();
  I3ExoticPropagator(const I3ExoticPropagator&);
  I3ExoticPropagator& operator=(const I3ExoticPropagator&);

  std::string inputTreeName_;
  std::string outputTreeName_;
  std::string infoName_;

  double meanFreePath_;
  double sigma0_;
  bool scaleEnergy_;
  double energyScaleFactor_;
  bool useCorrectDecay_;
  bool keepPrimary_;

  double maxdistfromcenter_;
  bool calcEn_;
  double stepSize_;
  double minlength_;
  double maxlength_;
  bool profiling_;
  bool checkParticle_;

  SET_LOGGER("I3ExoticPropagator");

};

#endif //EXOTIC_PROPAGATOR_I3EXOTICPROPAGATOR_H_INCLUDED
