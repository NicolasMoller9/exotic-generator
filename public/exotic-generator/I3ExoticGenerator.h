#ifndef EXOTIC_GENERATOR_I3EXOTICGENERATOR_H_INCLUDED
#define EXOTIC_GENERATOR_I3EXOTICGENERATOR_H_INCLUDED
/**
 * A Monopole/Qball/ShadowCharge Generator Module
 * Copyright (c) 2004 - 2014 IceCube Collaboration
 * Version $Id$
 *
 * @file I3ExoticGenerator.h
 * @date $Date$
 * @author jacobi
 * @author bchristy
 * @author olivas
 * @author flauber
 * @brief A module to generate exotic particles (Monopoles, Qballs and Shadow-charges)
 */

#include <string>
#include <cmath>
#include "icetray/I3Module.h"
#include "icetray/I3Frame.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "exotic-generator/I3ExoticGeneratorUtils.h"



class I3ExoticGenerator : public I3Module
{
 public:

  I3ExoticGenerator(const I3Context& ctx);
  ~I3ExoticGenerator();

  void Configure();
  void DAQ(I3FramePtr frame);
  void Finish();

 private:
  I3ExoticGenerator();
  I3ExoticGenerator(const I3ExoticGenerator&);
  I3ExoticGenerator& operator=(const I3ExoticGenerator&);

  std::string treeName_;
  std::string infoName_;
  int Nevents_;
  int part_id_;
  double mass_;
  int charge_;
  double beta_;
  double gamma_;
  double speed_;
  double energy_;
  bool useBetaRange_;
  double diskDist_;
  double diskRad_;
  double radOnDisk_;
  double aziOnDisk_;
  double startTime_;
  double weight_;
  double powerLawIndex_;
  double length_;
  double totalweight_;
  I3MapStringDoublePtr exinfo_dict_config_;
  I3Particle::ParticleType particleType_;
  /**
   * Zenith and Azimuth follow convention defined in I3Direction.h
   * Angle points to origin of particle, hence zenith=0 is down-going
   */
  std::vector<double> betaRange_;
  std::vector<double> zenithRange_;
  std::vector<double> azimuthRange_;
  std::vector<double> shiftCenter_;
  std::vector<double> precalculated_betas_;
  bool need_precalculate_;

  SET_LOGGER("I3ExoticGenerator");

};
#endif //EXOTIC_GENERATOR_I3EXOTICGENERATOR_H_INCLUDED
