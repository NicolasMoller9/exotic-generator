/**
 * class: I3ExoticGeneratorUtils.cxx
 * Copyright (c) 2004 - 2014 IceCube Collaboration
 * Version $Id$
 *
 * @date $Date$
 * @author Brian Christy <bchristy@icecube.umd.edu>
 * @author Alex Olivas <olivas@icecube.umd.edu>
 * @brief Utility Namespace for the I3ExoticGenerator Module
 */

#include "I3ExoticGeneratorUtils.h"
#include <cmath>

static const std::map<I3Particle::ParticleType, const std::pair<const double, const double> > 
             mass_validity_ranges = { { I3Particle::Monopole, {1e5 * I3Units::GeV, 1e17 * I3Units::GeV}}, 
                                      { I3Particle::Qball, {1e17 * I3Units::GeV, 1e29 * I3Units::GeV}}, 
                                      { I3Particle::ShadowCharge, {0 * I3Units::GeV, 1e100 * I3Units::GeV}},  // mass not constrained yet
                                    };

double
I3ExoticGeneratorUtils::CalcSolidAngle(double zmin, double zmax,
                                         double amin, double amax) {
  return (amax / I3Units::radian - amin / I3Units::radian) * \
    (cos(zmin / I3Units::radian) - cos(zmax / I3Units::radian));
}

double I3ExoticGeneratorUtils::beta2gamma(double beta) { //helper function to translate gamma <-> beta
  return 1.0 / (sqrt(1.0 - pow(beta, 2)));
}

double I3ExoticGeneratorUtils::gamma2beta(double gamma) {
  return sqrt(1.0 - pow(gamma, -2));
}

void
I3ExoticGeneratorUtils::CheckMass(const double &mass, std::vector<double> &betaRange, const I3Particle::ParticleType &particleType) {
  if (!std::isnan(mass)) {
    if (mass < std::get<0>(mass_validity_ranges.at(particleType)) || mass > std::get<1>(mass_validity_ranges.at(particleType))) {
      log_warn("%s Mass (%f GeV) out of range.", boost::lexical_cast<std::string>(particleType ).c_str(), mass / I3Units::GeV);
    }
  }
}

void
I3ExoticGeneratorUtils::CheckAndSetBeta(std::vector<double> &betaRange,
                                               const double &powerLawIndex) {

  //remove all nan values from the beta range
  betaRange.erase(std::remove_if(betaRange.begin(), betaRange.end(), [](const double & o) { return std::isnan(o);}), betaRange.end());

  if (betaRange.size() != 1 && betaRange.size() != 2) {
    log_fatal("Beta needs to be one or two non nan elements!");
  }

  if (betaRange.size() == 1){
    betaRange.push_back(betaRange[0]);
  }

  if (betaRange[0]>betaRange[1]){
    log_fatal("Negative beta range: minimum beta is bigger than maximum beta. Aborting!");
  }

  if (betaRange[0] == betaRange[1]) {
    log_trace("Single beta usecase");
    if (betaRange[0] < 1e-6 || betaRange[0] >= 1.) {
      log_fatal("Beta (%f) out of range.", betaRange[0]);
    }
    if (!std::isnan(powerLawIndex) && powerLawIndex != 0) {
      log_fatal("No range for beta is defined, cannot apply a power law.");
    }
  } else {
    if (betaRange[0] < 1e-6 || betaRange[0] >= 1.) {
      log_fatal("BetaMin (%f) out of range.", betaRange[0]);
    }
    if (betaRange[1] < 1e-6 || betaRange[1] >= 1.) {
      log_fatal("BetaMax (%f) out of range.", betaRange[1]);
    }
  }
}

void
I3ExoticGeneratorUtils::CheckDisk(const double &diskDist, const double &diskRad) {
  // stuff regarding the trajectory of the exotic
  if (diskDist < 0 || diskDist > 20. * I3Units::km) {
    log_fatal("Disk_dist (%f m) out of range.", diskDist / I3Units::m);
  }
  if (diskRad < 0) {
    log_fatal("Disk_Rad (%f m) out of range.", diskRad / I3Units::m);
  }
}

void
I3ExoticGeneratorUtils::CheckAndSetLength(double &length, const double &defaultValue) {
  if (length ==
      -1.) {   // The default length is 2* diskDist_ but it can be configured to any valid value, including NaN
    length = defaultValue;
  }

  if (length < 0) {
    log_fatal("Can't have a negative length besides -1 to indicate default value.");
  }
}

void
I3ExoticGeneratorUtils::CheckAndSetAngle(const std::string &name, const std::vector<double> &range, const double &validMin, const double &validMax) {
  if (range.size() != 2) {
    log_fatal("Please configure %sRange with a list containing *only* the minimum and the maximum %s.", name.c_str(),
              name.c_str());
  }

  if (range[0] / I3Units::degree < validMin || range[0] / I3Units::degree > validMax) {
    log_fatal("%sMin (%f rad) out of range.", name.c_str(), range[0] / I3Units::radian);
  }
  if (range[1] / I3Units::degree < validMin || range[1] / I3Units::degree > validMax) {
    log_fatal("%sMax (%f rad) out of range.", name.c_str(), range[1] / I3Units::radian);
  }
  if (range[0] > range[1]) {
    log_error("%s: (%f rad, %f rad)", name.c_str(), range[0] / I3Units::radian, range[1] / I3Units::radian);
    log_fatal("Can't have lower boundary > upper boundary.");
  }
  if (range[0] / I3Units::degree > validMin || range[1] / I3Units::degree < validMax) {
    log_warn("The direction of the exotic is restricted!\n%s range: [%.2f, %.2f] degree.", name.c_str(),
             range[0] / I3Units::degree, range[1] / I3Units::degree);
  }
}

void
I3ExoticGeneratorUtils::CheckAndSetZenith(const std::vector<double> &zenithRange) {
  CheckAndSetAngle("zenith", zenithRange, 0., 180.);
}

void
I3ExoticGeneratorUtils::CheckAndSetAzimuth(const std::vector<double> &azimuthRange) {
  CheckAndSetAngle("azimuth", azimuthRange, 0., 360.);
}

void
I3ExoticGeneratorUtils::CheckShiftCenter(const std::vector<double> &shiftCenter) {
  if (shiftCenter.size() != 3) {
    log_fatal("If you want to shift the center, please provide an array containing x,y and z coordinates.");
  }
  if (shiftCenter[0] != 0. || shiftCenter[1] != 0. || shiftCenter[2] != 0.) {
    log_warn("Shifting positions from center of IceCube by (%.2f,%.2f,%.2f).", shiftCenter[0], shiftCenter[1],
             shiftCenter[2]);
  }
}

void
I3ExoticGeneratorUtils::CheckDiskPos(const double &radOnDisk, const double &aziOnDisk,
                                       const double &diskRad) {
   if (std::isnan(diskRad)) {
    log_fatal("DiskRad cannot be NAN!");
   }
   if (!std::isnan(aziOnDisk)) {
      if (aziOnDisk / I3Units::degree < 0. || aziOnDisk / I3Units::degree > 360.) {
        log_fatal("Azimuth on generation disk (%f rad) out of range.", aziOnDisk / I3Units::radian);
      }
      log_warn("The azimuth on the generation disk is NOT randomized!");
   }
   if (!std::isnan(radOnDisk)) {
      if (radOnDisk < 0 || radOnDisk > diskRad) {
        log_fatal(
              "Radius on generation disk (%f m) out of range. Make sure it's not larger than the disk radius.",
              radOnDisk / I3Units::m);
      }
      log_warn("The radius on the generation disk is NOT randomized!");
   }
}



double I3ExoticGeneratorUtils::RandomUniformSampled(I3RandomServicePtr random, double const  min, double const max){
  return random->Uniform(min, max);
}

double I3ExoticGeneratorUtils::RandomCosinusSampled(I3RandomServicePtr random, double const min_angle, double const max_angle){
  return acos(random->Uniform(cos(max_angle), cos(min_angle)));
}

double I3ExoticGeneratorUtils::RandomCircularSampled(I3RandomServicePtr random, double const min_radius, double const max_radius)
{
  return sqrt(random->Uniform(pow(min_radius, 2), pow(max_radius, 2)));
}


double I3ExoticGeneratorUtils::RandomPowerLawSampled(I3RandomServicePtr random, double const min, double const max, double const powerLawIndex){
  if (std::isnan(powerLawIndex) || powerLawIndex == 0) {
    return random->Uniform(min, max);
  }
  else {
    double absIndex = std::abs(powerLawIndex);
    double x_prime;
    if (absIndex == 1.) {
      double border = log(max / min);
      double uniform = random->Uniform(0, border);
      x_prime = min * exp(uniform);
    }
    else {
      double oneMinusIndex = 1. - absIndex;
      double border = 1.0 / oneMinusIndex * (pow(max, oneMinusIndex) - pow(min, oneMinusIndex));
      double uniform = random->Uniform(0, border);
      x_prime = pow(oneMinusIndex * uniform + pow(min, oneMinusIndex), 1.0 / oneMinusIndex);
    }

    if (powerLawIndex < 0) {
      return max - (x_prime - min);
    } else {
      return x_prime;
    }
  }
}

