#ifndef EXOTIC_GENERATOR_I3EXOTICGENERATORUTILS_H_INCLUDED
#define EXOTIC_GENERATOR_I3EXOTICGENERATORUTILS_H_INCLUDED

/**
 * A Monopole/Qball/ShadowCharge Generator Module
 * Copyright (c) 2004 - 2014 IceCube Collaboration
 * Version $Id$
 *
 * @file I3ExoticGeneratorUtils.h
 * @date $Date$
 * @author bchristy
 * @author olivas
 * @brief Utility namespace for the I3ExoticGenerator module
 */

#include "dataclasses/I3Constants.h"
#include "dataclasses/physics/I3Particle.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Direction.h"
#include "phys-services/I3RandomService.h"

/**
 *Namespace used to generate the position of the particle
 */
namespace I3ExoticGeneratorUtils {
    /**
    *@brief Calculate the solid angle for the given configuration
    *@param zmin minimum zenith angle
    *@param zmax maximum zenith angle
    *@param amin minimum azimuth angle
    *@param amax maximum azimuth angle
    *@return solid angle
     */
    double CalcSolidAngle(double, double, double, double);

    /**
    *@brief Calculate corresponding gamma for given beta
    *@param beta beta
    *@return gamma
     */
    double beta2gamma(double);

    /**
    *@brief Calculate corresponding beta for given gamma
    *@param gamma gamma
    *@return beta
     */
    double gamma2beta(double);

    /**
    * @brief Sanity checks for the mass parameter
    * @param mass
    * @param betaRange
    */
    void CheckMass(const double &, std::vector<double> &, const I3Particle::ParticleType &);

    /**
    *@brief Check which usecase (gamma only, beta range, only one beta), Set other variables accordingly and apply sanity checks (0 < beta < 1, 1 < gamma < 1000, etc.)
    *@param betaRange vector containing the minimum and maximum beta value, Any NAN will be ignored. Single value or borders of beta range need to be set.
    *@param powerLawIndex power law factor used if a beta range is defined. If NAN or 0, an uniform distribution will be applied.
    */
    void CheckAndSetBeta(std::vector<double> &, const double &);

      /**
      *@brief Sanity checks for the Disk Distance parameters
      *@param diskDist
      *@param diskRad
      */
    void CheckDisk(const double &, const double &);

    /**
    *@brief Sanity checks for the length parameter and apply default value if length=-1
    *@param length
    *@param defaultValue value of length if length was set to -1
    */
    void CheckAndSetLength(double &, const double &);

    /**
    *@brief Only internally used to abstract CheckAndSetZenith and CheckAndSetAzimuth to a common function
    *@param name name used in logs, azimuth or zenith
    *@param range vector of size 2 with minimum and maximum angle
    *@param validMin minimum of range for sanity check
    *@param validMax maximum of range for sanity check
    */
    void CheckAndSetAngle(const std::string &, const std::vector<double> &, const double &, const double &);

    /**
    *@brief Applies sanity checks to the zenith and extracts the values from the vector the min/max variables
    *@param zenithRange vector of size 2 with minimum and maximum angle
    */
    void CheckAndSetZenith(const std::vector<double> &);

    /**
    *@brief Applies sanity checks to the azimuth and extracts the values from the vector the min/max variables
    *@param azimuthRange vector of size 2 with minimum and maximum angle
    */
    void CheckAndSetAzimuth(const std::vector<double> &);

    /**
    *@brief Sanity checks for the shiftCenter parameter
    *@param shiftCenter
    */
    void CheckShiftCenter(const std::vector<double> &);

    /**
    *@brief Applies sanity checks to position on generation disk parameters depending on a set position or random position
    *@param radOnDisk radius in disk
    *@param aziOnDisk azimuth on disk
    *@param diskRad radius of disk
    */
    void CheckDiskPos(const double &, const double &, const double &);

    /**
    *@brief Helper function for uniformal random numbers
    *@param random I3RandomServicePtr to randomness service
    *@param min minimal value returned
    *@param max maximal value returned
    */

    double RandomUniformSampled(I3RandomServicePtr, double const, double const);
    /**
    *@brief Helper function for cosinus weighted random numbers
    *@param random I3RandomServicePtr to randomness service
    *@param min_angle minimal value returned
    *@param max_angle maximal value returned
    */

    double RandomCosinusSampled(I3RandomServicePtr, double const, double const);
    /**
    *@brief Helper function for circular weighted random numbers (i.e. getting properly weighted radius of a circle)
    *@param random I3RandomServicePtr to randomness service
    *@param min_radius minimal value returned
    *@param max_radius maximal value returned
    */

    double RandomCircularSampled(I3RandomServicePtr, double const, double const);
    /**
    *@brief Helper function for powerlaw weighted random numbers
    *@param random I3RandomServicePtr to randomness service
    *@param min minimal value returned
    *@param max maximal value returned
    *@param powerLawIndex powerlaw to be weighted by, only positive values are supported but interpreted as negative
    */
    double RandomPowerLawSampled(I3RandomServicePtr, double const, double const, double const);

}


#endif //EXOTIC_GENERATOR_I3EXOTICGENERATORUTILS_H_INCLUDED
