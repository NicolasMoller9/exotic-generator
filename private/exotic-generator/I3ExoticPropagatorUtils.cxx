/**
 * I3ExoticPropagatorUtils.cxx
 * Propagating exotics
 * @version $Revision: $
 * @date $Date: $
 * @author $$
 * @brief Propagate exotics
*/

#include "I3ExoticPropagatorUtils.h"
#include "I3ExoticEnergyLoss.h"

#include <dataclasses/I3Constants.h>
#include <map>
#include "icetray/I3Frame.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

// Particles allowed for propagation:
static const std::set<I3Particle::ParticleType> AllowedPrimaries = {
    I3Particle::Monopole,
    I3Particle::Qball,
    I3Particle::ShadowCharge
};

// Allowed mass interval for these particles:
static const std::map<I3Particle::ParticleType, const std::pair<const double, const double> > 
             mass_validity_ranges = { { I3Particle::Monopole, {1e5 * I3Units::GeV, 1e17 * I3Units::GeV}}, 
                                      { I3Particle::Qball, {1e17 * I3Units::GeV, 1e29 * I3Units::GeV}}, 
                                      { I3Particle::ShadowCharge, {0 * I3Units::GeV, 1e100 * I3Units::GeV}},  // mass not constrained yet
                                    };

/** 
 *@param frame                Instance of current frame
 *@param mctree               Original mctree created by generator
 *@param prop_tree            Copy of mctree where new particles will be added
 *@param random               Pointer to random number service
 *@param meanFreePath_        Mean free path of the catalysis process
 *@param useCorrectDecay_     Whether to use the physical decay process or a simplification (ex: monopole proton decay --> 2 opposite cascades of 500 MeV VS one 1 GeV cascade)
 *@param deltaEnergy_         Energy of the catalyzed process 
 *@param maxdistfromcenter_   Propagation limit (stop propagating when limit reached)
 *@param keepPrimary_         Wether to keep the primary exotic track (allowing it to produce light), or set it dark 
 *@param calcEn_              (if keepPrimary_) Whether to calculate energy loss. If true, the speed will be decreased from a track-segment to the next
 *@param mass_                (if calcEn_) Mass of the particle
 *@param charge_              (if calcEn_) Charge of the particle. For monopoles, it is the n factor in the definition of the magnetic charge (g=n*gD)
 *@param profiling_           (if calcEn_): If the speed profile should be written into frame
 *@param stepSize_            (if splitPrimary) Length of the sub-tracks 
 *@param minlength_           (if calcEn_): If no step size set, this represents smallest segment value allowed
 *@param maxlength_           (if calcEn_): If no step size set, this represents largest segment value allowed
 *@param checkParticle_       Whether to include exhaustive checking
 */

void I3ExoticPropagatorUtils::PropagateExotic(I3FramePtr frame,
                                                I3MCTreeConstPtr mctree,
                                                I3MCTreePtr prop_tree,
                                                I3RandomServicePtr random,
                                                // params
                                                double meanFreePath_,
                                                bool useCorrectDecay_,
                                                double deltaEnergy_,
                                                double maxdistfromcenter_,
                                                bool keepPrimary_,
                                                bool calcEn_,
                                                double mass_,
                                                double charge_,
                                                bool profiling_,
                                                double stepSize_,
                                                double minlength_,
                                                double maxlength_,
                                                bool checkParticle_) {
    
    log_debug("Entering PropagateExotic");

    if (calcEn_ && !keepPrimary_) {
        log_warn("The primary track cannot be split if the option KeepPrimary is false! calcEn_ turned off.");
    }


    for (auto const &exo:*mctree) {  // loop through the exotic particles in the MCTree
        if (AllowedPrimaries.find(exo.GetType()) == AllowedPrimaries.end()) continue;  // if the particle is not an exotic allowed by the code, ignore it
        if (std::isnan(exo.GetLength())) {  // Check the particle has a defined track length
            log_error("The exotic has a length of NaN. I don't know where to put the track. I'M NOT PROPAGATING THIS!");}
        if (calcEn_ && exo.GetType() == I3Particle::ShadowCharge) {
        log_warn("You are trying to calculate the energy loss of a shadow charge, although shadow charges don't lose any energy. Be sure of what you are doing!");
        }
        // Get time and position coordinate of the particle
        const double X_exo(exo.GetX());
        const double Y_exo(exo.GetY());
        const double Z_exo(exo.GetZ());
        const double T_exo(exo.GetTime());
        const double THETA_exo = I3Constants::pi - exo.GetZenith();
        const double PHI_exo = exo.GetAzimuth() - I3Constants::pi;
        const double beta_exo = (exo.GetSpeed() * I3Units::ns / I3Units::m) / I3Constants::c;


        // CATALYSIS PROCESS (if a mean free path is given)
        int nint; // number of catalyzed process occurring along the track
        if (std::isnan(meanFreePath_)){
            nint = 0;
            log_debug("The mean free path was Nan, so no catalysis process were simulated");
        } else {
            nint = static_cast<int>(random->Poisson(exo.GetLength() / meanFreePath_));
            log_debug("%f catalysis processes will be simulated along the track", nint);
            if (exo.GetType() == I3Particle::Monopole && beta_exo>0.09){
                log_warn("The speed of the simulated monopole is too high to catalyze proton decay: %f > 0.09. Proton decay turned off!", beta_exo);
                nint = 0;
            };
        }
        for (int i(0); i < nint; ++i) {  // loop through the catalyzed events along the track
            double l(random->Uniform(0, exo.GetLength())); // distance between beginning of the track and the catalyzed event
            // Starting position of the daughter particle(s):
            double x = X_exo + l * sin(THETA_exo) * cos(PHI_exo);
            double y = Y_exo + l * sin(THETA_exo) * sin(PHI_exo);
            double z = Z_exo + l * cos(THETA_exo);
            // random direction of the daughter particle
            double phi(random->Uniform(0, 2. * I3Constants::pi));
            double theta(acos(random->Uniform(-1., 1)));

            // MONOPOLE: PROTON DECAY
            if (exo.GetType() == I3Particle::Monopole) {
                if (!useCorrectDecay_) {  // (default): gives a smaller MCTree and the light yield at the end is the same
                    I3Particle eplus; // Just one positron carrying all energy of the decay
                    eplus.SetType(I3Particle::EPlus);
                    eplus.SetLocationType(I3Particle::InIce);
                    eplus.SetEnergy(deltaEnergy_);
                    eplus.SetThetaPhi(theta, phi);
                    eplus.SetPos(x, y, z);
                    eplus.SetTime(T_exo + l / exo.GetSpeed());
                    prop_tree->append_child(exo, eplus); // add the positron in the new MCTree as a child of the monopole
                    log_debug("Added positron with position: (%f,%f,%f) and energy %f GeV",
                        eplus.GetX(), eplus.GetY(), eplus.GetZ(),
                        eplus.GetEnergy() / I3Units::GeV);
                }
                else {   // eplus + pi_0, the most dominant proton decay channel
                    I3Particle eplus;
                    eplus.SetType(I3Particle::EPlus);
                    eplus.SetLocationType(I3Particle::InIce);
                    eplus.SetEnergy(460 * I3Units::MeV);
                    eplus.SetThetaPhi(theta, phi);
                    eplus.SetPos(x, y, z);
                    eplus.SetTime(T_exo + l / exo.GetSpeed());
                    prop_tree->append_child(exo, eplus);
                    log_debug("Added positron with position: (%f,%f,%f) and energy %f GeV",
                        eplus.GetX(), eplus.GetY(), eplus.GetZ(),
                        eplus.GetEnergy() / I3Units::GeV);
        
                    I3Particle hadron;
                    hadron.SetType(I3Particle::Pi0);   // Note: pi_0 and eplus are treated equally with photonics
                    hadron.SetLocationType(I3Particle::InIce);
                    hadron.SetEnergy(480 * I3Units::MeV);
                    double phi_corr = phi + I3Constants::pi;
                    if (phi_corr > 2 * I3Constants::pi) {phi_corr = phi_corr - 2 * I3Constants::pi;}
                    hadron.SetThetaPhi(I3Constants::pi - theta, phi_corr);
                    hadron.SetPos(x, y, z);
                    hadron.SetTime(T_exo + l / exo.GetSpeed());
                    prop_tree->append_child(exo, hadron);
                    log_debug("Added hadron with position: (%f,%f,%f) and energy %f GeV",
                        hadron.GetX(), hadron.GetY(), hadron.GetZ(),
                        adron.GetEnergy() / I3Units::GeV);
                }
            } // end of monopole

            // OTHER EXOTIC EXEMPLE:
            //if (exo.GetType() == I3Particle::OtherExoticName) {
            //    ...
            //}
        } // end of the for loop through the catalyzed events along the primary track

        // Handle the primary exotic particle:
        auto initial_exotic_iter = prop_tree->find(exo);
        if (!keepPrimary_) {
            initial_exotic_iter->SetLength(nan(""));  // remove primary track
            continue;
        }
        
        // split the track if a stepsize or a step-length range is given
        if (!std::isnan(stepSize_) || (!std::isnan(minlength_) && !std::isnan(maxlength_))) {
            if (calcEn_) {
                if (mass_ < std::get<0>(mass_validity_ranges.at(exo.GetType())) || mass_ > std::get<1>(mass_validity_ranges.at(exo.GetType()))) {
                    log_fatal("Mass (%f GeV) is out of range", mass_ / I3Units::GeV);
                }
                if ((exo.GetType()==I3Particle::ShadowCharge) && charge_ > 137) {
                    log_fatal("Charge (%f) is out of range", charge_); 
                }
            }
            initial_exotic_iter->SetLength(I3ExoticPropagatorUtils::CalculateNextLength(stepSize_, initial_exotic_iter->GetEnergy(), mass_, maxlength_, minlength_));
            I3Particle initial_exotic = *initial_exotic_iter;

            //Check if particle is within the max distance its supposed to propagate
            //If not, particle start position will be used as max distance
            
            if ((initial_exotic.GetPos().GetR()) > maxdistfromcenter_) {
                log_warn("Start of particle further away from Detector Center than max distance its supposed to propagate; "
                        "Propagation will continue until it reaches same distance as original start (%f); " 
                        "Currently, MaxDistanceFromCenter variable set to %f, and WILL NOT BE USED; " 
                        "Please make sure this is what you intended (perhaps just didn't set MaxDistanceFromCenter?)",
                        initial_exotic.GetPos().GetR(), maxdistfromcenter_);
            }

            log_debug("-----Starting Propagate on Exotic----");
            log_debug("Position: (%f,%f,%f)", initial_exotic.GetX(), initial_exotic.GetY(), initial_exotic.GetZ());
            log_debug("Energy: %f GeV", initial_exotic.GetEnergy() / I3Units::GeV);
            log_debug("Length: %f m", initial_exotic.GetLength() / I3Units::m);

            //Each loop generates a new track segment and adds to the tree
            int count = 0;
            I3Particle temp0 = initial_exotic;
            I3ParticlePtr temp1;
            while (true) {
                ++count;

                temp1 = I3ExoticPropagatorUtils::AddExotic(temp0,
                                                mass_,
                                                charge_,
                                                calcEn_,
                                                stepSize_,
                                                minlength_,
                                                maxlength_, 
                                                checkParticle_);
                if (!temp1) {
                    log_fatal("New particle did not get filled!");
                    return; //silence temp could be null messages from IDEs //fhl
                }
                log_debug("New particle set with Energy %f(GeV)", temp1->GetEnergy() / I3Units::GeV);
                prop_tree->append_child(temp0.GetID(), *temp1); // Add new exotic to the output tree as a child of the previous monopole
                
                //End the while loop when particle went too far away
                if ((temp1->GetPos().GetR() > initial_exotic.GetPos().GetR()) && (temp1->GetPos().GetR() >= maxdistfromcenter_)) {
                    log_info("The exotic has moved through detector beyond %f(m) of detector center", temp1->GetPos().GetR() / I3Units::m);
                    break;
                }

                //catch just to ensure no infinite loops
                if (count > 100000) {
                    log_warn("Propagator Method Exited due to too many iterations.");
                    break;
                }

                //generated monopole becomes the previous monopole for the next loop iteration!
                temp0 = *temp1;
            }// End while loop

        } // of of track splitter

        else {
            if (calcEn_) {
                log_warn("Cannot calculate energy loss if no fixed step length or min+max step length are given! calcEn_ turned off.\n");
            }
            if (profiling_) {
                log_warn("Cannot write the speed profile if no fixed step length or min+max step length are given!\n");
            }
        }
                
    } // end of the for loop through primary particles
    
    // Add a profile (type `I3VectorDouble`) of the exotic speed for each track segment to the frame
    if (profiling_) {
        I3VectorDoublePtr profile(new I3VectorDouble()); // create a vector of doubles
        I3ExoticPropagatorUtils::ExoSpeedProfile(prop_tree, profile);
        frame->Put("ExoSpeedProfile", profile); // add to the frame a ExoSpeedProfile containing the speed profile
    }
}   




/**
 *Generates a new exotic segment at end of given segment
 *@param exo            The previous exotic particle segment
 *@param exoticMass     The mass of the exotic particle
 *@param charge         Charge of the particle. For monopoles, it is the n factor in the definition of the magnetic charge (g=n*gD)
 *@param calcEn         Whether to calculate energy loss or just segment the track
 *@param stepSize       Length of the segments when splitting the track (ex: 10m). 
                        If nan, the segments length will be calculated so that Energy loss along the segment is roughly 0.1% of total KE 
 *@param minlength      (if calcEn): If no step size set, this represents smallest segment value allowed
 *@param maxlength      (if calcEn): If no step size set, this represents largest segment value allowed
 *@param checkParticle  Whether to include exhaustive checking
 *@return The next exotic at the end of the previous exotic segment
 */
I3ParticlePtr I3ExoticPropagatorUtils::AddExotic(I3Particle exo,
                                                    double exoticMass,
                                                    double charge,
                                                    bool calcEn,
                                                    double stepSize,
                                                    double minlength,
                                                    double maxlength,
                                                    bool checkParticle) {
    log_debug("Entering AddExotic");
    double exoMass = exoticMass / I3Units::GeV;
    double startEnergy = exo.GetEnergy() / I3Units::GeV;
    double startBeta = I3ExoticPropagatorUtils::EnergyToBeta(startEnergy, exoMass);
    double startTime = exo.GetTime() / I3Units::ns;
    double trackLength = exo.GetLength() / I3Units::m;
    double timeNext = startTime + trackLength / (exo.GetSpeed() * I3Units::ns / I3Units::m);
    double energyNext = 0.0;
    if (calcEn) {
        double energyLoss = 0;
        if (exo.GetType() == I3Particle::Monopole) {
            energyLoss = I3ExoticEnergyLoss::CalculateMonopoleEnergyLoss(startBeta, trackLength, charge);
        }
        energyNext = startEnergy - energyLoss;
        log_debug("Energy Loss calculated as %f", energyLoss);
    } else {
        energyNext = startEnergy;
    }
    double const lengthNext = I3ExoticPropagatorUtils::CalculateNextLength(stepSize, energyNext, exoMass, maxlength, minlength);

    //Create new particle
    I3ParticlePtr exoNext(new I3Particle);
    exoNext->SetType(exo.GetType());
    exoNext->SetLocationType(I3Particle::InIce);
    exoNext->SetDir(exo.GetDir());
    exoNext->SetPos(exo.GetPos() + exo.GetDir() * trackLength * I3Units::m);
    exoNext->SetTime(timeNext * I3Units::ns);
    exoNext->SetEnergy(energyNext * I3Units::GeV);
    exoNext->SetSpeed(I3ExoticPropagatorUtils::EnergyToBeta(energyNext, exoMass) * I3Constants::c);
    exoNext->SetLength(lengthNext * I3Units::m);
    log_debug("New Particle Position (%f,%f,%f)", exoNext->GetX(), exoNext->GetY(), exoNext->GetZ());
    log_debug("New Particle Direction: %f*PI zenith, %f*PI azimuth", exoNext->GetZenith() / I3Constants::pi, exoNext->GetAzimuth() / I3Constants::pi);
    log_debug("New time %f (ns)", timeNext * I3Units::ns);
    log_debug("New Energy set to %f GeV", energyNext);
    log_debug("New Speed set to %f ", exoNext->GetSpeed() / I3Constants::c);
    log_debug("New Length set to %f", lengthNext);
 
    if (checkParticle) I3ExoticPropagatorUtils::CheckParticle(*exoNext, exoMass);
 
    return exoNext;
}

 

/**
 * Determines length to use for newly generated exotic segments
 * by estimating how far it could travel before losing 0.1%
 * of its kinetic energy, or by setting it manually
 *
 * @param stepSize      If given a value, the function will simply return this value
 * @param nextEnergy    Starting energy of particle
 * @param particleMass  Mass of particle
 * @param maxLength     Largest segment value allowed
 * @param minLength     Smallest segment value allowed
 *
 * @return A length over which the exotic would lose roughly 0.1% of its KE
 *  within user defined boundaries
 */

double I3ExoticPropagatorUtils::CalculateNextLength(double stepSize,
                                                double nextEnergy,
                                                double particleMass,
                                                double maxLength,
                                                double minLength) {
    if (!__builtin_isnan(stepSize)) {  // if stepSize is not Nan, i.e., if a stepSize is given, 
        return stepSize;               // the function CalculateNextLength will return the stepSize
    }

    double newLength = 0.0;
    //Pick length so Energy loss is roughly 0.1% of total KE
    const double AVG_ENERGY_LOSS = 1200.0;  //in GeV/m  ??? where does it come from?
    newLength = (0.001) * (nextEnergy - particleMass) / AVG_ENERGY_LOSS;

    if (__builtin_isnan(newLength)) {
        log_fatal("new length was NaN and was not set");
    }

    if (newLength > maxLength) {
        newLength = maxLength;
    } else if (newLength < minLength) {
        newLength = minLength;
    }

    if (newLength == 0.0) {
        log_fatal("new length was 0 and was not set");
    }
    return newLength;
}



/**
 * Function to convert kinetic energy to a beta value
 * @param energy Energy of particle
 * @param mass   Mass of particle
 * @return beta of particle
 */
double I3ExoticPropagatorUtils::EnergyToBeta(double energy, double mass) {
    double gamma = energy / mass;
    double beta = sqrt(1 - pow(gamma, -2));
    return beta;
}



/**
 * Responsible for performing extensive sanity checks on particle
 * result to ensure nothing went horribly wrong
 *
 * @param particle  Particle to check
 * @param checkmass Mass of particle to check
 */

 void I3ExoticPropagatorUtils::CheckParticle(I3Particle &particle, double checkmass) {
    if (AllowedPrimaries.find(particle.GetType()) == AllowedPrimaries.end())
        log_fatal("Particle type is not allowed! D'oh! (%d)", particle.GetType());
    if (__builtin_isnan(particle.GetX()) || __builtin_isnan(particle.GetY()) || __builtin_isnan(particle.GetZ()))
        log_fatal("Particle position is not set");
    if ((particle.GetPos().GetR() / I3Units::m) > 5000)
        log_fatal("Particle more than 5 km away...why am I looking at this?");
    if (__builtin_isnan(particle.GetZenith()) || __builtin_isnan(particle.GetAzimuth()))
        log_fatal("Particle direction is not set");
    if (__builtin_isnan(particle.GetEnergy()) || particle.GetEnergy() == 0.0)
        log_fatal("Particle has no energy");
    if (__builtin_isnan(particle.GetTime()))
        log_fatal("Particle has no time");
    if (__builtin_isnan(particle.GetLength()) || (particle.GetLength() <= 0.0))
        log_fatal("Something's wrong with the length: %f", particle.GetLength());
    double checkbeta = I3ExoticPropagatorUtils::EnergyToBeta(particle.GetEnergy(), checkmass);
    if ((particle.GetType() == I3Particle::Monopole) && (checkbeta > 0.99995))
        log_fatal("Pole out of Control!!! - speed of %f too much for assumptions of this module", checkbeta);
}



/**
 * Extract the speed profile of the exotic from the MCTree.
 * @param mctree      The mctree
 * @param speed_prof  Pointer to the vector storing the profile
 */

void I3ExoticPropagatorUtils::ExoSpeedProfile(I3MCTreeConstPtr mctree, I3VectorDoublePtr speed_prof) {
    speed_prof->clear();
 
    // Find the particle of highest energy in the MCTree
    I3Particle maxEnergeticPrimary = I3Particle();
    maxEnergeticPrimary.SetEnergy(0);
    for (auto const &particle:mctree->get_heads()) {
        if (maxEnergeticPrimary.GetEnergy() < particle.GetEnergy()) {
        maxEnergeticPrimary = particle;
        }
    }
 
    I3Particle tmp = maxEnergeticPrimary;
    while (mctree->number_of_children(tmp.GetID())) {
        if (mctree->number_of_children(tmp.GetID()) > 1) log_fatal("wasn't expecting more than one child.");
        if (AllowedPrimaries.find(tmp.GetType()) != AllowedPrimaries.end()) {
            speed_prof->push_back(tmp.GetSpeed() / I3Constants::c);
            tmp = mctree->first_child(tmp.GetID());
        } else {
            break;
        }
    }
}