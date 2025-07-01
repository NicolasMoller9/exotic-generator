/**
 * A Monopole/Qball/ShadowCharge Generator Module
 * Copyright (c) 2004-2014 IceCube Collaboration
 * Version $Id$
 *
 * @file I3ExoticGenerator.cxx
 * @date $Date$
 * @author jacobi
 * @author bchristy
 * @author olivas
 * @author flauber
 * @brief Implementation of a module to generate exotic particles (Monopoles, Qballs and Shadow charges)
 */

#include "exotic-generator/I3ExoticGenerator.h"

I3_MODULE(I3ExoticGenerator);

I3ExoticGenerator::I3ExoticGenerator(const I3Context &ctx) :
        I3Module(ctx),
        treeName_("I3MCTree"),
        infoName_("EXOInfoDict"),
        Nevents_(0),
        part_id_(-41), // monopole by default
        mass_(NAN),   // see validity limits on the mass in the CheckMass function in I3ExoticGeneratorUtils.cxx
        charge_(1),
        beta_(NAN),
        gamma_(NAN),
        speed_(NAN),
        energy_(NAN),
        useBetaRange_(NAN),
        diskDist_(1000. * I3Units::m),
        diskRad_(800. * I3Units::m),
        radOnDisk_(NAN),
        aziOnDisk_(NAN),
        startTime_(0.0),
        weight_(1.),
        powerLawIndex_(0),
        length_(-1.), //NAN is a valid value, using -1 to signal default value (2* diskDist_)
        totalweight_(0),
        exinfo_dict_config_(new I3MapStringDouble),
        betaRange_{NAN, NAN},
        zenithRange_{0.0 * I3Units::degree, 180.0 * I3Units::degree}, //min, max zenith
        azimuthRange_{0.0 * I3Units::degree, 360.0 * I3Units::degree}, //min, max azimuth
        shiftCenter_{0.0 * I3Units::meter, 0.0 * I3Units::meter, 0.0 * I3Units::meter},
        precalculated_betas_{},
        need_precalculate_(false)

{
  log_debug("Constructor I3ExoticGenerator");
  AddParameter("TreeName", "Name of the `I3MCTree` to write the generated particle to.", treeName_);
  AddParameter("InfoName", "Name of the particle info dictionary. The generator writes generation parameters into this dictionary.", infoName_);
  AddParameter("Nevents", "Number of events to be generated for the simulation.", Nevents_);
  AddParameter("PartID", "Particle identifier: -41 for Monopoles, -42 for Shadow charges, -43 for qballs.", part_id_);
  AddParameter("Disk_dist", "Distance of the generation disk from the center of IceCube. (Default: `1000 * I3Units::m`)", diskDist_);
  AddParameter("Disk_rad", "Radius of Generation Disk. (Default: `800 * I3Units::m`)", diskRad_);
  AddParameter("BetaRange", "Velocity range of the exotic as array of lower and upper boundary. If a fixed velocify is desired, set lower and upper boundary to the same velocity. Express the velocities as ratio of the speed of light, c. For example: `[0.1, 0.4]`.", betaRange_);
  AddParameter("PowerLawIndex", "If a power-law index is given, the velocities of the simulated exotics will follow a power-law distribution with the given index. If the given index is negative, the distribution will be inverted (i.e more statistics at high speed). If no power-law index is given, the velocities will be uniformly distributed (default).", powerLawIndex_);
  AddParameter("Mass", "Mass of the exotic particle. For slow monopoles and shadow charges, this parameter is optional. For fast monopoles, this parameter is needed to calculate energy losses along the trajectory. Example: `1e7 * I3Units::GeV`.", mass_);
  AddParameter("Charge", "Charge of the exotics. If monopole, this parameter is the n in g=n*gD, where gD is the Dirac charge and g the magnetic charge. Example: 2", charge_);
  AddParameter("Length", "Length of the exotic track. Can be `NaN` or any length. The default value is calculated to twice the disk distance. Set to `-1` to indicate that the default value, `2 * Disk_dist`, should be used.", length_);
  AddParameter("ZenithRange", "List of lower and upper zenith bound. Use this parameter to restrict the direction of the exotic. Example: `[0.0 * I3Units::deg, 180.0 * I3Units::deg]`", zenithRange_);
  AddParameter("AzimuthRange", "List of lower and upper azimuth bound. Use this parameter to restrict the direction of the exotic. Example: `[0.0 * I3Units::deg, 360.0 * I3Units::deg]`", azimuthRange_);
  AddParameter("Rad_on_disk", "Set the radius coordinate of the starting position on the generation disk. Randomized if `NaN`. Example: `5. * I3Units::m`", radOnDisk_);
  AddParameter("Azi_on_disk", "Set the azimuth coordinate of the starting position on the generation disk. Randomized if `NaN`. Example: `45. * I3Units::deg`", aziOnDisk_);
  AddParameter("ShiftCenter", "Shifts the exotic. This is useful to explore different geometries. To shift according to the center of DeepCore (IC86-I SLOP trigger only acts on DC), configure with `ShiftCenter = ([46.0 * icetray.I3Units.m, -34.5 * icetray.I3Units.m, -330.0 * icetray.I3Units.m])`.", shiftCenter_);
  AddParameter("StartTime", "The time, measured from the beginning of the event, the exotic particle should be started. Example: `0. * I3Units::s`", startTime_);
  AddOutBox("OutBox");

}

I3ExoticGenerator::~I3ExoticGenerator() {}

void I3ExoticGenerator::Configure() {
  GetParameter("NEvents", Nevents_);
  GetParameter("PartID", part_id_);
  GetParameter("TreeName", treeName_);
  GetParameter("Mass", mass_);
  GetParameter("Charge", charge_);
  GetParameter("Disk_dist", diskDist_);
  GetParameter("Disk_rad", diskRad_);
  GetParameter("Rad_on_disk", radOnDisk_);
  GetParameter("Azi_on_disk", aziOnDisk_);
  GetParameter("shiftCenter", shiftCenter_);
  GetParameter("InfoName", infoName_);
  GetParameter("Length", length_);
  GetParameter("BetaRange", betaRange_);
  GetParameter("ZenithRange", zenithRange_);
  GetParameter("AzimuthRange", azimuthRange_);
  GetParameter("PowerLawIndex", powerLawIndex_);
  GetParameter("StartTime", startTime_);
  if (Nevents_<= 0){
    log_fatal("Requested less than 1 event to be generated. Abort!");
  }

  /**
   * check that input parameters are sane and set some more variables
   */

  if (part_id_ == -41) {
        particleType_ = I3Particle::Monopole;
    } else if (part_id_ == -42) {
        particleType_ = I3Particle::ShadowCharge;
    } else if (part_id_ == -43) {
        particleType_ = I3Particle::Qball;
    } else {
        log_fatal("Invalid PartID: %d. Supported values are -41 (Monopole), -42 (ShadowCharge) and -43 (Qball).", part_id_);
    }

  I3ExoticGeneratorUtils::CheckMass(mass_, betaRange_, particleType_);
  I3ExoticGeneratorUtils::CheckAndSetBeta(betaRange_, powerLawIndex_);
  I3ExoticGeneratorUtils::CheckDisk(diskDist_, diskRad_);
  I3ExoticGeneratorUtils::CheckAndSetLength(length_, 2 * diskDist_);
  I3ExoticGeneratorUtils::CheckAndSetZenith(zenithRange_);
  I3ExoticGeneratorUtils::CheckAndSetAzimuth(azimuthRange_);
  I3ExoticGeneratorUtils::CheckShiftCenter(shiftCenter_);
  I3ExoticGeneratorUtils::CheckDiskPos(radOnDisk_, aziOnDisk_, diskRad_);

  if(betaRange_[0] == betaRange_[1]) {
    //cache result so we do not have to add it multiple times
    useBetaRange_ = false;
    beta_ = betaRange_[0];
    gamma_ = I3ExoticGeneratorUtils::beta2gamma(beta_);
    speed_ = beta_ * I3Constants::c;
    energy_ = mass_ * gamma_;
    weight_ = 1.;
    totalweight_ = Nevents_;
    (*exinfo_dict_config_)["Gamma"] = gamma_;
    (*exinfo_dict_config_)["Beta"] = beta_;
    (*exinfo_dict_config_)["LogEnergy"] = log10(energy_);
    (*exinfo_dict_config_)["Weight"] = weight_;
    (*exinfo_dict_config_)["OneWeight"] = weight_ / totalweight_;

  }else{
    useBetaRange_ = true;
    if (std::isnan(powerLawIndex_) || powerLawIndex_==0 ){
      totalweight_ = Nevents_;
      (*exinfo_dict_config_)["Weight"] = weight_;
      (*exinfo_dict_config_)["OneWeight"] = weight_ / totalweight_;
    }else{
      need_precalculate_=true;
    }
  }

  // create exotic-dict and save it
  (*exinfo_dict_config_)["PartID"] = part_id_;
  (*exinfo_dict_config_)["Mass"] = mass_;
  (*exinfo_dict_config_)["Charge"] = charge_;
  (*exinfo_dict_config_)["DiskRadius"] = diskRad_;
  (*exinfo_dict_config_)["DiskDistance"] = diskDist_;
  (*exinfo_dict_config_)["ZenithMin"] = zenithRange_[0];
  (*exinfo_dict_config_)["ZenithMax"] = zenithRange_[1];
  (*exinfo_dict_config_)["AzimuthMin"] = azimuthRange_[0];
  (*exinfo_dict_config_)["AzimuthMax"] = azimuthRange_[1];
  //(*exinfo_dict_config_)["SolidAngle"] = I3ExoticGeneratorUtils::CalcSolidAngle(zenithRange_[0], zenithRange_[1], azimuthRange_[0], azimuthRange_[1]);
  (*exinfo_dict_config_)["BetaMin"] = betaRange_[0];
  (*exinfo_dict_config_)["BetaMax"] = betaRange_[1];
  (*exinfo_dict_config_)["powerLawIndex"] = powerLawIndex_;
}


void I3ExoticGenerator::DAQ(I3FramePtr frame) {
  I3RandomServicePtr random = GetService<I3RandomServicePtr>();
  //in case of powerlaw set, we need to calculate all betas in advance to properly calculate OneWeight
  if(need_precalculate_) {
    need_precalculate_ = false;
    for (int i = 0; i < Nevents_; ++i) {
      beta_ = I3ExoticGeneratorUtils::RandomPowerLawSampled(random, betaRange_[0], betaRange_[1], powerLawIndex_);
      if (powerLawIndex_ > 0) {
        totalweight_ += pow(beta_, powerLawIndex_);
      }
      else {
        totalweight_ += pow((betaRange_[0]+betaRange_[1])-beta_, -powerLawIndex_);
      }
      precalculated_betas_.push_back(beta_);
    }
  }
  --Nevents_;
  I3MapStringDoublePtr exinfo(new I3MapStringDouble(*exinfo_dict_config_));
  if (!random) {
    log_fatal("Failed to Get Random Service.");
  }

  const double zenith = zenithRange_[0] == zenithRange_[1] ? zenithRange_[0] : I3ExoticGeneratorUtils::RandomCosinusSampled(random, zenithRange_[0] / I3Units::radian, zenithRange_[1] / I3Units::radian);
  const double azimuth = azimuthRange_[0] == azimuthRange_[1] ? azimuthRange_[0] : I3ExoticGeneratorUtils::RandomUniformSampled(random, azimuthRange_[0] / I3Units::radian, azimuthRange_[1] / I3Units::radian);
  I3Direction mpDir(zenith, azimuth);

  (*exinfo)["Zenith"] = zenith;
  (*exinfo)["Azimuth"] = azimuth;

  //r and theta on creation disk
  const double r = !std::isnan(radOnDisk_) ? radOnDisk_ : I3ExoticGeneratorUtils::RandomCircularSampled(random, 0 / I3Units::m, diskRad_ / I3Units::m);
  const double theta = !std::isnan(aziOnDisk_) ? aziOnDisk_ : I3ExoticGeneratorUtils::RandomUniformSampled(random, 0, 2 * I3Constants::pi);
  I3Position mpPos(r * cos(theta) * I3Units::m, r * sin(theta) * I3Units::m, diskDist_  * I3Units::m);

  (*exinfo)["OnDiskRadius"] = r;
  (*exinfo)["OnDiskAzimuth"] = theta;


  //rotate to get the proper arrival direction
  mpPos.RotateY(mpDir.GetZenith());
  mpPos.RotateZ(mpDir.GetAzimuth());

  //apply the shift
  mpPos.SetX(mpPos.GetX() + shiftCenter_[0]);
  mpPos.SetY(mpPos.GetY() + shiftCenter_[1]);
  mpPos.SetZ(mpPos.GetZ() + shiftCenter_[2]);

  // Randomize speed
  if (useBetaRange_) {
    if((std::isnan(powerLawIndex_) || powerLawIndex_==0 )){
      beta_ = I3ExoticGeneratorUtils::RandomUniformSampled(random, betaRange_[0], betaRange_[1]);
    }else{
      beta_ = precalculated_betas_[Nevents_];
      if (powerLawIndex_ > 0) {
        weight_ = pow(beta_, powerLawIndex_); //probability for this beta to show up
      }
      else {
        weight_ = pow((betaRange_[0]+betaRange_[1])-beta_, -powerLawIndex_);
      }
      (*exinfo)["Weight"] = weight_;
      (*exinfo)["OneWeight"] = weight_ / totalweight_;
    }
    gamma_ = I3ExoticGeneratorUtils::beta2gamma(beta_);
    speed_ = beta_ *I3Constants::c;
    energy_ = mass_ * gamma_;
    (*exinfo)["Gamma"] = gamma_;
    (*exinfo)["Beta"] = beta_;
    (*exinfo)["LogEnergy"] = log10(energy_);
  }

  I3Particle particle;
  
  particle.SetType(particleType_);
  particle.SetLocationType(I3Particle::InIce);
  particle.SetPos(mpPos);
  particle.SetDir(mpDir);
  particle.SetEnergy(energy_);
  particle.SetSpeed(speed_);
  particle.SetTime(startTime_);
  particle.SetLength(length_);

  I3MCTreePtr mptree(new I3MCTree);
  I3MCTreeUtils::AddPrimary(*mptree, particle);

  //add tree and info dict to frame
  frame->Put(treeName_, mptree);
  frame->Put(infoName_, exinfo);
  PushFrame(frame);
  if(Nevents_ <= 0){
    RequestSuspension();
  }
}

void I3ExoticGenerator::Finish() {
}


