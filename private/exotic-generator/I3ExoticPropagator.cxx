/**
 * class: I3ExoticPropagator.cxx
 * Copyright (c) 2008 IceCube Collaboration
 * Version $Id: I3ExoticPropagator.cxx 141809 2016-02-13 17:31:29Z anna.obertacke $
 *
 * Date 06 Feb 2008
 * @version $Revision: 141809 $
 * @date $Date: 2016-02-13 11:31:29 -0600 (Sa, 13. Feb 2016) $
 * @author Brian Christy <bchristy@icecube.umd.edu>
 * @author Alex Olivas <olivas@icecube.umd.edu>
 *
 * @brief A module to convert exotic particles into a chain of propagated
 * @brief particles through the ice
 *
 */

#include "exotic-generator/I3ExoticPropagator.h"
#include "I3ExoticPropagatorUtils.h"

#include "icetray/I3Frame.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"


I3_MODULE(I3ExoticPropagator);

I3ExoticPropagator::I3ExoticPropagator(const I3Context &ctx) :
        I3Module(ctx),
        inputTreeName_("I3MCTree"),      
        outputTreeName_("I3MCTree"),     
        infoName_("EXOInfoDict"),             

        // if catalyse cascades
        meanFreePath_(-1),           // default value = -1 => calculated from cross section 
                                     // => if cross section is NAN, mean free path will be NAN and no catalysis process will occur
                                     // to switch on catalysis, give a value to sigma0 or a positive value to meanFreePath_
        sigma0_(NAN),                // exemple: sigma=1e-31 m^2 (upper value for proton decay catalysis by monopoles)       
        scaleEnergy_(false),            
        energyScaleFactor_(1.0),        
        useCorrectDecay_(false),
        keepPrimary_(true),

        // for primary track
        maxdistfromcenter_(1300 * I3Units::m),  // only used when the primary track is split in smaller steps
        calcEn_(false),                               
        stepSize_(NAN),                  
        minlength_(NAN),    // ex: 0.001 * I3Units::m       
        maxlength_(NAN),     // ex: 10.0 * I3Units::m
        profiling_(false),          
        checkParticle_(true){

    log_debug("Constructor I3ExoticPropagator");

    AddParameter("InputTreeName", "Name of the `I3MCTree` containing the exotic generated by the exotic generator.", inputTreeName_);
    AddParameter("OutputTreeName", "Name of the `I3MCTree` to write the exotic into after propagation.", outputTreeName_);
    AddParameter("InfoName", "Name of the exotic info dictionary, containing all necessary information about the generation parameters.", infoName_);

    // Parameters for exotics producing cascades along their track:    
    AddParameter("MeanFreePath", "Mean free path (lambda) between catalyzed process. Example: `1 * I3Units::m`", meanFreePath_);
    AddParameter("CrossSection", "Cross section (sigma) of the catalyzed process (in the case of monopoles, it is sigma0, not sigma_cat). Example: `1e-31 * I3Units::m * I3Units::m`", sigma0_);
    AddParameter("ScaleEnergy", "Whether to set the mean free path (lambda) to 1 meter and scale up the energy by 1/lambda. The overall light output stays comparable, while the number of secondary particles in the `I3MCTree` is reduced. This saves computing resources especially for short mean free paths.", scaleEnergy_);
    AddParameter("EnergyScaleFactor", "Scale down the cascade energy in order to test the influence of other decay channels.", energyScaleFactor_);
    AddParameter("UseCorrectDecay", "For monopoles: Whether to simulate back-to-back positrons (460 MeV) and neutral pions (480 MeV) instead of just one positron. The overall light output is similar, but correct decay has twice as much secondary particles in the `I3MCTree`. This option cannot be used together with `ScaleEnergy` or `EnergyScaleFactor`, since the energies are hard coded.", useCorrectDecay_);
    AddParameter("KeepPrimary", "Wether to keep the exotic track in the detector, allowing it to produce light if possible (ex: monopole luminescence). If `false`, the track-length will be set to 0.", keepPrimary_);

    // Parameters for the primary exotic particle track:
    AddParameter("MaxDistanceFromCenter", "How far beyond the detector to add small track segments when propagating a track-like exotic. If the start of the exotic is further from the detector than this value, the propagator will IGNORE the parameter and propagate until it reaches the same distance away on far side of detector. Example: `800 * I3Units::m`", maxdistfromcenter_);
    AddParameter("CalculateEnergy", "If set to `true`, the energy loss of the exotic is calculated during the propagation, and its speed updated along the track.", calcEn_);
    AddParameter("StepSize", "Length of exotic track segments. If set this will override `MinLength` and `MaxLength`. Otherwise, `MinLength` and `MaxLength` are used to set the lower and upper bounds on the track segment lengths. Example: `1 * I3Units::m`", stepSize_);
    AddParameter("MinLength", "Assuming `StepSize` is `NaN`, this represents the smallest segment the propagator will generate. Example: `0.001 * I3Units::m`", minlength_);
    AddParameter("MaxLength", "Assuming stepsize is `NaN`, this represents the largest segment the propagator will generate. Example: `10 * I3Units::m`", maxlength_);
    AddParameter("Profiling", "If `true`, adds a profile (type `I3VectorDouble`) of the exotic speed for each track segment to the frame.", profiling_);
    AddParameter("ParticleCheck", "If `true`, each segment of the track will undergo some sanity check.", checkParticle_);

    AddOutBox("OutBox");
}

I3ExoticPropagator::~I3ExoticPropagator() { }


void I3ExoticPropagator::Configure() {
    log_debug("Configure I3ExoticPropagator");

    GetParameter("InputTreeName", inputTreeName_);
    GetParameter("OutputTreeName", outputTreeName_);
    GetParameter("InfoName", infoName_);

    GetParameter("MeanFreePath", meanFreePath_);
    GetParameter("CrossSection", sigma0_);
    GetParameter("ScaleEnergy", scaleEnergy_);
    GetParameter("EnergyScaleFactor", energyScaleFactor_);
    GetParameter("UseCorrectDecay", useCorrectDecay_);
    GetParameter("KeepPrimary", keepPrimary_);

    GetParameter("MaxDistanceFromCenter", maxdistfromcenter_);
    GetParameter("CalculateEnergy", calcEn_);
    GetParameter("StepSize", stepSize_);
    GetParameter("MinLength", minlength_);
    GetParameter("MaxLength", maxlength_);
    GetParameter("Profiling", profiling_);
    GetParameter("ParticleCheck", checkParticle_);

    /**
     * Check that input parameters are sane
     */
    if (minlength_ > maxlength_) {
        log_fatal("Oops, MaxLength<MinLength.  Did you switch them?");
    }
    if (std::isnan(meanFreePath_) || ((meanFreePath_ <= 0.0) && (meanFreePath_ != -1))) {
        log_fatal("Mean Free Path must be a positive number or -1 (fixed cross section)");
    }
    if (energyScaleFactor_ > 1.0 || energyScaleFactor_ <= 0) {
        log_fatal("The vaild range for the energy scale factor is between > 0 and <= 1.0, since then the cascade energy cannot be zero or negative and not be greater than the proton mass (rest energy). ");
    }
    if (useCorrectDecay_ && scaleEnergy_) {
        log_fatal("Correct decay and scale energy cannot be used together, since correct decay depends on a hard coded energy.");
    }
    if (useCorrectDecay_ && energyScaleFactor_ != 1.0) {
        log_fatal("Correct decay and energy scale factor cannot be used together, since correct decay depends on a hard coded energy.");
    }

    /**
    TODO rest
    */
}

void I3ExoticPropagator::DAQ(I3FramePtr frame) {
    log_debug("Entering DAQ I3ExoticPropagator");

    // get generated tree (ouput of I3ExoticGenerator)
    I3MCTreeConstPtr mctree = frame->Get<I3MCTreeConstPtr>(inputTreeName_);
    if(!mctree){
      log_fatal("Tree not found in frame, Aborting!");
    }
    // create new tree for propagated exotic
    I3MCTreePtr prop_tree(new I3MCTree(*mctree));  // create an independent copy of mctree

    // Extract some parameters from the info dictionnary created by I3ExoticGenerator
    if (!frame->Has(infoName_)) {
        std::cerr << "Error: Dictionnary '" << infoName_ << "' not found in the frame." << std::endl;
        return;
    }
    const I3MapStringDouble& exinfo = frame->Get<I3MapStringDouble>(infoName_);  // find the parameters used for the exotic generation in the EXOInfoDict
    double beta = exinfo.at("Beta");
    double charge_ = static_cast<double>(exinfo.at("Charge")); 
    double mass_ = exinfo.at("Mass");
    int pID = static_cast<int>(exinfo.at("PartID"));
    if (beta <= 0.0 or beta >= 1.0 or std::isnan(beta)) {
        log_fatal("Got beta from %s frame, but value is not sane: %g", infoName_.c_str(), beta);
    }

    // Get random service
    I3RandomServicePtr random = context_.Get<I3RandomServicePtr>();
    if (!random) {
        log_fatal("Failed to Get Random Service");
    }

    double Mean_Free_Path = meanFreePath_; 
    if (meanFreePath_==-1) {  // the mean free path value will be calculated from the cross section
        Mean_Free_Path = CalculateMeanFreePath(pID,sigma0_,beta);  // in m
        log_debug("Calculated Lambda=%f m from sigma=%f m^2 and beta=%f", Mean_Free_Path, sigma0_, beta);
    }

    double deltaEnergy_; // Energy freed by the catalysis process (if any)
    if (std::isnan(Mean_Free_Path)) {
        log_debug("No catalysis process will be simulated")
    } else {
        if (pID == -41) { // for monopoles
            deltaEnergy_ = 0.94 * I3Units::GeV * energyScaleFactor_;  // = the proton total energy, except if using the useCorrectDecay option
        }

        if (scaleEnergy_) {
            if (Mean_Free_Path >= 1.0 * I3Units::m) {
                log_error("Cannot scale energy down! Mean free path (%f) has to be smaller than 1m\n", Mean_Free_Path);
                log_error("Scale Energy turned off");
            }
            else { // scaling energy to 1 m mean free path
                log_info("Scaling Energy/MFP up");
                deltaEnergy_ = ((1.0 * I3Units::m) / Mean_Free_Path) * deltaEnergy_;
                Mean_Free_Path = 1.0 * I3Units::m;
            }
        }
        log_debug("Energy of catalyzed cascades %f (GeV)", deltaEnergy_);
        log_debug("Mean Free Path for catalysis process %f (m)", Mean_Free_Path);
    }
        
    // Update the mean free path value in the dictionary
    I3MapStringDouble exinfo_new = exinfo;
    exinfo_new["MeanFreePath"] = Mean_Free_Path;
    exinfo_new["Sigma0"] = sigma0_;
    frame->Delete(infoName_);
    frame->Put(infoName_, boost::make_shared<I3MapStringDouble>(exinfo_new));

    I3ExoticPropagatorUtils::PropagateExotic(frame,
                                            mctree,
                                            prop_tree,
                                            random,
                                            // params
                                            Mean_Free_Path,
                                            useCorrectDecay_,
                                            deltaEnergy_,
                                            maxdistfromcenter_,
                                            keepPrimary_,
                                            calcEn_,        
                                            mass_,        
                                            charge_,
                                            profiling_,
                                            stepSize_,
                                            minlength_,
                                            maxlength_,
                                            checkParticle_);


    //Remove old tree and replace with new, propagated one
    if (inputTreeName_ == outputTreeName_) {
        frame->Delete(inputTreeName_);
    }
    frame->Put(outputTreeName_, prop_tree);
    PushFrame(frame, "OutBox");
}//End DAQ


// Function to calculate the value of the mean free path of a catalysis process depending on the cross section
double I3ExoticPropagator::CalculateMeanFreePath(int partID, double sig0, double b) { 
    double meanFreePath = std::numeric_limits<double>::quiet_NaN();  // Défaut = NaN
    
    // MONOPOLE
    if (partID==-41) {
        // Constants
        const double NA = 6.022140857e23;   // mol^-1
        const double density = 9.17e5;      // g/m^3
        const double Aice = 18.015;         // g/mol (molar mass of ice)
        const double n_nucl = 18 * NA * density / Aice;  // m-3
        const double Renu_H = -0.5;         // from [1]
        const double Renu_O = 3.123/2.0;    // from [1]
        const double beta0_H = 0.173341;    // from [1]
        const double beta0_O = 0.004334;    // from [1]
        double fH;
        double fO;

        if (b >= beta0_H) {
            fH = 1.0;
            fO = 1.0;
        } else if (beta0_O <= b && b < beta0_H) {
            fH = pow(b / beta0_H, 2 * Renu_H);
            fO = 1.0;
        } else { // beta < beta0_O
            fH = pow(b / beta0_H, 2 * Renu_H);
            fO = pow(b / beta0_O, 2 * Renu_O);
        }
    
        double F = (2.0 / 18.0) * fH + (16.0 / 18.0) * fO;
        double sigma_cat = (sig0 / b) * F;
        meanFreePath = 1.0 / (sigma_cat * n_nucl);
    }
    if (std::isnan(meanFreePath) && !std::isnan(sig0)) {
        std::cerr << "Fatal error: meanFreePath is NaN! sigma0 was: " << sig0 << " and beta = " << b << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return meanFreePath;
}



void I3ExoticPropagator::Finish() {
}


// References:
// [1]: J. Arafune &  M. Fukugita, "Velocity-Dependent Factors for the Rubakov Process for Slowly Moving Magnetic Monopoles in Matter", PHYSICAL REVIEW LETTERS VOL50 NBR24 (1983)