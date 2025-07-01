
/**
 * class: I3ExoticEnergyLoss.cxx
 * Copyright (c) 2008 IceCube Collaboration
 * Version $Id: I3ExoticPropagator.cxx 141809 2016-02-13 17:31:29Z anna.obertacke $
 *
 * Date 06 Feb 2008
 * @version $Revision: 141809 $
 * @date $Date: 2016-02-13 11:31:29 -0600 (Sa, 13. Feb 2016) $
 * @author Brian Christy <bchristy@icecube.umd.edu>
 * @author Alex Olivas <olivas@icecube.umd.edu>
 *
 * @brief A module to calculate energy loss of exotic particles in the ice
 *
 */

#include "I3ExoticEnergyLoss.h"

#include <cmath>                 // for sqrt, pow, log, etc.
#include <dataclasses/I3Constants.h>
#include "icetray/I3Frame.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"


/**
 * Responsible for taking a given segment of the Monopole and calculating
 * the energy loss due to Ionization (with Density Correction)
 * @param sBeta               Starting beta value of monopole
 * @param tLength             Length of segment to calculate energy loss over
 * @param n                   The n factor in the definition of the magnetic charge (g=n*gD)
 * @return the energy loss for the given segment and starting speed
 */
double I3ExoticEnergyLoss::CalculateMonopoleEnergyLoss(double sBeta,
                                                        double tLength,
                                                        double n) {
    log_debug("Entering CalculateMonopoleEnergyLoss");
    double dEdX=0;
    if (sBeta>=0.99995) {
        log_warn("So far, energy calculation only does Ionization Loss; For a monopole of beta %f," 
                " stochastic's become important; Make sure you know what you're doing", sBeta);
        return 0;
    } else if (sBeta >= 0.05) {    // ionisation
        dEdX = BetheBlock(sBeta, n, true);
    } else if (sBeta <= 0.01) {    // electronic excitation
        double mag_corr = 0.25 * pow(n,2);
        dEdX = SlowIonInIce(sBeta, 1)*mag_corr;
    } else if ((0.01<sBeta) && (sBeta<0.05)) {  // linear interpolation
        double mag_corr = 0.25 * pow(n,2);
        double dEdX_005 = BetheBlock(0.05, n, true);
        double dEdX_001 = SlowIonInIce(0.01, 1)*mag_corr;
        double f = (log10(dEdX_005) - log10(dEdX_001)) / (log10(0.05) - log10(0.01));
        dEdX = dEdX_001*pow((sBeta/0.01),f);
        log_debug("dE/dX (interpolation): %f GeV/m", dEdX);
    }
    return dEdX * tLength;
}



/**
 *  -For electric charges, used formula from PDG 
 *  -For magnetic charges, used formula Eq 16 of Ahlen Phys Rev D Vol 17 # 1
 *   "Stopping Power formula for magnetic monopoles"
 *   and used values of K and B found in: J. Derkaoui et al. - Astroparticle Physics 9 (1998) 173-183
 *
 *  For Density Correction constants:
 *  Used values from Table C.2 in Dima's Thesis on mmc, which in turn came from
 *  https://web.archive.org/web/20070713230354/http://pdg.lbl.gov/AtomicNuclearProperties/substances/276.html
*/
/**
 * Responsible for taking a given segment of the exotic particle and calculating
 * the energy loss due to Ionization (with Density Correction)
 * @param sBeta               Starting beta value of the exotic
 * @param z                   The electric charge, or the n factor in the definition of the magnetic charge (g=n*gD)
 * @param magnetic            True to calculate the energy loss of a magnetic charge instead of an electric one
 * @return the energy loss for the given segment and starting speed
 */
double I3ExoticEnergyLoss::BetheBlock(double sBeta, double z, bool magnetic) {
    // Physical constants
    const double MASS_E = 5.11e-4;     // GeV/c^2
    const double ION_LOSS = 7.5e-8;    // GeV
    const double DENSITY = 9.17e5;     // g/m^3
    const double Zice = 10;
    const double Aice = 18.015;        // g/mol
    const double NA = 6.02214076e23;   //I3Constants::NA; // mol^-1
    const double Ne = NA * DENSITY * Zice / Aice; //m^-3
    const double CHARGE_E_SQRD = 1.439964459e-18; // GeV*m

    // Take care of magnetic vs electronic case
    double log_BB = 0;
    double QED_CORR = 0;
    double BLOCH_CORR = 0;
    if (magnetic==true) {
        const double n = z;
        const double gD = 137./2.; 
        z = n*gD*sBeta;

        // logarithmic term
        log_BB = log((2 * MASS_E * pow(sBeta, 2)) / (ION_LOSS * (1 - pow(sBeta, 2)))) - 0.5;

        // QED and Bloch corrections
        if (n == 1) {
            QED_CORR = 0.406; BLOCH_CORR = 0.248;
        } else if (n == 2) {
            QED_CORR = 0.346; BLOCH_CORR = 0.672;
        } else if (n == 3 || n == 4 || n == 5) {
            QED_CORR = 0.3; BLOCH_CORR = 1.022;
        } else if (n == 6 || n == 7 || n == 8) {
            QED_CORR = 0.3; BLOCH_CORR = 1.685;
        } else if (n >= 9) {
            QED_CORR = 0.3; BLOCH_CORR = 2.085;
        }
    } else {
        const double Tmax = 2 * MASS_E * pow(sBeta, 2) / (1 - pow(sBeta, 2)); // in the approximation me/M ~ 0
        log_BB = 0.5*log( (2 * MASS_E * pow(sBeta, 2) * Tmax) / (pow(ION_LOSS, 2) * (1 - pow(sBeta, 2))) ) - pow(sBeta, 2);
    }

    // BetheBloch coefficient
    const double COEFF1 = (4 * I3Constants::pi * Ne) / pow(sBeta,2); // m^-3
    const double COEFF2 = pow(z * CHARGE_E_SQRD, 2); // GeV^2·m^2
    const double COEFF = COEFF1 * COEFF2 / MASS_E; // GeV/m

    // Density correction
    double denCorr = 0.0;
    double sGamma = 1 / (sqrt(1 - pow(sBeta, 2)));
    double X = log10(sBeta * sGamma);
    const double A_DC = 0.09116;
    const double X0 = 0.240;
    const double X1 = 2.8004;
    const double M_DC = 3.477;
    const double C_BAR = 3.5017;
    if (X > X0 && X < X1) {
        denCorr = log(pow(sBeta * sGamma, 2)) - C_BAR + A_DC * pow((X1 - X), M_DC);
    } else if (X > X1) {
        denCorr = log(pow(sBeta * sGamma, 2)) - C_BAR;
    }

    // Final formula
    double dEdX = COEFF * ( log_BB - (denCorr/2.0) + (QED_CORR/2.0) - BLOCH_CORR );

    log_debug("Density Correction: %f", denCorr);
    log_debug("dE/dX (BetheBlock_MM): %f GeV/m", dEdX);
    return dEdX;
}


/**
 * Responsible for taking a given segment of the exotic particle and calculating
 * the energy loss due to excitation
 * @param sBeta               Starting beta value of the exotic
 * @param z                   The electric charge
 * @return the energy loss for the given segment and starting speed
 */
double I3ExoticEnergyLoss::ELEC_LOSS(double sBeta, double z) {
    // ELECTRIC ENERGY LOSS
    const double A0 = 5.29177210544e-11;  // m 
    const double CHARGE_E_SQRD = 1.439964459e-18; // GeV*m
    const double ALPHA_MINUS1 = 137; 
    const double DENSITY = 9.17e5;    // g/m^3
    const double Zice = 10.0/3.0;
    const double Aice = 18.015/3.0;   // g/mol
    const double NA = 6.02214076e23;   //I3Constants::NA; // mol^-1
    const double BETA0 = 7e-4;
    const double EL_DENSITY = NA*DENSITY*Zice/Aice; // m^-3 

    double COEFF1=8*I3Constants::pi*A0*CHARGE_E_SQRD*sBeta*ALPHA_MINUS1;     // GeV*m2 
    double COEFF2=pow(z,7.0/6.0)*EL_DENSITY/pow(pow(z,2.0/3.0)+pow(Zice,2.0/3.0),3.0/2.0); // m-3
    double corr=1-exp(-pow(sBeta/BETA0,2)); 

    double dEdX = COEFF1*COEFF2*corr;  // GeV/m

    log_debug("dE/dX (electric): %f GeV/m", dEdX);
    return dEdX;
}



/**
 * Responsible for taking a given segment of the exotic particle and calculating
 * the energy loss due to nuclear recoil
 * @param sBeta               Starting beta value of the exotic
 * @param z                   The electric charge
 * @param Z_NUC               The atomic number of the nucleus
 * @param M_NUC               The mass of the nucleus in GeV
 * @return the energy loss for the given segment and starting speed
 */
double I3ExoticEnergyLoss::NUCL_LOSS(double sBeta, double z, int Z_NUC, double M_NUC) {
    const double A0 = 5.29177210544e-11;  // m 
    const double CHARGE_E_SQRD = 1.439964459e-18; // GeV*m
    const double DENSITY = 9.17e5;    // g/m^3
    const double Aice = 18.015;        // g/mol
    const double NA = 6.02214076e23;   //I3Constants::NA; // mol^-1
    const double NUCL_DENSITY = NA*DENSITY/Aice; // m^-3 

    double sGamma = 1 / (sqrt(1 - pow(sBeta, 2)));
    double Tm = 4 * M_NUC * (sGamma - 1);              // GeV
    double a = 0.8854*A0 / (pow(z,0.23) + pow(Z_NUC,0.23));  // m
    double epsilon = ((sGamma-1)*a*M_NUC) / (z*Z_NUC*CHARGE_E_SQRD); // no unit
    double MAIN = I3Constants::pi * pow(a,2) * Tm * NUCL_DENSITY / epsilon;   // GeV/m
    double S_eps;
    if (epsilon<=30){
        S_eps = 0.5*log(1+1.1383*epsilon) / (epsilon+(0.01321*pow(epsilon,0.21226))+(0.19593*pow(epsilon,0.5))); // no units
    }else{
        S_eps = 0.5*log(epsilon) / epsilon;
    }
    double dEdX = MAIN * S_eps;

    log_debug("dE/dX (nuclear Z=%d): %f GeV/m", Z_NUC, dEdX);
    return dEdX;
}


/**
 * Responsible for taking a given segment of the exotic particle and calculating
 * the total energy loss at low speed
 * @param sBeta               Starting beta value of the exotic
 * @param z                   The electric charge
 * @return the energy loss for the given segment and starting speed
 */
double I3ExoticEnergyLoss::SlowIonInIce(double sBeta, double z) {
    const double MH = 0.93878;    // GeV/c2
    const double MO = 14.902973;  // GeV/c2

    double elec_loss = ELEC_LOSS(sBeta, z);
    double nucl_loss = (2*NUCL_LOSS(sBeta, z, 1, MH) + NUCL_LOSS(sBeta, z, 8, MO)) / 3;

    double dEdX = elec_loss + nucl_loss;

    log_debug("dE/dX (slow ion): %f GeV/m", dEdX);
    return dEdX;
}