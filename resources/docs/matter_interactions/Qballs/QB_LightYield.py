# Authors: Nicolas Moller & Hiroki Toribuchi (Chiba U.)

import numpy as np
import math
import matplotlib.pyplot as plt


# 0) USEFUL CONSTANTS AND FUNCTIONS ##################################################################################################################################

# ShadowCharge parameters
beta_min, beta_max = [1e-4,1e-2]         # Qball beta expected range [3]

# Conversion factors:
MeV_to_m = 1.97326e-13                   # m*MeV
MeV_to_g = 1.78266192e-27                # g/(MeV/c2)
J_to_MeV = 6241509647120.4               # MeV/J (conversion factor)
statC_to_MeV = np.sqrt(1e-9 * J_to_MeV)  # (MeV*m)^1/2 / statC

# Useful constants:
pi = np.pi
e_cgs = 4.80320425*1e-10                 # statC = g^1/2 * cm^3/2 * s^-1 (Elementary charge)
e = e_cgs*statC_to_MeV                   # (MeV*m)^1/2 (Elementary charge)
NA = 6.022*1e23                          # mol^-1 (Avogadro number)
alpha = 1/137                            # (fine structure constant)
a0 = 5.29177210544e-9                    # cm (Bohr radius)
c = 2.99792458e8                         # m/s (speed of light)
m_e = 0.511                              # MeV (electron rest energy)
vBohr = 2.188*1e6                        # m/s (Bohr velocity: speed of an electron in the first orbit of a hydrogen atom)

# Ice properties:
rho_ice = 0.917                          # g/cm^3 (Ice density. From [1])
ZH = 1                                   # (atomic number of hydrogen)
AH = 1.008                               # g/mol (molar mass of hydrogen)
MH = 938.78                              # MeV/c2 (hydrogen mass)                
ZO = 8                                   # (atomic number of oxygen)
AO = 15.999                              # g/mol (molar mass of oxygen)
MO = 14902.973                           # MeV/c2 (oxygen mass) 
#Z_H2O = 7.42                            # Effective number of electron per ice atom (from [2])
Z_H2O = 10                               # (number of electron per H2O)
A_H2O = 18.015                           # g/mol (molar mass of H2O)
n_H2O = rho_ice/(A_H2O/NA)               # cm^-3 (number density of H2O molecule)
ne_H2O = Z_H2O*n_H2O                     # cm-3 (electron density)

# Useful functions
def gamma_(beta):                        # gamma factor
    return 1./np.sqrt(1. - beta**2)           

# (OPTIONAL) Get Qball parameters (Qcharge, mass, radius) depending on the Qball type (gauge-mediation, new-type)
def MRQ(type):
    elips = 1             # ellipticity of the field orbit (1=>circular orbit, 0.1=>oblate case) (from [3])
    mgrav = 1             # MeV (gravitino mass. 1 is a chosen value, it ranges from keV to GeV or even smaller [4])
    dzeta = 1             # O(1) parameter from [3]
    g = 1                 # gauge coupling of the standard model (assumption: g=1 as in Fig.2.3 of [5])
    Kabs = 0.1            # coefficient of the one-loop correction (assumed value from Fig.2.5 in [5])
    #MP = 2.4e21          # MeV (reduced Planck mass)
    MF = 4e7 * np.sqrt(g) # MeV (scale of gauge-mediated SUSY breaking (lower limit) [4])
    TRH = 1e6             # MeV (reheating temperature of the universe after inflation (minimum value for the thermally produced gravitinos to be dark matter, maximum value = 1e7 MeV [4])
    betaG = 6e-4*elips    # dimensionless parameter related to the ellipticity of the field orbit [3]
    phi = 5.8e15 * (dzeta/2.5)**(-1/2) * (betaG/6e-4)**(1/8) * (TRH/1000)**(-1/2) # MeV (amplitude of the field at the beginning of the oscillation [4])

    if type=='gauge':
        Q = betaG*(phi/MF)**4                       # Qball Q-charge [4]
        M = (4*np.sqrt(2)*pi/3) * dzeta*MF*Q**(3/4) # MeV (Qball mass [4])
        R = Q**(1/4)/(np.sqrt(2)*dzeta*MF)          # MeV-1 (never exceding the Bohr radius (5e-11 m) in the intresting range [1])
    elif type =='new':
        Q = 0.02*(phi/mgrav)**2                     # [4]
        M = mgrav*Q                                 # [4]
        R = Kabs**(-1/2) / mgrav                    # [4]
    return M, R*MeV_to_m, Q  # MeV, m, /


# 1) ENERGY LOSSES ##################################################################################################################################

def SlowIon_ElecLoss(beta,z,Z=Z_H2O/3,ne=ne_H2O):  # energy loss of a non-relativistic ion when interacting with atomic electrons
    beta_max = min(vBohr * z**(2/3) / c,1e-2)
    if beta <= beta_max:  # [6]*
        main = (8*pi*a0*e_cgs**2)*(beta/alpha)*(z**(7/6)/(z**(2/3)+Z**(2/3))**(3/2))*ne  # in g * cm * s^-2 (cgs system) [7]
        main *= 1e-5 # in kg * m * s^-2 = J/m
        main *= J_to_MeV # in MeV/m
        beta0 = 7e-4        # Fit parameter from [8] (indicates a reduced efficiency for electronic excitation when beta is lower than beta0)
        correction_factor = 1 - math.exp(-(beta/beta0)**2)  # [8]
        return main*correction_factor # in MeV/m
    else:
        return 0
# *when v<v0 (ion slower than orbital electron), the energy tranfer to the electron is low, which means it can be more simply described by the Lindhard formalism
# and Lindhard formalism includes a screening effect based on the Thomas-Fermi model of the atom, but when v>v0, the screening is different and can't be described by this model anymore


def nuclear_loss(beta,z,Z_i,M_i):  # energy loss of a non-relativistic ion when interacting with atomic nuclei (from [9])
    gam_ = gamma_(beta)
    a = 0.8854*a0*1e-2 / (z**0.23 + Z_i**0.23)  # m
    epsilon = ((gam_-1)*a*M_i) / (z*Z_i*e**2)   # no unit
    Tm = 4 * M_i * (gam_ - 1)                   # MeV
    main = pi * a**2 * Tm * n_H2O*1e6 / epsilon # MeV/m
    if epsilon<=30:
        S_eps = 0.5*np.log(1+1.1383*epsilon) / (epsilon+(0.01321*epsilon**0.21226)+(0.19593*epsilon**0.5))    # no units (from table 2 of [9])
    else:
        S_eps = 0.5*np.log(epsilon) / epsilon
    return main * S_eps                         # MeV/m
def nuclear_loss_full(beta,z):   # total nuclear energy loss in ice
    return (2*nuclear_loss(beta,z,ZH,MH)+nuclear_loss(beta,z,ZO,MO))/3


def thermal_shock(beta): # [?]
    R = a0            # m
    A = pi * R**2     # m^2
    main = A * rho_ice*1e3 * (beta*c)**2  # J/m 
    main *= J_to_MeV  # MeV/m
    return main       # MeV/m


def KKST(beta,z,type):  # (eq.18 in [10])
    Rq = MRQ(type)[1]        # m
    Rq_min = z*e_cgs**2*1e-6/(0.5*m_e*beta**2*c**2) # m (minimal impact parameter for Coulomb scattering: E_{coulomb}(Rq_min)!=E_{kin} )
    if Rq >= Rq_min:
        dzet = rho_ice*1000  # MeV
        sigma = pi*Rq**2     # m^2
        return sigma*(3*n_H2O*1e6)*dzet             # MeV/m  
    else:    # coulomb barrier too strong for KKST process to happen
        return 0


# 3) LIGHT YIELD ########################################################################################################################################################################

def qball_light(beta):
    dN_dt = 4.38e5  # s^-1
    dN_dx = dN_dt/(beta*c) # m^-1
    return dN_dx



### PLOTS #####################################################################################################################################
# a)
def plot_Eloss(charge=1):
    beta_list = np.logspace(-2, -4, num=33)
    dE_el = [SlowIon_ElecLoss(b,charge) for b in beta_list]
    dE_nu = [nuclear_loss_full(b,charge) for b in beta_list]
    dE_th = [thermal_shock(b) for b in beta_list]
    dE_SENS_g = [KKST(b,0,'gauge') for b in beta_list]
    dE_SENS_n = [KKST(b,0,'new') for b in beta_list]
    plt.plot(beta_list,dE_el,label=f'SECS electronic loss (z={charge})', linewidth=3, alpha=0.6)
    plt.plot(beta_list,dE_nu,label=f'SECS nuclear loss (z={charge})')
    plt.plot(beta_list,np.array(dE_el)+np.array(dE_nu),label=f'SECS total loss (z={charge})')
    plt.plot(beta_list,dE_SENS_g,label=r'SENS KKST gauge-mediation ($\sigma =1.3$e-27 m$^2$)',linestyle='--')
    plt.plot(beta_list,dE_SENS_n,label=r'SENS KKST new-type ($\sigma =1.2$e-24 m$^2$)',linestyle='--')
    plt.plot(beta_list,dE_th,label=r'Thermal shock: $A \rho v^2$ (new type, m$_{3/2} = 1$ TeV)',linestyle='-.', color='yellow')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(fontsize = 8)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\frac{dE}{dx}$ [$\frac{MeV}{m}$]', fontsize=16)
    plt.title('Energy loss of SECS & SENS (gauge-mediation & new-type)')
    plt.show()
#plot_Eloss()


#### REFERENCES #############################################################################################################################################################################################################################################################
#[1]: https://wiki.icecube.wisc.edu/index.php/Ice_Density
#[2]: MURTY, R. Effective Atomic Numbers of Heterogeneous Materials. Nature 207, 398–399 (1965).
#     & Master thesis : "Signatures of Q-balls in the IceCube" - Sarah Pieper
#[3]: Shinta Kasuya and Masahiro Kawasaki, "Baryogenesis from the Gauge-mediation type Q ball and the New type of Q ball as dark matter",  arXiv:1402.4546v2 [hep-ph] (2014)
#[4]: Shinta Kasuya, Masahiro Kawasaki, and Tsutomu T. Yanagida, "IceCube potential for detecting Q-ball dark matter in gauge mediation",  arXiv:1502.00715v2 [hep-ph] (2015)
#[5]: Sarah Pieper's Master Thesis, "Signatures of Q-balls in the IceCube Detector" (2020)
#[6]: H. H. Andersen and J.F. Ziegler, "Hydrogen stopping power and ranges in all elements", Pergamon Press (1977)
#[7]: J. LINDHARD, "Energy Dissipation by Ions in the kev Region", PHYSICAL REVIEW VOL124 NBR1 (1961)
#[8]: D. J. Ficenec, S. P. Ahlen, and A. A. Marin, "Observation of electronic excitation by extremely slow protons with applications to the detection of supermassive charged particles", PHYSICAL REVIEW D VOL36 NBR1 (1987)
#[9]: W. D. Wilson and L. G. Haggmark, "Calculations of nuclear stopping, ranges, and straggling in the low-energy region", PHYSICAL REVIEW B VOL15 NBR5 (1977), https://journals.aps.org/prb/pdf/10.1103/PhysRevB.15.2458
#[10]: D. Bakari et al., " Energy Losses of Q-balls", arXiv:hep-ex/0003003v1 (2000)
