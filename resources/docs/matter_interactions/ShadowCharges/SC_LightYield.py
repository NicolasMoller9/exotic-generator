# Author: Nicolas Moller - PhD student (Chiba U.)

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
from NuclearCrossSection import differential_NCC


# 0) USEFUL CONSTANTS AND FUNCTIONS ##################################################################################################################################

# ShadowCharge parameters
SC_beta = 1e-3
SC_gamma = 1/np.sqrt(1-SC_beta**2)

# Conversion factors:
J_to_MeV = 6241509647120.4               # MeV/J
statC_to_MeV = np.sqrt(1e-9 * J_to_MeV)  # (MeV*m)^1/2 / StatC

# Useful constants:
pi = np.pi
e_cgs = 4.80320425*1e-10                 # statC = g^1/2 * cm^3/2 * s^-1 (Elementary charge)
e = e_cgs*statC_to_MeV                   # (MeV*m)^1/2 (Elementary charge)
NA = 6.022*1e23                          # mol^-1 (Avogadro number)
alpha = 1/137                            # (fine structure constant)
a0 = 5.29177210544e-9                    # cm (Bohr radius)
c = 2.99792458e8                         # m/s (speed of light)
m_e = 0.511                              # MeV (electron rest energy)
m_p = 938.272                            # MeV/c2 (proton rest energy)
vBohr = 2.188*1e6                        # m/s (Bohr velocity: speed of an electron in the first orbit of a hydrogen atom)

# Ice properties:
rho_ice = 0.917                          # g/cm^3 (Ice density. From [1])
I = 7.5e-5                               # MeV (mean ionization potential of the ice) (value used in ppc)
ZH = 1                                   # (atomic number of hydrogen)
AH = 1.008                               # g/mol (molar mass of hydrogen)
MH = 938.78                              # MeV/c2 (hydrogen mass)                
ZO = 8                                   # (atomic number of oxygen)
AO = 15.999                              # g/mol (molar mass of oxygen)
MO = 14902.973                           # MeV/c2 (oxygen mass) 
#Z_H2O = 7.42                            # Effective number of electron per ice atom (from [12])
Z_H2O = 10                               # (number of electron per H2O)      
A_H2O = 18.015                           # g/mol (molar mass of H2O)
n_H2O = rho_ice/(A_H2O/NA)               # cm^-3 (number density of H2O molecule)
ne_H2O = Z_H2O*n_H2O                     # cm-3 (electron density)
atomic_separation = 2 * (4*pi*n_H2O*1e6/3)**(-1/3) # m

# Useful functions
def gamma_(beta):                        # gamma factor
    return 1./np.sqrt(1. - beta**2)              

params = [0.240, 2.8004, 0.09116, 3.477, 3.5017]   # [X0, X1, A_DC, M_DC, C_BAR]  (values from [2])
def delta(beta, params):  #  density effect correction in the Bethe Bloch formula [3]
    bg = beta*gamma_(beta)
    X = math.log(bg,10)
    X0, X1, A_DC, M_DC, C_BAR = params
    if X0<X<X1:
       return 4.6052*X - C_BAR + A_DC*math.pow((X1-X),M_DC)
    elif X>X1:
       return 4.6052*X - C_BAR
    else:
       return 0

def E_proton(beta, M=m_p):    # (total energy of a proton = gamma*M*c^2)
    return gamma_(beta) * M   # MeV
  


# 1) ENERGY LOSSES ##################################################################################################################################

def SlowIon_ElecLoss(beta,z,Z=Z_H2O/3,ne=ne_H2O):  # energy loss of a non-relativistic ion when interacting with atomic electrons
    beta_max = min(vBohr * z**(2/3) / c,1e-2)
    if beta <= beta_max:  # [4]*
        main = (8*pi*a0*e_cgs**2)*(beta/alpha)*(z**(7/6)/(z**(2/3)+Z**(2/3))**(3/2))*ne  # in g * cm * s^-2 (cgs system) [5]
        main *= 1e-5 # in kg * m * s^-2 = J/m
        main *= J_to_MeV # in MeV/m
        beta0 = 7e-4        # Fit parameter from [6] (indicates a reduced efficiency for electronic excitation when beta is lower than beta0)
        correction_factor = 1 - math.exp(-(beta/beta0)**2)  # [6]
        return main*correction_factor # in MeV/m
    else:
        return 0
# *when v<v0 (ion slower than orbital electron), the energy tranfer to the electron is low, which means it can be more simply described by the Lindhard formalism
# and Lindhard formalism includes a screening effect based on the Thomas-Fermi model of the atom, but when v>v0, the screening is different and can't be described by this model anymore


def nuclear_loss(beta,z,Z_i,M_i):  # energy loss of a non-relativistic ion when interacting with atomic nuclei (from [7])
    gam_ = gamma_(beta)
    a = 0.8854*a0*1e-2 / (z**0.23 + Z_i**0.23)  # m
    epsilon = ((gam_-1)*a*M_i) / (z*Z_i*e**2)   # no unit
    Tm = 4 * M_i * (gam_ - 1)                   # MeV
    main = pi * a**2 * Tm * n_H2O*1e6 / epsilon # MeV / m
    if epsilon<=30:
        S_eps = 0.5*np.log(1+1.1383*epsilon) / (epsilon+(0.01321*epsilon**0.21226)+(0.19593*epsilon**0.5))    # no units (from table 2 of [7])
    else:
        S_eps = 0.5*np.log(epsilon) / epsilon
    return main * S_eps                         # MeV/m
def nuclear_loss_full(beta,z):   # total nuclear energy loss in ice
    return (2*nuclear_loss(beta,z,ZH,MH)+nuclear_loss(beta,z,ZO,MO))/3


def FastIon_BetheBloch(beta,z): # standard bethe bloch formula (from PDG)
    gam = gamma_(beta)
    bg = beta*gam
    if 0.1<=bg<=1000:  # from PDG
        Tm_el = 2 * m_e * bg**2   # MeV
        A = 4*pi * (ne_H2O*1e6) * z**2 * e**4 / m_e / beta**2   # MeV/m
        B = 0.5 * math.log(2 * m_e * bg**2 * Tm_el / I**2) - beta**2 - delta(beta, params)/2
        return A*B  # MeV/m
    else:
        return 0


"""
def SlowToFast_Varelas_interpolation(beta, z): # Varelas and Biersack interpolation [4] (used in [8])
    beta_min = min(vBohr*z**(2/3) / c, 1e-2)
    beta_max = np.sqrt(0.1**2 / (1+0.1**2))
    if beta_min<beta<beta_max:
        yL = SlowIon_ElecLoss(beta_min,z)   # MeV/m
        yH = FastIon_BetheBloch(beta_max,z) # MeV/m
        E = E_proton(beta)                  # MeV, energy of a proton with velocity beta  
        E_min = E_proton(beta_min)
        E_max = E_proton(beta_max)
        A1 = yL / E_min**(0.45)             # MeV^0.55 / m
        A2 = (243 - 0.375 * Z_H2O) * Z_H2O  # (what units)???  should be put in MeV^2 / m
        A4 = 4 * m_e / (I * m_p)            # MeV^-1
        #print(yH*EH/A2)
        A3 = E_max * (np.exp(yH*E_max/A2) - 1 - A4*E_max)  # MeV
        print(A3)
        SL = A1 * np.sqrt(E)
        SH = (A2 / E) * np.log(1 + (A3 / E) + E * A4)
        return SL * SH / (SL + SH)
    else:
        return 0
"""


def dEdx_Ion(beta,z):
    betaL = min(vBohr*z**(2/3) / c, 1e-2)
    betaH = np.sqrt(0.1**2 / (1+0.1**2))
    if beta>=betaH:          # fast ions can ionize nearby atoms due to the electric field they produce when moving
        dE_dx = FastIon_BetheBloch(beta,z)
    elif beta<=betaL:        # excitation
        dE_dx = SlowIon_ElecLoss(beta,z) + nuclear_loss_full(beta,z)
    elif betaL<beta<betaH:   # linear interpolation
        yL = SlowIon_ElecLoss(betaL,z) + nuclear_loss_full(beta,z) # MeV/m
        yH = FastIon_BetheBloch(betaH,z) # MeV/m
        f = (np.log10(yH) - np.log10(yL)) / (np.log10(betaH) - np.log10(betaL))
        dE_dx = yL*(beta/betaL)**f
    else:
        print(f"Error: invalid beta value: {beta}")
    return dE_dx

# 2) QUENCHING ##########################################################################################################################################################################

def dNdE(dEdx,quenching=True,S_shift=0,k_shift=0): # From [9]
    #T0 = 0               # minimum energy of a delta ray to be able to produce an exciton in the halo (region further away from the track)
    #F = 0.5 * (np.log(2*mec2*1e6*SC_beta**2*SC_gamma**2/T0)-SC_beta**2)/(np.log(2*mec2*1e6*SC_beta**2*SC_gamma**2/I)-SC_beta**2) # fraction of energy producing excitons in the halo (through delta-ray production)
    F = 0                 
    k = (1e-4 + k_shift)/ rho_ice      # m/MeV (from literature [10])
    S = 10.96478 + S_shift             # photon/MeV  (Lab experiment in Wuppertal [11])
    if quenching:
        Sbis = S*((1-F)/(1 + k*(1-F)*dEdx) + F)
        return Sbis  # photons/MeV
    else:
        return S     # photons/MeV


# 3) LIGHT YIELD ########################################################################################################################################################################

# 3.1) Luminescence light from the shadow charge itself
def luminescence_yield(beta,z,quenching=True,S_shift=0,k_shift=0):
    dedx = SlowIon_ElecLoss(beta,z)  # MeV/m
    dnde = dNdE(dedx,quenching,S_shift,k_shift)  # photons/MeV
    return dnde*dedx  # photons/m


# 3.2) Secondary yield: luminescence from recoiled nucleus
#a1. photons produced by a recoiled nucleus (luminescence) per energy deposit (total energy loss)
def recoiled_nucl_instantaneous_yield(E, Zrec, Mrec, nucl_quench=True,S_shift=0,k_shift=0):  # E in MeV and Mrec in MeV/c2
    beta = np.sqrt(2*E/Mrec)
    dNdx = luminescence_yield(beta,Zrec,quenching=nucl_quench, S_shift=S_shift, k_shift=k_shift)  # photons/m
    dEdx = dEdx_Ion(beta, Zrec)  # MeV/m
    dNdE = dNdx/dEdx             # photons/MeV
    return dNdE
#a2. inverse of energy loss of recoil nucleus (integrand for range calculation) 
def recoiled_nucl_range_integrand(E, Zrec, Mrec):  # E in MeV and Mrec in MeV/c2
    beta = np.sqrt(2*E/Mrec)
    dEdx = dEdx_Ion(beta, Zrec)  # MeV/m
    dxdE = 1/dEdx # m/MeV
    return dxdE
#b1. integral of the nucleus light yield over its path length (from initial energy T to 0)
def recoiled_nucl_tot_yield(E, Zrec, Mrec, nucl_quench=True,S_shift=0,k_shift=0): # E in MeV and Mrec in MeV/c2
    result = quad(recoiled_nucl_instantaneous_yield, 0, E, args=(Zrec, Mrec, nucl_quench, S_shift, k_shift))[0] # photons
    return result
#b2. distance traveled by the nucleus until it stops (E=0), in unit of atomic separation
def recoiled_nucl_range(E, Zrec, Mrec): # E in MeV and Mrec in MeV/c2
    result = quad(recoiled_nucl_range_integrand, 0, E, args=(Zrec, Mrec))[0] # photons
    return result / atomic_separation
#c1. number of photons produced by the nuclei of energy T multiplied by probability of energy transfer T from the shadow charge (cross section)
def secondary_yield_integrand(E, beta, z, Zrec, Mrec, nucl_quench=True,S_shift=0,k_shift=0,f0=0.16):
    L = recoiled_nucl_tot_yield(E, Zrec, Mrec, nucl_quench, S_shift, k_shift)  # photons
    dsigma = differential_NCC(E,beta,z,Zrec,Mrec,f0)[-1]                       # m^2 / MeV
    return L*dsigma  # photons / m / MeV
#c2. transfered energy T multiplied by probability of energy transfer T (cross section)
def secondary_Eloss_integrand(E, beta, z, Zrec, Mrec):
    dsigma = differential_NCC(E,beta,z,Zrec,Mrec)[-1] # m^2 / MeV
    return E*dsigma # m^2
#d1. integral of c1. over all possible energy transfer, for H and O
def secondary_yield(beta,z,nucl_quench=True,S_shift=0,k_shift=0,f0=0.16):
    gam_ = gamma_(beta)
    TmaxH = 4 * MH * (gam_ - 1)  # MeV
    TmaxO = 4 * MO * (gam_ - 1)  # MeV
    dNdx_H = n_H2O*1e6*quad(secondary_yield_integrand, 0, TmaxH, args=(beta, z, ZH, MH, nucl_quench, S_shift, k_shift, f0))[0] # photons / m
    dNdx_O = n_H2O*1e6*quad(secondary_yield_integrand, 0, TmaxO, args=(beta, z, ZO, MO, nucl_quench, S_shift, k_shift, f0))[0] # photons / m
    dNdx = (2*dNdx_H + dNdx_O)/3    # photons / m
    return dNdx
#d2. integral of c2. over all possible energy transfer, for H and O
def nuclear_loss_bis(beta,z):
    gam_ = gamma_(beta)
    TmaxH = 4 * MH * (gam_ - 1)  # MeV
    TmaxO = 4 * MO * (gam_ - 1)  # MeV
    dEdx_H = n_H2O*1e6*quad(secondary_Eloss_integrand, 0, TmaxH, args=(beta, z, ZH, MH))[0] # MeV / m
    dEdx_O = n_H2O*1e6*quad(secondary_Eloss_integrand, 0, TmaxO, args=(beta, z, ZO, MO))[0] # MeV / m
    dEdx = (2*dEdx_H + dEdx_O)/3  # MeV/m
    return dEdx


### PLOTS #####################################################################################################################################
# a)
def plot_eloss(zSC=10, option=1):
    beta_list = np.logspace(0, -5, num=310)[1:]
    if option==1:
        zlist=[zSC]
    else:
        zlist=[1, 2, 8, 14, 26]
    for z in zlist:
        dE_dx_elec = [SlowIon_ElecLoss(b, z) for b in beta_list]
        dE_dx_nucl = [nuclear_loss_full(b, z) for b in beta_list]
        dE_dx_nucl2 = [nuclear_loss_bis(b, z) for b in beta_list]
        dE_dx_bethe = [FastIon_BetheBloch(b, z) for b in beta_list]
        dE_dx_tot = [dEdx_Ion(b, z) for b in beta_list]

        x_elec_filtered, dE_dx_elec_filtered = zip(*[(xi, yi) for xi, yi in zip(beta_list, dE_dx_elec) if yi != 0])
        x_bethe_filtered, dE_dx_bethe_filtered = zip(*[(xi, yi) for xi, yi in zip(beta_list, dE_dx_bethe) if yi != 0])

        if option==1:
            plt.plot(x_elec_filtered,dE_dx_elec_filtered,label='electronic loss', color='red')
            plt.plot(beta_list,dE_dx_nucl,label='nuclear loss', color='blue')
            plt.plot(beta_list,dE_dx_nucl2, '--',label=r'$n_{nucl}\int(T*\frac{d\sigma}{dT} dT)$', color='blue')
            plt.plot(x_bethe_filtered,dE_dx_bethe_filtered,label='bethe-bloch', color='orange')
        else:
            plt.plot(beta_list,dE_dx_tot,label=fr'total loss ($z_{{SC}}={z}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'dE/dx / (MeV m$^{-1})$')
    plt.legend()
    plt.title("Energy loss of shadow charges in South Pole ice")
    if option==1: plt.text(0.05, 0.50, fr'$z_{{SC}} = {zSC}$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    #plt.show()
    if option==1: plt.savefig(f'dE_dx_z{zSC}.png')
    else: plt.savefig('dE_dx_tot.png')  # compare with fig. 5 of [7]
    plt.close()
#plot_eloss(option=1)


# b)
def plot_recoiled_nucleus_dNdE(beta):
    gam_ = gamma_(beta)
    plt.figure(figsize=(8, 5))
    for Zrec, Mrec in [[ZH,MH], [ZO,MO]]:
        Tmax = 4 * Mrec * (gam_ - 1)  # MeV
        E_values = np.logspace(np.log10(1e-5 * Tmax), np.log10(Tmax), 500)
        y_values = [recoiled_nucl_instantaneous_yield(E, Zrec, Mrec) for E in E_values]
        plt.plot(E_values, y_values, label=f'Z={Zrec}')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$E_{nuc}(t)$ /MeV", fontsize=16)
    plt.ylabel(r"$\frac{dN}{dE}$ / MeV$^{-1}$", fontsize=16)
    plt.title(r"Scintillation efficiency of a nuclei with charge $Z$ and kinetic energy $E_{nuc}$")
    plt.grid(True, which="both", ls="--", lw=0.5)
    plt.legend(fontsize=14)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig('RecoilNuclei_dNdE.png')
    plt.close()
#plot_recoiled_nucleus_dNdE(SC_beta)


# c)
def plot_nuclei_TotalYield(beta):
    gam_ = gamma_(beta)
    plt.figure(figsize=(8, 5))
    for Zrec, Mrec in [[ZH,MH], [ZO,MO]]:
        Tmax = 4 * Mrec * (gam_ - 1)  # MeV
        E_values = np.logspace(np.log10(1e-5 * Tmax), np.log10(Tmax), 500)
        y_values = [recoiled_nucl_tot_yield(E, Zrec, Mrec) for E in E_values]
        plt.plot(E_values, y_values, label=f'Z={Zrec}')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$E_{nuc}(t=0)$ /MeV", fontsize=16)
    plt.ylabel(r"$L_p$ (= total photons produced)", fontsize=16)
    plt.title(r"Number of photons produced by a recoil nuclei of initial energy $E_{nuc}(t_0)$")
    plt.grid(True, which="both", ls="--", lw=0.5)
    plt.legend(fontsize=14)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig('RecoilNuclei_tot_photons.png')
    plt.close()
#plot_nuclei_TotalYield(SC_beta)


# d)
def plot_nuclei_range(beta):
    gam_ = gamma_(beta)
    plt.figure(figsize=(8, 5))
    for Zrec, Mrec in [[ZH,MH], [ZO,MO]]:
        Tmax = 4 * Mrec * (gam_ - 1)  # MeV
        E_values = np.logspace(np.log10(1e-5 * Tmax), np.log10(Tmax), 500)
        y_values = [recoiled_nucl_range(E, Zrec, Mrec) for E in E_values]
        quenching_limit = [(deltaE/1e3) * (1/atomic_separation) for deltaE in E_values]
        plt.plot(E_values, y_values, label=f'Z={Zrec}')
    plt.plot(E_values, quenching_limit, color='grey', linestyle='--', label=r'$\langle \frac{dE}{dx} \rangle = 10^3$ MeV/m')
    plt.fill_between(E_values, 1e-10, quenching_limit, color='grey', alpha=0.3)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$E_{nuc}(t=0)$ /MeV", fontsize=16)
    plt.ylabel(r"Range / $d_0$", fontsize=16)
    plt.ylim(1e-2, 1e5)
    plt.title(r"Distance traveled by recoil nuclei of initial energy $E_{nuc}(t_0)$")
    plt.legend(fontsize=14)
    plt.text(0.60, 0.20, 'QUENCHING', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig('RecoilNuclei_range.png')
    plt.close()
#plot_nuclei_range(SC_beta)


# e)
def plot_secondary_yield_integrand(beta,zSC):
    gam_ = gamma_(beta)
    plt.figure(figsize=(8, 5))
    for Zrec, Mrec in [[ZH,MH], [ZO,MO]]:
        Tmax = 4 * Mrec * (gam_ - 1)  # MeV
        E_values = np.logspace(np.log10(1e-5 * Tmax), np.log10(Tmax), 500)
        y_values = [n_H2O*1e6*secondary_yield_integrand(E, beta, zSC, Zrec, Mrec) for E in E_values]
        plt.plot(E_values, y_values, label=f'Z={Zrec}')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("T (= transfered energy) / MeV", fontsize=16)
    plt.ylabel(r"$\frac{dN}{dxdE}$ / photons * m$^{-1}$ MeV$^{-1}$", fontsize=16)
    plt.title(r"Secondary Yield Integrand: $n_{nuc}\times L_p\times d\sigma/dT$")
    plt.grid(True, which="both", ls="--", lw=0.5)
    plt.legend(fontsize=14)
    plt.text(0.10, 0.80, fr'$\beta = {beta}$' + '\n' + fr'$z_{{SC}} = {zSC}$' + '\n' + fr'$n_{{nuc}} = {n_H2O*1e6:.2e}$ /m$^3$', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
    plt.tight_layout()
    #plt.show()
    plt.savefig('secondary_yield_integrand.png')
    plt.close()
#plot_secondary_yield_integrand(SC_beta, 10)


# f)
def plot_secondary_yield(beta):
    z_list = list(range(1, 138))
    dn_list = []
    de_list = []
    for z in z_list:
        dn = secondary_yield(beta, z)
        de = nuclear_loss_bis(beta, z)
        dn_list.append(dn)
        de_list.append(de)

    fig, ax1 = plt.subplots()

    # Axe gauche
    ax1.plot(z_list, dn_list, color='tab:blue', label=r'$n_{nucl}\int L_p(T)*\frac{d\sigma}{dT} dT$')
    ax1.set_yscale('log')
    ax1.set_xlabel(r"z$_{SC}$", fontsize=16)
    ax1.set_ylabel(r"$\frac{dN}{dx}$ / m$^{-1}$", fontsize=16, color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.grid(True, which="both", ls="--", lw=0.5)
    # Axe droit
    ax2 = ax1.twinx()
    ax2.plot(z_list, de_list, color='tab:red', label=r'$n_{nucl}\int(T*\frac{d\sigma}{dT} dT)$')
    ax2.set_ylabel(r'$\frac{dE}{dx}$ / MeV m$^{-1}$', fontsize=16, color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=16)
    ax2.set_yscale('log')

    yticks = ax1.get_yticks()
    ax2.set_yticks(yticks)
    ax2.set_ylim(ax1.get_ylim())

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=18, loc='lower right')
    plt.title(r"Secondary Light Yield")
    fig.tight_layout()
    #plt.show()
    plt.savefig('secondary_light_yield.png')
    plt.close()
#plot_secondary_yield(SC_beta)



# g)
def plot_SecondaryYield_f0_dependence(beta):
    gam_ = gamma_(beta)
    Tmax = 4 * MH * (gam_ - 1)  # MeV
    f0_list = [0.001,0.01,0.1,0.16,0.5]
    z_list = list(range(1, 138))

    plt.figure(figsize=(8, 5))
    for f0 in f0_list:
        y = [secondary_yield(beta, z, f0) for z in z_list]
        plt.plot(z_list, y, label=fr'$f_0={f0}$')
    plt.yscale('log')
    plt.xlabel(r"z$_{SC}$", fontsize=16)
    plt.ylabel(r"$\frac{dN}{dx}$ / m$^{-1}$", fontsize=16)
    plt.title(r"Secondary light yield for different assumptions of $f_0$")
    plt.grid(True, which="both", ls="--", lw=0.5)
    plt.legend(fontsize=14)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig('SecondaryYield_f0_dependence.png')
    plt.close()
#plot_SecondaryYield_f0_dependence(SC_beta)


# h)
def plot_yield(beta):
    z_list = list(range(1, 138))  # Schwinger limit
    dE_dx_z = [SlowIon_ElecLoss(beta, z) for z in z_list]
    dE_dx_z_sec = [nuclear_loss_full(beta,z) for z in z_list]
    dN_dx_z = [luminescence_yield(beta, z) for z in z_list]
    dN_dx_z_upper = [luminescence_yield(beta, z, S_shift=8,k_shift=-2e-5) for z in z_list]
    dN_dx_z_lower = [luminescence_yield(beta, z, S_shift=-5,k_shift=2e-5) for z in z_list]
    dN_dx_z_sec = [secondary_yield(beta,z, nucl_quench=False) for z in z_list]
    dN_dx_z_sec_quenched = [secondary_yield(beta,z,nucl_quench=True) for z in z_list]
    dN_dx_z_sec_quenched_upper = [secondary_yield(beta,z,nucl_quench=True,S_shift=8,k_shift=-2e-5) for z in z_list]
    dN_dx_z_sec_quenched_lower = [secondary_yield(beta,z,nucl_quench=True,S_shift=-5,k_shift=2e-5) for z in z_list]
    muon_08_light = 7840.1196086055  # m^-1

    # Create figure and primary axis
    fig, ax1 = plt.subplots()

    # Plot light yield on left y-axis
    ax1.plot(z_list, dN_dx_z, label='primary light yield (with quenching)', color='tab:blue')
    ax1.plot(z_list, dN_dx_z_upper, linewidth=0, color='tab:blue')
    ax1.plot(z_list, dN_dx_z_lower, linewidth=0, color='tab:blue')
    ax1.fill_between(z_list, dN_dx_z_lower, dN_dx_z_upper, color='tab:blue', alpha=0.2)

    ax1.plot(z_list, dN_dx_z_sec, '--', label='secondary Light yield (no quenching)', color='tab:orange')
    ax1.plot(z_list, dN_dx_z_sec_quenched, label='secondary Light yield (with quenching)', color='tab:orange')
    ax1.plot(z_list, dN_dx_z_sec_quenched_upper, linewidth=0, color='tab:orange')
    ax1.plot(z_list, dN_dx_z_sec_quenched_lower, linewidth=0, color='tab:orange')
    ax1.fill_between(z_list, dN_dx_z_sec_quenched_lower, dN_dx_z_sec_quenched_upper, color='tab:orange', alpha=0.2)

    ax1.axhline(y=muon_08_light, color='r', linestyle='--', linewidth=1, label=r"muon light yield ($\beta=0.8$)")
    # Set labels and scales
    ax1.set_xlabel(r'$z_{SC}$', fontsize=18)
    ax1.set_ylabel(r'$\frac{dN}{dx}$ / m$^{-1}$', fontsize=18)
    ax1.set_yscale('log')
    ax1.tick_params(axis='both', labelsize=16)

    # Create second y-axis and plot energy loss
    ax2 = ax1.twinx()
    ax2.plot(z_list, dE_dx_z, ':', label='electronic energy loss', color='tab:blue')
    ax2.plot(z_list, dE_dx_z_sec, ':', label='nuclear energy loss', color='tab:orange')
    ax2.set_ylabel(r'Energy loss / (MeV m$^{-1}$)', fontsize=18)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.set_yscale('log')

    # Synchronize Y-axis limits
    ymin = 10**(3.5) #min(min(dN_dx_z), min(dE_dx_z))
    ymax = 5*10**(6)
    ax1.set_ylim(ymin, ymax)
    ax2.set_ylim(ymin, ymax)
    
    # Combine and show legend
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(
    lines_1 + lines_2,
    labels_1 + labels_2,
    fontsize=8,
    loc='center',                     # position de référence dans bbox
    bbox_to_anchor=(0.52, 0.28),      # coordonnées dans l’axe
    bbox_transform=ax1.transAxes,     # ← clef : coordonnées relatives à ax1
    frameon=True)

    
    ax1.text(0.08, 0.94, fr'$\beta_{{SC}}={beta}$', transform=ax1.transAxes, fontsize=12, verticalalignment='top')
    plt.suptitle(r'Light yield of shadow charges in South Pole ice ($\rho = 0.917$ g/cm$^3$)')
    plt.tight_layout()
    #plt.show()
    plt.savefig('dN_dx.png')
    plt.close()
#plot_yield(SC_beta)


#### REFERENCES #############################################################################################################################################################################################################################################################
#[1]: https://wiki.icecube.wisc.edu/index.php/Ice_Density
#[2]: Dima's Thesis table C.2, "Cosmic Ray Energy Spectrum Measurement with the Antarctic Muon and Neutrino Detector Array (AMANDA)" (2003), https://user-web.icecube.wisc.edu/~dima/work/BKP/DCS/THESIS/main.pdf
#     or https://web.archive.org/web/20070713230354/http://pdg.lbl.gov/AtomicNuclearProperties/substances/276.html
#[3]: R. M. Sternheimer, S.M. Seltzer and M. J.Berger, "Density effect for the ionization loss of charged particles in various substances", PHYSICAL REVIEW B VOL26 NBR11 (1982), https://journals.aps.org/prb/pdf/10.1103/PhysRevB.26.6067
#[4]: H. H. Andersen and J.F. Ziegler, "Hydrogen stopping power and ranges in all elements", Pergamon Press (1977)
#[5]: J. LINDHARD, "Energy Dissipation by Ions in the kev Region", PHYSICAL REVIEW VOL124 NBR1 (1961)
#[6]: D. J. Ficenec, S. P. Ahlen, and A. A. Marin, "Observation of electronic excitation by extremely slow protons with applications to the detection of supermassive charged particles", PHYSICAL REVIEW D VOL36 NBR1 (1987)
#[7]: W. D. Wilson and L. G. Haggmark, "Calculations of nuclear stopping, ranges, and straggling in the low-energy region", PHYSICAL REVIEW B VOL15 NBR5 (1977), https://journals.aps.org/prb/pdf/10.1103/PhysRevB.15.2458
#[8]: Mohamed Ouchrif, "Energy losses of Q-balls in Matter, Earth and Detectors", arXiv:hep-ex/0004030v1 (2000)
#[9]: M.H. SALAMON and S.P. AHLEN, "PLASTIC SCINTILLATOR RESPONSE TO RELATIVISTIC Ne, Ar, Fe IONS", Nuclear Instruments and Methods 195 (1982) 557-568 
#[10]: J. B. Birks, "The Theory and Practice of Scintillation Counting", Pergamon Press, Oxford (1964)
#[11]: A. Pollmann PoS ICRC2021 (2021) 1093.
#[12]: MURTY, R. Effective Atomic Numbers of Heterogeneous Materials. Nature 207, 398–399 (1965).
#      & Master thesis : "Signatures of Q-balls in the IceCube" - Sarah Pieper