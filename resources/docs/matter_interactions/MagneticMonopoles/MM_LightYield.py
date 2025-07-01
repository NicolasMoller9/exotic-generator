# Author: Nicolas Moller - PhD student (Chiba U.)

import numpy as np
import math
import matplotlib.pyplot as plt
from PDecay_beta_to_MeanFreePath import beta_to_lamda
from scipy.integrate import quad


# 0) USEFUL CONSTANTS AND FUNCTIONS ##################################################################################################################################

# Conversion factors:
J_to_MeV = 6241509647120.4               # MeV/J
statC_to_MeV = np.sqrt(1e-9 * J_to_MeV)  # (MeV*m)^1/2 / StatC
MeV_to_kg = 1.78266192e-30               # kg/(MeV/c^2)

# Useful constants:
pi = np.pi
e_cgs = 4.80320425*1e-10                 # statC = g^1/2 * cm^3/2 * s^-1 (Elementary charge)
e = e_cgs*statC_to_MeV                   # (MeV*m)^1/2 (Elementary charge)  =3.7946863625337996e-08
gD = 137*e/2                             # (MeV*m)^1/2 (Dirac charge)
NA = 6.02214076e23                       # mol^-1 (Avogadro number)
mec2 = 0.511                             # MeV (electron rest energy)
M_mono = 1e14                            # MeV (monopole mass)
alpha = 1/137                            # (fine structure constant)
a0 = 5.29177210544e-9                    # cm (Bohr radius)
c = 2.99792458e8                         # m/s (speed of light)
L0, L1 = [300 * 10**-7,  600 * 10**-7]   # cm (The lower and upper integration limits of the Cherenkov spectrum [1])
N_PHOT_AMP = 2*pi*alpha * (1/L0 - 1/L1)*100 # m^-1 (integral of Frank-Tamm formula (PDG) = Amplitude of the number of photons emitted per unit path length)

# Ice properties:
rho_ice = 0.917                          # g/cm^3 (Ice density. From [2])
I = 7.5e-5                               # MeV (mean ionization potential of the ice) (value used in ppc)
ZH = 1                                   # (atomic number of hydrogen)
MH = 938.78                              # MeV/c2 (hydrogen mass)                
ZO = 8                                   # (atomic number of oxygen)
MO = 14902.973                           # MeV/c2 (oxygen mass) 
#Z_H2O = 7.42/3                          # Effective number of electron per ice atom (from [3])
Z_H2O = 10                               # Mean number of electron per ice atom
A_H2O = 18.015                           # g/mol (mean molecular weight of water atoms)
n_H2O = rho_ice/(A_H2O/NA)               # cm^-3 (number density of H2O molecule)
ne_H2O = Z_H2O*n_H2O                     # cm-3 (electron density)
REFR_index = 1.3195                     # refractive index (for wavelength of 400nm [4])
T_0 = mec2 * ((1/np.sqrt(1 - 1/REFR_index**2)) - 1)   # MeV (minimum kinetic energy required for an electron to produce Cherenkov radiation [5])
N_DELTA_AMP = 2*pi * ne_H2O*1e6 * gD**2 * e**2 / mec2  # MeV/m (delta-ray production amplitude [6])
X0_ICE = 36.08                           # g/cm^2 (radiation length, from PDG)

# Useful functions
def gamma_(beta):                        # gamma factor
    return 1./np.sqrt(1. - beta**2)     

def Cher_angle(beta):                    # Cherenkov angle in radian
    cos_th = 1/(beta*REFR_index)
    return np.acos(cos_th)

def k(n):                                # QED correction from [7]
    if n==1:
        return 0.406
    elif n==2:
        return 0.346
    elif n>=3:
        return 0.3
    else:
        return 0
    
def Bm(n):                               # Bloch correction [7] 
    if n==1:
        return 0.248
    elif n==2:
        return 0.672
    elif n==3 or n==4 or n==5:
        return 1.022
    elif n==6 or n==7 or n==8:
        return 1.685
    elif n==9:
        return 2.085
    else:
        return 0 

params = [0.240, 2.8004, 0.09116, 3.477, 3.5017]   # [X0, X1, A_DC, M_DC, C_BAR]  (values from [8])
def delta(beta, params):  #  density effect correction in the Bethe Bloch formula [9]
    bg = beta*gamma_(beta)
    X = math.log(bg,10)
    X0, X1, A_DC, M_DC, C_BAR = params
    if X0<X<X1:
       return 4.6052*X - C_BAR + A_DC*math.pow((X1-X),M_DC)
    elif X>X1:
       return 4.6052*X - C_BAR
    else:
       return 0


# 1) ENERGY LOSSES ##################################################################################################################################

def BetheBloch_mag(beta,n):  # bethe bloch formula for a magnetic charge [10]
    gam = gamma_(beta)
    bg = beta*gam
    if 0.1<=bg<=1000:  # from PDG
        A1 = 4*pi*(ne_H2O*1e6)/mec2  # MeV^-1 m^-3
        A2 = (n*gD)**2 * e**2        # MeV^2 m^2
        A = A1*A2                    # MeV/m
        B = np.log(2*mec2*bg**2/I) - 0.5 + k(n)/2 - delta(beta, params)/2 - Bm(n)
        return A*B                   # MeV/m
    else:
        return 0


def SlowMopo_ElecLoss(beta,n,Z=Z_H2O/3,ne=ne_H2O):  # energy loss of a non-relativistic proton when interacting with atomic electrons
    beta_max = 1e-2       # [11]
    if beta <= beta_max:  
        main = (8*pi*a0*e_cgs**2)*(beta/alpha)*(1**(7/6)/(1**(2/3)+Z**(2/3))**(3/2))*ne  # g * cm * s^-2 (cgs system) [12]
        main *= 1e-5      # kg * m * s^-2 = J/m
        main *= J_to_MeV  # MeV/m
        beta0 = 7e-4      # Fit parameter from [13] (indicates a reduced efficiency for electronic excitation when beta is lower than beta0)
        correction_factor = 1 - math.exp(-(beta/beta0)**2)  # [13]
        proton_to_monopole = 0.25 * n**2  # eq.9 in [11]
        return main*correction_factor*proton_to_monopole # MeV/m
    else:
        return 0


def nuclear_loss(beta,n,Z_i,M_i):  # energy loss of a non-relativistic proton when interacting with atomic nuclei (from [14])
    gam_ = gamma_(beta)
    a = 0.8854*a0*1e-2 / (1**0.23 + Z_i**0.23)  # m
    epsilon = ((gam_-1)*a*M_i) / (1*Z_i*e**2)   # no unit
    Tm = 4 * M_i * (gam_ - 1)                   # MeV
    main = pi * a**2 * Tm * n_H2O*1e6 / epsilon # MeV / m
    if epsilon<=30:
        S_eps = 0.5*np.log(1+1.1383*epsilon) / (epsilon+(0.01321*epsilon**0.21226)+(0.19593*epsilon**0.5))    # no units (from table 2 of [14])
    else:
        S_eps = 0.5*np.log(epsilon) / epsilon
    proton_to_monopole = 0.25 * n**2            # eq.9 in [11]  
    return main*S_eps*proton_to_monopole        # MeV/m
def nuclear_loss_full(beta,n):   # total nuclear energy loss in ice
    return (2*nuclear_loss(beta,n,ZH,MH)+nuclear_loss(beta,n,ZO,MO))/3


def dEdx_monopole(beta, n):
    betaL = 1e-2
    betaH = np.sqrt(0.1**2 / (1+0.1**2))
    if beta>=betaH:          # fast monopoles can ionize nearby atoms due to the electric field they produce when moving
        dE_dx = BetheBloch_mag(beta,n)
    elif beta<=betaL:        # excitation
        dE_dx = SlowMopo_ElecLoss(beta,n) + nuclear_loss_full(beta,n)   # *1.37 if contribution from the electron magnetic moment? [10]
    elif betaL<beta<betaH:   # linear interpolation
        yL = SlowMopo_ElecLoss(betaL,n) + nuclear_loss_full(betaL,n)
        yH = BetheBloch_mag(betaH,n)
        f = (np.log10(yH) - np.log10(yL)) / (np.log10(betaH) - np.log10(betaL))
        dE_dx = yL*(beta/0.01)**f 
    else:
        print(f"Error: invalid beta value: {beta}")
    return dE_dx


def dEdx_FastElectron(beta): # normal BetheBloch & Bremsstrahlung
    gam = gamma_(beta)

    # normal bethe block (PDG)
    Tm_el = mec2 * beta**2 * gam**2 / (1 + gam)       # MeV
    A = 4*pi * (ne_H2O*1e6) * e**4 / mec2 / beta**2   # MeV/m
    B = 0.5 * math.log(2 * mec2 * beta**2 * gam**2 * Tm_el / I**2) - beta**2 - delta(beta, params)/2

    # Bremsstrahlung  (PDG)
    T = gam * mec2 * beta**2           # MeV
    brems = 100 * T / (X0_ICE/rho_ice) # MeV/m
    return A*B + brems   # MeV/m



# 2) QUENCHING ##########################################################################################################################################################################

def dNdE(dEdx,quenching=True,S_shift=0,k_shift=0): # From [15]
    #T0 = 0               # minimum energy of a delta ray to be able to produce an exciton in the halo (region further away from the track)
    #F = 0.5 * (np.log(2*mec2*1e6*SC_beta**2*SC_gamma**2/T0)-SC_beta**2)/(np.log(2*mec2*1e6*SC_beta**2*SC_gamma**2/I)-SC_beta**2) # fraction of energy producing excitons in the halo (through delta-ray production)
    F = 0                 
    k = (1e-4 + k_shift)/ rho_ice      # m/MeV (from literature [16])
    S = 10.96478 + S_shift             # photon/MeV  (Lab experiment in Wuppertal [17])
    if quenching:
        Sbis = S*((1-F)/(1 + k*(1-F)*dEdx) + F)
        return Sbis  # photons/MeV
    else:
        return S     # photons/MeV



# 3) LIGHT YIELD ########################################################################################################################################################################

# 3.1) Luminescence light from the monopole itself
# 3.1.1) Luminescence
def luminescence_yield(beta,n,quenching=True,S_shift=0,k_shift=0):
    dedx = SlowMopo_ElecLoss(beta,n)  # MeV/m
    dnde = dNdE(dedx,quenching,S_shift,k_shift)  # photons/MeV
    return dnde*dedx  # photons/m

# 3.1.2) Cherenkov light
def monopole_DirectCherenkov_light(beta, n):  # (from [18])
    if (beta*REFR_index) > 1:
        return N_PHOT_AMP * (n*gD*REFR_index/e)**2 * (1 - 1/(beta**2 * REFR_index**2))  # m^-1
    else:
       return 0
#print(N_PHOT_AMP * (1*gD*REFR_index/e)**2)


# 3.2) Light emission from secondary particles
# 3.2.1) Proton decay catalysis
def EM_cascade_yield(E):  # Light yield of an electromagnetic cascade of energy E/MeV (from [19])
    dx_de = 5.21e-3  # m/MeV (dx/dE for water)
    ratio = 0.9216   # ratio of ice/water density
    eff_length = dx_de * (0.924/ratio) * E # m (effective length of photon emission of the electromagnetic cascades)
    beta = np.sqrt(E**2 - mec2**2)/E       # speed of the positron
    thetaC = Cher_angle(beta)              # Cherenkov angle
    dndx = 49000 * np.sin(thetaC)          # m^-1 (from Frank-Tamm formula)
    gamma_per_cascade = eff_length * dndx  # gammas / cascade
    return gamma_per_cascade
#print(f"A 940 MeV positron cascade produces {EM_cascade_yield(940)} photons!")
def ProtonDecay_light(beta,sigma0=1e-31,E=940.): # sigma0 in m^2, E in MeV
    if beta<0.1:  # [20]
        lamda = beta_to_lamda(beta, sigma0)[0] # m
        positron_per_meter = 1/lamda           # m^-1
        gamma_per_positron = EM_cascade_yield(E)
        return gamma_per_positron * positron_per_meter  # m^-1
    else:
        return 0

# 3.2.2) Indirect Cherenkov
def electron_DirectCherenkov_light(beta): # Number of cherenkov photons produced by a single electron
    if (beta*REFR_index) > 1:
       return N_PHOT_AMP * (1 - 1/(beta**2 * REFR_index**2)) # m^-1
    else:
       return 0
def int_int_kyg(T):  # dn/de: number of photons per energy deposited by delta-electron (integrand of eq. 3.58 in [21])
    gam = 1 + T/mec2
    beta = np.sqrt(1 - 1/gam**2)
    return electron_DirectCherenkov_light(beta)/dEdx_FastElectron(beta)   # (m^-1) / (MeV/m) = MeV^-1
# KYG form factor constants [22]
p0 =  1.00001
p1 = -0.000186945
p2 =  3.01349e-05
p3 = -1.51855e-06
p4 =  2.93075e-08
p5 = -1.73281e-10
p6 =  3.20185e-13
def int_kyg(T, T_m, beta, n):  # total cherenkov light from a delta-electron being stopped in the medium (integrand of eq 3.59 in [21])
                               # T is the electron energy in MeV, beta is the monopole speed
    if T_0<T<T_m:
        gam = gamma_(beta)
        E = gam * mec2    # MeV
        p = E * beta      # MeV
        px = -1.0 * (T + mec2 - gam*E)/(gam * beta)  # MeV
        x = math.degrees(math.acos(px/p))
        F_kyg = p0 + p1*x + p2*x**2 + p3*x**3 + p4*x**4 + p5*x**5 + p6*x**6   # KYG form factor (3.23?)
        DeltaElectronDistribution_KYG = beta**2 * REFR_index * (N_DELTA_AMP*n**2) * F_kyg / T**2  # eq 3.52 in [21]
        return quad(int_int_kyg, T_0, T)[0] * DeltaElectronDistribution_KYG   # in m^-1 * MeV^-1
    else:
        return 0
def monopole_IndirectCherenkov_light(beta, n):  # eq 3.59 in [21]
    gam = gamma_(beta)
    m_ratio = mec2/M_mono
    denom = 1 + 2*gam*m_ratio + m_ratio**2
    T_m = 2 * mec2 * beta**2 * gam**2 / denom    # MeV (maximum kinetic energy of recoiled electrons (see PDG p.442, Eq. 33.4). me/M << 1 => k~1)
    return quad(int_kyg, T_0, T_m, args=(T_m,beta,n))[0]  # in m^-1


### PLOTS #####################################################################################################################################
# a)
def plot_dedx():
    beta_list = np.logspace(0, -4, num=310)[1:]
    dE_dx_1 = [dEdx_monopole(b,1) for b in beta_list]
    dE_dx_2 = [dEdx_monopole(b,2) for b in beta_list]
    dE_dx_3 = [dEdx_monopole(b,3) for b in beta_list]
    dE_dx_4 = [dEdx_FastElectron(b) for b in beta_list]
    plt.plot(beta_list,dE_dx_1,label='g=gD')
    plt.plot(beta_list,dE_dx_2,label='g=2*gD')
    plt.plot(beta_list,dE_dx_3,label='g=3*gD')
    plt.plot(beta_list,dE_dx_4,label='electron (BB+Brem)')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel(r'$\frac{dE}{dx}$ [$\frac{MeV}{m}$]', fontsize=16)
    plt.xlabel(r'$\beta_{monopole}$')
    plt.legend()
    plt.title(r'Energy loss of monopoles in South Pole ice ($\rho = 0.917$ g/cm$^3$)')
    
    #plt.show()
    plt.savefig('Monopole_dEdx.png')
    plt.close()
#plot_dedx()


# b)
def plot_dndx():
    beta_list_slow = np.logspace(-1, -4, num=210)[1:]
    dN_dx_1_slow = [luminescence_yield(b,1) for b in beta_list_slow]
    dN_dx_2_slow = [ProtonDecay_light(b,1e-31) for b in beta_list_slow]
    #dN_dx_3_slow = np.array(dN_dx_1_slow)+np.array(dN_dx_2_slow)
    dN_dx_4_slow = [ProtonDecay_light(b,1e-35) for b in beta_list_slow]
    
    beta_list_fast = np.linspace(0.1, 0.999, num=300)
    dN_dx_2_fast = [electron_DirectCherenkov_light(b) for b in beta_list_fast]
    dN_dx_3_fast = [monopole_DirectCherenkov_light(b,1) for b in beta_list_fast]
    dN_dx_4_fast = [monopole_IndirectCherenkov_light(b,1) for b in beta_list_fast]
    #dN_dx_5_fast = np.array(dN_dx_3_fast)+np.array(dN_dx_4_fast)

    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(12, 5))

    # slow
    mask = np.array(beta_list_slow) <= 1e-2
    x_filtered = np.array(beta_list_slow)[mask]
    y_filtered = np.array(dN_dx_1_slow)[mask]
    axes[0].plot(x_filtered, y_filtered,label='Luminescence (g=g$_D$)', color='green')
    axes[0].plot(beta_list_slow,dN_dx_2_slow,label=r'Proton Decay ($\sigma_0 = 10^{-31}$ m$^2$)', color='darkblue')
    #axes[0].plot(beta_list_slow,dN_dx_3_slow,label='TOTAL Monopole light', color='gray', linewidth=5, alpha=0.4)
    axes[0].plot(beta_list_slow,dN_dx_4_slow,'--',label=r'Proton Decay ($\sigma_0 = 10^{-35}$ m$^2$)', color='darkblue')
    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[0].set_ylabel(r'$\frac{dN}{dx}$ [m$^{-1}$]', fontsize=16)
    axes[0].legend(loc='lower left', bbox_to_anchor=(0.15, 0.01))
    axes[0].grid(axis='y', linestyle=(0, (1, 7)), color='gray', alpha=0.6)
    axes[0].set_xlim(10**(-4.2), beta_list_slow.max())
    # Create the lambda axis
    beta_ticks_slow = [min(beta_list_slow, key=lambda x: abs(x - v)) for v in [0.001, 10**(-2.5), 0.01, 10**(-1.5)]]
    lambda_ticks_slow = np.array([beta_to_lamda(b, 1e-31)[0] for b in beta_ticks_slow])
    ax2 = axes[0].twiny()  # Création d'un deuxième axe x pour 'lambda'
    ax2.set_xscale('log')
    ax2.set_xlim(axes[0].get_xlim())
    ax2.set_xticks(beta_ticks_slow)
    ax2.set_xticklabels([f"{val:.4g}" for val in lambda_ticks_slow])
    ax2.set_xlabel(r'$\lambda$ [m]', x=0.25, fontsize=13)
    
    # fast
    axes[1].plot(beta_list_fast,dN_dx_3_fast, label='Monopole Direct Cherenkov (g=g$_D$)', color='blue', linestyle='-.')
    axes[1].plot(beta_list_fast,dN_dx_4_fast, label='Monopole Indirect Cherenkov (g=g$_D$)', color='red', linestyle='-.')
    #axes[1].plot(beta_list_fast,dN_dx_5_fast, color='gray', linewidth=5, alpha=0.4)
    axes[1].plot(beta_list_fast,dN_dx_2_fast, label='Muon Direct Cherenkov', color='black', linewidth=3, linestyle='dotted')
    axes[1].legend(loc='upper left')
    axes[1].grid(axis='y', linestyle=(0, (1, 7)), color='gray', alpha=0.6)
    axes[1].yaxis.set_label_position("right")
    axes[1].yaxis.tick_right()
    axes[1].tick_params(labelright=True)
    axes[1].set_xlim(beta_list_fast.min(), 1.05)

    fig.supxlabel(r'$\beta$')
    plt.suptitle(r'Light yield of monopoles in South Pole ice ($\rho = 0.917$ g/cm$^3$)')
    plt.subplots_adjust(wspace=0.)
    #plt.show()
    plt.savefig('Monopole_dNdx.png')
    plt.close()
#plot_dndx()


#### REFERENCES #############################################################################################################################################################################################################################################################
#[1]: Aartsen et al., "IceCube Sensitivity for Low-Energy Neutrinos from Nearby Supernovae" (2011)
#[2]: https://wiki.icecube.wisc.edu/index.php/Ice_Density
#[3]: MURTY, R. Effective Atomic Numbers of Heterogeneous Materials. Nature 207, 398–399 (1965).
#     & Master thesis : "Signatures of Q-balls in the IceCube" - Sarah Pieper
#[4]: Stephen G. Warren, "Optical constants of ice from the ultraviolet to the microwave" (1984)
#[5]: Esther Ciarrocchi and NicolaBelcari, "Cerenkov luminescence imaging: physics principles and potential applications in biomedical sciences" (2017)
#[6]: Jackson, "Classical Electrodynamics"
#[7]: J. Derkaoui et al. - Astroparticle Physics 9 (1998) 173-183
#[8]: Dima's Thesis table C.2, "Cosmic Ray Energy Spectrum Measurement with the Antarctic Muon and Neutrino Detector Array (AMANDA)" (2003), https://user-web.icecube.wisc.edu/~dima/work/BKP/DCS/THESIS/main.pdf
#     or https://web.archive.org/web/20070713230354/http://pdg.lbl.gov/AtomicNuclearProperties/substances/276.html
#[9]: R. M. Sternheimer, S.M. Seltzer and M. J.Berger, "Density effect for the ionization loss of charged particles in various substances", PHYSICAL REVIEW B VOL26 NBR11 (1982), https://journals.aps.org/prb/pdf/10.1103/PhysRevB.26.6067
#[10]: L. Patrizii1 and M. Spurio - Annu. Rev. Nucl. Part. Sci. 2015. 65:279–302
#[11]: J. Derkaoui et al. - Astroparticle Physics 10 (1999) 339-352
#[12]: J. LINDHARD, "Energy Dissipation by Ions in the kev Region", PHYSICAL REVIEW VOL124 NBR1 (1961)
#[13]: D. J. Ficenec, S. P. Ahlen, and A. A. Marin, "Observation of electronic excitation by extremely slow protons with applications to the detection of supermassive charged particles", PHYSICAL REVIEW D VOL36 NBR1 (1987)
#[14]: W. D. Wilson and L. G. Haggmark, "Calculations of nuclear stopping, ranges, and straggling in the low-energy region", PHYSICAL REVIEW B VOL15 NBR5 (1977), https://journals.aps.org/prb/pdf/10.1103/PhysRevB.15.2458
#[15]: M.H. SALAMON and S.P. AHLEN, "PLASTIC SCINTILLATOR RESPONSE TO RELATIVISTIC Ne, Ar, Fe IONS", Nuclear Instruments and Methods 195 (1982) 557-568 
#[16]: J. B. Birks, "The Theory and Practice of Scintillation Counting", Pergamon Press, Oxford (1964)
#[17]: A. Pollmann PoS ICRC2021 (2021) 1093.
#[18]: D.R. Tompkins., "Total Energy Loss and Cherenkov Emission from Monopoles", Phys. Rev. 138, B248, (1965)
#      or Anna Obertacke, "Searches for Relativistic Magnetic Monopoles in IceCube", arXiv:1511.01350v2 [astro-ph.HE] (2015)
#[19]: IceCube collaboration, "Measurement of South Pole ice transparency with the IceCube LED calibration system", https://arxiv.org/pdf/1301.5361 (2013)
#[20]: Frederik Lauber (IceCube), "Ongoin magnetic monopole searches with IceCube", EPJ Web of Conferences 182, 02071 (2018)
#[21]: Anna Pollmann, "Search for mildly relativistic Magnetic Monopoles with IceCube", Thesis (2015)
#[22]: KYG form factors??????? Fig. 3.7 Anna's Thesis
