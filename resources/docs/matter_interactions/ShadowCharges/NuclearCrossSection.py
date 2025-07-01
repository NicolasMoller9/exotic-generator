# Author: Nicolas Moller - PhD student (Chiba U.)

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# 0) USEFUL CONSTANTS AND FUNCTIONS ##################################################################################################################################

# Conversion factors:
J_to_MeV = 6241509647120.4               # MeV/J
statC_to_MeV = np.sqrt(1e-9 * J_to_MeV)  # (MeV*m)^1/2 / StatC

a0 = 5.29177210544e-9                    # cm (Bohr radius)
e_cgs = 4.80320425*1e-10                 # statC = g^1/2 * cm^3/2 * s^-1 (Elementary charge)
e = e_cgs*statC_to_MeV                   # (MeV*m)^1/2 (Elementary charge)
NA = 6.022*1e23                          # mol^-1 (Avogadro number)
c = 2.99792458e8                         # m/s (speed of light)

# Ice properties:
rho_ice = 0.917                          # g/cm^3 (Ice density. From [1])
I = 7.5e-5                               # MeV (mean ionization potential of the ice) (value used in ppc)
ZH = 1                                   # (atomic number of hydrogen)
AH = 1.008                               # g/mol (molar mass of hydrogen)
MH = 938.78                              # MeV/c2 (hydrogen mass)                
ZO = 8                                   # (atomic number of oxygen)
AO = 15.999                              # g/mol (molar mass of oxygen)
MO = 14902.973                           # MeV/c2 (oxygen mass) 
Z_H2O = 10                               # (number of electron per H2O)      
A_H2O = 18.015                           # g/mol (molar mass of H2O)
n_H2O = rho_ice/(A_H2O/NA)               # cm^-3 (number density of H2O molecule)

# Useful functions
def gamma_(beta):                        # gamma factor
    return 1./np.sqrt(1. - beta**2)  



# 1) Recover function f(t^1/2) from fig.1 in [2] #################################################################################################################
data = np.loadtxt('ft_table.txt', skiprows=1)
x_plot = data[:, 0]
y_plot = data[:, 1]
spline = UnivariateSpline(x_plot, y_plot, s=0) # 's' is a smoothing parameter
xmin_spline, xmax_spline = 2e-3, 10  # range of the function f plotted in figure 1 of [1]
# Above xmax_spline, the function f can be approximated by the Rutherford formalism
# Rutherford cross section is proportional to T^-2 => f(x) must be proportional to 1/x
def Rutherford(x): # proportional to 1/x
    A = xmax_spline * spline(xmax_spline)
    return A/x
# For very small values of x (<xmin_spline) we approximate the function f(x) to be constant!
# Last paragraph of section "Nuclear stopping and scattering cross section" in [2]:
#   "It must be emphasized though, that at extremely low epsilon-values, epsilon < 10-2, 
#   the nuclear scattering and stopping becomes somewhat uncertain, because the 
#   Thomas-Fermi treatment is a crude approximation when the ion and the atom 
#   do not come close to each other"
def f_universal(x, f0=spline(xmin_spline)):
    if x>xmax_spline:
        return Rutherford(x)
    elif xmin_spline<=x<=xmax_spline:
        return spline(x)
    else:
        return f0


# 2) Calculate Nuclear Cross Section (NCC) ###################################################################################################
def differential_NCC(T,beta,Z1,Z2,M2,f0=spline(xmin_spline)):  # T in MeV and M2 in MeV/c2 (from [2])
    # epsilon(beta):
    gam_ = gamma_(beta)
    a = a0*1e-2 * 0.8853 * (Z1**(2/3) + Z2**(2/3))**(-1/2)  # m
    epsilon = ((gam_-1)*a*M2) / (Z1*Z2*e**2)  # no unit
    # t(epsilon):
    Tmax = 4*M2*(gam_-1) # MeV
    t12 = epsilon * np.sqrt(T/Tmax) # assuming elastic collision
    # f(t):
    ft12 = f_universal(t12,f0)
    # sigma(t):
    dtdT = epsilon**2 / Tmax  # MeV^-1
    dsigma_ = np.pi * a**2 * dtdT * ft12 / (2 * t12**3)  # m^2 / MeV
    return epsilon, t12, ft12, dsigma_



### PLOTS #####################################################################################################################################
# a)
def plot_f_spline():
    plt.figure(figsize=(8, 5))
    x_dense = np.logspace(np.log10(np.min(x_plot)), np.log10(np.max(x_plot)), 500)
    y_dense = spline(x_dense)
    plt.plot(x_plot, y_plot, 'o', label='points extracted from fig.1')
    plt.plot(x_dense, y_dense, '-', label='Spline fit')
    plt.xlabel(r'$t^{1/2}$', fontsize=18)
    plt.ylabel(r'f($t^{1/2}$)', fontsize=18)
    plt.xscale('log')
    plt.legend(fontsize = 14)
    plt.grid(True)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig('f_t12_splinefit.png')
#plot_f_spline()


# b)
def plot_Ruth():
    x_spline = np.linspace(xmin_spline, xmax_spline, 500)
    x_ruth = np.linspace(1, 100, 200)
    y_spline = spline(x_spline)
    y_ruth = Rutherford(x_ruth)

    plt.figure(figsize=(8, 5))
    plt.plot(x_spline, y_spline, label=r"Spline: $f(x)$", lw=2)
    plt.plot(x_ruth, y_ruth, '--', color='black', label=r"Rutherford extension: $f_R(x)=\frac{0.48335}{x}$", lw=2)

    plt.xlabel(r'$x = \sqrt{t}$', fontsize=14)
    plt.ylabel(r'$f(x)$', fontsize=14)
    plt.legend(fontsize=14)
    plt.grid(True, which='both', ls=':')
    plt.xscale('log')
    plt.tight_layout()
    #plt.show()
    plt.savefig('f_Rutherford.png')
#plot_Ruth()


# c) 
def plot_f_full():
    plt.figure(figsize=(8, 5))
    x_list = np.logspace(-6, 2, 500)
    y = [f_universal(x) for x in x_list]
    plt.plot(x_list, y)
    plt.xlabel(r'$x=t^{1/2}$', fontsize=18)
    plt.ylabel(r'$f(x)$', fontsize=18)
    plt.xscale('log')
    #plt.legend(fontsize = 14)
    plt.grid(True)
    plt.tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig('f_t12_full.png')
#plot_f_full()


# d)
def plot_NCC(beta, z):
    gam_ = gamma_(beta)
    epsilonH = differential_NCC(1,beta,z,ZH,MH)[0]  # no unit
    epsilonO = differential_NCC(1,beta,z,ZO,MO)[0]  # no unit
    TmaxH = 4 * MH * (gam_ - 1)  # MeV
    TmaxO = 4 * MO * (gam_ - 1)  # MeV
    Tmin_plotH = xmin_spline**2 / epsilonH**2  # lowest value of T/Tm for hydrogen plotted in [1]
    Tmin_plotO = xmin_spline**2 / epsilonO**2  # same for oxygen
    Tmax_plotH = xmax_spline**2 / epsilonH**2  # highest
    Tmax_plotO = xmax_spline**2 / epsilonO**2

    xlistH = np.logspace(np.log10(Tmin_plotH), min(0,np.log10(Tmax_plotH)), 200)
    xlistO = np.logspace(np.log10(Tmin_plotO), min(0,np.log10(Tmax_plotO)), 200)
    y1_vals = []
    y2_vals = []
    y3_vals = []
    y4_vals = []
    for x in xlistH:
        epsilon, t12, ft12, sigH = differential_NCC(x*TmaxH, beta, z, ZH, MH)
        y1_vals.append(t12)
        y2_vals.append(ft12)
        y3_vals.append(sigH)
    for x in xlistO:
        sigO = differential_NCC(x*TmaxO, beta, z, ZO, MO)[-1]
        y4_vals.append(sigO)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax2 = ax.twinx()  # axe Y de droite

    ax.plot(xlistH, y3_vals, color='red', label=r'$\frac{d\sigma_{H}}{dT}$')
    ax.plot(xlistO, y4_vals, color='blue', label=r'$\frac{d\sigma_{O}}{dT}$')
    ax.set_ylabel(r'$\frac{d\sigma}{dT}$ / m$^2$MeV$^{-1}$', fontsize=16)
    ax.set_yscale('log')

    ax2.plot(xlistH, y1_vals, '-.', color='black', label=r'$t^{1/2} \sim \epsilon*\sqrt{\frac{T}{T_{\mathrm{max}}}}$')
    ax2.plot(xlistH, y2_vals, ':', color='black', label=r'$f(t^{1/2})$')
    ax2.set_ylabel('dimensionless variables', fontsize=16)

    ax.set_xlabel(r'T/T$_{max}$', fontsize=16)
    ax.set_xscale('log')

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=12, loc='upper right', bbox_to_anchor=(0.7, 1))

    textstr = '\n'.join((
        fr'$\beta$ = {beta:.3f}',
        fr'$z_{{SC}}$ = {z}',
        fr'$\epsilon_H$ = {epsilonH:.4f}',
        fr'$\epsilon_O$ = {epsilonO:.4f}',
        fr'$T_{{max,H}}$ = {TmaxH:.6f} MeV',
        fr'$T_{{max,O}}$ = {TmaxO:.6f} MeV'))
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.8)
    ax.text(0.30, 0.40, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)

    plt.tight_layout()
    #plt.show()
    plt.savefig(f'NCC_{z}.png')
#plot_NCC(1e-3, 10, spline(xmin_spline))


# e)
def plot_NCC_charge(beta, z1, z2, z3):
    gam_ = gamma_(beta)
    epsilonH = differential_NCC(1,beta,z1,ZH,MH)[0]  # no unit
    epsilonO = differential_NCC(1,beta,z1,ZO,MO)[0]  # no unit
    TmaxH = 4 * MH * (gam_ - 1)  # MeV
    TmaxO = 4 * MO * (gam_ - 1)  # MeV
    Tmin_plotH = xmin_spline**2 / epsilonH**2  # lowest value of T/Tm for hydrogen plotted in [1]
    Tmin_plotO = xmin_spline**2 / epsilonO**2
    Tmax_plotH = xmax_spline**2 / epsilonH**2
    Tmax_plotO = xmax_spline**2 / epsilonO**2

    xlistH = np.logspace(np.log10(Tmin_plotH), min(0,np.log10(Tmax_plotH)), 200)
    xlistO = np.logspace(np.log10(Tmin_plotO), min(0,np.log10(Tmax_plotO)), 200)

    y1_vals = []
    y2_vals = []
    y3_vals = []
    y4_vals = []

    for x in xlistH:
        sigH1 = differential_NCC(x*TmaxH, beta, z1, ZH, MH)[-1]
        y2_vals.append(sigH1)
        sigH2 = differential_NCC(x*TmaxH, beta, z2, ZH, MH)[-1]
        y3_vals.append(sigH2)
        sigH3 = differential_NCC(x*TmaxH, beta, z3, ZH, MH)[-1]
        y4_vals.append(sigH3)
    for x in xlistO:
        sigO = differential_NCC(x*TmaxO, beta, z1, ZO, MO)[-1]
        y1_vals.append(sigO)

    plt.figure(figsize=(8, 5))

    plt.plot(xlistH, y2_vals, color='red', label=fr'$\frac{{d\sigma_H}}{{dT}}$(z={z1})')
    plt.plot(xlistH, y3_vals, 'r--', label=fr'$\frac{{d\sigma_H}}{{dT}}$(z={z2})')
    plt.plot(xlistH, y4_vals, 'r-.', label=fr'$\frac{{d\sigma_H}}{{dT}}$(z={z3})')
    plt.plot(xlistO, y1_vals, color='blue', label=fr'$\frac{{d\sigma_O}}{{dT}}$(z={z1})')
    plt.ylabel(r'$\frac{d\sigma}{dT}$ / m$^2$MeV$^{-1}$', fontsize=16)
    plt.yscale('log')
    
    plt.xlabel(r'T/T$_{max}$', fontsize=16)
    plt.xscale('log')
    plt.legend(fontsize=12, loc='upper right')

    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.8)
    plt.text(0.30, 0.40, r'$\beta = 0.001$', transform=plt.gca().transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    #plt.show()
    plt.savefig(f'NCC_{z1}_{z2}_{z3}.png')
#plot_NCC_charge(1e-3, 1, 10, 100)


# f)
def plot_f0_dependence(beta,z):
    gam_ = gamma_(beta)
    TmaxH = 4 * MH * (gam_ - 1)  # MeV
    xlistH = np.logspace(0,-10, 200)
    y1_vals = []
    y2_vals = []
    y3_vals = []
    y4_vals = []
    for x in xlistH:
        sigH1 = differential_NCC(x*TmaxH, beta, z, ZH, MH, f0=0.0001)[-1]
        y1_vals.append(sigH1)
        sigH2 = differential_NCC(x*TmaxH, beta, z, ZH, MH, f0=0.01)[-1]
        y2_vals.append(sigH2)
        sigH3 = differential_NCC(x*TmaxH, beta, z, ZH, MH, f0=0.1)[-1]
        y3_vals.append(sigH3)
        sigH4 = differential_NCC(x*TmaxH, beta, z, ZH, MH, f0=0.5)[-1]
        y4_vals.append(sigH4)


    plt.figure(figsize=(8, 5))
    plt.plot(xlistH, y1_vals, label=r'$f_0 = 0.0001$')
    plt.plot(xlistH, y2_vals, label=r'$f_0 = 0.01$')
    plt.plot(xlistH, y3_vals, label=r'$f_0 = 0.1$')
    plt.plot(xlistH, y4_vals, label=r'$f_0 = 0.5$')
    plt.ylabel(r'$\frac{d\sigma_H}{dT}$ / m$^2$MeV$^{-1}$', fontsize=16)
    plt.yscale('log')
    
    plt.xlabel(r'T/T$_{max}$', fontsize=16)
    plt.xscale('log')
    plt.legend(fontsize=12, loc='upper right')

    textstr = '\n'.join((
        fr'$\beta$ = {beta:.3f}',
        fr'$z_{{SC}}$ = {z}'))
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.8)
    plt.text(0.30, 0.20, textstr, transform=plt.gca().transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    #plt.show()
    plt.savefig(f'NCC_f0_dependence.png')
#plot_f0_dependence(1e-3,10)



# Notes: For modern values of stopping power and cross section, look at SRIM data [3], which contains many experimental data and is updated every 6 months
#### REFERENCES #############################################################################################################################################################################################################################################################
#[1]: https://wiki.icecube.wisc.edu/index.php/Ice_Density
#[2]: J. LINDHARD, M. SCHARFF(t) AND H. E. SCHIØTT, "RANGE CONCEPTS AND HEAVY ION RANGES (NOTES ON ATOMIC COLLISIONS, II)"
#[3]: Ziegler et al., The Stopping and Range of Ions in Solids (SRIM framework)