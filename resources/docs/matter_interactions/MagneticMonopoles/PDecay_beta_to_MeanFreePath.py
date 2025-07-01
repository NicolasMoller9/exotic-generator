# Author: Nicolas Moller - PhD student (Chiba U.)

import matplotlib.pyplot as plt
import numpy as np


# 0) USEFUL CONSTANTS AND FUNCTIONS ##################################################################################################################################

# Useful constants:
NA = 6.022*1e23                          # mol^-1 (Avogadro number)

# Ice properties:
rho_ice = 9.17e5                         # g/m^3 (Ice density. From [1])
AH = 1.008                               # g/mol (molar mass of hydrogen)
AO = 15.999                              # g/mol (molar mass of oxygen)  
A_H2O = 18.015                           # g/mol (molar mass of H2O)
A_He = 4.0026                            # g/mol (molar mass of Helium)
n_H2O = rho_ice/(A_H2O/NA)               # m^-3 (number density of H2O molecule)
n_nucl = 18 * n_H2O                      # m^-3 (number density of nucleons)

# 1) Mean Free Path calculation ########################################################################################################################################

# There is a suppression factor for the monopole catalysis due to the distortion of the wave function:
# F ~ (beta/beta_0)^(2*nu) 
def beta0(A): 
    beta0_He = 0.0275 # from [2]
    return beta0_He * (round(A_He)**(1/3) * A_He) / (round(A)**(1/3) * A) 
beta0_H = beta0(AH)
beta0_O = beta0(AO)
Renu_H = -0.5         # from [2]
Renu_O = 3.123/2      # from [2]


def beta_to_lamda(beta, sig0=1e-31):
    if beta >= beta0_H:
        fH = 1
        fO = 1
    elif beta0_O <= beta < beta0_H:
        fH = (beta / beta0_H) ** (2 * Renu_H)
        fO = 1
    else:  # beta < beta0_O
        fH = (beta / beta0_H)**(2 * Renu_H)
        fO = (beta / beta0_O)**(2 * Renu_O)
    F = (2 / 18) * fH + (16 / 18) * fO
    sigma_cat = (sig0 / beta) * F            # from [2]
    meanFreePath = 1 / (sigma_cat * n_nucl)
    meanFreePath_bis = beta / (sig0 * n_nucl)
    return meanFreePath, meanFreePath_bis, F


### PLOTS #####################################################################################################################################
# a)
def plotter():
    beta_list = np.logspace(np.log10(0.001), np.log10(0.9), num=30)
    sigma0 = 1e-31  # m^2 (from [3])

    lam_list1, lam_list2, F_list = zip(*[beta_to_lamda(b, sigma0) for b in beta_list])

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(beta_list, lam_list1, label='sig_cat = sig0/beta * F', color='orange')
    ax1.plot(beta_list, lam_list2, label='sig_cat = sig0/beta', color='red')
    ax2.plot(beta_list, F_list, color='blue', linestyle='dashed', linewidth=0.2, label='Suppression factor (F)')

    ax1.axvline(x=beta0_H, color='black', linestyle='--', label=f'beta0_H = {beta0_H:.2e}')
    ax1.axvline(x=beta0_O, color='black', linestyle='dotted', label=f'beta0_O = {beta0_O:.2e}')

    ax1.set_xlabel('beta')
    ax1.set_ylabel('mean free path [m]', color='black')
    ax2.set_ylabel('F', color='blue')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('linear')

    plt.title(fr'$sigma_0={sigma0}$ m${{^2}}$')
    ax1.legend(loc='upper left')
    ax2.legend(loc='lower right', bbox_to_anchor=(0.7,0.2))

    #plt.show()
    plt.savefig('plots/suppression_factor.png')
#plotter()


# b)
def plotter2():
    beta_list = np.logspace(np.log10(0.0001), np.log10(0.1), num=300)
    sigma0_list = [1e-31, 1e-32, 1e-33, 1e-34, 1e-35]  # m2, from [2]

    fig, ax1 = plt.subplots()
    for sig in sigma0_list:
        lam_list = [beta_to_lamda(b, sig)[0] for b in beta_list]
        plt.plot(beta_list, lam_list, label=fr'$\sigma_0 = {sig} m^2$')

    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\lambda$ / m')

    plt.xscale('log')
    plt.yscale('log')

    plt.title("Mean free path of proton decay catalysis")
    ax1.legend(loc='upper left')

    #plt.show()
    plt.savefig('plots/beta_lamda.png')
#plotter2()


#### REFERENCES #############################################################################################################################################################################################################################################################
#[1]: https://wiki.icecube.wisc.edu/index.php/Ice_Density
#[2]: J. Arafune &  M. Fukugita, "Velocity-Dependent Factors for the Rubakov Process for Slowly Moving Magnetic Monopoles in Matter", PHYSICAL REVIEW LETTERS VOL50 NBR24 (1983)
#[3]: Kolb, E.W. (1984). Monopole Catalyzed Nucleon Decay: The Astrophysical Connection. In: Stone, J.L. (eds) Monopole ’83. NATO ASI Series, vol 111. Springer, Boston, MA. https://doi.org/10.1007/978-1-4757-0375-7_26
