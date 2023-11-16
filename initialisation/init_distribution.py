###############################################################
# This python script initialise the distribution function that
# will be read by GYSELA
###############################################################
### modules ###
import numpy as np
import matplotlib.pyplot as plt

########### TEMPORARY: will be changed to use YAML file ######
### The grid ###
Ntor1 = 16       # number of cells in the radial direction
Ntor2 = 16       # number of cells in the poloidal direction
Ntor3 = 16       # number of cells in the toroidal direction
Nvpar = 128     # number of cells in the parallel velocity direction
Nmu   = 64      # number of cells in the adiabatic invariant direction


tor1_min = 0.0
tor1_max = 1.0
tor2_min = 0.0
tor2_max = 2*np.pi
tor3_min = 0.0
tor3_max = 2*np.pi
vpar_max = 5.0
mu_max   = 12.0


grid_tor1 = np.linspace( tor1_min, tor1_max, Ntor1 )
grid_tor2 = np.linspace( tor2_min, tor2_max, Ntor2 )
grid_tor3 = np.linspace( tor3_min, tor3_max, Ntor3 )
grid_vpar = np.linspace( -vpar_max, vpar_max, Nvpar )
grid_mu   = np.linspace( 0.0, mu_max, Nmu )


## Normalisation
#R_ref = 1.0    # Major radius of reference  [m]
#N_ref = 1e19   # Electron density in R_ref [m^-3]
#T_ref = 1e3    # Electron temperature in R_ref [eV]
#B_ref = 1.0    # Magnetic field on R_ref [T]


# species caracteristics
As = 1
Zs = 1


# Construct a Maxellian distribution function
def Maxwellian_func( As_loc, N_loc, Upar_loc, T_loc, B_loc, vpar_loc, mu_loc):

#   Compute  Energy = (0.5 * (vpar_loc - Upar_loc)**2 + mu_loc * B_loc) / T_loc in a vectorize way
    Energy = np.add.outer(0.5 * (vpar_loc - Upar_loc)**2,  mu_loc * B_loc) / T_loc

    return N_loc * (As / (2.0 * np.pi * T_loc))**1.5 * np.exp( - Energy)


# Construct the density, mean parallel velocity and temperature profiles in 1D, assuming a parabolic radial dependance:
N_vec    = np.zeros(Ntor1)
Upar_vec = np.zeros(Ntor1)
T_vec    = np.zeros(Ntor1)

N_min = 0.5  #In normalized units, N_max=N_ref is on the axis
T_min = 0.2  #In normalized units, T_max=T_ref is on the axis
    
N_vec = 1.0  - (1.0 - N_min) * grid_tor1**2
T_vec = 1.0  - (1.0 - T_min) * grid_tor1**2


# Construction of the 5D distribution function
F_distribution_5D = np.zeros( (Ntor1, Ntor2, Ntor3) + (Nvpar, Nmu) )

for ir in range(Ntor1):
    N_loc    = N_vec[ir]
    Upar_loc = Upar_vec[ir]
    T_loc    = T_vec[ir]

    Maxwellian_loc = Maxwellian_func( As, N_loc, Upar_loc, T_loc, 1.0, grid_vpar, grid_mu)

    for itheta in range(Ntor2):
     
        for iphi in range(Ntor3):


            F_distribution_5D[ir, itheta, iphi, :, :] = Maxwellian_loc



plotting = input("Do you want to plot the figures? [y/n] (default = n)")
if plotting == "y":
    
    
    # Figure 0: Distribution function on the magnetic axis

    Maxwellian_loc[:, :] = F_distribution_5D[ 0, 0, 0, :, :]

    fig0 = plt.figure(0,figsize=(10, 10))
    ax0 = fig0.add_subplot(111)
    ax0.pcolor(grid_vpar,grid_mu,np.transpose(Maxwellian_loc))
    ax0.set_xlabel("$v_{\parallel}$")
    ax0.set_ylabel("$\mu$")
 

    # Figure 1: Profiles
    fig1 = plt.figure(1,figsize=(10, 10))

    ax11 = fig1.add_subplot(311)
    ax11.plot(grid_tor1, N_vec)
    ax11.set_ylabel("$N_e / N_e^{ref}$", fontsize = 20)

    ax12 = fig1.add_subplot(312)
    ax12.plot(grid_tor1, Upar_vec)
    ax12.set_ylabel("$U_{\parallel, e} / v_{Te}^{ref}$", fontsize = 20)

    ax13 = fig1.add_subplot(313)
    ax13.plot(grid_tor1, T_vec)
    ax13.set_ylabel("$T_e / T_e^{ref}$", fontsize = 20)

    plt.show()
