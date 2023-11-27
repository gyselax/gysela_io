import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# Construct a Maxellian distribution function
def Maxwellian_func( As_loc, N_loc, Upar_loc, T_loc, B_loc, vpar_loc, mu_loc):

    # Compute  Energy = (0.5 * (vpar_loc - Upar_loc)**2 + mu_loc * B_loc) / T_loc in a vectorize way
    Energy = np.add.outer(0.5 * (vpar_loc - Upar_loc)**2,  mu_loc * B_loc) / T_loc

    return N_loc * (As_loc / (2.0 * np.pi * T_loc))**1.5 * np.exp( - Energy)


def plot_profiles( ds_gysela, figure_name ):

    # Figure : Profiles
    fig1 = plt.figure(1,figsize=(10, 10))

    ax11 = fig1.add_subplot(311)
    ax11.plot(ds_gysela.grid_tor1.data, ds_gysela.densityTorCS.data[:, 0])
    ax11.set_ylabel("$N_e / N_e^{ref}$", fontsize = 20)

    ax12 = fig1.add_subplot(312)
    ax12.plot(ds_gysela.grid_tor1.data, ds_gysela.UparTorCS.data[:,0])
    ax12.set_ylabel("$U_{\parallel, e} / v_{Te}^{ref}$", fontsize = 20)

    ax13 = fig1.add_subplot(313)
    ax13.plot(ds_gysela.grid_tor1.data, ds_gysela.temperatureTorCS.data[:,0])
    ax13.set_ylabel("$T_e / T_e^{ref}$", fontsize = 20)

    print('--> Saving of the figure {}'.format(figure_name))
    plt.savefig(figure_name)
    
