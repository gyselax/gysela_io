import h5py
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# Construct a Maxellian distribution function
def Maxwellian_func( As_loc, N_loc, Upar_loc, T_loc, B_loc, vpar_loc, mu_loc):
    """
    Construct a Maxellian distribution function
    """

    # Compute  Energy = (0.5 * (vpar_loc - Upar_loc)**2 + mu_loc * B_loc) / T_loc in a vectorize way
    Energy = np.add.outer(0.5 * (vpar_loc - Upar_loc)**2,  mu_loc * B_loc) / T_loc

    return N_loc * (As_loc / (2.0 * np.pi * T_loc))**1.5 * np.exp( - Energy)


# Plot the cross-section of the 3 fluid moments: density, Upar and temperature
def plot_profiles( ds_gysela, figure_name ):
    """
    Plot the cross-section of the 3 fluid moments: density, Upar and temperature
    """

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


# Create an Xarray from the reading of an HDF5 GyselaX restart file
def create_Xarray_from_hdf5_restartfile( hdf5_restart_file ):
  """ 
  Create an Xarray from the reading of an HDF5 GyselaX restart file
  """

  if hdf5_restart_file.suffix == '.nc':
    ds_gysela = xr.open_dataset(hdf5_restart_file)
  elif hdf5_restart_file.suffix == '.h5':
    with h5py.File(hdf5_restart_file, "r") as h5f:
      if 'species_name' in h5f.keys():
        species_name = h5f['species_name'][()]
        for i in range(0,len(species_name)):
          species_name[i]=species_name[i].decode('utf-8')
      ds_gysela = xr.Dataset(
          data_vars=dict(
          densityTorCS=(['tor2','tor1'], h5f['densityTorCS'][()]),
          UparTorCS=(['tor2','tor1'], h5f['UparTorCS'][()]),
          temperatureTorCS=(['tor2','tor1'], h5f['temperatureTorCS'][()]),
          fdistribu=(['species', 'tor3', 'tor2', 'tor1', 'vpar', 'mu'], h5f['fdistribu'][()] )
          ),
          coords=dict(
            tor1=h5f['grid_tor1'][()],
            tor2=h5f['grid_tor2'][()],
            tor3=h5f['grid_tor3'][()],
            vpar=h5f['grid_vpar'][()],
            mu=h5f['grid_mu'][()],
            species=h5f['species'][()],
          ),
          attrs=dict(description="Mesh and initial profiles of GyselaX"),
      )

  return ds_gysela
