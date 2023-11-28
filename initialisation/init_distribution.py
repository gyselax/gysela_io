# SPDX-License-Identifier: MIT

"""
Created on 2023-11-17

 This python script initialise the distribution function that
 will be read by GYSELA

@author: P. Donnel & V. Grandgirard
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import xarray as xr
import yaml

from argparse import ArgumentParser
from pathlib import Path
from yaml.loader import SafeLoader

import utils.function_utils as fut
import utils.plot_utils as put
import utils.math_utils as mut

if __name__ == '__main__':
    parser = ArgumentParser(description="Initialisation of the distribution function for GyselaX")
    parser.add_argument('-i','--input_file',
                        action='store',
                        nargs='?',
                        default=Path('input_params_ref.yaml'),
                        type=Path,
                        help='input YAML file (default: input_params_ref.yaml)')
    parser.add_argument('-o','--output_file',
                        action='store',
                        nargs='?',
                        default=Path('GyselaX_restart_00000.nc'),
                        type=Path,
                        help='output NetCDF or HDF5 file (default: GyselaX_restart_00000.nc)')
    args = parser.parse_args()

    start = time.time()

    # Open the input YAML file and load the file
    input_file = args.input_file
    input_filename = input_file.name
    print('Input file={}'.format(input_filename))
    output_file = args.output_file

    with open(input_file) as f:
        d_params = yaml.load(f, Loader=SafeLoader)
    print(d_params)

    time_saved = 0.
    ncell_tor1 = d_params['Mesh']['ncell_tor1']
    ncell_tor2 = d_params['Mesh']['ncell_tor2']
    ncell_tor3 = d_params['Mesh']['ncell_tor3']
    ncell_vpar = d_params['Mesh']['ncell_vpar']
    ncell_mu = d_params['Mesh']['ncell_mu']
    
    min_tor1 = d_params['Mesh']['min_tor1']
    max_tor1 = d_params['Mesh']['max_tor1']
    min_tor2 = d_params['Mesh']['min_tor2']
    max_tor2 = d_params['Mesh']['max_tor2']
    min_tor3 = d_params['Mesh']['min_tor3']
    max_tor3 = d_params['Mesh']['max_tor3']
    max_vpar = d_params['Mesh']['max_vpar']
    max_mu = d_params['Mesh']['max_mu']

    grid_tor1 = np.linspace( min_tor1, max_tor1, ncell_tor1 )
    grid_tor2 = np.linspace( min_tor2, max_tor2, ncell_tor2, endpoint = False )
    grid_tor3 = np.linspace( min_tor3, max_tor3, ncell_tor3, endpoint = False )
    grid_vpar = np.linspace( -max_vpar, max_vpar, ncell_vpar )
    grid_mu   = np.linspace( 0.0, max_mu, ncell_mu )

#VG##VG#grid_grev_tor1 = mut.get_greville_points(grid_tor1, periodic= False, spline_degree = 3)
#VG##VG#grid_grev_tor2 = mut.get_greville_points(grid_tor2, periodic= True, spline_degree = 3)
#VG##VG#print("grid_tor2=",grid_tor2)
#VG##VG#print("grid_grev_tor2=",grid_grev_tor2)

    # species caracteristics (filled with dummy values now)
    nspecies = len(d_params['SpeciesInfo'])
    As = np.ones(nspecies)
    Zs = np.ones(nspecies)
    species_name = [];
    ispecies = np.arange(nspecies);
    for ispec in range(nspecies):
      species_name.append(d_params['SpeciesInfo'][ispec]['name'])
      As[ispec] = d_params['SpeciesInfo'][ispec]['mass']
      Zs[ispec] = d_params['SpeciesInfo'][ispec]['charge']

    # Construct the density, mean parallel velocity and temperature profiles in 1D, 
    #  assuming a parabolic radial dependance:
    N_vec    = np.zeros((ncell_tor2, ncell_tor1), dtype=float)
    Upar_vec = np.zeros((ncell_tor2, ncell_tor1), dtype=float)
    T_vec    = np.zeros((ncell_tor2, ncell_tor1), dtype=float)

    N_min = 0.5  #In normalized units, N_max=N_ref is on the axis
    T_min = 0.2  #In normalized units, T_max=T_ref is on the axis
    
    for ir in range(ncell_tor1):
      N_vec[:, ir] = 1.0  - (1.0 - N_min) * grid_tor1[ir]**2
      T_vec[:, ir] = 1.0  - (1.0 - T_min) * grid_tor1[ir]**2

    # Construction of the 5D distribution function
    F_distribution_5D = np.zeros( (nspecies, ncell_tor3, ncell_tor2, 
      ncell_tor1, ncell_vpar, ncell_mu), dtype=float )

    for ispec in range(nspecies):
      As_loc = As[ispec]
      for ir in range(ncell_tor1):
          for itheta in range(ncell_tor2):
              N_loc    = N_vec[itheta, ir]
              Upar_loc = Upar_vec[itheta, ir]
              T_loc    = T_vec[itheta, ir]

              Maxwellian_loc = fut.Maxwellian_func( As_loc, N_loc, Upar_loc, T_loc, 1.0, grid_vpar, grid_mu)

              for iphi in range(ncell_tor3):
                  F_distribution_5D[ispec, iphi, itheta, ir, :, :] = Maxwellian_loc

    # Create the associated xarray
    ds_gysela = xr.Dataset(
        data_vars=dict(
        grid_tor1=(['tor1'], grid_tor1),
        grid_tor2=(['tor2'], grid_tor2),
        grid_tor3=(['tor3'], grid_tor3),
        grid_vpar=(['vpar'], grid_vpar),
        grid_mu=(['mu'], grid_mu),
        densityTorCS=(['tor2', 'tor1'], N_vec),
        UparTorCS=(['tor2', 'tor1'], Upar_vec),
        temperatureTorCS=(['tor2', 'tor1'], T_vec),
        species_name=(['species'], species_name),
        masses=(['species'], As),
        charges=(['species'], Zs),
        fdistribu=(['species', 'tor3', 'tor2', 'tor1', 'vpar', 'mu'], F_distribution_5D )
        ),
        coords=dict(
          tor1=grid_tor1,
          tor2=grid_tor2,
          tor3=grid_tor3,
          vpar=grid_vpar,
          mu=grid_mu,
          species=ispecies,
          time_saved=time_saved
        ),
        attrs=dict(description="Mesh and initial profiles of GyselaX"),
    )

    if output_file.suffix == '.nc':
      #--> Saving in NetCDF file 
      print('--> Save NetCDF file {}'.format(output_file.name))
      ds_gysela.to_netcdf( path=output_file, mode='w', engine='h5netcdf' );
      ds_gysela.close();
    elif output_file.suffix == '.h5':
      #--> Saving in HDF5 files
      print('--> Save HDF5 file {}'.format(output_file.name))
      with h5.File(output_file.name, 'w') as h5file:
        h5file.create_dataset('grid_tor1', data=grid_tor1)
        h5file.create_dataset('grid_tor2', data=grid_tor2)
        h5file.create_dataset('grid_tor3', data=grid_tor3)
        h5file.create_dataset('grid_vpar', data=grid_vpar)
        h5file.create_dataset('grid_mu', data=grid_mu)
        h5file.create_dataset('species', data=ispecies)
        h5file.create_dataset('charges', data=As)
        h5file.create_dataset('masses', data=Zs)
        h5file.create_dataset('species_name', data=species_name)
        h5file.create_dataset('time_saved', data=time_saved)
        h5file.create_dataset('densityTorCS', data=N_vec)
        h5file.create_dataset('UparTorCS', data=Upar_vec)
        h5file.create_dataset('temperatureTorCS', data=T_vec)
        h5file.create_dataset('fdistribu', data=F_distribution_5D)

    end = time.time()
    print('The time of execution of above program is : {} s'.format(end-start))

