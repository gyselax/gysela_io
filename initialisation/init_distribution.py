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

import sys
import os

# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))

# Getting the parent directory name
# where the current directory is present.
parent = os.path.dirname(current)

# adding the parent directory to
# the sys.path.
sys.path.append(parent)

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

    # Geometry parameters
    R0 = d_params['Geometry']['R0']   # major radius
    a  = d_params['Geometry']['a']    # minor radius
    Ea = d_params['Geometry']['Ea']   # E(a) for Culham equilibrium
    Ta = d_params['Geometry']['Ta']   # T(a) for Culham equilibrium
    q0 = d_params['Geometry']['q0']   # safety factor on the magnetic axis (parabolic profile)
    qa = d_params['Geometry']['qa']   # safety factor at the edge (parabolic profile)


    # species caracteristics (filled with dummy values now)
    nspecies = len(d_params['SpeciesInfo'])
    As = np.ones(nspecies)
    Zs = np.ones(nspecies)
    Ns_min = np.ones(nspecies)
    Ts_min = np.ones(nspecies)
    species_name = [];
    ispecies = np.arange(nspecies);
    for ispec in range(nspecies):
      species_name.append(d_params['SpeciesInfo'][ispec]['name'])
      As[ispec] = d_params['SpeciesInfo'][ispec]['mass']
      Zs[ispec] = d_params['SpeciesInfo'][ispec]['charge']
      Ns_min[ispec] = d_params['SpeciesInfo'][ispec]['N_min']
      Ts_min[ispec] = d_params['SpeciesInfo'][ispec]['T_min']

    # Construct the density, mean parallel velocity and temperature profiles in 1D, 
    N_vec    = np.zeros((nspecies, ncell_tor2, ncell_tor1), dtype=float)
    Upar_vec = np.zeros((nspecies, ncell_tor2, ncell_tor1), dtype=float)
    T_vec    = np.zeros((nspecies, ncell_tor2, ncell_tor1), dtype=float)

    for ispec in range(nspecies):    
        for ir in range(ncell_tor1):
            N_vec[ispec, :, ir] = fut.parabolic_prof( 1.0, Ns_min[ispec], grid_tor1[ir] )
            T_vec[ispec, :, ir] = fut.parabolic_prof( 1.0, Ts_min[ispec], grid_tor1[ir] )

    # Construction of the 5D distribution function
    F_distribution_5D = np.zeros( (nspecies, ncell_tor3, ncell_tor2, 
      ncell_tor1, ncell_vpar, ncell_mu), dtype=float )

    for ispec in range(nspecies):
      As_loc = As[ispec]
      for ir in range(ncell_tor1):
          for itheta in range(ncell_tor2):
              N_loc    = N_vec[ispec, itheta, ir]
              Upar_loc = Upar_vec[ispec, itheta, ir]
              T_loc    = T_vec[ispec, itheta, ir]

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
        densityTorCS=(['species','tor2', 'tor1'], N_vec),
        UparTorCS=(['species','tor2', 'tor1'], Upar_vec),
        temperatureTorCS=(['species','tor2', 'tor1'], T_vec),
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
        h5file.create_dataset('masses', data=As)
        h5file.create_dataset('charges', data=Zs)
        h5file.create_dataset('species_name', data=species_name)
        h5file.create_dataset('time_saved', data=time_saved)
        h5file.create_dataset('densityTorCS', data=N_vec)
        h5file.create_dataset('UparTorCS', data=Upar_vec)
        h5file.create_dataset('temperatureTorCS', data=T_vec)
        h5file.create_dataset('fdistribu', data=F_distribution_5D)

    end = time.time()
    print('The time of execution of above program is : {} s'.format(end-start))

