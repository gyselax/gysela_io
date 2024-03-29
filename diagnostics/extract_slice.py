# SPDX-License-Identifier: MIT

"""
Created on 2023-11-23

 This python script initialise the distribution function that
 will be read by GYSELA

@author: P. Donnel & V. Grandgirard
"""

import h5py
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import xarray as xr
import yaml

from argparse import ArgumentParser
from pathlib import Path
from yaml.loader import SafeLoader

import utils.function_utils as fut

if __name__ == '__main__':
    parser = ArgumentParser(description="Extract slice from 5D distribution function")
    parser.add_argument('-i','--input_file',
                        action='store',
                        nargs='?',
                        default=Path('GyselaX_restart_00000.nc'),
                        type=Path,
                        required=False,
                        help='NetCDF or HDF5 file containing the distribution function')
    parser.add_argument('--keys', action='store_true',help='List the dataset contained in the input file')
    parser.add_argument('-d','--dataset',
                        action='store',
                        nargs='?',
                        default='fdistribu',
                        type=str,
                        help='Name of the dataset to plot')
    parser.add_argument('-s','--select',
                        action='store',
                        nargs='?',
                        default='tor1=0,tor2=0,tor3=0,species=0',
                        type=str,
                        help="Position of the slice to extract, by default 'tor1=0,tor2=0,tor3=0,species=0'")
    parser.add_argument('-o','--output',
                        action='store',
                        nargs='?',
                        default=None,
                        type=str,
                        help='output PNG filename for the plot')

    args = parser.parse_args()

    # Read the input file
    input_file = args.input_file
    dataset_name = args.dataset
    str_select = args.select.replace(" ","")

    if args.keys:
      os.system('h5ls {}'.format(input_file))
      sys.exit()

    # Create the Xarray
    ds_gysela = fut.create_Xarray_from_hdf5_restartfile( input_file )
    
    list_select_OK = []
    str_select_OK = ''
    for i in str_select.split(','):
      if i.split('=')[0] in ds_gysela[dataset_name].coords:
        list_select_OK.append(i)
        str_select_OK = str_select_OK+'_i'+str(i).replace('=','eq')

    dict_select = {key: int(value) for key, value in (pair.split('=') for pair in list_select_OK)}
    ds_select = ds_gysela[dataset_name].isel(dict_select)
    nb_dim_select = len(ds_select.dims)
    str_dim_select = ''
    for idim in ds_select.dims:
      str_dim_select = str_dim_select + '_' + idim

    # Read the output file name
    if args.output is None:
      figname = '{}{}{}.png'.format(dataset_name,str_dim_select,str_select_OK)
      figname = figname.replace('species','sp')
      figname = figname.replace('fdistribu','f')
    else:
      figname = args.output

    print(ds_select)
    if nb_dim_select==0 or nb_dim_select>2:
      print('--> no plot for {}D slice is saved'.format(nb_dim_select))
    else:
      ax = plt.axes()
      ds_select.plot()
      print('--> Figure saved in: {}'.format(figname))
      plt.savefig(figname)
