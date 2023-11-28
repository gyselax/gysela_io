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
    parser.add_argument('-i','--input_files',
                        action='store',
                        nargs='?',
                        default='GyselaX_restart_00000.h5 GyselaX_restart_00001.h5',
                        type=str,
                        required=False,
                        help='two NetCDF or HDF5 restart files (default: GyselaX_restart_00000.h5 GyselaX_restart_00001.h5)')
    parser.add_argument('--h5diff', action='store_true',help='Apply h5diff to the two input files')
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

    # Read the two input files
    str_input_files = args.input_files
    list_input_files = str_input_files.split(' ')
    if len(list_input_files) != 2:
      print('Error: Two input files required')
      sys.exit()
    input_file1 = Path(list_input_files[0])
    input_file2 = Path(list_input_files[1])
    if not input_file1.exists():
      print('The file {} does not exist'.format(input_file1))
      sys.exit()
    elif not input_file2.exists():
      print('The file {} does not exist'.format(input_file2))
      sys.exit()

    if args.h5diff:
      os.system('h5diff {} {}'.format(input_file1, input_file2))
      sys.exit()
    
    # Create the two Xarray
    ds_gysela1 = fut.create_Xarray_from_hdf5_restartfile( input_file1 )
    ds_gysela2 = fut.create_Xarray_from_hdf5_restartfile( input_file2 )
    
    # Read the dataset 
    dataset_name = args.dataset
    if (dataset_name not in ds_gysela1.keys()) or (dataset_name not in ds_gysela2.keys()):
      print('The dataset {} is not in one of the two files'.format(dataset_name))
      sys.exit()

    # Read the selection for the slice and create the corresponding dictionary
    str_select = args.select.replace(" ","")
    list_select_OK = []
    str_select_OK = ''
    for i in str_select.split(','):
      if i.split('=')[0] in ds_gysela1[dataset_name].coords:
        list_select_OK.append(i)
        str_select_OK = str_select_OK+'_i'+str(i).replace('=','eq')
    dict_select = {key: int(value) for key, value in (pair.split('=') for pair in list_select_OK)}

    data1_select = ds_gysela1[dataset_name].isel(dict_select)
    data2_select = ds_gysela2[dataset_name].isel(dict_select)
    diff_select = data1_select - data2_select
    nb_dim_select = len(diff_select.dims)
    str_dim_select = ''
    for idim in data1_select.dims:
      str_dim_select = str_dim_select + '_' + idim

    # Read the output file name
    if args.output is None:
      figname = 'diff_{}{}{}.png'.format(dataset_name,str_dim_select,str_select_OK)
      figname = figname.replace('species','sp')
      figname = figname.replace('fdistribu','f')
    else:
      figname = args.output

    print(diff_select)
    if nb_dim_select==0 or nb_dim_select>2:
      print('--> no plot for {}D slice is saved'.format(nb_dim_select))
    else:
      ax = plt.axes()
      diff_select.plot()
      print('--> Figure saved in: {}'.format(figname))
      plt.savefig(figname)
