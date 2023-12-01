init_distribution.py: Create an initial restart file for GyselaX

`python init_distribution.py --help` for help
or
`python init_distribution.py`
which is equivalent to:
`python init_distribution.py -i input_params.yaml`
or
`python init_distribution.py -i input_params.yaml -o GyselaX_restart_00000.nc`

By default the results is saved in NetCDF format. If you want to save the results in HDF5 format you have to specify an '.h5' extension for the output file name as for instance:
`python init_distribution.py -i <YAMLinputFileName> -o GyselaX_restart_00000.h5`

