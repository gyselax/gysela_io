# SPDX-License-Identifier: MIT

'''
General functions used for plot
'''

# pylint: disable=import-error

from io import BytesIO
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import imageio

def plot_field2d(data, titlename, figure_name, **opt_args):
    """
    Plots an image of a 2D field.

    Parameters
    ----------
    figure_name: pathlib.Path or str
        path of the file to generate

    data: xarray.DataArray
        2D data to plot
    """

    if data.ndim != 2:
        raise ValueError("Expected 2D data, got {}D, {}".format(
            data.ndim, data.coords))

    figsize = opt_args.get('figsize', (10, 10))
    fontsize = opt_args.get('fontsize', figsize[0]*2)
    title = opt_args.get('title', data.name+' '+titlename)
    axis_font = opt_args.get(
        'axis_font', {'size': fontsize})
    title_font = opt_args.get('title_font', {'size': fontsize,
                                             'verticalalignment': 'bottom'})
    scale = opt_args.get('scale', '')
    data2D = data.values
    if scale == 'log':
        print('log')
        data2D = abs(np.log10(data2D))
    vmin = opt_args.get('vmin', np.nanmin(data2D))
    vmax = opt_args.get('vmax', np.nanmax(data2D))

    # Plotting with naive matplotlib
    mpl.style.use('classic')
    plt.rc('xtick', labelsize=fontsize*.8)
    plt.rc('ytick', labelsize=fontsize*.8)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.pcolormesh(data.coords[data.dims[1]],
                       data.coords[data.dims[0]],
                       data2D,
                       vmin=vmin,
                       vmax=vmax,
                       shading='auto',
                       cmap='jet',
                       rasterized=True)
    ax.set_title(title, **title_font)
    ax.set_xlabel(r'${}$'.format(data.dims[1]), **axis_font)
    ax.set_ylabel(r'${}$'.format(data.dims[0]), **axis_font)
    ax.axis('tight')
    fig.colorbar(im, ax=ax)

    print('--> Save figure {}.png'.format(figure_name))
    fig.savefig(figure_name)
    plt.close(fig)



def anim_field3d_append(data, out_writer, anim_dim='time', **opt_args):
    """
    Plots an animation of a 3D field over multiple images.

    The function can be used to add images to an already partially
    plotted animation.

    Parameters
    ----------
    out_writer: imageio.Writer
        image to append

    data: xarray.DataArray
        3D data to plot, the 1st

    anim_dim: str
        the dimension over which to generate the animation
    """
    if data.ndim != 3:
        raise ValueError("Expected 3D data, got {}D, {}".format(
            data.ndim, data.coords))

    for it in range(data.sizes[anim_dim]):
        it_val = float(data.coords[anim_dim][it])
        data_it = data.isel(indexers={anim_dim: it}, drop=True)
        with BytesIO() as imgdata:
            plot_field2d(data_it, figure_name=imgdata,
                         titlename=r'${} = {:02.2f}$'.format(
                         anim_dim, it_val),
                         **opt_args)
            imgdata.seek(0)
            out_writer.append_data(imageio.imread(imgdata))


def anim_field3d(data, out_file, **opt_args):
    """
    Plots an animation of a 3D field over multiple images.

    The function can be used to add images to an already partially
    plotted animation.

    Parameters
    ----------
    out_file: pathlib.Path or str
        path of the gif output file

    data: xarray.DataArray
        3D data to plot, the 1st

    anim_dim: str
        the dimension over which to generate the animation
    """
    if out_file.is_file():
        out_file.unlink()
    with imageio.get_writer(out_file, mode='I') as writer:
        anim_field3d_append(data, out_writer=writer, **opt_args)
