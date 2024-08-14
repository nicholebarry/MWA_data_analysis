import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pyuvdata import UVData, UVBeam


def make_polar_contour_plot(
    ax,
    plot_vals,
    az_vals,
    za_vals,
    vmin=-1,
    vmax=1,
    cyclic_colorbar=False,
    ncontours=11,
):

    if cyclic_colorbar:
        use_cmap = matplotlib.cm.get_cmap("twilight_shifted").copy()
    else:
        if vmin >= 0:
            use_cmap = matplotlib.cm.get_cmap("inferno").copy()
        else:
            use_cmap = matplotlib.cm.get_cmap("Spectral").copy()
    use_cmap.set_bad(color="whitesmoke")


    # Fill in plotting gap by copying az=0 values to az=2Pi
    az_zeros = np.where(az_vals == 0.0)
    az_vals = np.concatenate(
        (az_vals, (az_vals[az_zeros[0], az_zeros[1]] + 2 * np.pi)[np.newaxis, :]),
        axis=0,
    )
    za_vals = np.concatenate(
        (za_vals, (za_vals[az_zeros[0], az_zeros[1]])[np.newaxis, :]), axis=0
    )
    plot_vals = np.concatenate(
        (plot_vals, (plot_vals[az_zeros[0], az_zeros[1]])[np.newaxis, :]), axis=0
    )


    # Set contour levels
    levels = np.linspace(vmin, vmax, num=ncontours)


    contourplot = ax.contourf(
        az_vals,
        za_vals,
        plot_vals,
        levels,
        vmin=vmin,
        vmax=vmax,
        cmap=use_cmap,
    )
    contourplot.set_clim(vmin=vmin, vmax=vmax)
    return contourplot


def plot_beam(
    beam,  # pyuvdata beam object
    plot_freq=50.0e6,  # Frequency in Hz, must be included in the beam obj
    real_part=True,
    plot_amplitude=False,
    plot_pols=[0, 1],
    vmin=-1,
    vmax=1,
    horizon_cut=True,  # Truncate the beam model at the horizon
    az_za_basis=True,  # Label plots as Az/ZA
    savepath=None,
):

    #Naxes_vec, Nfeeds or Npols, Nfreqs, Npixels
    use_beam = beam.select(frequencies=[plot_freq], inplace=False)
    az_axis = np.degrees(beam.axis1_array)
    za_axis = np.degrees(beam.axis2_array)
    if horizon_cut:
        use_za_inds = np.where(za_axis < 90.)[0]
        use_beam.select(axis2_inds=use_za_inds, inplace=True)
        za_axis = za_axis[use_za_inds]


    if plot_amplitude:
        plot_jones_vals = np.sqrt(
            np.abs(use_beam.data_array[0, :, :, :, :]) ** 2.0
            + np.abs(use_beam.data_array[1, :, :, :, :]) ** 2.0
        )
        # Normalize
        plot_jones_vals /= np.max(plot_jones_vals)
        title = f"Beam Amplitude, {plot_freq/1e6} MHz"
    elif real_part:
        plot_jones_vals = np.real(use_beam.data_array)
        title = f"Jones Matrix Components at {plot_freq/1e6} MHz, Real Part"
    else:
        plot_jones_vals = np.imag(use_beam.data_array)
        title = f"Jones Matrix Components at {plot_freq/1e6} MHz, Imaginary Part"


    use_cmap = matplotlib.cm.get_cmap("Spectral").copy()
    use_cmap.set_bad(color="whitesmoke")
    za_vals, az_vals = np.meshgrid(za_axis, az_axis)

    #Naxes_vec, Nfeeds or Npols, Nfreqs, Npixels

    feed_names = ["X", "Y"]
    sky_pol_names = ["Az", "ZA"]
    if plot_amplitude:
        fig, ax = plt.subplots(
            nrows=1, ncols=2, subplot_kw=dict(projection="polar"), figsize=(9, 6)
        )
        
        for pol in plot_pols:
            contourplot = make_polar_contour_plot(
                ax[pol],
                (plot_jones_vals[pol, 0, :, :]).T,
                np.radians(az_vals),
                za_vals,
                vmin=vmin,
                vmax=vmax,
            )
            fig.colorbar(contourplot, ax=ax[pol])
            ax[pol].set_title(f"Pol {feed_names[pol]}")
    else:
        fig, ax = plt.subplots(
            nrows=2, ncols=2, subplot_kw=dict(projection="polar"), figsize=(9, 9)
        )
        for pol1 in plot_pols:
            for pol2 in [0, 1]:
                contourplot = make_polar_contour_plot(
                    ax[pol1, pol2],
                    (plot_jones_vals[pol1, pol2, 0, :, :]).T,
                    np.radians(az_vals),
                    za_vals,
                    vmin=vmin,
                    vmax=vmax,
                )
                fig.colorbar(contourplot, ax=ax[pol1, pol2])
                subtitle = f"Feed {feed_names[pol2]}"
                if az_za_basis:
                    subtitle += f", {sky_pol_names[pol1]}"
                ax[pol1, pol2].set_title(subtitle)
    fig.suptitle(title)
    fig.tight_layout()
    if savepath is None:
        plt.show()
        plt.close()
    else:
        plt.savefig(savepath, dpi=300)
        plt.close()


def main():
    delays = np.zeros((2, 16), dtype='int')
    # delays[0,:] = [0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6] #-2
    # delays[1,:] = [0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6]
    # delays[0,:] = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3] #-1
    # delays[1,:] = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
    delays[0,:] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # 0
    delays[1,:] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    #3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0
    #delays[:, 0] = 32
    #beam = UVBeam.from_file('/fred/oz048/MWA/CODE/FHD/mwa_full_embedded_element_pattern.h5', delays=delays, use_future_array_shapes=True)
    beam = UVBeam.from_file('/fred/oz048/MWA/CODE/FHD/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5', delays=delays, use_future_array_shapes=True)
    plot_beam(beam,plot_freq=1.8176e+08,savepath='/fred/oz048/MWA/CODE/FHD/test.png',plot_amplitude=True,vmin=0,vmax=1)

    
    #beam = UVBeam.from_file('/fred/oz048/MWA/CODE/FHD/mwa_full_embedded_element_pattern.h5', delays=delays, use_future_array_shapes=True)
    #plot_beam(beam,plot_freq=1.8176e+08,savepath='/fred/oz048/MWA/CODE/FHD/test.png',plot_amplitude=True,vmin=0,vmax=1)

#********************************

if __name__ == '__main__':
        main()