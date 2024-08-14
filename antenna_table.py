from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

NUM_FREQS = 16
NUM_TIMES = 56
NUM_BASELINES = 8128

def get_data(file_prepend, sim, num_bands):

    all_xx = np.empty((NUM_BASELINES*NUM_TIMES, NUM_FREQS*num_bands), dtype=complex)
    all_uu = np.empty((NUM_BASELINES*NUM_TIMES, num_bands))
    all_vv = np.empty((NUM_BASELINES*NUM_TIMES, num_bands))

    for band in range(1, num_bands + 1):
        with fits.open(f"{file_prepend}{band:02d}.uvfits") as hdu:

            if sim == 'FHD':
                this_xx = hdu[0].data.data[:,0,0,0,:,0,0] + 1j*hdu[0].data.data[:,0,0,0,:,0,1]

            elif sim == 'WODEN':
                this_xx = hdu[0].data.data[:,0,0,:,0,0] + 1j*hdu[0].data.data[:,0,0,:,0,1]

            # print(this_xx.shape, band*NUM_FREQS, (band+1)*NUM_FREQS)

            all_xx[:,(band-1)*NUM_FREQS:band*NUM_FREQS] = this_xx
            all_uu[:,(band-1)] = hdu[0].data['UU']
            all_vv[:,(band-1)] = hdu[0].data['VV']

    return all_xx, all_uu, all_vv

def plot_uv_coords(wod_uu, wod_vv, fhd_uu, fhd_vv, num_bands):

    for band in range(num_bands):
        fig, axs = plt.subplots(2,2,figsize=(12,12))

        wod_uu_plus_conj = np.append(wod_uu, -wod_uu)
        wod_vv_plus_conj = np.append(wod_vv, -wod_vv)

        fhd_uu_plus_conj = np.append(fhd_uu, -fhd_uu)
        fhd_vv_plus_conj = np.append(fhd_vv, -fhd_vv)

        print("THIS_MANY", len(np.where(fhd_uu == 0)[0]))
        print("THIS_MANY", len(np.where(fhd_vv == 0)[0]))


        print(np.where(fhd_vv == 0)[0])

        # axs[0,0].plot(wod_uu[:, band], wod_vv[:, band],marker='.',linestyle='none',
        #             label='WODEN',ms=0.5, color='C0')
        # axs[0,0].plot(-wod_uu[:, band], -wod_vv[:, band],marker='.',linestyle='none',
        #             ms=0.5, color='C0')

        axs[0,0].plot(wod_uu_plus_conj, wod_vv_plus_conj, marker='.',linestyle='none',
                    label='WODEN',ms=0.5, color='C0')

        axs[0,1].plot(fhd_uu_plus_conj, fhd_vv_plus_conj,marker='.',linestyle='none',
                    label='FHD',ms=0.5, color='C0')

        bins = np.linspace(-1e-6, 1e-6, 76)



        axs[1,0].hist(wod_uu_plus_conj, histtype='step',label='WODEN UU',
                 bins=bins, lw=1.5)
        axs[1,0].hist(wod_vv_plus_conj, histtype='step',label='WODEN VV',
                 bins=bins, lw=1.5)

        axs[1,1].hist(fhd_uu_plus_conj, histtype='step',label='FHD UU',
                 bins=bins, lw=1.5)
        axs[1,1].hist(fhd_vv_plus_conj, histtype='step',label='FHD VV',
                 bins=bins, lw=1.5)

        axs[1,0].axvline(0.0,color='k',linestyle='--', alpha=0.4)
        axs[1,1].axvline(0.0,color='k',linestyle='--', alpha=0.4)

        for ax in axs[0,:]:
            ax.set_xlabel('$u$ (seconds)')
            ax.set_ylabel('$v$ (seconds)')

        for ax in axs.flatten():
            ax.legend()


        plt.tight_layout()
        fig.savefig(f'uvw_comp_band{band+1:02d}.png',bbox_inches='tight')

        # plt.show()
        plt.close()


if __name__ == '__main__':
    num_bands = 2
    wod_xx, wod_uu, wod_vv = get_data(f"./data/MWA_analy_pumav3_band", 'WODEN', num_bands)
    fhd_xx, fhd_uu, fhd_vv = get_data(f"./data/uv_model_", 'FHD', num_bands)

    plot_uv_coords(wod_uu, wod_vv, fhd_uu, fhd_vv, num_bands)

