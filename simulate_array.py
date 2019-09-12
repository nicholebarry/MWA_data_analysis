#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import csv
from pyuvdata import UVData
from pyuvdata import uvutils


def create_hex_array(side_length, antenna_spacing, output_path,
                     save_uvfits=True, plot_array=True):

    # get the number of antennas in an array with this side length
    antennas = ((side_length-1)**2 + (side_length-1))*3 + 1

    # initialize antenna locations with the center antenna
    antenna_xlocs = []
    antenna_ylocs = []

    for row in range(0, side_length):
        yloc = row*antenna_spacing*(3.**.5)/2
        if row % 2 == 0:
            new_xlocs = [
                i*antenna_spacing for i in range(
                    -side_length+row/2+1,
                    side_length-row/2
                )
            ]
            antenna_xlocs.extend(new_xlocs)
            antenna_ylocs.extend([yloc]*len(new_xlocs))
            if row != 0:
                antenna_xlocs.extend(new_xlocs)
                antenna_ylocs.extend([-yloc]*len(new_xlocs))
        else:
            new_xlocs = [
                         (i+.5)*antenna_spacing for i in range(
                            int(-side_length+row/2.+.5),
                            int(side_length-row/2.-.5)
                          )
                         ]
            antenna_xlocs.extend(new_xlocs*2)
            antenna_ylocs.extend([yloc]*len(new_xlocs))
            antenna_ylocs.extend([-yloc]*len(new_xlocs))
    antennas = len(antenna_xlocs)

    print antennas

    if plot_array:
        plt.figure()
        plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1.)
        # plt.xlim(ra_plot_range[0],ra_plot_range[1])
        # plt.ylim(dec_plot_range[0],dec_plot_range[1])
        plt.xlabel('East/West Location (m)')
        plt.ylabel('North/South Location (m)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(
            '{}/antenna_hex_{}.png'.format(output_path, int(antennas))
        )
        plt.close()

    # find the array's radial distribution
    radial_dist = [
        (antenna_xlocs[i]**2+antenna_ylocs[i]**2)**.5 for i in range(
            len(antenna_xlocs)
        )
    ]
    radial_hist, bin_edges = np.histogram(
        radial_dist, bins=int(side_length*.75)
    )
    bin_centers = [
        (bin_edges[i]+bin_edges[i+1])/2 for i in range(len(bin_edges)-1)
    ]
    plt.figure()
    plt.scatter(bin_centers, radial_hist)
    plt.savefig(
        '{}/antenna_hex_dist_{}.png'.format(output_path, antennas)
    )
    plt.close()

    print antennas

    if save_uvfits:
        create_uvfits(antenna_xlocs,
                      antenna_ylocs,
                      '{}/hex_array_sim_{}.uvfits'.format(output_path, int(antennas))
                      )
        # Save antenna locations to a csv
        csv_outfile = open(
            '{}/hex_array_sim_{}_antenna_locs.csv'.format(output_path, int(antennas)),
            'w'
        )
        outfile_writer = csv.writer(csv_outfile)
        outfile_writer.writerow(['Antenna Number', 'E-W Location (m)',
                                 'N-S Location (m)', 'Altitude (m)'])
        for i in range(antennas):
            outfile_writer.writerow(
                [i, antenna_xlocs[i], antenna_ylocs[i], 0.]
            )
        csv_outfile.close()

    return antennas, radial_hist, bin_centers


def create_random_array(side_length, antenna_spacing, output_path,
                        plot_array=True, save_uvfits=True):

    array_numbers = [1]
    number_of_arrays = len(array_numbers)

    antenna_size = antenna_spacing/3.  # Minimum antenna spacing
    antennas, radial_hist, bin_centers = create_hex_array(side_length,
                                                          antenna_spacing,
        						  output_path,
                                                          save_uvfits=False)
    print antennas
    radial_vals = np.arange(bin_centers[0], bin_centers[-1], .1)
    radial_distribution = np.interp(radial_vals, bin_centers, radial_hist)
    norm_factor = sum(radial_distribution)
    radial_distribution = [i/norm_factor for i in radial_distribution]

    for array in array_numbers:

        antenna_xlocs = []
        antenna_ylocs = []
        for ant in range(antennas):

            random = np.random.uniform()
            dist_sum = 0.
            for i, dist_val in enumerate(radial_distribution):
                dist_sum += dist_val
                if random < dist_sum:
                    break
            radius = radial_vals[i]

            check_overlap = False
            while not check_overlap:
                azimuth = np.random.uniform()*2*math.pi
                xloc = radius*math.cos(azimuth)
                yloc = radius*math.sin(azimuth)

                # Check that the antenna does not overlap others
                check_overlap = True
                for i in range(len(antenna_xlocs)):
                    if (
                        abs(antenna_xlocs[i]-xloc) < antenna_size
                        and abs(antenna_ylocs[i]-yloc) < antenna_size
                    ):
                        check_overlap = False
                if check_overlap:
                    antenna_xlocs.append(xloc)
                    antenna_ylocs.append(yloc)

        if plot_array:
            plt.figure()
            plt.scatter(antenna_xlocs, antenna_ylocs, marker='o', s=1.)
            # plt.xlim(ra_plot_range[0],ra_plot_range[1])
            # plt.ylim(dec_plot_range[0],dec_plot_range[1])
            plt.xlabel('East/West Location (m)')
            plt.ylabel('North/South Location (m)')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.savefig('{}/antenna_random{}_{}.png'.format(
                            output_path,
                            array,
                            int(antennas))
                        )
            plt.close()

        if save_uvfits:
            print 'creating uvfits'
            create_uvfits(antenna_xlocs, antenna_ylocs,
                          '{}/random{}_array_sim_{}.uvfits'.format(
                              output_path,
                              array,
                              int(antennas))
                          )
            # Save antenna locations to a csv
            csv_outfile = open(
                '{}/random{}_array_sim_{}_antenna_locs.csv'.format(
                    output_path,
                    array,
                    int(antennas)
                ),
                'w'
            )
            outfile_writer = csv.writer(csv_outfile)
            outfile_writer.writerow(['Antenna Number', 'E-W Location (m)',
                                     'N-S Location (m)', 'Altitude (m)'])
            for i in range(antennas):
                outfile_writer.writerow(
                    [i, antenna_xlocs[i], antenna_ylocs[i], 0.]
                )
            csv_outfile.close()


def create_hera_array(side_length, antenna_spacing, output_path,
                      plot_array=True, save_uvfits=True):

    antenna_spacing = float(antenna_spacing)

    a1 = np.array([antenna_spacing/2., antenna_spacing*np.sqrt(3)/2.])
    a2 = np.array([antenna_spacing/2., -antenna_spacing*np.sqrt(3)/2.])
    a3 = -a1-a2
    d0 = np.array([antenna_spacing*2./3., 0])
    d1 = np.array([-antenna_spacing/3., antenna_spacing/3.*np.sqrt(3)])
    d2 = np.array([-antenna_spacing/3., -antenna_spacing/3.*np.sqrt(3)])
    pos = []

    n = 0
    for ii in range(side_length-1):
        for jj in range(side_length-1):
            pos.append(d0+ii*a1+jj*a2)
            #pos[n, :] = d0+ii*a1+jj*a2
            n += 1
    for ii in range(side_length):
        for jj in range(side_length-1):
            pos.append(d1+ii*a3+jj*a1)
            #pos[n, :] = d1+ii*a3+jj*a1
            n += 1
    for ii in range(side_length):
        for jj in range(side_length):
            pos.append(d2+ii*a2+jj*a3)
            #pos[n, :] = d2+ii*a2+jj*a3
            n += 1

    antennas = len(pos)
    pos = np.array(pos)

    if plot_array:
        plt.figure()
        plt.scatter(pos[:,0], pos[:,1], marker='o', s=1.)
        plt.xlabel('East/West Location (m)')
        plt.ylabel('North/South Location (m)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(
            '{}/antenna_split_hex_{}.png'.format(output_path, int(antennas))
        )
        plt.close()

    if save_uvfits:
        print 'creating uvfits'
        create_uvfits(pos[:,0], pos[:,1],
                      '{}/split_hex_array_sim_{}.uvfits'.format(
                          output_path,
                          int(antennas))
                      )
        # Save antenna locations to a csv
        csv_outfile = open(
            '{}/split_hex_array_sim_{}_antenna_locs.csv'.format(output_path, int(antennas)),
            'w'
        )
        outfile_writer = csv.writer(csv_outfile)
        outfile_writer.writerow(['Antenna Number', 'E-W Location (m)',
                                 'N-S Location (m)', 'Altitude (m)'])
        for i in range(antennas):
            outfile_writer.writerow(
                [i, pos[i,0], pos[i,1], 0.]
            )
        csv_outfile.close()


def create_uvfits(antenna_xlocs, antenna_ylocs, save_filename):

    Nants = len(antenna_xlocs)
    antenna_locs_ENU = np.zeros((Nants, 3))
    antenna_locs_ENU[:, 0] = antenna_xlocs
    antenna_locs_ENU[:, 1] = antenna_ylocs

    filename = '/home/ubuntu/1061316296.uvfits'
    UV = UVData()
    UV.read_uvfits(filename)

    UV.Nants_data = Nants
    UV.Nants_telescope = Nants
    UV.Nbls = (Nants**2 + Nants)/2
    UV.Nblts = UV.Nbls*UV.Ntimes
    UV.antenna_numbers = np.array(range(Nants), dtype=int)
    UV.antenna_names = [str(ant) for ant in UV.antenna_numbers]
    ant_1_array = np.zeros(UV.Nbls, dtype=int)
    ant_2_array = np.zeros(UV.Nbls, dtype=int)
    baseline_array = np.zeros(UV.Nbls, dtype=int)
    index = 0
    for ant1 in range(Nants):
        for ant2 in range(ant1, Nants):
            ant_1_array[index] = ant1
            ant_2_array[index] = ant2
            baseline_array[index] =  2048 * (ant1+1) + (ant2+1) + 2**16
            index += 1
    UV.ant_1_array = np.tile(ant_1_array, UV.Ntimes)
    UV.ant_2_array = np.tile(ant_2_array, UV.Ntimes)
    UV.baseline_array = np.tile(baseline_array, UV.Ntimes)
    old_time_array = np.copy(UV.time_array)
    UV.time_array = np.repeat(np.array(list(set(UV.time_array))), UV.Nbls)
    UV.lst_array = np.array(
        [UV.lst_array[
                      np.where(old_time_array == time)[0][0]
                      ] for time in UV.time_array]
    )
    # Add dummy data
    UV.data_array = np.full((UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols),
                            UV.data_array[0, 0, 0, 0],
                            dtype=complex)
    UV.nsample_array = np.full((UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols),
                            1.,
                            dtype=float)
    # Unflag all
    UV.flag_array = np.full((UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols),
                            False,
                            dtype=bool)
    # Calculate UVWs
    UV.uvw_array = np.zeros((UV.Nblts, 3), dtype=float)
    antenna_locs_ECEF = uvutils.ECEF_from_ENU(
        antenna_locs_ENU.T, *UV.telescope_location_lat_lon_alt
        ).T
    UV.antenna_positions = antenna_locs_ECEF - UV.telescope_location
    UV.set_uvws_from_antenna_positions(
        allow_phasing=True, orig_phase_frame='gcrs', output_phase_frame='icrs'
    )
    UV.phase_center_frame = 'icrs'

    print 'Saving uvfits to {}'.format(save_filename)
    UV.write_uvfits(save_filename, spoof_nonessential=True)


if __name__ == '__main__':
    output_path = '/home/ubuntu'
    create_hera_array(11, 15., output_path)
