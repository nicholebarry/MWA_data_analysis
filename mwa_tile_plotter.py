### 
#
# Plot the tile locations of the MWA given a metafits file
#
###

import matplotlib
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import argparse
import glob


#********************************
def main():

	# Parse the command line inputs. 
	parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=
		"Plot the tile locations of the MWA given an input metafits file \n" +
		"Example call: \n"
		"python mwa_tile_plotter.py --metafits=/path/to/metafits/<gps_time>.metafits --out_dir=/path/to/save/directory")
	parser.add_argument("-m", "--metafits", default=None, type=str,
		help="Absolute path to the metafits file")
	parser.add_argument("-o", "--out_dir", default='./', type=str,
		help="Path to the plot output directory")
	args = parser.parse_args()

	# Check for a metafits file in the current dir if a path is not provided
	if not args.metafits:
		try:
			# Default to the first one found
			metafits = glob.glob('*.metafits')
			metafits = metafits[0]
		except:
			exit("Provide a path to a metafits file")
	else:
		metafits = args.metafits

	# Parse the tile locations from the metafits
	tile_names, north, east = get_tile_locations(metafits=metafits)

	# Plot the tile locations using matplotlib
	out = plot_tile_locations(args.out_dir, tile_names, north, east)

#********************************
def get_tile_locations(metafits):

	tile_names = []
	north = []
	east = []
	hdu = fits.open(metafits)

	# Get the tile names from the metafits, which in older versions
	# are under the key "Tile" and newer versions are under the key 
	# "TileName"
	try:
		tile_names_all_pols = hdu[1].data["TileName"]
	except:
		tile_names_all_pols = hdu[1].data["Tile"]

	north_all_pols = hdu[1].data["North"]
	east_all_pols = hdu[1].data["East"]
	pols = hdu[1].data["Pol"]

	# Get locations from only the X polarization to prevent duplication
	# X locations are the same as Y locations
	for tile_i in range(len(tile_names_all_pols)):
		if "X" in pols[tile_i]:
			tile_names.append(tile_names_all_pols[tile_i])
			north.append(north_all_pols[tile_i])
			east.append(east_all_pols[tile_i])

	# Return the name of the tile, its northing, and its easting
	return tile_names, north, east


#********************************
def plot_tile_locations(out_dir, tile_names, north, east):

	# Calculate the mean North and East value for plot aethestics 
	east = east - np.mean(east)
	north = north - np.mean(north)
	
	fig, ax = plt.subplots(figsize=(3.6, 3.6))
	ax.scatter(east, north, marker="s", s=10.0, color="orange", label="MWA")

	# Calculate the maximum location in m to find the plot limits
	max_east = np.max(np.abs(east))
	max_north = np.max(np.abs(north))
	max_val = np.max([max_east,max_north])

	ax.set_aspect("equal")
	ax.set_ylim(-max_val, max_val)
	ax.set_xlim(-max_val, max_val)
	ax.set_ylabel("North [m]")
	ax.set_xlabel("East [m]")

	plt.savefig(f"{out_dir}/mwa_map.pdf", transparent=True, bbox_inches="tight")
	print(f"MWA tile location plot saved to {out_dir}")

#********************************

#********************************

if __name__ == '__main__':
        main()

