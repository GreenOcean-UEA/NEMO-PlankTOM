#!/usr/bin/env python3

import numpy as np
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import sys
import pathlib

# Script for creating vertical depth plots for given variables
# Will only be run at the end of the model run process

tmod    = sys.argv[1]
baseDir = sys.argv[2]
year    = sys.argv[3]

saveDir = f'/{baseDir}/monitor/{tmod}/'
pathlib.Path(saveDir).mkdir(parents=True, exist_ok=True)

## Define functions
def verticalDepth(infile, var_list, var_calc, mask, metres, m2_glob):

	var_out = []

	for i in range(0,len(var_list)):
		
		varName = var_list[i]
	
		# Read in the variable
		var_in     = infile.variables[varName][:]
		var_glob   = np.where(mask == 1, var_in, np.nan)
		var_glob_m = var_glob * metres

		# Calculate average over all depths
		var_glob_a1 = np.apply_over_axes(np.nansum, var_glob_m, [2,3]) # sum over x and y
		var_glob_a2 = np.apply_over_axes(np.nanmean, var_glob_a1, [0]) # average over t

		var_glob_tot = np.zeros( (31) );

		# Use a loop to avoid dividing by zero
		for k in range(0,31):
			if m2_glob[k] == 0:
				var_glob_tot[k] = np.nan # this should only occur at the very lowest depth
			else:
				if var_calc[i] == "ave":
					var_glob_tot[k] = var_glob_a2[0,k,0,0] / m2_glob[k] # ave : average over the ocean
				else:
					var_glob_tot[k] = var_glob_a2[0,k,0,0] # sum : sum over the ocean
	
		var_out.append(var_glob_tot)

	return var_out


## Import mask files and creates metres variable
meshmask  = Dataset(f"/gpfs/data/greenocean/software/resources/breakdown/mesh_mask3_6_low_res.nc", "r")
basinmask = Dataset(f"/gpfs/data/greenocean/software/resources/breakdown/basin_mask.nc", "r")

metres = meshmask.variables['e1t'][:] * meshmask.variables['e2t'][:]
vol    = basinmask.variables['VOLUME'][:]

vol_mask    = np.where(vol > 0, 1, np.nan)
m2_glob     = np.where(vol_mask == 1, metres, np.nan)
m2_glob_tot = np.apply_over_axes(np.nansum, m2_glob, [1,2])


## Import model data
in_file_diad = Dataset(f"{baseDir}/{tmod}/ORCA2_1m_{year}0101_{year}1231_diad_T.nc", "r")
in_file_grid = Dataset(f"{baseDir}/{tmod}/ORCA2_1m_{year}0101_{year}1231_grid_T.nc", "r")
in_file_ptrc = Dataset(f"{baseDir}/{tmod}/ORCA2_1m_{year}0101_{year}1231_ptrc_T.nc", "r")

depthsLR = np.array((5, 15.00029, 25.00176, 35.00541, 45.01332, 55.0295, 65.06181, 75.12551, 85.25037, 95.49429, 105.9699, 116.8962, 128.6979, 142.1953, 158.9606, 181.9628,
                     216.6479, 272.4767, 364.303, 511.5348, 732.2009, 1033.217, 1405.698, 1830.885, 2289.768, 2768.242, 3257.479, 3752.442, 4250.401, 4749.913, 5250.227))

cols = ['b', 'r', 'g', 'm', 'y', 'c', 'b', 'r', 'g', 'm', 'y', 'c']


## Create ecosystem plot
var_list_eco = ["TChl", "PPT", "EXP", "ExpCO3", "sinksil"]
var_unit_eco = ["ug Chl/L", "umol/m3/s", "umol/m2/s", "umol/m2/s", "umol/m2/s"]
var_conv_eco = [1e6, 1e6, 1e6, 1e6, 1e6]           # units conversion multiplier
var_calc_eco = ["ave", "ave", "ave", "ave", "ave"] # whether to sum or average over time

var_out_eco = verticalDepth(in_file_diad, var_list_eco, var_calc_eco, vol_mask, metres, m2_glob_tot)

# Create plot figure
ratio = 2
fig = plt.figure()
fig.set_figheight(4*ratio)
fig.set_figwidth(9*ratio)

ax1 = plt.subplot2grid(shape=(4,9), loc=(0,0), rowspan=2, colspan=3)
ax2 = plt.subplot2grid(shape=(4,9), loc=(0,3), rowspan=2, colspan=3)
ax3 = plt.subplot2grid(shape=(4,9), loc=(2,0), rowspan=2, colspan=3)
ax4 = plt.subplot2grid(shape=(4,9), loc=(2,3), rowspan=2, colspan=3)
ax5 = plt.subplot2grid(shape=(4,9), loc=(2,6), rowspan=2, colspan=3)

axs = [ax1, ax2, ax3, ax4, ax5]

axs[0].plot(var_out_eco[0][0:20] * var_conv_eco[0], depthsLR[0:20], cols[0]); axs[0].set_title(f'TChl vertical depth profile [ug Chl/L]');
axs[1].plot(var_out_eco[1][0:20] * var_conv_eco[1], depthsLR[0:20], cols[1]); axs[1].set_title(f'PPT vertical depth profile [umol/m3/s]');
axs[2].plot(var_out_eco[2][1:30] * var_conv_eco[2], depthsLR[1:30], cols[2]); axs[2].set_title(f'EXP vertical depth profile [umol/m2/s]');
axs[3].plot(var_out_eco[3][1:30] * var_conv_eco[3], depthsLR[1:30], cols[3]); axs[3].set_title(f'CaCO3 flux vertical depth profile [umol/m2/s]');
axs[4].plot(var_out_eco[4][1:30] * var_conv_eco[4], depthsLR[1:30], cols[4]); axs[4].set_title(f'Si flux vertical depth profile [umol/m2/s]');

for ax in axs:
	ax.invert_yaxis()
	ax.grid(linestyle='--', linewidth=0.25)

axs[0].set_ylabel('depth [m]')
axs[2].set_ylabel('depth [m]')

plt.tight_layout()
fig.savefig(f'{saveDir}/{tmod}_vertical_depth_ecosystem.png')
print(f'Created {tmod} vertical depth figure for ecosystem')


## Create PFTs plot
var_list_pfts = ["PIC", "FIX", "COC", "MIX", "DIA", "PHA", "BAC", "PRO", "PTE", "MES", "GEL", "CRU"]
var_unit_pfts = ["umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L", "umol/L"]
var_conv_pfts = [1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6]                         # units conversion multiplier
var_calc_pfts = ["ave", "ave", "ave", "ave", "ave", "ave", "ave", "ave", "ave", "ave", "ave", "ave"] # whether to sum or average over time

var_out_pfts = verticalDepth(in_file_ptrc, var_list_pfts, var_calc_pfts, vol_mask, metres, m2_glob_tot)

# Create plot figure
ratio = 2
fig = plt.figure()
fig.set_figheight(8*ratio)
fig.set_figwidth(12*ratio)

ax1  = plt.subplot2grid(shape=(8,12), loc=(0,0), rowspan=2, colspan=3)
ax2  = plt.subplot2grid(shape=(8,12), loc=(0,3), rowspan=2, colspan=3)
ax3  = plt.subplot2grid(shape=(8,12), loc=(0,6), rowspan=2, colspan=3)
ax4  = plt.subplot2grid(shape=(8,12), loc=(0,9), rowspan=2, colspan=3)
ax5  = plt.subplot2grid(shape=(8,12), loc=(2,0), rowspan=2, colspan=3)
ax6  = plt.subplot2grid(shape=(8,12), loc=(2,3), rowspan=2, colspan=3)

ax7  = plt.subplot2grid(shape=(8,12), loc=(4,0), rowspan=2, colspan=3)
ax8  = plt.subplot2grid(shape=(8,12), loc=(4,3), rowspan=2, colspan=3)
ax9  = plt.subplot2grid(shape=(8,12), loc=(4,6), rowspan=2, colspan=3)
ax10 = plt.subplot2grid(shape=(8,12), loc=(4,9), rowspan=2, colspan=3)
ax11 = plt.subplot2grid(shape=(8,12), loc=(6,0), rowspan=2, colspan=3)
ax12 = plt.subplot2grid(shape=(8,12), loc=(6,3), rowspan=2, colspan=3)

axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]

for i in range(0,len(var_out_pfts)):
	v_l = var_list_pfts[i]
	v_u = var_unit_pfts[i]	

	axs[i].plot(var_out_pfts[i][0:20] * var_conv_pfts[i], depthsLR[0:20], cols[i]); axs[i].set_title(f'{v_l} vertical depth profile [{v_u}]');

for ax in axs:
	ax.invert_yaxis()
	ax.grid(linestyle='--', linewidth=0.25)

axs[0].set_ylabel('depth [m]')
axs[4].set_ylabel('depth [m]')
axs[6].set_ylabel('depth [m]')
axs[10].set_ylabel('depth [m]')

plt.tight_layout()
fig.savefig(f'{saveDir}/{tmod}_vertical_depth_pfts.png')
print(f'Created {tmod} vertical depth figure for pfts')


## Create nutrient plot
var_list_nutr = ["PO4", "NO3", "Fer", "Si", "O2", "DIC", "Alkalini"]
var_unit_nutr = ["umol/L", "umol/L", "nmol/L", "umol/L", "umol/L", "umol/L", "ueq/L"]
var_conv_nutr = [1e6 / 122, 1e6, 1e9, 1e6, 1e6, 1e6, 1e6]         # units conversion multiplier : for PO4, we need to divide by 122 for the Redfield ratio
var_calc_nutr = ["ave", "ave", "ave", "ave", "ave", "ave", "ave"] # whether to sum or average over time

var_out_nutr = verticalDepth(in_file_ptrc, var_list_nutr, var_calc_nutr, vol_mask, metres, m2_glob_tot)

# Observation data for nutrients
obs_po4 = np.array((0.52987817, 0.53723771, 0.55302459, 0.57752131, 0.60796679, 0.64315427, 0.68543891, 0.72999259, 0.77588365, 0.82197094, 0.87047604, 0.91896358, 0.97032602, 1.02294042, 1.0885159,  1.17099291, 1.28269611, 1.43557452, 1.65979373, 1.95546142, 2.26561946, 2.43762347, 2.43400088, 2.35771974, 2.2984275, 2.27786267, 2.26746146, 2.25801942, 2.25585054, 2.25091653, np.nan))

obs_no3 = np.array((5.15596731, 5.29250436, 5.49719698, 5.83949275, 6.27647375, 6.81155932, 7.47264204, 8.16599841, 8.86565344, 9.58907617, 10.35772892, 11.17176228, 12.02842699, 12.89016145, 13.95589097, 15.2886822, 17.1511473, 19.48546483, 23.0522826, 27.56186998, 32.22938658, 34.66204542, 34.68560257, 33.77784652, 33.0720577, 32.86223635, 32.77188099, 32.69169026, 32.65494948, 32.66317508, np.nan))

obs_fer = np.array((0.20410695, 0.20392409, 0.20346693, 0.20286527, 0.20364223, 0.20517888, 0.20117731, 0.20842939, 0.22054434, 0.23136719, 0.24208713, 0.27503688, 0.28285428, 0.31288732, 0.3288177, 0.36131601, 0.39949994, 0.43956002, 0.48749718, 0.55256646, 0.63194108, 0.68715239, 0.71547418, 0.72451647, 0.74744286, 0.7667307, 0.70055715, 0.64873491, 0.6147887, 0.60248194, np.nan))

obs_si  = np.array((7.47860313, 7.54843832, 7.67942662, 7.91786445, 8.27133267, 8.70980563, 9.27420487, 9.83478652, 10.45451495, 11.08533897, 11.81268545, 12.63385613, 13.5323516, 14.5139029, 15.71062878, 17.23699873, 19.50530787, 22.25874185, 27.10183926, 35.72050873, 50.64933101, 70.29555107, 88.06209122, 100.54722904, 108.07068799, 113.42978688, 117.28819568, 120.03367631, 121.75350525, 122.50957118, np.nan))

obs_o2  = np.array((250.97528313, 250.72531186, 249.9767884, 248.60442148, 246.6196827, 243.7979899, 239.89543791, 235.72111309, 231.17703255, 226.34979905, 221.37365453, 216.38875869, 211.1814609, 205.99116117, 199.94491128, 193.06273313, 184.94413176, 176.13303087, 166.42425088, 158.40550163, 148.72185209, 143.69675657, 150.12807541, 163.77241594, 175.33226097, 183.1712375, 189.71122936, 196.45112296, 202.29371039, 204.28782027, np.nan))

obs_alk = np.array((2352.92, 2355.41, 2359.23, 2362.62, 2365.28, 2367.67, 2369.90, 2372.04, 2373.23, 2374.55, 2375.51, 2376.17, 2376.68, 2376.78, 2376.57, 2375.86, 2374.47, 2372.36, 2370.82, 2372.23, 2382.88, 2403.02, 2423.12, 2436.43, 2444.32, 2449.16, 2451.55, 2451.82, 2450.59, 2450.58, np.nan))

# Create plot figure
ratio = 2
fig = plt.figure()
fig.set_figheight(6*ratio)
fig.set_figwidth(9*ratio)

ax1 = plt.subplot2grid(shape=(6,9), loc=(0,0), rowspan=2, colspan=3)
ax2 = plt.subplot2grid(shape=(6,9), loc=(0,3), rowspan=2, colspan=3)
ax3 = plt.subplot2grid(shape=(6,9), loc=(0,6), rowspan=2, colspan=3)
ax4 = plt.subplot2grid(shape=(6,9), loc=(2,0), rowspan=2, colspan=3)
ax5 = plt.subplot2grid(shape=(6,9), loc=(2,3), rowspan=2, colspan=3)
ax6 = plt.subplot2grid(shape=(6,9), loc=(4,0), rowspan=2, colspan=3)
ax7 = plt.subplot2grid(shape=(6,9), loc=(4,3), rowspan=2, colspan=3)

axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

for i in range(0,len(var_out_nutr)):
	v_l = var_list_nutr[i]
	v_u = var_unit_nutr[i]	

	axs[i].plot(var_out_nutr[i][0:30] * var_conv_nutr[i], depthsLR[0:30], cols[i]); axs[i].set_title(f'{v_l} vertical depth profile [{v_u}]');

axs[0].plot(obs_po4[0:30], depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);
axs[1].plot(obs_no3[0:30], depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);
axs[2].plot(obs_fer[0:30], depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);
axs[3].plot(obs_si[0:30],  depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);
axs[4].plot(obs_o2[0:30],  depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);
axs[6].plot(obs_alk[0:30], depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);

for ax in axs:
	ax.invert_yaxis()
	ax.grid(linestyle='--', linewidth=0.25)

axs[0].set_ylabel('depth [m]')
axs[3].set_ylabel('depth [m]')
axs[5].set_ylabel('depth [m]')

plt.tight_layout()
fig.savefig(f'{saveDir}/{tmod}_vertical_depth_nutrients.png')
print(f'Created {tmod} vertical depth figure for nutrients')


## Create physics plot
var_list_phys = ["votemper", "vosaline"]
var_unit_phys = ["degC", "1e-3"]
var_conv_phys = [1, 1]         # units conversion multiplier
var_calc_phys = ["ave", "ave"] # whether to sum or average over time

var_out_phys = verticalDepth(in_file_grid, var_list_phys, var_calc_phys, vol_mask, metres, m2_glob_tot)

# Observartion data for physics
obs_vot = np.array((18.46968323, 18.35156491, 18.16075393, 17.91191125, 17.55027907, 17.14119307, 16.70302896, 16.26867056, 15.83042865, 15.4261324, 14.99602221, 14.56798708, 14.12363121, 13.6793829, 13.13801873, 12.47545093, 11.60176232, 10.50469098, 9.05916524, 7.28409231, 5.55997472, 4.07699053, 3.07618598, 2.38382614, 1.91813595, 1.59916395, 1.36290165, 1.14330104, 0.96715039, 0.90128792, np.nan))

obs_vos = np.array((34.64622841, 34.69757408, 34.74632434, 34.81183937, 34.86226256, 34.90520628, 34.93940252, 34.96984264, 34.99678594, 35.02844861, 35.04024812, 35.04804427, 35.05274072, 35.04937885, 35.03927575, 35.01633672, 34.96953203, 34.89715145, 34.79069656, 34.66331274, 34.59205452, 34.61075754, 34.68039212, 34.72562831, 34.74371412, 34.74098658, 34.73739271, 34.73210276, 34.72774482, 34.72516131, np.nan))

# Create plot figure
ratio = 2
fig = plt.figure()
fig.set_figheight(2*ratio)
fig.set_figwidth(6*ratio)

ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), rowspan=2, colspan=3)
ax2 = plt.subplot2grid(shape=(2,6), loc=(0,3), rowspan=2, colspan=3)

axs = [ax1, ax2]

for i in range(0,len(var_out_phys)):
	v_l = var_list_phys[i]
	v_u = var_unit_phys[i]

	axs[i].plot(var_out_phys[i][0:30] * var_conv_phys[i], depthsLR[0:30], cols[i]); axs[i].set_title(f'{v_l} vertical depth profile [{v_u}]');

axs[0].plot(obs_vot[0:30], depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);
axs[1].plot(obs_vos[0:30], depthsLR[0:30], color='k', linestyle='dashed', alpha=0.75);

for ax in axs:
	ax.invert_yaxis()
	ax.grid(linestyle='--', linewidth=0.25)

axs[0].set_ylabel('depth [m]')
axs[1].set_xlim(34.4,35.2) # adjust the x scale for the salinity plot

plt.tight_layout()
fig.savefig(f'{saveDir}/{tmod}_vertical_depth_physics.png')
print(f'Created {tmod} vertical depth figure for physics')


## Create OC plot
var_list_oc = ["DOC", "POC", "GOC"]
var_unit_oc = ["umol/L", "umol/L", "umol/L"]
var_conv_oc = [1e6, 1e6, 1e6]       # units conversion multiplier
var_calc_oc = ["ave", "ave", "ave"] # whether to sum or average over time

var_out_oc = verticalDepth(in_file_ptrc, var_list_oc, var_calc_oc, vol_mask, metres, m2_glob_tot)

# Create plot figure
ratio = 2
fig = plt.figure()
fig.set_figheight(2*ratio)
fig.set_figwidth(9*ratio)

ax1 = plt.subplot2grid(shape=(2,9), loc=(0,0), rowspan=2, colspan=3)
ax2 = plt.subplot2grid(shape=(2,9), loc=(0,3), rowspan=2, colspan=3)
ax3 = plt.subplot2grid(shape=(2,9), loc=(0,6), rowspan=2, colspan=3)

axs = [ax1, ax2, ax3]

for i in range(0,len(var_out_oc)):
	v_l = var_list_oc[i]
	v_u = var_unit_oc[i]

	axs[i].plot(var_out_oc[i][0:30] * var_conv_oc[i], depthsLR[0:30], cols[i]); axs[i].set_title(f'{v_l} vertical depth profile [{v_u}]');

for ax in axs:
	ax.invert_yaxis()
	ax.grid(linestyle='--', linewidth=0.25)

axs[0].set_ylabel('depth [m]')

plt.tight_layout()
fig.savefig(f'{saveDir}/{tmod}_vertical_depth_orgcarbon.png')
print(f'Created {tmod} vertical depth figure for organic carbon')

