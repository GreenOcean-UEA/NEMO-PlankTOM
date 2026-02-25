#!/usr/bin/env python
import sys
import shutil
import os.path
import datetime
from dateutil import relativedelta
import numpy as np
from netCDF4 import Dataset
import math
import logging
import glob
from decimal import Decimal
from breakdown_functions import surfaceData, volumeData, intergrateData, volumeDataAverage, observationData, levelData
from breakdown_functions import bloom, trophic, regrid
from breakdown_observations import observationDatasets

# ---------- 1 UTILITIES ----------
# Get command line arguments
if len(sys.argv) != 4:
	sys.exit("Stopped - Incorrect arguements. Use: breakdown.py <parameter file> <year from> <year to>")

parmFile = sys.argv[1]
yearFrom = int(sys.argv[2])
yearTo = int(sys.argv[3])

# set up logging to file 
logging.basicConfig(level=logging.DEBUG,
			format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
			datefmt='%m-%d %H:%M',
			filename='breakdown.'+str(yearFrom)+'_'+str(yearTo)+'.log',
			filemode='w')

# define a Handler which writes INFO messages or higher
console = logging.StreamHandler()
console.setLevel(logging.INFO)

# set a format which is simpler for console use
formatter = logging.Formatter('%(levelname)-5s %(message)s')

# tell the handler to use this format
console.setFormatter(formatter)

# add the handler to the root logger
logging.getLogger('').addHandler(console)
log = logging.getLogger("Run")

# Get command line arguments
if len(sys.argv) != 4:
	log.error("Command line arguments not correct")
	log.error("Use: breakdown.py <parameter file> <year from> <year to>")
	sys.exit("Stopped")

log.info("Processing for parameters: "+str(parmFile))
log.info("Processing for years: "+str(yearFrom)+" to "+str(yearTo))


# ---------- 2 OBSERVATION DATA ----------
# get from exterior file
ListOfObservations = observationDatasets()

for l in ListOfObservations:
	log.info(l['name']+" available from "+str(l['origin']))


# ---------- 3 UNITS ----------
# constants
missingVal = 1E20
secondsInYear = 3600.*24.*365.  # seconds in a year (non leap)
peta = 1e-15
terra = 1e-12
giga = 1e-9
carbon = 12
litre = 1000
micro = 1e6
deptharea = 200.*3.62e+14 # depth in meters * surface ocean area

ListOfUnits = {
		"PgCarbonPerYr" : peta*carbon*secondsInYear,
		"TgCarbonPerYr" : terra*carbon*secondsInYear,
		"TmolPerYr" : terra*secondsInYear,
		"TgPerYr" : terra*secondsInYear,
		"PgPerYr" : peta*secondsInYear,
		"PerYr" : 1,
		"giga" : 1/giga,
		"1/giga" : 1/giga,
		"PgCarbon" : peta*carbon,
		"Conc->PgCarbon" : peta*carbon*litre, # for converting concentraion (mol/l) to total mass
		"u/L" : micro, # for converting gChl/L to ugChl/L OR molC/L to umolC/L: Becci
		"PgCarbonPer200m" : peta*carbon*litre*deptharea # for converting PFT mol/l to PgC for ocean top 200m: Becci
		}

log.info("UNITS:")
for key, value in ListOfUnits.items():
	log.info(str(str(key)+" "+str(value)))

# set area of each pixel given 1 degree grid
Er = 6.3781E6 # meters
Ec = 2*math.pi*Er
Ea = 4*math.pi*Er*Er

area = np.zeros((180,360))
for y in range(0,180):
	ang = np.radians(y-90)
	area[y,:] = Ec*(1/360.)*math.cos(ang) * Ec*(1/360.)


# ---------- 4 GET PARAMETERS ----------
# Define variable lists
varSurface = []
varLevel = []
varVolume = []
varInt = []
varTotalAve = []
obsComparisons = []
properties = []
varMap = []
nc_woa_id = -1

# Get parameters from file
f_parmFile = open(parmFile, "r")
lines = f_parmFile.readlines()

for l in lines:
	l = l.strip('\n')
	if len(l) > 0 and l[0] != '#':
		words = l.split(':')
		type = words[0]
		value = words[1]

		if type == 'BasinMask': basinFileName = value
		if type == 'WOAMask': WOAFileName = value
		if type == 'RegionMask': RegionFileName = value
		if type == 'RECCAPmask': RECCAPFileName = value
		if type == 'Meshmask': MeshMaskFileName = value
		if type == 'AncillaryData': AncilMaskFileName = value

		lonLim = [-90,90]
		latLim = [-180,180]
		reg = -1

		# the last -1 value is a place holder for the value calculated
		if type == 'Surface': 
			variable = value.split(',')[0]
			units = value.split(',')[1] 
			if value.split(',')[2] == 'Region':
				reg = int(value.split(',')[3])
			else:
				lonLim = value.split(',')[2]
				latLim = value.split(',')[3]
			key = value.split(',')[4]

			varSurface.append([variable, units, lonLim, latLim, key, reg, []])
			log.info("Surface: "+variable+" "+units)

		if type == 'Level': 
			variable = value.split(',')[0] 
			level = int(value.split(',')[1])
			units = value.split(',')[2] 
			if value.split(',')[3] == 'Region':
				reg = int(value.split(',')[4])
			else:
				lonLim = value.split(',')[3]
				latLim = value.split(',')[4]
			key = value.split(',')[5]

			varLevel.append([variable, level, units, lonLim, latLim, key, reg, []])
			log.info("Level: "+variable+" "+units)

		if type == 'Volume': 
			variable = value.split(',')[0] 
			units = value.split(',')[1] 
			if value.split(',')[2] == 'Region':
				reg = int(value.split(',')[3])
			else:
				lonLim = value.split(',')[2]
				latLim = value.split(',')[3]
			key = value.split(',')[4]

			varVolume.append([variable, units, lonLim, latLim, key, reg, []])
			log.info("Volume: "+variable+" "+units)

		if type == 'Integration':
			variable = value.split(',')[0] 
			startDepth= int(value.split(',')[1])
			endDepth = int(value.split(',')[2])
			units = value.split(',')[3] 
			if value.split(',')[4] == 'Region':
				reg = int(value.split(',')[5])
			else:
				lonLim = value.split(',')[4]
				latLim = value.split(',')[5]
			key = value.split(',')[6]

			varInt.append([variable,startDepth,endDepth,units,lonLim, latLim, key, reg, []])
			log.info("Integration: "+variable+" "+str(startDepth)+" "+str(endDepth)+" "+units)

		if type == 'Avg':
			variables = value.split(',')[0] 
			startDepth= int(value.split(',')[1])
			endDepth = int(value.split(',')[2])
			units = value.split(',')[3] 
			if value.split(',')[4] == 'Region':
				reg = int(value.split(',')[5])
			else:
				lonLim = value.split(',')[4]
				latLim = value.split(',')[5]
			key = value.split(',')[6]

			varTotalAve.append([variables,startDepth,endDepth,units,lonLim, latLim, key, reg, []])
			log.info("Avg: "+str(variables)+" "+units)

		if type == 'Observations':
			obsData = value.split(',')[0]
			obsVar = value.split(',')[1]       
			# handle multiple entries later
			var = value.split(',')[2]
			depthObs = int(value.split(',')[3])
			depthVar = int(value.split(',')[4])
			if value.split(',')[5] == 'T':
				gamFlag = True
			else:
				gamFlag = False

			if value.split(',')[6] == 'Region':
				reg = int(value.split(',')[7])
			else:
				lonLim = value.split(',')[6]
				latLim = value.split(',')[7]
			key = value.split(',')[8] 

			obsComparisons.append([obsData,obsVar,var, depthObs, depthVar, gamFlag, lonLim, latLim, key, reg, [] ])
			log.info("Observations: "+str(obsData)+" "+var) 

		if type == 'Property':
            
			propName = value.split(',')[0]
			if propName == "Bloom":
				propVar = value.split(',')[1]
				depthFrom = value.split(',')[2]
				depthTo = value.split(',')[3]
				lonLim = value.split(',')[4] 
				latLim = value.split(',')[5] 
				key = value.split(',')[6] 

				properties.append([propName, propVar, depthFrom, depthTo, lonLim, latLim, key, [] ])
				log.info("Property: "+str(propName)+" from depth level "+depthFrom+" to "+depthTo)

			if propName == "Trophic":
				propVar_1 = value.split(',')[1]
				propVar_2 = value.split(',')[2]
				propVar_3 = value.split(',')[3]

				depthFrom = value.split(',')[4]
				depthTo = value.split(',')[5]
				lonLim = value.split(',')[6] 
				latLim = value.split(',')[7] 
				key = value.split(',')[8] 

				propVar = [propVar_1, propVar_2, propVar_3]

				properties.append([propName, propVar, depthFrom, depthTo, lonLim, latLim, key, [] ])
				log.info("Property: "+str(propName)+" from depth level "+depthFrom+" to "+depthTo)

		if type == "Map":
			variable = value.split(',')[0] 
			level = value.split(',')[1]
			varMap.append([variable,level])
			log.info("Map: "+variable)
			if nc_woa_id == -1:
				# open WOA mask
				WOAFile = glob.glob(WOAFileName)
				nc_woa_id = Dataset(WOAFile[0], 'r' )
				woamask = nc_woa_id.variables["mask"][:]
				# shift mask as not greenwich centered, not the WOA way
				# woamask = np.roll(woamask, int(woamask.shape[1]/2.), axis = 1) 


# ---------- 5 OPEN NETCDF ----------
# open Ancillary Data, WOA masks with depth and areas
AncilMaskFile = glob.glob(AncilMaskFileName)
nc_ancil_id = Dataset(AncilMaskFile[0], 'r' )
vmasked = nc_ancil_id.variables["VOLUME_MASKED"][:,:,:].data 

# open mesh mask.nc
meshMaskFile = glob.glob(MeshMaskFileName)
nc_mesh_id = Dataset(meshMaskFile[0], 'r' )
tmeshmask = nc_mesh_id.variables["tmask"][0,:,:,:].data 

# open basin_mask.nc
basinFile = glob.glob(basinFileName)
nc_basin_id = Dataset(basinFile[0], 'r' )
mask_area = nc_basin_id.variables["AREA"][:].data #m^2
mask_vol = nc_basin_id.variables["VOLUME"][:].data #m^3
landMask = np.copy(mask_area)
landMask[ landMask > 0 ] = 1
landMask[ landMask == 0 ] = np.nan
volMask = np.copy(mask_vol)
volMask[ volMask > 0 ] = 1
volMask[ volMask == 0 ] = np.nan

# open regions file 
regFile = glob.glob(RegionFileName)
nc_reg_id = Dataset(regFile[0], 'r' )
regions = []
regions.append(nc_reg_id.variables['ARCTIC'][:].data) 
regions.append(nc_reg_id.variables['A1'][:].data) 
regions.append(nc_reg_id.variables['P1'][:].data) 
regions.append(nc_reg_id.variables['A2'][:].data + nc_reg_id.variables['P2'][:].data ) 
regions.append(nc_reg_id.variables['A3'][:].data + nc_reg_id.variables['P3'][:].data + nc_reg_id.variables['I3'][:].data ) 
regions.append(nc_reg_id.variables['A4'][:].data + nc_reg_id.variables['P4'][:].data + nc_reg_id.variables['I4'][:].data ) 
regions.append(nc_reg_id.variables['A5'][:].data) 
regions.append(nc_reg_id.variables['P5'][:].data) 
regions.append(nc_reg_id.variables['I5'][:].data) 

# open RECCAP regions file
RECCAPFile = glob.glob(RECCAPFileName)
nc_reccap_id = Dataset(RECCAPFile[0], 'r' )
regions.append(nc_reccap_id.variables['open_ocean_0'][:].data) #9
regions.append(nc_reccap_id.variables['open_ocean_1'][:].data) #10
regions.append(nc_reccap_id.variables['open_ocean_2'][:].data) #11
regions.append(nc_reccap_id.variables['open_ocean_3'][:].data) #12
regions.append(nc_reccap_id.variables['open_ocean_4'][:].data) #13
regions.append(nc_reccap_id.variables['atlantic_0'][:].data) #14
regions.append(nc_reccap_id.variables['atlantic_1'][:].data) #15
regions.append(nc_reccap_id.variables['atlantic_2'][:].data) #16
regions.append(nc_reccap_id.variables['atlantic_3'][:].data) #17
regions.append(nc_reccap_id.variables['atlantic_4'][:].data) #18
regions.append(nc_reccap_id.variables['atlantic_5'][:].data) #19
regions.append(nc_reccap_id.variables['pacific_0'][:].data) #20
regions.append(nc_reccap_id.variables['pacific_1'][:].data) #21
regions.append(nc_reccap_id.variables['pacific_2'][:].data) #22
regions.append(nc_reccap_id.variables['pacific_3'][:].data) #23
regions.append(nc_reccap_id.variables['pacific_4'][:].data) #24
regions.append(nc_reccap_id.variables['pacific_5'][:].data) #25
regions.append(nc_reccap_id.variables['indian_0'][:].data) #26
regions.append(nc_reccap_id.variables['indian_1'][:].data) #27
regions.append(nc_reccap_id.variables['arctic_0'][:].data) #28
regions.append(nc_reccap_id.variables['arctic_1'][:].data) #29
regions.append(nc_reccap_id.variables['arctic_2'][:].data) #30
regions.append(nc_reccap_id.variables['arctic_3'][:].data) #31
regions.append(nc_reccap_id.variables['southern_0'][:].data) #32
regions.append(nc_reccap_id.variables['southern_1'][:].data) #33
regions.append(nc_reccap_id.variables['southern_2'][:].data) #34
regions.append(nc_reccap_id.variables['seamask'][:].data) #35
regions.append(nc_reccap_id.variables['coast'][:].data) #36

#  additional regions (the above is a dashboard output), for tuning we need to seperate the regions
regions.append(nc_reg_id.variables['ATL'][:].data)  # 37
regions.append(nc_reg_id.variables['PAC'][:].data)  # 38
regions.append(nc_reg_id.variables['IND'][:].data)  # 39
regions.append(nc_reg_id.variables['SO'][:].data)  # 40
regions.append(nc_reg_id.variables['ARCTIC'][:].data)  # 41
regions.append(nc_reg_id.variables['P1'][:].data)  # 42
regions.append(nc_reg_id.variables['P2'][:].data)  # 43
regions.append(nc_reg_id.variables['P3'][:].data)  # 44
regions.append(nc_reg_id.variables['P4'][:].data)  # 45
regions.append(nc_reg_id.variables['P5'][:].data)  # 46
regions.append(nc_reg_id.variables['A1'][:].data)  # 47
regions.append(nc_reg_id.variables['A2'][:].data)  # 48
regions.append(nc_reg_id.variables['A3'][:].data)  # 49
regions.append(nc_reg_id.variables['A4'][:].data)  # 50
regions.append(nc_reg_id.variables['A5'][:].data)  # 51
regions.append(nc_reg_id.variables['I3'][:].data)  # 52
regions.append(nc_reg_id.variables['I4'][:].data)  # 53
regions.append(nc_reg_id.variables['I5'][:].data)  # 54

for region in regions:
	region[ region == 0 ] = np.nan

log.info("Data read from: "+RegionFileName)

# format and open netcdf file
runFileNames=[]
years = []

for year in range(yearFrom,yearTo+1):
	gridFileName = ''
	ptrcFileName = ''
	diadFileName = ''
	diadFileName2d = ''
	diadFileName3d = ''

	if len(glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_grid_T.nc")) > 0:
		gridFileName = glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_grid_T.nc")[0]

	if len(glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_ptrc_T.nc")) > 0:
		ptrcFileName = glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_ptrc_T.nc")[0]

	if len(glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_diad_T.nc")) > 0:
		diadFileName = glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_diad_T.nc")[0]

	if len(glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_dia2d_T.nc")) > 0:
		diadFileName2d = glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_dia2d_T.nc")[0]

	if len(glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_dia3d_T.nc")) > 0:
		diadFileName3d = glob.glob("ORCA2_1m_"+str(year)+"0101_"+str(year)+"1231_dia3d_T.nc")[0]

	runFileNames.append((gridFileName,ptrcFileName,diadFileName, diadFileName2d, diadFileName3d))
	log.info("Run data for year "+str(year)+": "+gridFileName+" "+ptrcFileName+" "+diadFileName+" "+diadFileName2d+" "+diadFileName3d)
	years.append(year)

# open files and get ids
nc_run_ids = []
nc_runFileNames=[]

for runID in runFileNames:  # each runID is the names of each output file for a specific year
	nc_ids = np.zeros((5)) - 1 # set to -1 as null value
	nc_avail = []
	nc_runFileNames_avail = []
	for f in range(0,len(runID)): # loop over the output files for that year
		if runID[f] != '':
			nc_id = Dataset(runID[f], 'r')
			nc_avail.append(nc_id)
			nc_runFileNames_avail.append(runID[f])

	nc_run_ids.append(nc_avail) # nc_run_ids is a list per year, each item being the nc_id for each output file
	nc_runFileNames.append(nc_runFileNames_avail)


# ---------- 6 PROCESS NETCDF OUTPUTS ----------
# Null outputs
nullAnnual = -1
nullMonthly = np.array( [ np.zeros((6))-1 for r in range(0,12)] )

# ----- 6.1 PROCESS SURFACE VARIABLES -----
for var in varSurface:
	for n in range(0,len(nc_run_ids)):
		# try each file to find variables
		varName = var[0]
		units = var[1]
		reg = var[5]
		latLim = [-90,90]
		lonLim = [-180,180]
		if reg == -1:
			lonLim = var[2].split(';')
			latLim = var[3].split(';')

		found = False
		for i in range(0,len(nc_run_ids[n])):
			try:
				data = nc_run_ids[n][i].variables[varName][:].data
				if len(data.shape) == 4:
					log.info(varName+" is volume data, taking surface values" )
					data = data[:,0,:,:]
				val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
				val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
				found = True
				log.info(varName+" found in "+str(nc_runFileNames[n][i]) )
			except KeyError as e:
				t = 1 

		try: 
			unitsToUse = ListOfUnits[units]
		except KeyError as e:
			log.info("Unit: "+str(units)+" not found, using raw data")
			unitsToUse = 1

		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionLandMask = np.copy(landMask)
		if reg != -1:
			regionLandMask = regionLandMask * regions[reg]

		if found == False:
			log.info(varName+" not found")
			var[-1].append((nullAnnual, nullMonthly))
		else:
			outputTotal= surfaceData(data, val_lons, val_lats, unitsToUse, mask_area, regionLandMask, volMask, missingVal, lonLim, latLim)
			var[-1].append(outputTotal) 

# ----- 6.2 PROCESS LEVEL VARIABLES -----
for var in varLevel:
	for n in range(0,len(nc_run_ids)):

		# try each file to find variables
		varName = var[0]
		level = var[1]
		units = var[2]
		reg = var[6]
		latLim = [-90,90]
		lonLim = [-180,180]
		if reg == -1:
			lonLim = var[3].split(';')
			latLim = var[4].split(';')

		found = False
		for i in range(0,len(nc_run_ids[n])):
			try:
				if len(nc_run_ids[n][i].variables[varName].dimensions) == 3:
					log.info(varName+" in "+str(nc_runFileNames[n][i])+ " is 2D data, try using Surface")
				if len(nc_run_ids[n][i].variables[varName].dimensions) == 4:
					data = nc_run_ids[n][i].variables[varName][:].data
					val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
					val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
					found = True
					log.info(varName+" found in "+str(nc_runFileNames[n][i]) )
			except KeyError as e:
				t = 1

		try: 
			unitsToUse = ListOfUnits[units]
		except KeyError as e:
			log.info("Unit: "+str(units)+" not found, using raw data")
			unitsToUse = 1

		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionVolMask = np.copy(volMask)
		if reg != -1:
			for z in range(0,regionVolMask.shape[0]):
				regionVolMask[z,:,:] = regionVolMask[z,:,:] * regions[reg]

		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionLandMask = np.copy(landMask)
		if reg != -1:
			regionLandMask = regionLandMask * regions[reg]

		if found == False:
			log.info(varName+" not found")
			var[-1].append((nullAnnual, nullMonthly))
		else:
			outputTotal = levelData(data, val_lons, val_lats, unitsToUse, mask_area, regionLandMask, regionVolMask, missingVal, lonLim, latLim, level)
			var[-1].append(outputTotal) 

# ----- 6.3 PROCESS VOLUME VARIABLES -----
for var in varVolume:
	for n in range(0,len(nc_run_ids)): # loop over years
		# try each file to find variables
		varName = var[0]
		units = var[1]
		reg = var[5]
		latLim = [-90,90]
		lonLim = [-180,180]
		if reg == -1:
			lonLim = var[2].split(';')
			latLim = var[3].split(';')

		found = False
		for i in range(0,len(nc_run_ids[n])):
			try:
				if len(nc_run_ids[n][i].variables[varName].dimensions) == 3:
					log.info(varName+" in "+str(nc_runFileNames[n][i])+ " is 2D data, try using Surface")
				if len(nc_run_ids[n][i].variables[varName].dimensions) == 4:
					data = nc_run_ids[n][i].variables[varName][:].data
					val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
					val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
					found = True
					log.info(varName+" found in "+str(nc_runFileNames[n][i]) )
			except KeyError as e:
				t = 1

		try: 
			unitsToUse = ListOfUnits[units]
		except KeyError as e:
			log.info("Unit: "+str(units)+" not found, using raw data")
			unitsToUse = 1

		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionVolMask = np.copy(volMask)
		if reg != -1:
			for z in range(0,regionVolMask.shape[0]):
				regionVolMask[z,:,:] = regionVolMask[z,:,:] * regions[reg]

		if found == False:
			log.info(varName+" not found")
			var[-1].append((nullAnnual, nullMonthly))
		else:
			# do a check for values integrated in the model output
			if len(data.shape) == 4:
				outputTotal= volumeData(data, val_lons, val_lats, unitsToUse, mask_vol, landMask, regionVolMask, missingVal, lonLim, latLim)
				var[-1].append(outputTotal) 

# ----- 6.4 PROCESS INTEGRATION VARIABLES -----
for var in varInt:
	for n in range(0,len(nc_run_ids)):
		# try each file to find variables
		varName = var[0]
		depthFrom = var[1]
		depthTo = var[2]
		units = var[3]
		reg = var[7]
		latLim = [-90,90]
		lonLim = [-180,180]
		if reg == -1:
			lonLim = var[4].split(';')
			latLim = var[5].split(';')

		found = False
		for i in range(0,len(nc_run_ids[n])):
			try:
				if len(nc_run_ids[n][i].variables[varName].dimensions) == 3:
					log.info(varName+" in "+str(nc_runFileNames[n][i])+ " is 2D data, try using Surface")
				if len(nc_run_ids[n][i].variables[varName].dimensions) == 4:
					data = nc_run_ids[n][i].variables[varName][:].data
					val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
					val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
					found = True
					log.info(varName+" found in "+str(nc_runFileNames[n][i]) )
			except KeyError as e:
				t = 1

		try: 
			unitsToUse = ListOfUnits[units]
		except KeyError as e:
			log.info("Unit: "+str(units)+" not found, using raw data")
			unitsToUse = 1
        
 		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionVolMask = np.copy(volMask)
		if reg != -1:
			for z in range(0,regionVolMask.shape[0]):
				regionVolMask[z,:,:] = regionVolMask[z,:,:] * regions[reg]

		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionLandMask = np.copy(landMask)
		if reg != -1:
			regionLandMask = regionLandMask * regions[reg]

		if found == False:
			log.info(varName+" not found")
			var[-1].append((nullAnnual, nullMonthly))
		else:
			outputTotal = intergrateData(data, val_lons, val_lats, depthFrom, depthTo, unitsToUse, mask_vol, landMask, regionVolMask, missingVal, lonLim, latLim)
			var[-1].append(outputTotal) 

# ----- 6.5 PROCESS AVERAGE VARIABLES -----
for var in varTotalAve:
	for n in range(0,len(nc_run_ids)):
		# try each file to find variables
		varNames = var[0].split('+')
		nVars = len(varNames)
		depthFrom = var[1]
		depthTo = var[2]
		units = var[3]
		reg = var[7]
		latLim = [-90,90]
		lonLim = [-180,180]
		if reg == -1:
			lonLim = var[4].split(';')
			latLim = var[5].split(';')

		allData = []
		for varName in varNames:
			found = False
			for i in range(0,len(nc_run_ids[n])):
				try:
					data = nc_run_ids[n][i].variables[varName][:].data
					val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
					val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
                    
					if len(data.shape) == 3:
						# reshape to 4 dims, add depth dims
						newData = np.zeros((data.shape[0],volMask.shape[0],data.shape[1],data.shape[2]))
						for z in range(0,volMask.shape[0]):
							newData[:,z,:,:] = data
						data = newData
                    
					if len(data.shape) == 4:
						allData.append(data)
						found = True
						log.info(varName+" found in "+str(nc_runFileNames[n][i]) )
				except KeyError as e:
					t = 1

			if found == False:
				log.info("Not all of "+str(varName)+" are found in files, "+str(varName)+" not found")

		try:
			unitsToUse = ListOfUnits[units]
		except KeyError as e:
			log.info("Unit: "+str(units)+" not found, using raw data")
			unitsToUse = 1

		# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
		regionVolMask = np.copy(volMask)
		if reg != -1:
			for z in range(0,regionVolMask.shape[0]):
				regionVolMask[z,:,:] = regionVolMask[z,:,:] * regions[reg]

		if found == False:
			log.info(varName+" not found")
			var[-1].append((nullAnnual, nullMonthly))
		else:
			outputAllTotal = 0
			monthlyAllOutput = np.array( [ np.zeros((6)) for r in range(0,12)] )

			for d in range(0,len(allData)):
				output = volumeDataAverage(allData[d], val_lons, val_lats,  depthFrom, depthTo, unitsToUse, mask_vol, landMask, regionVolMask, missingVal, lonLim, latLim)
				outputAllTotal = outputAllTotal + output[0] 
				for m in range(0,12):
					for f in range(0,6):
						if len(output[1]) > m:
							monthlyAllOutput[m][f] = monthlyAllOutput[m][f]+ output[1][m][f]  

			outputTotal = (outputAllTotal, monthlyAllOutput)
			var[-1].append(outputTotal) 

# ----- 6.6 PROCESS OBSERVATION COMPARSION -----
for obs in obsComparisons:

	obsName = obs[0]
	obsVar = obs[1]
	varName = obs[2].split('+') # this is an array for obs to allow for totalling up (e.g. CHL's)
	depthObs = obs[3]
	depthVar = obs[4]
	gamFlag = obs[5]
	reg = obs[9]
	latLim = [-90,90]
	lonLim = [-180,180]
	if reg == -1:
		lonLim = obs[6].split(';') # this splits the value into 2 limits
		latLim = obs[7].split(';') # this splits the value into 2 limits

	sourceObs = None
	for l in ListOfObservations:
		if obs[0] == l['name']:
			sourceObs = l
	nc_obs_id = Dataset(sourceObs['path'], 'r')

	# adjust landmask for region definitions - this will not cover everything if you add them together due to inland seas etc...
	print(reg, volMask.shape)
	regionVolMask = np.copy(volMask)
	if reg != -1:
		for z in range(0,regionVolMask.shape[0]):
			regionVolMask[z,:,:] = regionVolMask[z,:,:] * regions[reg]

	# look for var in model data
	found_var = False
	for n in range(0,len(nc_run_ids)):
		var_data_arr = [] # to sum up variables like CHL
		for i in range(0,len(nc_run_ids[n])):
			for v in range(0,len(varName)):
				try:
					var_data = nc_run_ids[n][i].variables[varName[v]][:].data
					if len(nc_run_ids[n][i].variables[varName[v]][:].data.shape) == 4: # take surface data only for observations
						log.info(varName[v]+" has 4 dimensions, only taking data from depth "+str(depthVar) )
						var_data = nc_run_ids[n][i].variables[varName[v]][:,depthVar,:,:].data

					var_data_arr.append(var_data)
					val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
					val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
					found_var = True
					log.info(varName[v]+" found in "+str(nc_runFileNames[n][i]) )
				except KeyError as e:
					t = 1

		if np.array(var_data_arr).shape[0] == 0: # if no data found for this year skip
			found_var = False
		else:
			var_data = np.sum(np.array(var_data_arr),axis=0)

			# multiply by a factor
			var_data = var_data * sourceObs['factor']

			# look for var in obs data
			found_obs = False
			try:
				obs_lat = nc_obs_id.variables[sourceObs['latVar']][:].data
				obs_lon = nc_obs_id.variables[sourceObs['lonVar']][:].data
				obs_time = nc_obs_id.variables[sourceObs['timeVar']][:].data

				t_factor = 1
				if "days" in nc_obs_id.variables[sourceObs['timeVar']].units:
					# convert days to seconds
					t_factor = 60*60*24

				obs_time = obs_time * t_factor
				first_time = sourceObs['origin'] + datetime.timedelta(0,int(obs_time[0])) # adds seconds onto origin date

				# assume data is monthly same as model outputs
				year = years[n]
				timePoint = (relativedelta.relativedelta(datetime.datetime(year,1,1,0,0,0) , first_time)).years*12

				# if data is averaged then ignore the start point and take the first point
				if sourceObs['climatological'] == True:
					timePoint = 0

				if len(nc_obs_id.variables[obsVar].dimensions) == 3:
					obs_data = nc_obs_id.variables[obsVar][timePoint:timePoint+12,:,:].data
				if len(nc_obs_id.variables[obsVar].dimensions) == 4:
					log.info(obsVar+" has 4 dimensions, only taking data from depth "+str(depthObs) )
					obs_data = nc_obs_id.variables[obsVar][timePoint:timePoint+12,depthObs,:,:].data

				log.info("Timepoint found: "+str(timePoint))
				found_obs = True
				log.info(obsVar+" found in "+obsName)

				# incomplete data
				if obs_data.shape[0] != 12 or timePoint < 0:
					found_obs = False
					log.info(obsVar+" data incomplete or model data is pre obs data")
				else:
					# get conversion data, make sure we add missing val tidy up routines
					if sourceObs['conversion'] != None:
						log.info("Conversion File: "+sourceObs['conversion']+" for "+sourceObs['conversionName'])
						nc_conv_id = Dataset(sourceObs['conversion'], 'r')

					if len(nc_conv_id.variables[sourceObs['conversionName']].dimensions) == 3:
						conv_data = nc_conv_id.variables[sourceObs['conversionName']][timePoint:timePoint+12,:,:].data
					if len(nc_conv_id.variables[sourceObs['conversionName']].dimensions) == 4:
						conv_data = nc_conv_id.variables[sourceObs['conversionName']][timePoint:timePoint+12,depthObs,:,:].data

					obs_data = obs_data * conv_data
					obs_data[ obs_data < (-missingVal/100.) ] = missingVal
					obs_data[ obs_data > (missingVal/100.) ] = missingVal
			except KeyError as e:
				t = 1

		if found_var == True and found_obs == True:
			log.info("Comparing to obs data")
			outputTotal= observationData(obs_data, obs_lon, obs_lat, var_data, val_lons, val_lats, sourceObs['centered'], regionVolMask[depthVar,:,:] , missingVal, gamFlag, lonLim, latLim)
			obs[-1].append(outputTotal) 
		else:
			log.info(obs[2]+" or "+obsVar+" not found correctly")
			nullMonthly = np.array( [ -1 for r in range(0,12)] )
			obs[-1].append((nullAnnual, nullMonthly, [[],[]] ))

# ----- 6.7 PROCESS EMERGING PROPERTY VALUES -----
for prop in properties:

	propName = prop[0]

	if propName == "Bloom":
		propVar = prop[1].split('+')
		depthFrom = int(prop[2])
		depthTo = int(prop[3])
		lonLim = prop[4].split(';') # this splits the value into 2 limits
		latLim = prop[5].split(';') # this splits the value into 2 limits

		# look for var in model data
		found_var = False
		for n in range(0,len(nc_run_ids)):
			var_data_arr = []
			for i in range(0,len(nc_run_ids[n])):
				for v in range(0,len(propVar)):
					try:
						var_data = nc_run_ids[n][i].variables[propVar[v]][:].data

						if len(nc_run_ids[n][i].variables[propVar[v]][:].data.shape) == 4: # take surface data only for observations
							log.info(propVar[v]+" has 4 dimensions, summing data for depths "+str(depthFrom)+" "+ str(depthTo))
							var_data = nc_run_ids[n][i].variables[propVar[v]][:,depthFrom:depthTo+1,:,:].data
							var_data = np.sum(var_data,axis=1)

						var_data_arr.append(var_data)
						val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
						val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
						found_var = True
						log.info(propVar[v]+" found in "+str(nc_runFileNames[n][i]) )
					except KeyError as e:
						t = 1

			var_data = np.sum(np.array(var_data_arr),axis=0)

			if found_var == False:
				log.info(prop[1]+" not found")
				prop[-1].append((-1, -1, -1, -1, ["none", "none", "none", "none"], -1 ))
			else:
				outputProp = bloom(var_data, val_lons, val_lats, missingVal, lonLim, latLim)
				prop[-1].append(outputProp) 

	if propName == "Trophic":
		propVarList = [ prop[1][0].split('+'), prop[1][1].split('+'), prop[1][2].split('+') ]
		depthFrom = int(prop[2])
		depthTo = int(prop[3])
		lonLim = prop[4].split(';') # this splits the value into 2 limits
		latLim = prop[5].split(';') # this splits the value into 2 limits

		# Three sets of data defined by propVar   
		all_var_data = [ None, None, None ]
        
		# look for var in model data
		found_var = False
		for n in range(0,len(nc_run_ids)):
			for p in range(0,3):
				var_data_arr = []  # for adding up variables
				propVar = propVarList[p]
				for i in range(0,len(nc_run_ids[n])):
					try:
						for v in range(0,len(propVar)):
							var_data = nc_run_ids[n][i].variables[propVar[v]][:].data
							if len(nc_run_ids[n][i].variables[propVar[v]][:].data.shape) == 4: # take surface data only for observations
								log.info(propVar[v]+" has 4 dimensions, summing data for depths "+str(depthFrom)+" "+ str(depthTo))
								var_data = nc_run_ids[n][i].variables[propVar[v]][:,depthFrom:depthTo+1,:,:].data
								var_data = np.sum(var_data,axis=1)

							var_data_arr.append(var_data)
							val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
							val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
							found_var = True
							log.info(propVar[v]+" found in "+str(nc_runFileNames[n][i]) )
					except KeyError as e:
						t = 1

				all_var_data[p] = np.sum(np.array(var_data_arr),axis=0)

			if found_var == False:
				log.info(prop[1]+" none found")
				prop[-1].append((-1, -1, -1, -1, ["none", "none", "none", "none"], -1 ))
			else:
				outputProp = trophic(all_var_data, val_lons, val_lats, missingVal, lonLim, latLim)
				prop[-1].append(outputProp)


# ---------- 7 PROCESS NETCDF REGRIDDED MAP ----------
# get target mesh values
target_lon = np.arange(0.5,360.5,1) # triggers setting in regrid as < 2 dims
target_lat = np.arange(-89.5,90.5,1)

# process variables to map
for n in range(0,len(nc_run_ids)): # one file for each year
	# create output file
	if len(varMap) > 0:
		for var in varMap:

			# try each file to find variables
			varName = var[0]
			level = var[1]
			if level != 'all':
				lev = int(level)

			outputFileName = "WOA_"+varName+"_"+str(years[n])+"_"+level+".nc"

			nc_out_id = Dataset(outputFileName, 'w', format='NETCDF4_CLASSIC')
			nc_out_id.createDimension("lon", 360)
			nc_out_id.createDimension("lat", 180)

			times = nc_run_ids[n][0].variables['time_counter'][:].data
			nc_out_id.createDimension("time", None)
			if level == 'all':
				depths = nc_run_ids[n][0].variables['deptht'][:].data
				nc_out_id.createDimension("deptht", len(depths))
			else:
				depths = nc_run_ids[n][0].variables['deptht'][:].data
				nc_out_id.createDimension("deptht", 1)
            
			nc_out_id.createVariable("lon", 'f', ('lon'))
			nc_out_id.createVariable("lat", 'f', ('lat'))
			nc_out_id.createVariable("deptht", 'f', ('deptht'))
			nc_out_id.createVariable("time", 'f', ('time'))

			found = False

			for i in range(0,len(nc_run_ids[n])):
				try:
					data = nc_run_ids[n][i].variables[varName][:].data
					if len(data.shape) == 4:
						log.info(varName+" is volume data" )
						nc_v = nc_out_id.createVariable(varName, 'f', ('time','deptht','lat','lon'), fill_value=missingVal)
						if level == 'all':
							nc_out_id.variables['deptht'][:] = depths
						else:
							nc_out_id.variables['deptht'][:] = depths[lev]
					else:
						nc_v = nc_out_id.createVariable(varName, 'f', ('time','lat','lon'), fill_value=missingVal)

					for ncattr in nc_run_ids[n][i].variables[varName].ncattrs():
						if ncattr != '_FillValue':
							nc_v.setncattr(ncattr, nc_run_ids[n][i].variables[varName].getncattr(ncattr) )
					nc_v.setncattr("missing_value", np.array(missingVal,'f'))

					val_lats = nc_run_ids[n][i].variables['nav_lat'][:].data
					val_lons = nc_run_ids[n][i].variables['nav_lon'][:].data
					nc_out_id.variables['lon'][:] = target_lon
					nc_out_id.variables['lat'][:] = target_lat
					found = True
					log.info(varName+" found in "+str(nc_runFileNames[n][i]) )
				except KeyError as e:
					t = 1 

			if found == False:
				log.info(varName+" not found")
			else:
				woaind = np.where( woamask == 0 )
				if len(data.shape) == 3:
					regriddedData = regrid(data, val_lons, val_lats, target_lon, target_lat, tmeshmask[0,:,:], missingVal)
					for t in range(0,data.shape[0]):
						# tidy up
						regriddedData_slice = regriddedData[t,:,:]
						regriddedData_slice[ woaind ] = missingVal
						regriddedData_slice[ np.isnan(regriddedData_slice) ] = missingVal
						nc_out_id.variables[varName][t,:,:] = regriddedData_slice
                
						nc_out_id.variables['time'][t] = times[t]

				if len(data.shape) == 4:
					if level == 'all':
						regriddedData = np.zeros( (data.shape[0], data.shape[1], 180, 360) )
						for z in range(0,data.shape[1]):
							regriddedData[:,z,:,:] = regrid(data[:,z,:,:], val_lons, val_lats, target_lon, target_lat, tmeshmask[z,:,:], missingVal)
					else:
						print(tmeshmask.shape, lev)
						regriddedData = np.zeros( (data.shape[0], 1, 180, 360) )
						regriddedData[:,0,:,:] = regrid(data[:,lev,:,:], val_lons, val_lats, target_lon, target_lat, tmeshmask[lev,:,:], missingVal)

					# mask at depth using tmeshmask
					for t in range(0,regriddedData.shape[0]):
						for z in range(0,regriddedData.shape[1]):
							regriddedData_slice = regriddedData[t,z,:,:]
							if level == 'all':
								meshind = np.where( vmasked[z,:,:] == 0)
							else:
								meshind = np.where( vmasked[lev,:,:] == 0)

							regriddedData_slice[ meshind ] = missingVal
							regriddedData_slice[ np.isnan(regriddedData_slice) ] = missingVal
							nc_out_id.variables[varName][t,z,:,:] = regriddedData[t,z,:,:]
                            
						nc_out_id.variables['time'][t] = times[t]

		nc_out_id.close()


# ---------- 8 OUTPUT ANNUAL VALUES ----------

outputFileName_sur = "breakdown.sur.annual.dat"
outputFileName_lev = "breakdown.lev.annual.dat"
outputFileName_vol = "breakdown.vol.annual.dat"
outputFileName_int = "breakdown.int.annual.dat"
outputFileName_ave = "breakdown.ave.annual.dat"
outputFileName_obs_annual = "breakdown.obs.annual.dat"

# ----- 8.1 SURFACE FILE -----
if (os.path.exists(outputFileName_sur)):
	file_sur = open(outputFileName_sur, "a") 
else:
	file_sur = open(outputFileName_sur, "w") 

	# write column headers
	file_sur.write("year\t")
	for var in varSurface:
		file_sur.write(str(var[0])+"\t")
	file_sur.write("\n\t")
	for var in varSurface:
		file_sur.write(str(var[1])+"\t")
	file_sur.write("\n\t")
	for var in varSurface:
		file_sur.write(str(var[-3])+"\t")  # key
	file_sur.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	file_sur.write(str(year)+"\t")
	for var in varSurface:
		file_sur.write(str( format(var[-1][y][0],".4e") )+"\t")
	file_sur.write("\n")


# ----- 8.2 LEVEL FILE -----
if (os.path.exists(outputFileName_lev)):
	file_lev = open(outputFileName_lev, "a") 
else:
	file_lev = open(outputFileName_lev, "w") 

	# write column headers
	file_lev.write("year\t")
	for var in varLevel:
		file_lev.write(str(var[0])+"\t")
	file_lev.write("\n\t")
	for var in varLevel:
		file_lev.write(str(var[2])+"\t")
	file_lev.write("\n\t")
	for var in varLevel:
		file_lev.write(str(var[-3])+"\t")  # key
	file_lev.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	file_lev.write(str(year)+"\t")
	for var in varLevel:
		file_lev.write(str( format(var[-1][y][0],".4e") )+"\t")
	file_lev.write("\n")


# ----- 8.3 VOLUME FILE -----
if (os.path.exists(outputFileName_vol)):
	file_vol = open(outputFileName_vol, "a") 
else:
	file_vol = open(outputFileName_vol, "w") 

	# write column headers
	file_vol.write("year\t")
	for var in varVolume:
		file_vol.write(str(var[0])+"\t")
	file_vol.write("\n\t")
	for var in varVolume:
		file_vol.write(str(var[1])+"\t")
	file_vol.write("\n\t")
	for var in varVolume:
		file_vol.write(str(var[-3])+"\t")
	file_vol.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	file_vol.write(str(year)+"\t")
	for var in varVolume:
		file_vol.write(str( format(var[-1][y][0],".4e") )+"\t")
	file_vol.write("\n")


# ----- 8.4 INTEGRATION FILE -----
if (os.path.exists(outputFileName_int)):
	file_int = open(outputFileName_int, "a") 
else:
	file_int = open(outputFileName_int, "w") 

	# write column headers
	file_int.write("year\t")
	for var in varInt:
		file_int.write(str(var[0])+"\t")
	file_int.write("\n\t")
	for var in varInt:
		file_int.write(str(var[3])+"\t")
	file_int.write("\n\t")
	for var in varInt:
		file_int.write(str(var[-3])+"\t")
	file_int.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	file_int.write(str(year)+"\t")
	for var in varInt:
		file_int.write(str( format(var[-1][y][0],".4e") )+"\t")
	file_int.write("\n")


# ----- 8.5 AVERAGE FILE -----
if (os.path.exists(outputFileName_ave)):
	file_ave = open(outputFileName_ave, "a") 
else:
	file_ave = open(outputFileName_ave, "w") 

	# write column headers
	file_ave.write("year\t")
	for var in varTotalAve:
		name = var[0]
		file_ave.write(name+"\t")
	file_ave.write("\n\t")
	for var in varTotalAve:
		file_ave.write(str(var[3])+"\t")
	file_ave.write("\n\t")
	for var in varTotalAve:
		file_ave.write(str(var[-3])+"\t")
	file_ave.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	file_ave.write(str(year)+"\t")
	for var in varTotalAve:
		file_ave.write(str( format(var[-1][y][0],".4e") )+"\t")
	file_ave.write("\n")


# ----- 8.6 OBSERVATIONS FILE -----
if (os.path.exists(outputFileName_obs_annual)):
	file_obs = open(outputFileName_obs_annual, "a") 
else:
	file_obs = open(outputFileName_obs_annual, "w") 

	# write column headers
	file_obs.write("year\t")
	for obs in obsComparisons:
		file_obs.write(str(obs[2])+"("+str(obs[1])+")\t")
	file_obs.write("\n\t")

	# write column key
	for obs in obsComparisons:
		file_obs.write(str(obs[-2])+"\t")
	file_obs.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	file_obs.write(str(year)+"\t")
	for obs in obsComparisons:
		file_obs.write(str( format(obs[-1][y][0],".4e") )+"\t")
	file_obs.write("\n")


# ----- 8.7 PROPERTIES FILE -----
for prop in properties:

	outputFileName_prop = "breakdown."+prop[0]+"."+prop[-2]+".dat" # use key for output file

	if (os.path.exists(outputFileName_prop)):
		file_prop = open(outputFileName_prop, "a") 
	else:
		file_prop = open(outputFileName_prop, "w") 

		headings = prop[-1][0][-2] #last var is list of results, first result, last entry is headings
		# write column headers
		file_prop.write("year\t")
		for col in headings:
			file_prop.write(col+"\t")
		file_prop.write("\n\t")

		for col in headings:
			if prop[0] == "Trophic":
				header = prop[1][0]+";"+prop[1][1]+";"+prop[1][2]
			else:
				header = prop[1]
			file_prop.write(header+"\t")
		file_prop.write("\n")

	# write data
	for year in range(yearFrom,yearTo+1):
		y = year-yearFrom
		file_prop.write(str(year)+"\t")
		headings = prop[-1][0][-2] 
		for c in range(0,len(headings)):
			file_prop.write(str( format(prop[-1][y][c],".4e") )+"\t")
		file_prop.write("\n")


# ---------- 9 OUTPUT SPREAD VALUES ----------

outputFileName_sur = "breakdown.sur.spread.dat"
outputFileName_lev = "breakdown.lev.spread.dat"
outputFileName_vol = "breakdown.vol.spread.dat"
outputFileName_int = "breakdown.int.spread.dat"
outputFileName_ave = "breakdown.ave.spread.dat"

# ----- 9.1 SURFACE FILE -----
if (os.path.exists(outputFileName_sur)):
	file_sur = open(outputFileName_sur, "a") 
else:
	file_sur = open(outputFileName_sur, "w") 

	# write column headers
	file_sur.write("year\t")
	file_sur.write("month\t")
	for var in varSurface:
		for i in range(0,5):
			file_sur.write(str(var[0])+"\t")
	file_sur.write("\n\t\t")
    
	for var in varSurface:
		for i in range(0,5):
			file_sur.write(str(var[1])+"\t")
	file_sur.write("\n\t\t")

	for var in varSurface:
		for i in range(0,5):
			file_sur.write(str(var[-3])+"\t")
	file_sur.write("\n\t\t")

	for var in varSurface:
		file_sur.write("min(5pt)\t")
		file_sur.write("25pt\t")
		file_sur.write("median\t")
		file_sur.write("75pt\t")
		file_sur.write("max(95pt)\t")
	file_sur.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_sur.write(str(year)+"\t")
		file_sur.write(str(m)+"\t")

		for var in varSurface:
			# monthlyVals.append(var[-1][n][1][m][0]) # -1=list of results, n=year, 1=monthlyvals, m=month, 0=sum
			minVals = var[-1][y][1][m][1]
			firstPCVals = var[-1][y][1][m][2]
			medianVals = var[-1][y][1][m][3]
			thirdPCVals = var[-1][y][1][m][4]
			maxVals = var[-1][y][1][m][5]

			# file_sur.write(str( format( np.mean(monthlyVals),".4e" ) )+"\t") 
			file_sur.write(str( format( minVals,".4e" ) )+"\t") 
			file_sur.write(str( format( firstPCVals,".4e" ) )+"\t") 
			file_sur.write(str( format( medianVals,".4e" ) )+"\t") 
			file_sur.write(str( format( thirdPCVals,".4e" ) )+"\t") 
			file_sur.write(str( format( maxVals,".4e" ) )+"\t") 

		file_sur.write("\n")


# ----- 9.2 LEVEL FILE -----
if (os.path.exists(outputFileName_lev)):
	file_lev = open(outputFileName_lev, "a") 
else:
	file_lev = open(outputFileName_lev, "w") 

	# write column headers
	file_lev.write("year\t")
	file_lev.write("month\t")
	for var in varLevel:
		for i in range(0,5):
			file_lev.write(str(var[0])+"\t")
	file_lev.write("\n\t\t")
    
	for var in varLevel:
		for i in range(0,5):
			file_lev.write(str(var[2])+"\t")
	file_lev.write("\n\t\t")

	for var in varLevel:
		for i in range(0,5):
			file_lev.write(str(var[-3])+"\t")
	file_lev.write("\n\t\t")

	for var in varLevel:
		file_lev.write("min(5pt)\t")
		file_lev.write("25pt\t")
		file_lev.write("median\t")
		file_lev.write("75pt\t")
		file_lev.write("max(95pt)\t")
	file_lev.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_lev.write(str(year)+"\t")
		file_lev.write(str(m)+"\t")

		for var in varLevel:
			# monthlyVals.append(var[-1][n][1][m][0]) # -1=list of results, n=year, 1=monthlyvals, m=month, 0=sum
			minVals = var[-1][y][1][m][1]
			firstPCVals = var[-1][y][1][m][2]
			medianVals = var[-1][y][1][m][3]
			thirdPCVals = var[-1][y][1][m][4]
			maxVals = var[-1][y][1][m][5]

			# file_lev.write(str( format( np.mean(monthlyVals),".4e" ) )+"\t") 
			file_lev.write(str( format( minVals,".4e" ) )+"\t") 
			file_lev.write(str( format( firstPCVals,".4e" ) )+"\t") 
			file_lev.write(str( format( medianVals,".4e" ) )+"\t") 
			file_lev.write(str( format( thirdPCVals,".4e" ) )+"\t") 
			file_lev.write(str( format( maxVals,".4e" ) )+"\t") 

		file_lev.write("\n")


# ----- 9.3 VOLUME FILE -----
if (os.path.exists(outputFileName_vol)):
	file_vol = open(outputFileName_vol, "a") 
else:
	file_vol = open(outputFileName_vol, "w") 

	# write column headers
	file_vol.write("month\t")
	for var in varVolume:
		for i in range(0,5):
			file_vol.write(str(var[0])+"\t")
	file_vol.write("\n\t\t")

	for var in varVolume:
		for i in range(0,5):
			file_vol.write(str(var[1])+"\t")
	file_vol.write("\n\t\t")

	for var in varVolume:
		for i in range(0,5):
			file_vol.write(str(var[-3])+"\t")
	file_vol.write("\n\t\t")

	for var in varVolume:
		file_vol.write("min(5pt)\t")
		file_vol.write("25pt\t")
		file_vol.write("median\t")
		file_vol.write("75pt\t")
		file_vol.write("max(95pt)\t")
	file_vol.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_vol.write(str(year)+"\t")
		file_vol.write(str(m)+"\t")

		for var in varVolume:

			# monthlyVals.append(var[-1][n][1][m][0]) # -1=list of results, n=year, 1=monthlyvals, m=month, 0=sum
			minVals = var[-1][y][1][m][1]
			firstPCVals = var[-1][y][1][m][2]
			medianVals = var[-1][y][1][m][3]
			thirdPCVals = var[-1][y][1][m][4]
			maxVals = var[-1][y][1][m][5]

			# file_sur.write(str( format( np.mean(monthlyVals),".4e" ) )+"\t") 
			file_vol.write(str( format( minVals,".4e" ) )+"\t") 
			file_vol.write(str( format( firstPCVals,".4e" ) )+"\t") 
			file_vol.write(str( format( medianVals,".4e" ) )+"\t") 
			file_vol.write(str( format( thirdPCVals,".4e" ) )+"\t") 
			file_vol.write(str( format( maxVals,".4e" ) )+"\t") 
		file_vol.write("\n")


# ----- 9.4 INTEGRATION FILE -----
if (os.path.exists(outputFileName_int)):
	file_int = open(outputFileName_int, "a") 
else:
	file_int = open(outputFileName_int, "w") 

	# write column headers
	file_int.write("month\t")
	for var in varInt:
		for i in range(0,5):
			file_int.write(str(var[0])+"\t")
	file_int.write("\n\t\t")

	for var in varInt:
		for i in range(0,5):
			file_int.write(str(var[3])+"\t")
	file_int.write("\n\t\t")

	for var in varInt:
		for i in range(0,5):
			file_int.write(str(var[-3])+"\t")
	file_int.write("\n\t\t")

	for var in varVolume:
		file_int.write("min(5pt)\t")
		file_int.write("25pt\t")
		file_int.write("median\t")
		file_int.write("75pt\t")
		file_int.write("max(95pt)\t")
	file_int.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_int.write(str(year)+"\t")
		file_int.write(str(m)+"\t")

		for var in varInt:

			# monthlyVals.append(var[-1][n][1][m][0]) # -1=list of results, n=year, 1=monthlyvals, m=month, 0=sum
			minVals = var[-1][y][1][m][1]
			firstPCVals = var[-1][y][1][m][2]
			medianVals = var[-1][y][1][m][3]
			thirdPCVals = var[-1][y][1][m][4]
			maxVals = var[-1][y][1][m][5]

			# file_sur.write(str( format( np.mean(monthlyVals),".4e" ) )+"\t") 
			file_int.write(str( format( minVals,".4e" ) )+"\t") 
			file_int.write(str( format( firstPCVals,".4e" ) )+"\t") 
			file_int.write(str( format( medianVals,".4e" ) )+"\t") 
			file_int.write(str( format( thirdPCVals,".4e" ) )+"\t") 
			file_int.write(str( format( maxVals,".4e" ) )+"\t") 
		file_int.write("\n")


# ----- 9.5 AVERAGE FILE -----
if (os.path.exists(outputFileName_ave)):
	file_ave = open(outputFileName_ave, "a") 
else:
	file_ave = open(outputFileName_ave, "w") 

	# write column headers
	file_ave.write("month\t")
	for var in varTotalAve:
		for i in range(0,5):
			name = var[0]
			file_ave.write(name+"\t")
	file_ave.write("\n\t\t")

	for var in varTotalAve:
		for i in range(0,5):
			file_ave.write(str(var[1])+"\t")
	file_ave.write("\n\t\t")

	for var in varTotalAve:
		for i in range(0,5):
			file_ave.write(str(var[-3])+"\t")
	file_ave.write("\n\t\t")

	for var in varTotalAve:
		file_ave.write("min(5pt)\t")
		file_ave.write("25pt\t")
		file_ave.write("median\t")
		file_ave.write("75pt\t")
		file_ave.write("max(95pt)\t")
	file_ave.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_ave.write(str(year)+"\t")
		file_ave.write(str(m)+"\t")

		for var in varTotalAve:

			# monthlyVals.append(var[-1][n][1][m][0]) # -1=list of results, n=year, 1=monthlyvals, m=month, 0=sum
			minVals = var[-1][y][1][m][1]
			firstPCVals = var[-1][y][1][m][2]
			medianVals = var[-1][y][1][m][3]
			thirdPCVals = var[-1][y][1][m][4]
			maxVals = var[-1][y][1][m][5]

			# file_sur.write(str( format( np.mean(monthlyVals),".4e" ) )+"\t") 
			file_ave.write(str( format( minVals,".4e" ) )+"\t") 
			file_ave.write(str( format( firstPCVals,".4e" ) )+"\t") 
			file_ave.write(str( format( medianVals,".4e" ) )+"\t") 
			file_ave.write(str( format( thirdPCVals,".4e" ) )+"\t") 
			file_ave.write(str( format( maxVals,".4e" ) )+"\t") 
		file_ave.write("\n")


# ---------- 10 OUTPUT MONTHLY VALUES ----------

outputFileName_sur = "breakdown.sur.monthly.dat"
outputFileName_ave = "breakdown.ave.monthly.dat"

# ----- 10.1 SURFACE FILE -----
if (os.path.exists(outputFileName_sur)):
	file_sur = open(outputFileName_sur, "a") 
else:
	file_sur = open(outputFileName_sur, "w") 

	# write column headers
	file_sur.write("year\t")
	file_sur.write("month\t")
	for var in varSurface:
		file_sur.write(str(var[0])+"\t")
	file_sur.write("\n\t\t")
    
	for var in varSurface:
		file_sur.write(str(var[1])+"\t")
	file_sur.write("\n\t\t")

	for var in varSurface:
		file_sur.write(str(var[-3])+"\t")
	file_sur.write("\n\t\t")

	for var in varSurface:
		file_sur.write("total\t")
	file_sur.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_sur.write(str(year)+"\t")
		file_sur.write(str(m)+"\t")

		for var in varSurface:
			month_tot = var[-1][y][1][m][0]
			file_sur.write(str( format( month_tot,".4e" ) )+"\t") 

		file_sur.write("\n")


# ----- 10.2 LEVEL FILE -----
# ----- 10.3 VOLUME FILE -----
# ----- 10.4 INTEGRATION FILE -----

# ----- 10.5 AVERAGE FILE -----
if (os.path.exists(outputFileName_ave)):
	file_ave = open(outputFileName_ave, "a") 
else:
	file_ave = open(outputFileName_ave, "w") 

	# write column headers
	file_ave.write("year\t")
	file_ave.write("month\t")
	for var in varTotalAve:
		name = var[0]
		file_ave.write(name+"\t")
	file_ave.write("\n\t\t")

	for var in varTotalAve:
		file_ave.write(str(var[1])+"\t")
	file_ave.write("\n\t\t")

	for var in varTotalAve:
		file_ave.write(str(var[-3])+"\t")
	file_ave.write("\n\t\t")

	for var in varTotalAve:
		file_ave.write("total\t")
	file_ave.write("\n")

# write data
for year in range(yearFrom,yearTo+1):
	y = year-yearFrom
	for m in range(0,12):
		file_ave.write(str(year)+"\t")
		file_ave.write(str(m)+"\t")

		for var in varTotalAve:
			month_tot = var[-1][y][1][m][0]
			file_ave.write(str( format( month_tot,".4e" ) )+"\t") 
		file_ave.write("\n")

