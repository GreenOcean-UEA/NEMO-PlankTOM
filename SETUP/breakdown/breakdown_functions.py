import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset
import math
from scipy.ndimage import gaussian_filter

# ---------- 1 REGRID DATA ----------
def regrid(var, var_lons, var_lats, target_lon, target_lat, tmask, missingVal):

	tDim = var.shape[0]

	# duplicate edge pixels 
	for iy in range(0,var_lons.shape[0]):
		var[:,iy,0] = (var[:,iy,1]+var[:,iy,-2])/2
		var[:,iy,-1] =  var[:,iy,0]
	for ix in range(0,var_lons.shape[0]):
		var[:,-1,ix] = var[:,-2,ix] 

	# adjust nav_lon to WOA grid 
	var_lons[ var_lons < 0 ] = var_lons[ var_lons < 0 ] + 360

	list_lons = var_lons.ravel()
	list_lats = var_lats.ravel()
	points = np.column_stack([list_lats,list_lons])

	if( len(target_lon.shape) < 2):
		# create obs_lon and obs_lat
		# -179.5->179.5; -89.5->89.5, 1 deg res
		target_lon = np.arange(0.5,360.5,1)
		target_lat = np.arange(-89.5,90.5,1)

	# re-interpolate onto WOA grid, centred on greenwich meridian
	data_out = np.zeros((var.shape[0],target_lat.shape[0],target_lon.shape[0]),dtype=float) + np.nan

	for t in range(0,tDim):
		# preprocesss input data to remove mask
		vals = []
		maskVals = []
		for iy in range(0,var_lons.shape[0]):
			for ix in range(0,var_lons.shape[1]):
				maskVals.append(tmask[iy,ix])
				vals.append(var[t,iy,ix])

		vals = np.array(vals)

		valsFilt = []
		pointsFilt = []
		for p in range(0,points.shape[0]):
			if maskVals[p] == 1:
				valsFilt.append(vals[p])
				pointsFilt.append( ( points[p,0], points[p,1] ) )

		pointsFilt = np.array(pointsFilt)
		valsFilt = np.array(valsFilt)

		if pointsFilt.shape[0] > 0:
			grid_lon,grid_lat = np.meshgrid(target_lon,target_lat)
			data_out[t,:,:] = griddata(pointsFilt,valsFilt,(grid_lat,grid_lon), method='linear')

			data_out_near = griddata(pointsFilt,valsFilt,(grid_lat[:,0:2],grid_lon[:,0:2]), method='nearest') # take these values as has errors at date line due to itnerpoloation
			data_out[t,:,0:2] = data_out_near
			data_out_near = griddata(pointsFilt,valsFilt,(grid_lat[:,-3:],grid_lon[:,-3:]), method='nearest') # take these values as has errors at date line due to itnerpoloation
			data_out[t,:,-3:] = data_out_near[:,-3:]

			# northern boundary
			data_out_near = griddata(pointsFilt,valsFilt,(grid_lat[-3:,:],grid_lon[-3:,:]), method='nearest') # take these values as has errors at date line due to itnerpoloation
			data_out[t,-3:,:] = data_out_near[-3:,:]

	# tidy up
	data_out[ data_out > missingVal/1000. ] = missingVal

	return data_out

# ---------- 2 SUB DOMAIN DATA ----------
def subDomain(lonLim, latLim, in_data):
	lonStart = int(lonLim[0])+180
	lonEnd = int(lonLim[1])+180
	latStart = int(latLim[0])+90
	latEnd = int(latLim[1])+90
	subDiff = np.zeros(in_data.shape) * np.nan
	subDiff[:, latStart:latEnd, lonStart:lonEnd] = in_data[:, latStart:latEnd, lonStart:lonEnd]
	return subDiff

# ---------- 3 SUB DOMAIN ORCA DATA ----------
def subDomainORCA(lonLim, latLim, var_lons, var_lats, in_data, landMask, volMask, missingVal):

	lonStart = int(lonLim[0])
	lonEnd = int(lonLim[1])
	latStart = int(latLim[0])
	latEnd = int(latLim[1])

	if len(in_data.shape) == 3:
		mask = np.zeros(in_data[0,:,:].shape) + 1
	if len(in_data.shape) == 4:
		mask = np.zeros(in_data[0,0,:,:].shape) + 1

	mask[var_lons < lonStart] = missingVal
	mask[var_lons > lonEnd] = missingVal
	mask[var_lats < latStart] = missingVal
	mask[var_lats > latEnd] = missingVal
    
	ind_mask = np.where( mask == missingVal )
	ind_land = np.isnan(landMask)
	ind_vol = np.isnan(volMask)
	print(ind_vol.shape, in_data.shape)

	if len(in_data.shape) == 3:
		for t in range(0,in_data.shape[0]):
			temparr = in_data[t,:,:]
			temparr[ ind_mask ] = missingVal
			temparr[ ind_land ] = missingVal
			in_data[t,:,:] = temparr

	if len(in_data.shape) == 4:
		for t in range(0,in_data.shape[0]):
			for z in range(0,in_data.shape[1]):
				temparr = in_data[t,z,:,:]
				temparr[ ind_mask ] = missingVal
				in_data[t,z,:,:] = temparr

			temparr = in_data[t,:,:,:]
			if ind_vol.shape[0] == in_data.shape[1]:
				temparr[ ind_vol ] = missingVal
			in_data[t,:,:,:] = temparr

	return in_data

# ---------- 4 SURFACE DATA (sum at surface) ----------
def surfaceData(var, var_lons, var_lats, units, area, landMask, volMask, missingVal, lonLim, latLim):

	tDim = var.shape[0]

	var = subDomainORCA(lonLim, latLim, var_lons, var_lats, var, landMask, volMask, missingVal)

	# filter data for missingVal 
	varNan = np.copy(var)
	varNan[ varNan > missingVal/10. ] = np.nan
	varNan[ varNan < -missingVal/10. ] = np.nan

	total = 0
	monthly = []
	for t in range(0,tDim):
		total = total + np.nansum(varNan[t,:,:] * area[:,:] * units / tDim)
		# min = np.nanmin(varNan[t,:,:] * area[:,:] * units)
		# max = np.nanmax(varNan[t,:,:] * area[:,:] * units)
		min = np.nanpercentile(varNan[t,:,:] * units, 5)
		max = np.nanpercentile(varNan[t,:,:] * units, 95)
		median = np.nanmedian(varNan[t,:,:] * units)
		first = np.nanpercentile(varNan[t,:,:] * units, 25)
		third = np.nanpercentile(varNan[t,:,:] * units, 75)

		monthly.append([ np.nansum(varNan[t,:,:] * area[:,:] * units / tDim), min, first, median, third, max ])

	return total, monthly

# ---------- 5 VOLUME DATA (sum over all depths) ----------
def volumeData(var, var_lons, var_lats, units, vol, landMask, volMask, missingVal, lonLim, latLim):

	var = subDomainORCA(lonLim, latLim, var_lons, var_lats, var, landMask, volMask, missingVal)

	tDim = var.shape[0]

	# filter data for missingVal
	varNan = np.copy(var)
	varNan[ varNan > missingVal/10. ] = np.nan
	varNan[ varNan < -missingVal/10. ] = np.nan
	
	total = 0
	monthly = []
	for t in range(0,tDim):
		total = total + np.nansum(varNan[t,:,:,:] * vol[:,:,:] * units / tDim)

		min = np.nanpercentile(varNan[t,:,:,:] * units, 5)
		max = np.nanpercentile(varNan[t,:,:,:] * units, 95)
		median = np.nanmedian(varNan[t,:,:,:] * units)
		first = np.nanpercentile(varNan[t,:,:,:] * units, 25)
		third = np.nanpercentile(varNan[t,:,:,:] * units, 75)

		monthly.append([ np.nansum(varNan[t,:,:,:] * vol[:,:,:] * units / tDim), min, first, median, third, max ]) 

	return total, monthly

# ---------- 6 LEVEL DATA (sum at certain depth) ----------
def levelData(var, var_lons, var_lats, units, area, landMask, volMask, missingVal, lonLim, latLim, level):

	var = subDomainORCA(lonLim, latLim, var_lons, var_lats, var, landMask, volMask, missingVal)

	tDim = var.shape[0]

	# filter data for missingVal
	varNan = np.copy(var)
	varNan[ varNan > missingVal/10. ] = np.nan
	varNan[ varNan < -missingVal/10. ] = np.nan

	total = 0
	monthly = []
	for t in range(0,tDim):
		total = total + np.nansum(varNan[t,level,:,:] * area[:,:] * units / tDim)

		min = np.nanpercentile(varNan[t,level:,:] * units, 5)
		max = np.nanpercentile(varNan[t,level,:,:] * units, 95)
		median = np.nanmedian(varNan[t,level,:,:] * units)
		first = np.nanpercentile(varNan[t,level,:,:] * units, 25)
		third = np.nanpercentile(varNan[t,level,:,:] * units, 75)

		monthly.append([ np.nansum(varNan[t,level,:,:] * area[:,:] * units / tDim), min, first, median, third, max ]) 

	return total, monthly

# ---------- 7 INTEGRATE DATA (sum over range of depths) ----------
def intergrateData(var, var_lons, var_lats, depthFrom, depthTo, units, vol, landMask, volMask, missingVal, lonLim, latLim):

	var = subDomainORCA(lonLim, latLim, var_lons, var_lats, var, landMask, volMask, missingVal)

	tDim = var.shape[0]

	# filter data for missingVal
	varNan = np.copy(var)
	varNan[ varNan > missingVal/10. ] = np.nan
	varNan[ varNan < -missingVal/10. ] = np.nan

	total = 0
	monthly = []
	for t in range(0,tDim):
		total = total + np.nansum(varNan[t,depthFrom:depthTo+1,:,:] * vol[depthFrom:depthTo+1,:,:] * units / tDim)

		min = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 5)
		max = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 95)
		median = np.nanmedian(varNan[t,depthFrom:depthTo+1,:,:] * units)
		first = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 25)
		third = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 75)

		monthly.append([ np.nansum(varNan[t,depthFrom:depthTo+1,:,:] * vol[depthFrom:depthTo+1,:,:] * units / tDim), min, first, median, third, max ]) 

	return total, monthly

# ---------- 8 AVERAGE DATA (average over range of depths) ----------
def volumeDataAverage(var, var_lons, var_lats, depthFrom, depthTo, units, vol, landMask, volMask, missingVal, lonLim, latLim):

	var = subDomainORCA(lonLim, latLim, var_lons, var_lats, var, landMask, volMask, missingVal)

	tDim = var.shape[0]

	# filter data for missingVal
	varNan = np.copy(var)
	varNan[ varNan > missingVal/10. ] = np.nan
	varNan[ varNan < -missingVal/10. ] = np.nan

	if vol.shape[0] == varNan.shape[1]:
		vol_masked = np.copy(vol)
	if varNan.shape[1] == 1:
		vol_masked = np.reshape( np.copy( vol[0,:,:] ), (1,varNan.shape[2], varNan.shape[3]) )
	ind_masked = np.isnan(varNan[0,:,:,:])
    
	vol_masked[ ind_masked ] = np.nan
	volNorm = np.nansum(vol_masked[depthFrom:depthTo+1,:,:])

	total = 0
	monthly = []
	for t in range(0,tDim):
		total = total + np.nansum(varNan[t,depthFrom:depthTo+1,:,:] * vol_masked[depthFrom:depthTo+1,:,:]/volNorm * units / tDim)

		# the following are not properly scaled
		min = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 5)
		max = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 95)
		# median = np.nansum(varNan[t,depthFrom:depthTo+1,:,:] * vol_masked[depthFrom:depthTo+1,:,:]/volNorm * units )
		median = np.nanmedian(varNan[t,depthFrom:depthTo+1,:,:] * units)
		first = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 25)
		third = np.nanpercentile(varNan[t,depthFrom:depthTo+1,:,:] * units, 75)

		# monthly.append([ np.nansum(varNan[t,depthFrom:depthTo+1,:,:] * vol_masked[depthFrom:depthTo+1,:,:]/volNorm * units / tDim), min, first, median, third, max ]) 
		monthly.append([ np.nansum(varNan[t,depthFrom:depthTo+1,:,:] * vol_masked[depthFrom:depthTo+1,:,:]/volNorm * units), min, first, median, third, max ]) 

	return total, monthly

# ---------- 9 OBS DATA (RMSE to class fo observations) ----------
def unit_vector(vector):
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	
	return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def observationData(obs, obs_lon, obs_lat, var, var_lons, var_lats, obsCenter, tmask, missingVal, gammaFlag, lonLim, latLim):

	# wrap in long if obs data centered differently (target is greenwich centering of array)
	if obsCenter == 'zero':
		obs = np.roll(obs,int(len(obs_lon)/2),axis=2)
		obs_lon = np.roll(obs_lon,int(len(obs_lon)/2),axis=0)

		for i in range(0,len(obs_lon)):
			if obs_lon[i] > 180:
				obs_lon[i] = obs_lon[i] - 360

	data_out = regrid(var, var_lons, var_lats, obs_lon, obs_lat, tmask, missingVal)

	obs[ obs > missingVal/1000. ] = missingVal

	# spread obs to account for DTA
	spread = 2 # degree std
    
	modified_weights = np.copy(obs)
	modified_obs = np.copy(obs)

	for t in range(0,12):
		slice = np.copy(obs[t,:,:])
		weight = np.copy(slice)
		slice[ slice <= -missingVal ] = 0
		slice[ slice != 0 ] = 1
		slice_spead = gaussian_filter(slice, sigma=spread)
		ind = np.where( slice_spead < 0.1)
		slice_spead[ ind ] = 0
		grid_lon,grid_lat = np.meshgrid(obs_lon,obs_lat)

		if np.max(slice) == 1: # not empty
			# source data
			ind_obs = np.where( slice == 1)
			list_lons = grid_lon[ ind_obs ].ravel()
			list_lats = grid_lat[ ind_obs ].ravel()
			points = np.column_stack([list_lats,list_lons])
			vals = []

			for p in range(0, len(ind_obs[0])):
				vals.append(obs[t, ind_obs[0][p], ind_obs[1][p] ])

			vals = np.array(vals)
			filled_obs = griddata(points,vals,(grid_lat,grid_lon), method='nearest')
			filled_obs[ ind ] = np.nan
			modified_weights[t,:,:] = slice_spead
			modified_obs[t,:,:] = filled_obs
		else:
			modified_weights[t,:,:] = np.zeros( slice.shape )
			modified_obs[t,:,:] = np.zeros( slice.shape )

	# Now, compare to obs data as its now on the same grid
	data_out[ data_out > missingVal/1000. ] = np.nan

	diff = (data_out-obs) * (data_out-obs)
	diff_mod = (data_out-modified_obs) * (data_out-modified_obs)
	diff[ diff < -missingVal/1000. ] = np.nan
	diff_mod[ diff_mod < -missingVal/1000. ] = np.nan

	diff[ diff >  1E6 ] = np.nan
	diff_mod[ diff_mod >  1E6 ] = np.nan
    
	# calc mean and std differences for filtering outliers
	mean = np.nanmean(diff)
	std = np.nanstd(diff)
	z_dat = (diff-mean)/std # Z statistic for outliers 
	diff[ z_dat > 10 ] = np.nan
	diff[ z_dat < -10 ] = np.nan
	mean = np.nanmean(diff_mod)
	std = np.nanstd(diff_mod)
	z_dat = (diff-mean)/std 
	diff_mod[ z_dat > 10 ] = np.nan
	diff_mod[ z_dat < -10 ] = np.nan

	# filter the diff by the limits to long and lat
	diff = subDomain(lonLim, latLim, diff)
	diff_mod = subDomain(lonLim, latLim, diff_mod)

	ind = np.isnan(diff_mod)
	modified_weights[ ind ] = 0

	# annual rmse
	total = np.nansum(diff[:,:,:])
	count = np.count_nonzero(~np.isnan(diff[:,:,:]))
	rmse_annual = np.sqrt( (total/count) )

	total_mod = np.nansum( diff_mod[:,:,:] * modified_weights[:,:,:] )
	count_mod = np.count_nonzero( modified_weights[:,:,:] )
	rmse_annual_mod = np.sqrt( (total_mod/count_mod) )

	# repeat on monthly basis using same z-dat criteria for outliers
	monthly_rms = []
	for t in range(0,12):
		total = np.nansum(diff[t,:,:])
		count = np.count_nonzero(~np.isnan(diff[t,:,:]))
		rmse = np.sqrt( (total/count) )

		if np.max(modified_weights[t,:,:]) != np.min(modified_weights[t,:,:]) :
			total_mod = np.nansum(diff_mod[t,:,:] * modified_weights[t,:,:])
			count_mod = np.count_nonzero(modified_weights[t,:,:])
			rmse_mod = np.sqrt( (total_mod/count_mod) )

			monthly_rms.append(rmse_mod)
		else:
			monthly_rms.append(-1)

	gammaResults = [[],[]]
	if gammaFlag == True:
		# create an array of vectors for each point on the surface
		vectors = np.zeros((3,180,360))
		latlon = np.zeros((2,180,360))
		for ix in range(0,360):
			for iy in range(0,180):
				vectors[0,iy,ix] = math.cos(np.radians(obs_lon[ix]))*math.cos(np.radians(obs_lat[iy]))
				vectors[1,iy,ix] = math.sin(np.radians(obs_lon[ix]))*math.cos(np.radians(obs_lat[iy]))
				vectors[2,iy,ix] = math.sin(np.radians(obs_lat[iy]))
				latlon[0,iy,ix] = obs_lat[iy]
				latlon[1,iy,ix] = obs_lon[ix]

		a_vals = [0.1,1,2,3,5,10] # this sampling of the relationship
		b_vals = [ [],[],[],[],[],[] ] #create empty lists to average
        
		d1 = np.absolute(data_out - obs)

		for ix in range(0,obs.shape[2]):
			for iy in range(0,obs.shape[1]):
				for t in range(0,obs.shape[0]):
					if obs[t,iy,ix] > -999 and obs[t,iy,ix] < missingVal/10 and data_out[t,iy,ix] < missingVal/10:
						# create dot product array
						angs = np.zeros((180,360))
						vect = vectors[:,iy,ix]

						for x in range(0,360):
							for y in range(0,180):
								angs[y,x] = angle_between(vect,vectors[:,y,x])

						for ai in range(0,len(a_vals)):
							a = a_vals[ai]

							prev_r = False

							for r in range(1,10):
								if prev_r == False:
									sub = np.copy(angs)
									rad = np.radians(r)
									sub[ sub > rad ] = missingVal
									sub_dist = np.copy(sub)
									sub[ sub < rad ] = 1

									d = np.copy(data_out[t,:,:])
									d = d * sub

									sub_diff = np.abs(d-obs[t,iy,ix])/obs[t,iy,ix] * 100 # % change
									min_diff = np.amin(sub_diff)

									c1 = sub_diff / a
									c2 = np.degrees(sub_dist) / r 
									gamma = np.sqrt( c1*c1 + c2*c2 )
									min_gamma = np.amin(gamma)
									if min_gamma < 1:
										b_vals[ai].append(r)
										print(ix,iy,obs[t,iy,ix], data_out[t,iy,ix],a,r,min_gamma)
										prev_r = True

		# output average b-vals for each a-val
		for ai in range(0,len(a_vals)):
			mean_b = np.mean(b_vals[ai])

			gammaResults[0].append(a_vals[ai])
			gammaResults[1].append(mean_b)

	return rmse_annual_mod, monthly_rms #, gammaResults

# ---------- 10 BLOOM FUNCTION ----------
def interpolate(x, x1, x2, y1, y2):
	return y1 + (x-x1)*(y2-y1)/(x2-x1)

def bloom(var, var_lons, var_lats, missingVal, lonLim, latLim):

	data_out = regrid(var, var_lons, var_lats, np.array([0]), np.array([0]), missingVal)

	bloomVal=np.zeros((5,data_out.shape[1],data_out.shape[2]))

	# fit a time series for each point
	for x in range(0,data_out.shape[2]):
		for y in range(0,data_out.shape[1]):
			# get time series 
			if np.sum(data_out[0,y,x]) < missingVal/10 :
				timeSeries = data_out[:,y,x]

				# extend for interpolation
				timeSeries = np.concatenate((timeSeries,timeSeries,timeSeries), axis=0)

				# set threshold and initialise indecis
				thresh = 1.05 * np.median(timeSeries)
				max_t = np.argmax(timeSeries) + 12
				t_low = max_t  
				t_high = max_t
				meanAtPoint = np.mean(timeSeries)
                
				# find start/end points
				while timeSeries[t_low] > thresh:
					t_low = t_low - 1
				t_low = interpolate(thresh, timeSeries[t_low], timeSeries[t_low+1], t_low, t_low+1)
				while timeSeries[t_high] > thresh:
					t_high = t_high + 1
				t_high = interpolate(thresh, timeSeries[t_high-1], timeSeries[t_high], t_high-1, t_high)

				# get the peak bloom value
				maxVal = timeSeries[max_t]

				# set the start/end points to month
				while t_low > 11:
					t_low = t_low - 12
				while t_high > 11:
					t_high = t_high - 12
				while max_t > 11:
					max_t = max_t - 12
				if t_low < 0:
					t_low = t_low + 12

				# get duration
				duration = 0
				if t_high > t_low:
					duration = t_high - t_low
				else:
					duration = (t_high +12) - t_low
                
				bloomVal[0,y,x] = maxVal/meanAtPoint
				bloomVal[1,y,x] = max_t
				bloomVal[2,y,x] = t_low
				bloomVal[3,y,x] = t_high
				bloomVal[4,y,x] = duration

	# filter the diff by the limits to long and lat
	bloomVal = subDomain(lonLim, latLim, bloomVal)

	# scale according to area and give stats
	# max peak, median peak, avg duration, s.d. duration

	maxPeak = np.nanmax(bloomVal[0,:,:])
	peaks = bloomVal[0,:,:]
	peaks[ peaks == 0 ] = np.nan
	medPeak = np.nanmedian(peaks)

	duration = bloomVal[4,:,:]
	duration[ duration == 0 ] = np.nan
	meanDur = np.nanmean(duration)
	sdDur = np.nanstd(duration)

	# add column headings we'll need for output (last one is monthly output holder)
	return maxPeak, medPeak, meanDur, sdDur, ['maxPeak','medPeak', 'meanDur', 'sdDur'], None

def getSlope(xs,ys):
	mx = np.mean(xs)
	my = np.mean(ys)
	mxy = np.mean(xs*ys)
	mxx = np.mean(xs*xs)
	m = ( mx*my - mxy) / ( mx*mx - mxx ) 
	c = my - m * mx

	return m, c

def trophic(var, var_lons, var_lats, missingVal, lonLim, latLim):

	trophicVals=np.zeros((4, 12, 180, 360 ))
 
	for p in range(0,3):
		data_out = regrid(var[p], var_lons, var_lats, np.array([0]), np.array([0]), missingVal)
		trophicVals[p,:,:,:] = data_out

	y_vals = [ 0, 1, 2 ]
	for x in range(0,trophicVals[0,:,:,:].shape[2]):
		for y in range(0,trophicVals[0,:,:,:].shape[1]):
			for t in range(0,trophicVals[0,:,:,:].shape[0]):
				x_vals = trophicVals[0:3,t,y,x]
				# fit
				if np.min(x_vals) != 0:
					if ~np.isnan(x_vals).any():
						m, c = getSlope(x_vals,y_vals)
						trophicVals[3,t,y,x] = m

	trophicVals[ trophicVals > missingVal/10. ] = np.nan

	trophicVals = subDomain(lonLim, latLim, trophicVals)

	mean0 = np.nanmean(trophicVals[0,:,:,:])
	mean1 = np.nanmean(trophicVals[1,:,:,:])
	mean2 = np.nanmean(trophicVals[2,:,:,:])
	mean_m = np.nanmean(trophicVals[3,:,:,:])
	mean_m_monthly = np.zeros((12))

	for t in range(0,12):
		mean_m_monthly[t]= np.nanmean(trophicVals[3,t,:,:])

	return mean0, mean1, mean2, mean_m, ['level1','level2', 'level3', 'meanSlope'], mean_m_monthly
