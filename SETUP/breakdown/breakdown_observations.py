import datetime

def observationDatasets():

	ListOfObservations = [
		{   
			'name' : 'SOCAT',
			'path' : '/gpfs/data/greenocean/software/resources/Observations/SOCATv2020_tracks_gridded_monthly.nc', 
			'lonVar': 'xlon', 'latVar': 'ylat', 'timeVar' : 'tmnth',
			'origin': datetime.datetime(1970,1,1,0,0,0),
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': False
		},

		{   
			'name' : 'GLODAP_talk',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/GLODAPv2.2020_Merged_Master_File_talk_avg.nc', 
			'lonVar': 'x', 'latVar': 'y', 'timeVar' : 'time_counter',
			'origin': datetime.datetime(1900,1,1,0,0,0),
			'conversion': '/gpfs/home/yzh17dvu/scratch/RunBreakdown/GLODAPv2.2020_density_avg.nc',
			'conversionName': 'DENSITY',
			'factor': 1E6,
			'centered': 'zero',
			'climatological': False
		},

		{   
			'name' : 'GLODAP_tco2',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/GLODAPv2.2020_Merged_Master_File_tco2_avg.nc', 
			'lonVar': 'x', 'latVar': 'y', 'timeVar' : 'time_counter',
			'origin': datetime.datetime(1900,1,1,0,0,0),
			'conversion': '/gpfs/home/yzh17dvu/scratch/RunBreakdown/GLODAPv2.2020_density_avg.nc',
			'conversionName': 'DENSITY',
			'factor': 1E6,
			'centered': 'zero',
			'climatological': False
		},

		{   
			'name' : 'GLODAP_fco2',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/GLODAPv2.2020_Merged_Master_File_fco2_avg.nc', 
			'lonVar': 'x', 'latVar': 'y', 'timeVar' : 'time_counter',
			'origin': datetime.datetime(1900,1,1,0,0,0),
			'conversion': '/gpfs/home/yzh17dvu/scratch/RunBreakdown/GLODAPv2.2020_density_avg.nc',
			'conversionName': 'DENSITY',
			'factor': 1E6,
			'centered': 'zero',
			'climatological': False
		},

		{   
			'name' : 'occci',
			'path' : '/gpfs/data/greenocean/software/resources/Observations/ESACCI-OC-L3S-CHLOR_A-MERGED-MONTHLY_1deg_GEO_OC4v6-199709-201207-fv1.0.nc', 
			'lonVar': 'LONGITUDE', 'latVar': 'LATITUDE', 'timeVar' : 'TIME',
			'origin': datetime.datetime(1970,1,1,0,0,0),
			'conversion': None,
			'conversionName': None,
			'factor': 1E6,
			'centered': 'greenwich',
			'climatological': False
		},

		{   
			'name' : 'cfc',   # need to test this as model outputs don't look complete
			'path' : './CFC11_monthly_171108.nc', 
			'lonVar': 'LONGITUDE', 'latVar': 'LATITUDE', 'timeVar' : 'TIME',
			'origin': datetime.datetime(1981,1,1,0,0,0),
			'conversion': None,
			'conversionName': None,
			'factor': 1E6,
			'centered': 'greenwich',
			'climatological': False
		},

		{
			'name' : 'WOA_T_1981_2010',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_1981_2010_ave_VOTEMP_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		},

		{
			'name' : 'WOA_T_2005_2017',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_2005_2017_ave_VOTEMP_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		},

		{
			'name' : 'WOA_T_1975_1984',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_1975_1984_ave_VOTEMP_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		},

		{
			'name' : 'WOA_T_1965_1974',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_1965_1974_ave_VOTEMP_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		},

		{
			'name' : 'WOA_T_1955_1964',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_1955_1964_ave_VOTEMP_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		},

		{
			'name' : 'WOA_S',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_1981_2010_ave_VOSAL_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		},
                        
		{
			'name' : 'WOA_MLD',
			'path' : '/gpfs/home/yzh17dvu/scratch/RunBreakdown/WOA/WOA_1981_2010_ave_MLD_v1.nc', 
			'lonVar': 'lon', 'latVar': 'lat', 'timeVar' : 'time',
			'origin': datetime.datetime(1900,1,1,0,0,0),  # ignored as climatological
			'conversion': None,
			'conversionName': None,
			'factor': 1,
			'centered': 'greenwich',
			'climatological': True
		}
	]

	return ListOfObservations
