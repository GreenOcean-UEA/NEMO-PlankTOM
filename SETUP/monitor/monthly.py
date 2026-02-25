#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pathlib
import sys
import calendar

# Script for creating monthly plots of variables using breakdowns created during the model run process
# Will only be run at the end of the model run process

tModel = sys.argv[1]
baseDir = sys.argv[2]

def makeSummaryFromBreakdowns(tmod, modBaseDir):
    
	saveDir = f'/{modBaseDir}/monitor/{tmod}/'
	pathlib.Path(saveDir).mkdir(parents=True, exist_ok=True)
    
	try:
		### Surface monthly
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.sur.monthly.dat', sep='\t')
		print('Read surface monthly')
        
		mnth = w.month[-12:].to_numpy().astype(int)
		Cflx_total = w.Cflx[-12:].to_numpy().astype(float)
		Cflx_reg1 = w['Cflx.1'][-12:].to_numpy().astype(float)
		Cflx_reg2 = w['Cflx.2'][-12:].to_numpy().astype(float)
		Cflx_reg3 = w['Cflx.3'][-12:].to_numpy().astype(float)
		Cflx_reg4 = w['Cflx.4'][-12:].to_numpy().astype(float)
		Cflx_reg5 = w['Cflx.5'][-12:].to_numpy().astype(float)

		### Average monthly
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.ave.monthly.dat', sep='\t')
		print('Read average monthly')
        
		mnth = w.month[-12:].to_numpy().astype(int)
		TChl_total = w.TChl[-12:].to_numpy().astype(float)
		TChl_reg1 = w['TChl.1'][-12:].to_numpy().astype(float)
		TChl_reg2 = w['TChl.2'][-12:].to_numpy().astype(float)
		TChl_reg3 = w['TChl.3'][-12:].to_numpy().astype(float)
		TChl_reg4 = w['TChl.4'][-12:].to_numpy().astype(float)
		TChl_reg5 = w['TChl.5'][-12:].to_numpy().astype(float)

		pCO2_total = w.pCO2[-12:].to_numpy().astype(float)
		pCO2_reg1 = w['pCO2.1'][-12:].to_numpy().astype(float)
		pCO2_reg2 = w['pCO2.2'][-12:].to_numpy().astype(float)
		pCO2_reg3 = w['pCO2.3'][-12:].to_numpy().astype(float)
		pCO2_reg4 = w['pCO2.4'][-12:].to_numpy().astype(float)
		pCO2_reg5 = w['pCO2.5'][-12:].to_numpy().astype(float)

		SST_total = w.tos[-12:].to_numpy().astype(float)
		SST_reg1 = w['tos.1'][-12:].to_numpy().astype(float)
		SST_reg2 = w['tos.2'][-12:].to_numpy().astype(float)
		SST_reg3 = w['tos.3'][-12:].to_numpy().astype(float)
		SST_reg4 = w['tos.4'][-12:].to_numpy().astype(float)
		SST_reg5 = w['tos.5'][-12:].to_numpy().astype(float)

		SSS_total = w.sos[-12:].to_numpy().astype(float)
		SSS_reg1 = w['sos.1'][-12:].to_numpy().astype(float)
		SSS_reg2 = w['sos.2'][-12:].to_numpy().astype(float)
		SSS_reg3 = w['sos.3'][-12:].to_numpy().astype(float)
		SSS_reg4 = w['sos.4'][-12:].to_numpy().astype(float)
		SSS_reg5 = w['sos.5'][-12:].to_numpy().astype(float)

		MLD_total = w.mldr10_1[-12:].to_numpy().astype(float)
		MLD_reg1 = w['mldr10_1.1'][-12:].to_numpy().astype(float)
		MLD_reg2 = w['mldr10_1.2'][-12:].to_numpy().astype(float)
		MLD_reg3 = w['mldr10_1.3'][-12:].to_numpy().astype(float)
		MLD_reg4 = w['mldr10_1.4'][-12:].to_numpy().astype(float)
		MLD_reg5 = w['mldr10_1.5'][-12:].to_numpy().astype(float)

		mnth_name = []
		for i in range (0,12):
			mnth_name.append(calendar.month_abbr[mnth[i]+1])

		# GCB 2022 data-products ensemble mean: average monthly over the period 2010-2019
		data_glob = np.array([374.5267, 376.9849, 378.6273, 377.6980, 374.4194, 372.0030, 372.8390, 373.1397, 373.7141, 374.7667, 375.3111, 376.1471])
		data_reg1 = np.array([367.6449, 374.6633, 376.1620, 367.9230, 348.2934, 327.9850, 319.9174, 315.5227, 320.6641, 336.6024, 353.5276, 366.8002])
		data_reg2 = np.array([360.3107, 359.7239, 360.9054, 364.4423, 372.4738, 384.9409, 398.0022, 403.2199, 398.5489, 386.5432, 374.1430, 366.2098])
		data_reg3 = np.array([399.8266, 401.4656, 403.5920, 404.5923, 404.1028, 402.5178, 401.7164, 401.2903, 401.1736, 401.0325, 400.8897, 401.5201])
		data_reg4 = np.array([383.5775, 385.4783, 381.7295, 373.6719, 366.9351, 362.8020, 361.2315, 360.9230, 361.6792, 363.8154, 368.7272, 377.6358])
		data_reg5 = np.array([360.1319, 361.5651, 368.3551, 376.4319, 381.8729, 387.2647, 391.9770, 394.6131, 394.8327, 390.6098, 380.8251, 368.3117])
		
		### Creating monthly region figures for CFLX
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(5*ratio)
		fig.set_figwidth(8*ratio)

		ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), rowspan=3, colspan=5)
		ax2 = plt.subplot2grid(shape=(5,8), loc=(3,0), rowspan=2, colspan=5)

		ax3 = plt.subplot2grid(shape=(5,8), loc=(0,5), colspan=3)
		ax4 = plt.subplot2grid(shape=(5,8), loc=(1,5), colspan=3, sharey=ax3)
		ax5 = plt.subplot2grid(shape=(5,8), loc=(2,5), colspan=3, sharey=ax3)
		ax6 = plt.subplot2grid(shape=(5,8), loc=(3,5), colspan=3, sharey=ax3)
		ax7 = plt.subplot2grid(shape=(5,8), loc=(4,5), colspan=3, sharey=ax3)

		for ax in [ax1, ax3, ax4, ax5, ax6, ax7]:
			ax.axhline(color='black', linewidth=0.5)

		l1, = ax1.plot(mnth_name, Cflx_total, "b", label='global total'); ax1.set_title('Surface Cflx (global) [PgC/yr]'); ax1.grid(linestyle='--', linewidth=0.25)
		ax1.plot(mnth_name, Cflx_reg1, "r", linewidth=0.5);
		ax1.plot(mnth_name, Cflx_reg2, "g", linewidth=0.5);
		ax1.plot(mnth_name, Cflx_reg3, "c", linewidth=0.5);
		ax1.plot(mnth_name, Cflx_reg4, "m", linewidth=0.5);
		ax1.plot(mnth_name, Cflx_reg5, "y", linewidth=0.5);

		ax2.plot(mnth_name, Cflx_total, "b", label='global total'); 
	
		l3, = ax3.plot(mnth_name, Cflx_reg1, "r", label='45N-90N'); ax3.set_title('Surface Cflx (regional) [PgC/yr]');
		l4, = ax4.plot(mnth_name, Cflx_reg2, "g", label='15N-45N'); 
		l5, = ax5.plot(mnth_name, Cflx_reg3, "c", label='15S-15N'); 
		l6, = ax6.plot(mnth_name, Cflx_reg4, "m", label='45S-15S'); 
		l7, = ax7.plot(mnth_name, Cflx_reg5, "y", label='90S-45S');

		lines = [l1, l3, l4, l5, l6, l7]
		labels = ['global total', '45N-90N', '15N-45N', '15S-15N', '45S-15S', '90S-45S']

		ax1.legend(lines, labels, loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
			ax.legend(loc='upper right')
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_monthly_Cflx.jpg')
		print(f'Created CFLX monthly summary figure for {tmod}')

		### Creating monthly region figures for TCHL
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(5*ratio)
		fig.set_figwidth(8*ratio)

		ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), rowspan=3, colspan=5)
		ax2 = plt.subplot2grid(shape=(5,8), loc=(3,0), rowspan=2, colspan=5)

		ax3 = plt.subplot2grid(shape=(5,8), loc=(0,5), colspan=3)
		ax4 = plt.subplot2grid(shape=(5,8), loc=(1,5), colspan=3, sharey=ax3)
		ax5 = plt.subplot2grid(shape=(5,8), loc=(2,5), colspan=3, sharey=ax3)
		ax6 = plt.subplot2grid(shape=(5,8), loc=(3,5), colspan=3, sharey=ax3)
		ax7 = plt.subplot2grid(shape=(5,8), loc=(4,5), colspan=3, sharey=ax3)

		for ax in [ax1, ax3, ax4, ax5, ax6, ax7]:
			ax.axhline(linestyle='--', linewidth=0.25)

		l1, = ax1.plot(mnth_name, TChl_total, "b", label='global average'); ax1.set_title('Average TChl over top 100m (global) [ug Chl/L]'); ax1.grid(linestyle='--', linewidth=0.25)
		ax1.plot(mnth_name, TChl_reg1, "r", linewidth=0.5);
		ax1.plot(mnth_name, TChl_reg2, "g", linewidth=0.5);
		ax1.plot(mnth_name, TChl_reg3, "c", linewidth=0.5);
		ax1.plot(mnth_name, TChl_reg4, "m", linewidth=0.5);
		ax1.plot(mnth_name, TChl_reg5, "y", linewidth=0.5);

		ax2.plot(mnth_name, TChl_total, "b", label='global average'); 
	
		l3, = ax3.plot(mnth_name, TChl_reg1, "r", label='45N-90N'); ax3.set_title('Average TChl (regional) [ug Chl/L]');
		l4, = ax4.plot(mnth_name, TChl_reg2, "g", label='15N-45N'); 
		l5, = ax5.plot(mnth_name, TChl_reg3, "c", label='15S-15N'); 
		l6, = ax6.plot(mnth_name, TChl_reg4, "m", label='45S-15S');
		l7, = ax7.plot(mnth_name, TChl_reg5, "y", label='90S-45S');

		lines = [l1, l3, l4, l5, l6, l7]
		labels = ['global average', '45N-90N', '15N-45N', '15S-15N', '45S-15S', '90S-45S']

		ax1.legend(lines, labels, loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
			ax.legend(loc='upper right')
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_monthly_TChl.jpg')
		print(f'Created TCHL monthly summary figure for {tmod}')

		### Creating monthly region figures for pCO2
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(5*ratio)
		fig.set_figwidth(8*ratio)

		ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), rowspan=3, colspan=5)
		ax2 = plt.subplot2grid(shape=(5,8), loc=(3,0), rowspan=2, colspan=5)

		ax3 = plt.subplot2grid(shape=(5,8), loc=(0,5), colspan=3)
		ax4 = plt.subplot2grid(shape=(5,8), loc=(1,5), colspan=3, sharey=ax3)
		ax5 = plt.subplot2grid(shape=(5,8), loc=(2,5), colspan=3, sharey=ax3)
		ax6 = plt.subplot2grid(shape=(5,8), loc=(3,5), colspan=3, sharey=ax3)
		ax7 = plt.subplot2grid(shape=(5,8), loc=(4,5), colspan=3, sharey=ax3)

		l1, = ax1.plot(mnth_name, pCO2_total, "b", label='global average'); ax1.set_title('Average Surface pCO2 (global) [ppm]'); ax1.grid(linestyle='--', linewidth=0.25)
		ax1.plot(mnth_name, pCO2_reg1, "r", linewidth=0.5);
		ax1.plot(mnth_name, pCO2_reg2, "g", linewidth=0.5);
		ax1.plot(mnth_name, pCO2_reg3, "c", linewidth=0.5);
		ax1.plot(mnth_name, pCO2_reg4, "m", linewidth=0.5);
		ax1.plot(mnth_name, pCO2_reg5, "y", linewidth=0.5);

		ax2.plot(mnth_name, pCO2_total, "b", label='global average'); ax2.plot(mnth_name, data_glob, "k", label='data-products (2010-2019 ave)');

		l3, = ax3.plot(mnth_name, pCO2_reg1, "r", label='45N-90N'); ax3.set_title('Average Surface pCO2 (regional) [ppm]'); ax3.plot(mnth_name, data_reg1, "k", label='data-products');
		l4, = ax4.plot(mnth_name, pCO2_reg2, "g", label='15N-45N'); ax4.plot(mnth_name, data_reg2, "k", label='data-products');
		l5, = ax5.plot(mnth_name, pCO2_reg3, "c", label='15S-15N'); ax5.plot(mnth_name, data_reg3, "k", label='data-products');
		l6, = ax6.plot(mnth_name, pCO2_reg4, "m", label='45S-15S'); ax6.plot(mnth_name, data_reg4, "k", label='data-products');
		l7, = ax7.plot(mnth_name, pCO2_reg5, "y", label='90S-45S'); ax7.plot(mnth_name, data_reg5, "k", label='data-products');

		lines = [l1, l3, l4, l5, l6, l7]
		labels = ['global average', '45N-90N', '15N-45N', '15S-15N', '45S-15S', '90S-45S']

		ax1.legend(lines, labels, loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
			ax.legend(loc='upper right')
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_monthly_pCO2.jpg')
		print(f'Created pCO2 monthly summary figure for {tmod}')
		
		### Creating monthly region figures for pCO2 (normalised)
		npCO2_total = pCO2_total - pCO2_total[0]
		npCO2_reg1 = pCO2_reg1 - pCO2_reg1[0]
		npCO2_reg2 = pCO2_reg2 - pCO2_reg2[0]
		npCO2_reg3 = pCO2_reg3 - pCO2_reg3[0]
		npCO2_reg4 = pCO2_reg4 - pCO2_reg4[0]
		npCO2_reg5 = pCO2_reg5 - pCO2_reg5[0]

		ndata_glob = data_glob - data_glob[0]
		ndata_reg1 = data_reg1 - data_reg1[0]
		ndata_reg2 = data_reg2 - data_reg2[0]
		ndata_reg3 = data_reg3 - data_reg3[0]
		ndata_reg4 = data_reg4 - data_reg4[0]
		ndata_reg5 = data_reg5 - data_reg5[0]

		ratio = 2
		fig = plt.figure()
		fig.set_figheight(5*ratio)
		fig.set_figwidth(8*ratio)

		ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), rowspan=3, colspan=5)
		ax2 = plt.subplot2grid(shape=(5,8), loc=(3,0), rowspan=2, colspan=5)

		ax3 = plt.subplot2grid(shape=(5,8), loc=(0,5), colspan=3)
		ax4 = plt.subplot2grid(shape=(5,8), loc=(1,5), colspan=3, sharey=ax3)
		ax5 = plt.subplot2grid(shape=(5,8), loc=(2,5), colspan=3, sharey=ax3)
		ax6 = plt.subplot2grid(shape=(5,8), loc=(3,5), colspan=3, sharey=ax3)
		ax7 = plt.subplot2grid(shape=(5,8), loc=(4,5), colspan=3, sharey=ax3)

		l1, = ax1.plot(mnth_name, npCO2_total, "b", label='global average'); ax1.set_title('Average Surface pCO2 (global) [ppm]'); ax1.grid(linestyle='--', linewidth=0.25)
		ax1.plot(mnth_name, npCO2_reg1, "r", linewidth=0.5);
		ax1.plot(mnth_name, npCO2_reg2, "g", linewidth=0.5);
		ax1.plot(mnth_name, npCO2_reg3, "c", linewidth=0.5);
		ax1.plot(mnth_name, npCO2_reg4, "m", linewidth=0.5);
		ax1.plot(mnth_name, npCO2_reg5, "y", linewidth=0.5);

		ax2.plot(mnth_name, npCO2_total, "b", label='global average'); ax2.plot(mnth_name, ndata_glob, "k", label='data-products (2010-2019 ave)');

		l3, = ax3.plot(mnth_name, npCO2_reg1, "r", label='45N-90N'); ax3.set_title('Average Surface pCO2 (regional) [ppm]'); ax3.plot(mnth_name, ndata_reg1, "k", label='data-products');
		l4, = ax4.plot(mnth_name, npCO2_reg2, "g", label='15N-45N'); ax4.plot(mnth_name, ndata_reg2, "k", label='data-products');
		l5, = ax5.plot(mnth_name, npCO2_reg3, "c", label='15S-15N'); ax5.plot(mnth_name, ndata_reg3, "k", label='data-products');
		l6, = ax6.plot(mnth_name, npCO2_reg4, "m", label='45S-15S'); ax6.plot(mnth_name, ndata_reg4, "k", label='data-products');
		l7, = ax7.plot(mnth_name, npCO2_reg5, "y", label='90S-45S'); ax7.plot(mnth_name, ndata_reg5, "k", label='data-products');

		lines = [l1, l3, l4, l5, l6, l7]
		labels = ['global average', '45N-90N', '15N-45N', '15S-15N', '45S-15S', '90S-45S']

		ax1.legend(lines, labels, loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
			ax.legend(loc='upper right')
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_monthly_pCO2_normalised.jpg')
		print(f'Created normalised pCO2 monthly summary figure for {tmod}')

		### Creating monthly region figures for SST
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(5*ratio)
		fig.set_figwidth(8*ratio)

		ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), rowspan=3, colspan=5)
		ax2 = plt.subplot2grid(shape=(5,8), loc=(3,0), rowspan=2, colspan=5)

		ax3 = plt.subplot2grid(shape=(5,8), loc=(0,5), colspan=3)
		ax4 = plt.subplot2grid(shape=(5,8), loc=(1,5), colspan=3)
		ax5 = plt.subplot2grid(shape=(5,8), loc=(2,5), colspan=3)
		ax6 = plt.subplot2grid(shape=(5,8), loc=(3,5), colspan=3)
		ax7 = plt.subplot2grid(shape=(5,8), loc=(4,5), colspan=3)

		ax1.axhline(linestyle='--', linewidth=0.25)

		l1, = ax1.plot(mnth_name, SST_total, "b"); ax1.set_title('Average Sea Surface Temperature (global) [degC]'); ax1.grid(linestyle='--', linewidth=0.25)
		ax1.plot(mnth_name, SST_reg1, "r", linewidth=0.5);
		ax1.plot(mnth_name, SST_reg2, "g", linewidth=0.5);
		ax1.plot(mnth_name, SST_reg3, "c", linewidth=0.5);
		ax1.plot(mnth_name, SST_reg4, "m", linewidth=0.5);
		ax1.plot(mnth_name, SST_reg5, "y", linewidth=0.5);

		ax2.plot(mnth_name, SST_total, "b", label='global average'); 
	
		l3, = ax3.plot(mnth_name, SST_reg1, "r", label='45N-90N'); ax3.set_title('Average Sea Surface Temperature (regional) [degC]');
		l4, = ax4.plot(mnth_name, SST_reg2, "g", label='15N-45N'); 
		l5, = ax5.plot(mnth_name, SST_reg3, "c", label='15S-15N'); 
		l6, = ax6.plot(mnth_name, SST_reg4, "m", label='45S-15S');
		l7, = ax7.plot(mnth_name, SST_reg5, "y", label='90S-45S');

		lines = [l1, l3, l4, l5, l6, l7]
		labels = ['global average', '45N-90N', '15N-45N', '15S-15N', '45S-15S', '90S-45S']

		ax1.legend(lines, labels, loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
			ax.legend(loc='upper right')
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_monthly_SST.jpg')
		print(f'Created SST monthly summary figure for {tmod}')

		### Creating monthly region figures for MLD
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(5*ratio)
		fig.set_figwidth(8*ratio)

		ax1 = plt.subplot2grid(shape=(5,8), loc=(0,0), rowspan=3, colspan=5)
		ax2 = plt.subplot2grid(shape=(5,8), loc=(3,0), rowspan=2, colspan=5)

		ax3 = plt.subplot2grid(shape=(5,8), loc=(0,5), colspan=3)
		ax4 = plt.subplot2grid(shape=(5,8), loc=(1,5), colspan=3, sharey=ax3)
		ax5 = plt.subplot2grid(shape=(5,8), loc=(2,5), colspan=3, sharey=ax3)
		ax6 = plt.subplot2grid(shape=(5,8), loc=(3,5), colspan=3, sharey=ax3)
		ax7 = plt.subplot2grid(shape=(5,8), loc=(4,5), colspan=3, sharey=ax3)

		ax1.axhline(linestyle='--', linewidth=0.25)

		l1, = ax1.plot(mnth_name, MLD_total, "b"); ax1.set_title('Average Mixed Layer Depth (global) [m]'); ax1.grid(linestyle='--', linewidth=0.25)
		ax1.plot(mnth_name, MLD_reg1, "r", linewidth=0.5);
		ax1.plot(mnth_name, MLD_reg2, "g", linewidth=0.5);
		ax1.plot(mnth_name, MLD_reg3, "c", linewidth=0.5);
		ax1.plot(mnth_name, MLD_reg4, "m", linewidth=0.5);
		ax1.plot(mnth_name, MLD_reg5, "y", linewidth=0.5);

		ax2.plot(mnth_name, MLD_total, "b", label='global average'); 
	
		l3, = ax3.plot(mnth_name, MLD_reg1, "r", label='45N-90N'); ax3.set_title('Average Mixed Layer Depth (regional) [m]');
		l4, = ax4.plot(mnth_name, MLD_reg2, "g", label='15N-45N'); 
		l5, = ax5.plot(mnth_name, MLD_reg3, "c", label='15S-15N'); 
		l6, = ax6.plot(mnth_name, MLD_reg4, "m", label='45S-15S');
		l7, = ax7.plot(mnth_name, MLD_reg5, "y", label='90S-45S');

		lines = [l1, l3, l4, l5, l6, l7]
		labels = ['global average', '45N-90N', '15N-45N', '15S-15N', '45S-15S', '90S-45S']

		ax1.legend(lines, labels, loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7]:
			ax.legend(loc='upper right')
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_monthly_MLD.jpg')
		print(f'Created MLD monthly summary figure for {tmod}')

	except:
		print('Oops, something went awry')

# Run summary on input model
makeSummaryFromBreakdowns(tModel, baseDir)
