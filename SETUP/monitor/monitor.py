#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pathlib
import sys

# Script for monitoring variables using breakdowns created during the model run process
# Creates figure summaries at the end of each year

tModel = sys.argv[1]
baseDir = sys.argv[2]

def makeSummaryFromBreakdowns(tmod, modBaseDir):
    
	saveDir = f'/{modBaseDir}/monitor/{tmod}/'
	pathlib.Path(saveDir).mkdir(parents=True, exist_ok=True)
    
	try:
		### Surface annual
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.sur.annual.dat', sep='\t')
		print('Read surface annual')
        
		year = w.year[2:].to_numpy().astype(float)
		Cflx = w.Cflx[2:].to_numpy().astype(float)

		### Volume annual
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.vol.annual.dat', sep='\t')
		print('Read volume annual')

		year = w.year[2:].to_numpy().astype(float)
		PPT = w.PPT[2:].to_numpy().astype(float)
		PROARA = w.proara[2:].to_numpy().astype(float)
		PROCOC = w.prococ[2:].to_numpy().astype(float)
		PROBSI = w.probsi[2:].to_numpy().astype(float)

		PROCACO3 = PROARA + PROCOC

		GRACRU = w.GRACRU[2:].to_numpy().astype(float)
		GRAGEL = w.GRAGEL[2:].to_numpy().astype(float)	
		GRAMES = w.GRAMES[2:].to_numpy().astype(float)
		GRAPRO = w.GRAPRO[2:].to_numpy().astype(float)
		GRAPTE = w.GRAPTE[2:].to_numpy().astype(float)

		SPT = GRACRU + GRAGEL + GRAMES + GRAPRO + GRAPTE

		### Level annual
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.lev.annual.dat', sep='\t')
		print('Read level annual')

		year = w.year[2:].to_numpy().astype(float)
		EXP = w.EXP[2:].to_numpy().astype(float)
		EXPCO3 = w.ExpCO3[2:].to_numpy().astype(float)
		SI_FLX = w.sinksil[2:].to_numpy().astype(float)

		### Average annual
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.ave.annual.dat', sep='\t')
		print('Read average annual')

		year = w.year[2:].to_numpy().astype(float)
		TChl = w.TChl[2:].to_numpy().astype(float)

		PO4 = w.PO4[2:].to_numpy().astype(float)
		NO3 = w.NO3[2:].to_numpy().astype(float)
		Fer = w.Fer[2:].to_numpy().astype(float)
		Si = w.Si[2:].to_numpy().astype(float)
		O2 = w.O2[2:].to_numpy().astype(float)

		nFer = Fer * 1000 # convert units from micromoles to nanomoles
		nPO4 = PO4 / 122 # adjustment for units

		SST  = w.tos[2:].to_numpy().astype(float)
		SSS  = w.sos[2:].to_numpy().astype(float)
		MLD  = w.mldr10_1[2:].to_numpy().astype(float)

		### Integration annual
		w = pd.read_csv(f'/{modBaseDir}/{tmod}/breakdown.int.annual.dat', sep='\t')
		print('Read integration annual')

		year = w.year[2:].to_numpy().astype(float)
		BAC = w.BAC[2:].to_numpy().astype(float)
		COC = w.COC[2:].to_numpy().astype(float)
		CRU = w.CRU[2:].to_numpy().astype(float)
		DIA = w.DIA[2:].to_numpy().astype(float)
		FIX = w.FIX[2:].to_numpy().astype(float)
		GEL = w.GEL[2:].to_numpy().astype(float)
		MES = w.MES[2:].to_numpy().astype(float)
		MIX = w.MIX[2:].to_numpy().astype(float)
		PHA = w.PHA[2:].to_numpy().astype(float)
		PIC = w.PIC[2:].to_numpy().astype(float)
		PRO = w.PRO[2:].to_numpy().astype(float)
		PTE = w.PTE[2:].to_numpy().astype(float)

		VIR = w.BAC[2:].to_numpy().astype(float)

		PHY = COC + DIA + FIX + MIX + PHA + PIC
		ZOO = BAC + CRU + GEL + MES + PRO + PTE
		TOT = PHY + ZOO

		### Creating global figures
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(6*ratio) # rows
		fig.set_figwidth(11*ratio) # columns

		ax1 = plt.subplot2grid(shape=(6,11), loc=(0,0), rowspan=2, colspan=3)
		ax2 = plt.subplot2grid(shape=(6,11), loc=(0,3), rowspan=2, colspan=3)
		ax3 = plt.subplot2grid(shape=(6,11), loc=(0,6), rowspan=2, colspan=3)

		ax4 = plt.subplot2grid(shape=(6,11), loc=(2,0), rowspan=2, colspan=3)
		ax5 = plt.subplot2grid(shape=(6,11), loc=(2,3), rowspan=2, colspan=3)
		ax6 = plt.subplot2grid(shape=(6,11), loc=(2,6), rowspan=2, colspan=3)

		ax7 = plt.subplot2grid(shape=(6,11), loc=(4,0), rowspan=2, colspan=3)
		ax8 = plt.subplot2grid(shape=(6,11), loc=(4,3), rowspan=2, colspan=3)
		ax9 = plt.subplot2grid(shape=(6,11), loc=(4,6), rowspan=2, colspan=3)

		ax1.plot(year, Cflx, "b"); ax1.set_title('Surface Cflx [PgCarbonPerYr]');
		ax2.plot(year, TChl, "r"); ax2.set_title('Average Surface TChl [ug Chl/L]');
		ax3.plot(year, SPT,  "y"); ax3.set_title('Total Volume of Secondary Production [PgCarbonPerYr]');

		ax4.plot(year, PPT,      "g"); ax4.set_title('Total Volume of PP C [PgCarbonPerYr]');
		ax5.plot(year, PROCACO3, "g"); ax5.set_title('Total Volume of PP CaCO3 [PgCarbonPerYr]');
		ax6.plot(year, PROBSI,   "g"); ax6.set_title('Total Volume of PP Si [PgCarbonPerYr]');

		ax7.plot(year, EXP,    "m"); ax7.set_title('OC flux (export) at 100m depth [PgCarbonPerYr]');
		ax8.plot(year, EXPCO3, "m"); ax8.set_title('CaCO3 flux at 100m depth [PgCarbonPerYr]');
		ax9.plot(year, SI_FLX, "m"); ax9.set_title('Si flux at 100m depth [PgCarbonPerYr]');

		for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
			ax.grid(linestyle='--', linewidth=0.25)

		# Add in shaded regions for observations
		y1 = 0.2921;         ax2.axhline(y1, color='k', linestyle='dashed', alpha=0.75) # TCHL
		y1 = 51;  y2 = 65;   ax4.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\') # PPT
		y1 = 7.8; y2 = 12.2; ax7.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\') # EXP
		y1 = 0.6; y2 = 1.1;  ax8.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\', label="observation range") # EXP CaCO3
		y1 = 3.4;            ax9.axhline(y1, color='k', linestyle='dashed', alpha=0.75, label="observation value") # SI flux

		fig.legend(bbox_to_anchor=(0.818,0.969+10/1163), loc='upper left', borderpad=1, fontsize=15)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_global.jpg')
		print(f'Created global summary figure for {tmod}')

		### Creating PFT figures
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(6*ratio)
		fig.set_figwidth(9*ratio)

		ax1 = plt.subplot2grid(shape=(6,9), loc=(0,0), rowspan=3, colspan=5)

		ax2 = plt.subplot2grid(shape=(6,9), loc=(0,5), colspan=2)
		ax3 = plt.subplot2grid(shape=(6,9), loc=(1,5), colspan=2)
		ax4 = plt.subplot2grid(shape=(6,9), loc=(2,5), colspan=2)
		ax5 = plt.subplot2grid(shape=(6,9), loc=(3,5), colspan=2)
		ax6 = plt.subplot2grid(shape=(6,9), loc=(4,5), colspan=2)
		ax7 = plt.subplot2grid(shape=(6,9), loc=(5,5), colspan=2)

		ax8  = plt.subplot2grid(shape=(6,9), loc=(0,7), colspan=2)
		ax9  = plt.subplot2grid(shape=(6,9), loc=(1,7), colspan=2)
		ax10 = plt.subplot2grid(shape=(6,9), loc=(2,7), colspan=2)
		ax11 = plt.subplot2grid(shape=(6,9), loc=(3,7), colspan=2)
		ax12 = plt.subplot2grid(shape=(6,9), loc=(4,7), colspan=2)
		ax13 = plt.subplot2grid(shape=(6,9), loc=(5,7), colspan=2)

		ax14 = plt.subplot2grid(shape=(6,9), loc=(4,3), colspan=2)

		ax1.plot(year, TOT, "g", label="Total biomass"); ax1.set_title('Total biomass integrated over top 200m [PgCarbon]');
		ax1.plot(year, PHY, "b", label="Phytoplankton biomass");
		ax1.plot(year, ZOO, "r", label="Zooplankton biomass");

		ax2.plot(year, PIC, "b", label="PIC"); y1 = 0.28;  y2 = 0.52;  ax2.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax3.plot(year, PHA, "b", label="PHA"); y1 = 0.11;  y2 = 0.69;  ax3.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax4.plot(year, MIX, "b", label="MIX");
		ax5.plot(year, DIA, "b", label="DIA"); y1 = 0.013; y2 = 0.75;  ax5.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax6.plot(year, COC, "b", label="COC"); y1 = 0.001; y2 = 0.032; ax6.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax7.plot(year, FIX, "b", label="FIX"); y1 = 0.008; y2 = 0.12;  ax7.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')

		ax8.plot(year, GEL, "r", label="GEL");  y1 = 0.10; y2 = 3.11; ax8.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax9.plot(year, PRO, "r", label="PRO");  y1 = 0.10; y2 = 0.37; ax9.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax10.plot(year, BAC, "r", label="BAC"); y1 = 0.25; y2 = 0.26; ax10.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax11.plot(year, CRU, "r", label="CRU"); y1 = 0.01; y2 = 0.64; ax11.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax12.plot(year, PTE, "r", label="PTE");
		ax13.plot(year, MES, "r", label="MES"); y1 = 0.21; y2 = 0.34; ax13.axhspan(y1, y2, color='black', alpha=0.25, fill=0, lw=0, hatch='\\\\')
		ax14.plot(year, VIR, "m", label="VIR");

		ax1.set_ylim(bottom=0)
		ax1.grid(linestyle='--', linewidth=0.25)
		ax1.legend(loc='upper right')

		for ax in [ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14]:
			ax.set_ylim(0,0.5)
			ax.grid(linestyle='--', linewidth=0.25)
			ax.legend(loc='upper right')

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_pfts.jpg')
		print(f'Created PFT summary figure for {tmod}')

		### Creating nutrient figures
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(4*ratio) # rows
		fig.set_figwidth(11*ratio) # columns

		ax1 = plt.subplot2grid(shape=(4,11), loc=(0,0), rowspan=2, colspan=3)
		ax2 = plt.subplot2grid(shape=(4,11), loc=(0,3), rowspan=2, colspan=3)
		ax3 = plt.subplot2grid(shape=(4,11), loc=(0,6), rowspan=2, colspan=3)
		ax4 = plt.subplot2grid(shape=(4,11), loc=(2,0), rowspan=2, colspan=3)
		ax5 = plt.subplot2grid(shape=(4,11), loc=(2,3), rowspan=2, colspan=3)

		ax1.plot(year, nPO4, "b"); ax1.set_title('Average Surface Phosphate [umol/L]');
		ax2.plot(year, NO3,  "r"); ax2.set_title('Average Surface Nitrates [umol/L]');
		ax3.plot(year, nFer, "g"); ax3.set_title('Average Surface Iron [nmol/L]');
		ax4.plot(year, Si,   "m"); ax4.set_title('Average Surface Silica [umol/L]');
		ax5.plot(year, O2,   "y"); ax5.set_title('Average Surface Oxygen [umol/L]');

		for ax in [ax1, ax2, ax3, ax4, ax5]:
			ax.grid(linestyle='--', linewidth=0.25)

		# Observation values
		y1 = 0.530; ax1.axhline(y1, color='k', linestyle='dashed', alpha=0.75) # PO4
		y1 = 5.152; ax2.axhline(y1, color='k', linestyle='dashed', alpha=0.75) # NO3
		y1 = 7.485; ax4.axhline(y1, color='k', linestyle='dashed', alpha=0.75) # Sil
		y1 = 251.1; ax5.axhline(y1, color='k', linestyle='dashed', alpha=0.75, label="observation value") # O2

		fig.legend(bbox_to_anchor=(0.818,0.9525+10/762), loc='upper left', borderpad=1, fontsize=15)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_nutrients.jpg')
		print(f'Created nutrients summary figure for {tmod}')

		### Creating global average figures for physics
		ratio = 2
		fig = plt.figure()
		fig.set_figheight(2*ratio) # rows
		fig.set_figwidth(9*ratio) # columns

		ax1 = plt.subplot2grid(shape=(2,9), loc=(0,0), rowspan=2, colspan=3)
		ax2 = plt.subplot2grid(shape=(2,9), loc=(0,3), rowspan=2, colspan=3)
		ax3 = plt.subplot2grid(shape=(2,9), loc=(0,6), rowspan=2, colspan=3)

		ax1.plot(year, SST, "b"); ax1.set_title('Average Sea Surface Temperature [degC]');
		ax2.plot(year, SSS, "r"); ax2.set_title('Average Sea Surface Salinity [1e-3]');
		ax3.plot(year, MLD, "g"); ax3.set_title('Average Mixed Layer Depth [m]');

		for ax in [ax1, ax2, ax3]:
			ax.grid(linestyle='--', linewidth=0.25)

		plt.tight_layout()
		fig.savefig(f'{saveDir}/{tmod}_summary_physics.jpg')
		print(f'Created physics summary figure for {tmod}')
		
	except:
		print('Oops, something went awry')

# Run summary on input model
makeSummaryFromBreakdowns(tModel, baseDir)
