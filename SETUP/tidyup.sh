#!/bin/sh

# parameters for tidying up
yearFrom=$1
yearTo=$2
version=$3
homeDir=$4

# Get parameters as specified in setUpData.dat file
read -r -a parms < tidy_parms 

spinupStart=${parms[0]}
spinupEnd=${parms[1]}
spinupRestartKeepFrequency=${parms[2]}
spinupOutputKeepFrequency=${parms[3]}
runRestartKeepFrequency=${parms[4]}
runOutputKeepFrequency=${parms[5]}
keepGrid_T=${parms[6]}
keepDiad=${parms[7]}
keepPtrc=${parms[8]}
keepIce=${parms[9]}
keepGrid_V=${parms[10]}
keepGrid_U=${parms[11]}
keepGrid_W=${parms[12]}
keepLimPhy=${parms[13]}

# Echo parameters to tidy.log
echo $spinupStart
echo $spinupEnd
echo $spinupRestartKeepFrequency
echo $spinupOutputKeepFrequency
echo $runRestartKeepFrequency
echo $runOutputKeepFrequency
echo $keepGrid_T
echo $keepDiad
echo $keepPtrc
echo $keepIce
echo $keepGrid_V
echo $keepGrid_U
echo $keepGrid_W
echo $keepLimPhy

pointsPerYear=5475
freq=1m

model=$(basename "$PWD")
id=${model: -4}
echo "Tidying up year: $1 to $2 for $model"

# Specify directory for copying files to central area
baseDir="/gpfs/afm/greenocean/software/runs/"

# Make directory for run in central area
if [ ! -d $baseDir$model ]; then
	echo "creating directory"
	mkdir $baseDir$model
fi

echo "Copying data centrally"
for (( y=$yearFrom; y<=$yearTo; y++ )); do

	# Python based totalling
	./breakdown.py breakdown_parms ${y} ${y}

	# Run monitor script
	./monitor.py $model $homeDir

	# Process output files
	if [[ $y < $spinupEnd ]]; then
		# years from start
		since=$((y-spinupStart))
		remainder=$(( since % spinupOutputKeepFrequency ))
	else
        	since=$((y-spinupEnd))
		remainder=$(( since % runOutputKeepFrequency ))
	fi
    
	echo $remainder $since

	if [[ "$remainder" -eq 0 ]]; then

		echo "Copying output $y"
		if [[ $keepGrid_T -eq 1 ]]; then cp ORCA2_${freq}_${y}0101_${y}1231_grid_T.nc $baseDir$model; fi
		if [[ $keepDiad -eq 1 ]];   then cp ORCA2_${freq}_${y}0101_${y}1231_diad_T.nc $baseDir$model; fi
		if [[ $keepPtrc -eq 1 ]];   then cp ORCA2_${freq}_${y}0101_${y}1231_ptrc_T.nc $baseDir$model; fi
		if [[ $keepIce -eq 1 ]];    then cp ORCA2_${freq}_${y}0101_${y}1231_icemod.nc $baseDir$model; fi
		if [[ $keepGrid_U -eq 1 ]]; then cp ORCA2_${freq}_${y}0101_${y}1231_grid_U.nc $baseDir$model; fi
		if [[ $keepGrid_V -eq 1 ]]; then cp ORCA2_${freq}_${y}0101_${y}1231_grid_V.nc $baseDir$model; fi
		if [[ $keepGrid_W -eq 1 ]]; then cp ORCA2_${freq}_${y}0101_${y}1231_grid_W.nc $baseDir$model; fi
		if [[ $keepLimPhy -eq 1 ]]; then cp ORCA2_${freq}_${y}0101_${y}1231_limphy.nc $baseDir$model; fi

		echo "Copying extra set up data and EMPave files"
		cp EMPave_${y}.dat $baseDir$model
		cp namelist* $baseDir$model
		cp *xml $baseDir$model
		cp setUpData*dat $baseDir$model
		cp ocean.output $baseDir$model
		cp opa $baseDir$model

		echo "Deleting local data if it exists centrally"
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_grid_T.nc || $keepGrid_T -eq 0 ]]; then rm -f ORCA2_${freq}_${y}0101_${y}1231_grid_T.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_diad_T.nc || $keepDiad -eq 0 ]];   then rm -f ORCA2_${freq}_${y}0101_${y}1231_diad_T.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_ptrc_T.nc || $keepPtrc -eq 0 ]];   then rm -f ORCA2_${freq}_${y}0101_${y}1231_ptrc_T.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_icemod.nc || $keepIce -eq 0 ]];    then rm -f ORCA2_${freq}_${y}0101_${y}1231_icemod.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_grid_U.nc || $keepGrid_U -eq 0 ]]; then rm -f ORCA2_${freq}_${y}0101_${y}1231_grid_U.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_grid_V.nc || $keepGrid_V -eq 0 ]]; then rm -f ORCA2_${freq}_${y}0101_${y}1231_grid_V.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_grid_W.nc || $keepGrid_W -eq 0 ]]; then rm -f ORCA2_${freq}_${y}0101_${y}1231_grid_W.nc; fi
		if [[ -f $baseDir$model/ORCA2_${freq}_${y}0101_${y}1231_limphy.nc || $keepLimPhy -eq 0 ]]; then rm -f ORCA2_${freq}_${y}0101_${y}1231_limphy.nc; fi

		echo "Creating symlinks for copied data"
		if [[ $keepGrid_T -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_grid_T.nc; fi
		if [[ $keepDiad -eq 1 ]];   then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_diad_T.nc; fi
		if [[ $keepPtrc -eq 1 ]];   then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_ptrc_T.nc; fi
		if [[ $keepIce -eq 1 ]];    then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_icemod.nc; fi
		if [[ $keepGrid_U -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_grid_U.nc; fi
		if [[ $keepGrid_V -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_grid_V.nc; fi
		if [[ $keepGrid_W -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_grid_W.nc; fi
		if [[ $keepLimPhy -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_${freq}_${y}0101_${y}1231_limphy.nc; fi
	else
		rm -f ORCA2_${freq}_${y}0101_${y}1231_grid_*.nc
		rm -f ORCA2_${freq}_${y}0101_${y}1231_diad_T.nc
		rm -f ORCA2_${freq}_${y}0101_${y}1231_ptrc_T.nc
		rm -f ORCA2_${freq}_${y}0101_${y}1231_icemod.nc
		rm -f ORCA2_${freq}_${y}0101_${y}1231_limphy.nc
	fi

	# Process restart files
	if [[ $y < $spinupEnd ]]; then
		# Years from start
		since=$(( y-spinupStart ))
		remainder=$(( since % spinupRestartKeepFrequency ))
	else
		since=$(( y-spinupEnd ))
		remainder=$(( since % runRestartKeepFrequency ))
	fi

	# Reset since to correctly count points for linking and deleting purposes
	yr=$yearStart
	since=$(( y-yr ))

	# This is concerned with the previous year's restart files, no restarts are deleted until the next year is complete    
	points=$(( pointsPerYear * since ))
	printf -v timestep "%08d" $points

	echo $points $timestep $remainder $since

	if [[ "$remainder" -eq 0 ]]; then
		echo "Copying restart $y"
		cp ORCA2_*${timestep}_restart_*.nc $baseDir$model
		rm -f ORCA2_*${timestep}_restart_*.nc
		ls -1 ${baseDir}${model}/ORCA2_*${timestep}_restart_*.nc | awk '{print "ln -s "$1 }' | bash
	else
		echo "Deleting local data"
		rm -f ORCA2_*${timestep}_restart_*.nc
	fi
done

# Commands to execute at the end of the simulation (final year only)
if [[ $yearTo -eq $yearEnd ]]; then

	# copy all breakdown files to base model directory
	echo "copying breakdown files"
	cp breakdown* $baseDir$model

	# run script to generate monthly regional plots
	./monthly.py $model $homeDir

	# run script to generate vertical depth plots
	./verticalDepth.py $model $homeDir $yearTo

	# run script to generate maps in FERRET for the final year
	./annualMaps.sh $model $yearTo $homeDir

	# run script to create html file
	./createHTML.sh $homeDir
fi

