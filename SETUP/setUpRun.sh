#!/bin/sh 

date
echo "To use: setUpRun <setUpData.dat> <Full Run ID>"

if [ "$#" -ne 2 ]; then
	exit 1
fi

# Input variables read from command line
setUpDatafile=$1
id=$2

# ----- Meta variables -----
version=$(echo $id | awk -F'_' '{print $1}')
initials=$(echo $id | awk -F'_' '{print $2}')
simulation=$(echo $id | awk -F'_' '{print $3}')

echo "Setting up:                  " $version $initials $simulation
echo "Reading run parameters from: " $setUpDatafile

# ----- Read setup data -----
dataFileFullPath=$(pwd)"/"$setUpDatafile

while IFS= read -r line; do
	if [ ! ${line:0:1} == "#" ]; then

		# Pre-link processing
		name=$(echo $line | awk -F':' '{print $1}')
		val=$(echo $line | awk -F':' '{print $2}')
        
		if [[ $name != *"."* && $name != "namelist"* && $name != "breakdown_parms" ]]; then
			if [[ $name == "yearStart" ]]; then yearStart=$val; fi
			if [[ $name == "yearEnd" ]]; then yearEnd=$val; fi
			if [[ $name == "CO2" ]]; then CO2=$val; fi
			if [[ $name == "forcing" ]]; then forcing=$val; fi
			if [[ $name == "basedir" ]]; then basedir=$val; fi
			if [[ $name == "EMPaveFile" ]]; then EMPaveFile=$val; fi
			if [[ $name == "model" ]]; then Model=$val; fi
			if [[ $name == "type" ]]; then type=$val; fi
			if [[ $name == "compilerKey" ]]; then compKey=$val; fi

			# Tidy up parameters
			if [[ $name == "spinupStart" ]]; then spinupStart=$val; fi
			if [[ $name == "spinupEnd" ]]; then spinupEnd=$val; fi
			if [[ $name == "spinupRestartKeepFrequency" ]]; then spinupRestartKeepFrequency=$val; fi
			if [[ $name == "spinupOutputKeepFrequency" ]]; then spinupOutputKeepFrequency=$val; fi
			if [[ $name == "runRestartKeepFrequency" ]]; then runRestartKeepFrequency=$val; fi
			if [[ $name == "runOutputKeepFrequency" ]]; then runOutputKeepFrequency=$val; fi
			if [[ $name == "keepGrid_T" ]]; then keepGrid_T=$val; fi
			if [[ $name == "keepDiad" ]]; then keepDiad=$val; fi
			if [[ $name == "keepPtrc" ]]; then keepPtrc=$val; fi
			if [[ $name == "keepIce" ]]; then keepIce=$val; fi
			if [[ $name == "keepGrid_U" ]]; then keepGrid_U=$val; fi
			if [[ $name == "keepGrid_V" ]]; then keepGrid_V=$val; fi
			if [[ $name == "keepGrid_W" ]]; then keepGrid_W=$val; fi
			if [[ $name == "keepLimPhy" ]]; then keepLimPhy=$val; fi
		fi
	fi
done < $dataFileFullPath

prevYear=$(($yearStart-1))

# ----- Move to or create model directory -----
# Adjust for a possible ~ expansion problem
if [ ${basedir:0:1} == "~" ]; then 
	homearea=$(readlink -f ~)
	basedir=$homearea${basedir:1:${#basedir}-1}
fi

modelDir=$basedir$id
if [ ! -d $modelDir ]; then
	mkdir $modelDir
fi

# Copy the setUpData file to the directory
cp $setUpDatafile $modelDir
cd $modelDir

# ----- Output tidy up parameters -----
echo $spinupStart $spinupEnd $spinupRestartKeepFrequency $spinupOutputKeepFrequency $runRestartKeepFrequency $runOutputKeepFrequency $keepGrid_T $keepDiad $keepPtrc $keepIce $keepGrid_V $keepGrid_U $keepGrid_W $keepLimPhy > tidy_parms

# ----- Create links -----
rm -f opa

while IFS= read -r line; do
	if [ ! ${line:0:1} == "#" ]; then
		
		# Pre-link processing
		name=$(echo $line | awk -F':' '{print $1}')
		val=$(echo $line | awk -F':' '{print $2}')

		# Make links for all .nc and xml files 
		if [[ $name == *"."* && $name != "namelist"* && $name != *".xml" ]]; then 

			if [[ $name == "restart"* ]]; then 
				# Only create link if restart file does not exist already (i.e. run already started in folder)
				if [ ! -f restart_0000.nc ]; then
					ln -fs $val $name
				fi
			else
				ln -fs $val $name
			fi
        	fi

		# Copy namelists in a way so changes can be made
		if [[ $name == "namelist"* ]]; then
			if [ -f $name ]; then
				echo $name " exists so no fresh copy made "
			else
				cp $val $name
			fi
		fi
		
		# Copy xml files in a way so changes can be made
		if [[ $name == *".xml" ]]; then
			if [ -f $name ]; then
				echo $name " exists so no fresh copy made "
			else
				cp $val $name
			fi
		fi

		# Copy breakdown_parms in a way so changes can be made
		if [[ $name == "breakdown_parms" ]]; then
			if [ -f $name ]; then
				echo $name " exists so no fresh copy made "
			else
				cp $val $name
			fi
		fi

		# Copy the executable over, good to keep these.
		if [[ $name == "opa"*$Model ]]; then
			echo "copying " $val $name " for " $Model
			cp $val $name
			ln -s $name opa
		fi
	fi
done < $dataFileFullPath

# Link EMP file
if [ ! -f EMPave_${prevYear}.dat ]; then
	rm -f EMPave_${prevYear}.dat
	ln -fs $EMPaveFile EMPave_${prevYear}.dat
	ln -fs EMPave_${prevYear}.dat EMPave_old.dat
else
	echo "EMPave exists, using existing file!!!"
	ln -fs EMPave_${prevYear}.dat EMPave_old.dat
fi

# ----- Check compiler keys -----
cp $compKey .

grep key_trc_piic $compKey > tmp
if [ -s tmp ]; then
	PIIC=piic
fi

grep key_c14b $compKey > tmp
if [ -s tmp ]; then
	C14=c14
fi

rm tmp
echo ${PIIC^^} ${C14^^}

# ----- Process flags -----
# CO2
rm -f atmco2.dat

echo $CO2
if [ $CO2 == "VARIABLE" ]; then
	ln -s atmco2.dat.variable atmco2.dat
else
	ln -s atmco2.dat.static atmco2.dat
fi

# Forcing 
rm -f namelist_ref

echo $forcing
if [ $forcing == "NCEP" ]; then
	ln -s namelist_ref_ncep_first_year namelist_ref_first_year
	ln -s namelist_ref_ncep namelist_ref_all_years 
	ln -s namelist_ref_ncep_looping namelist_ref_loop
elif [ $forcing == "ERA" ]; then
	ln -s namelist_ref_era_first_year namelist_ref_first_year
	ln -s namelist_ref_era namelist_ref_all_years 
	ln -s namelist_ref_era_looping namelist_ref_loop
else
	ln -s namelist_ref_jra_first_year namelist_ref_first_year
	ln -s namelist_ref_jra namelist_ref_all_years 
	ln -s namelist_ref_jra_looping namelist_ref_loop
fi

# Type
echo $type

if [ $type == "BIAS" ]; then
	rm -f namelist_ref_all_years
	ln -sf namelist_ref_loop namelist_ref_all_years
fi

# Initial year; if a CPU based restart file does not exist, then this is the first year
if [ ! -f restart_0000.nc ]; then
	ln -s namelist_ref_first_year namelist_ref
else
	ln -s namelist_ref_all_years namelist_ref
fi

# Check for variables set in the namelists
TR=$( grep "nn_sstr " namelist_ref | awk -F' ' '{print $3}' )
SR=$( grep "nn_sssr " namelist_ref | awk -F' ' '{print $3}' )
LP=$( grep "ln_lop" namelist_top_ref | awk -F' ' '{print $3}' )

if [ $TR = 1 ]; then
	echo "TEMPERATURE RESTORING"
else
	echo "NO TEMPERATURE RESTORING"
fi

if [ $SR = 1 ]; then
	echo "SALINITY RESTORING"
else
	echo "NO SALINITY RESTORING"
fi

# Check that files for LIMPHY are set correctly
IO=$( grep "^iodef.xml:" $setUpDatafile | awk -F'/' '{print $NF}' )
KP=$( grep "^keepLimPhy:" $setUpDatafile | awk -F':' '{print $NF}' )
err=0

if [ $LP = ".true." ]; then
	echo "LIMPHY TRUE"

	if [ $IO != iodef_tom12limphy.xml ]; then
		echo -e "\e[1;31mWARNING\e[0m : IODEF file does not match"
		err=1
	fi

	if [ $KP != 1 ]; then
		echo -e "\e[1;31mWARNING\e[0m : KEEP value for LimPhy not set to 1"
		err=1
	fi
else
	echo "LIMPHY FALSE"

	if [ $IO != iodef_tom12${PIIC}${C14}.xml ]; then
		echo -e "\e[1;31mWARNING\e[0m : IODEF file does not match"
		err=1
	fi

	if [ $KP != 0 ]; then
		echo -e "\e[1;31mWARNING\e[0m : KEEP value for LimPhy not set to 0"
		err=1
	fi
fi

if [ $err == 1 ]; then
	exit 2
fi

# Get code version
codeVersion=$( grep "opa_*$Model" $setUpDatafile | awk -F'/' '{print$(NF-5)}' )

# ----- Create copies of files used for run -----
# Get dgom file
if [ ! -f dgom ]; then
	cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/dgom_v2 dgom
fi

# Get tidying up scripts
cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/tidyup.sh .
cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/tidyupJob .

# Get breakdown scripts
cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/breakdown/breakdown*.py .

# If breakdown_parms does not exist (as specified in setUpData file) copy in default
if [ ! -f breakdown_parms ]; then
	echo "Copying in default script for breakdown_parms"
	cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/breakdown_parms breakdown_parms
fi

# Get monitor scripts and files
cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/monitor/* .

# Save parameters needed for creating html file
echo $id $codeVersion $(date '+%d-%b-%Y') $yearStart $yearEnd ${CO2,,} $forcing ${type,,} $TR $SR > html_parms

# Get setUpRun script
cp /gpfs/data/greenocean/software/CodeReview/KAR12/SETUP/setUpRun.sh .

# ----- Export parameters the dgom file will need -----
yearToRun=$yearStart

echo "Exporting " $yearToRun $yearEnd $basedir $modelDir $simulation $Model
export yearToRun yearStart yearEnd basedir modelDir simulation Model

read -p "Press any key to run it? (cntr+c otherwise)"

sbatch -J${simulation}${yearToRun} < dgom

echo "To check status of job 'squeue | grep <you user name>' "

# ----- Save model details -----
echo $2 "("$(date '+%a %d %b %T %Z %Y')")" >> /gpfs/afm/greenocean/software/runs/modelRuns.txt
echo $2 "("$(date '+%a %d %b %T %Z %Y')")" >> ${basedir}/modelRuns.txt
vi ${basedir}/modelRuns.txt

