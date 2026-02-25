#!/bin/sh

# Read in parameters for creating html file
read -r -a parms < html_parms

run=${parms[0]}
version=${parms[1]}
date=${parms[2]}
start=${parms[3]}
end=${parms[4]}
co2=${parms[5]^}
forcing=${parms[6]}
type=${parms[7]^}

if [[ ${parms[8]} -eq 0 ]]; then
	temperature="No"
else
	temperature="Yes"
fi

if [[ ${parms[9]} -eq 0 ]]; then
	salinity="No"
else
	salinity="Yes"
fi

# Set up save directory
baseDir=$1
monitor="monitor/"
saveDir="${baseDir}${monitor}${run}/"

# Move template to match run and replace parameters
mv template.html ${run}.html

sed -i "s/{{identifier}}/${run}/g" ${run}.html
sed -i "s/{{version}}/${version}/g" ${run}.html
sed -i "s/{{date}}/${date}/g" ${run}.html
sed -i "s/{{start_year}}/${start}/g" ${run}.html
sed -i "s/{{end_year}}/${end}/g" ${run}.html
sed -i "s/{{co2}}/${co2}/g" ${run}.html
sed -i "s/{{forcing}}/${forcing}/g" ${run}.html
sed -i "s/{{type}}/${type}/g" ${run}.html
sed -i "s/{{temperature_restoring}}/${temperature}/g" ${run}.html
sed -i "s/{{salinity_restoring}}/${salinity}/g" ${run}.html

# Move html file to save directory
mv ${run}.html ${saveDir}
