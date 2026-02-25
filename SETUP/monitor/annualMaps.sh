#!/bin/sh

module add ferret/7.6.0
module add ImageMagick

run=$1
year=$2
baseDir=$3

monitor="monitor/"
saveDir="${baseDir}${monitor}${run}/"

# Create FERRET maps
echo "ferret -gif -script maps.jnl" $baseDir $run $year | bash
echo "ferret -gif -script mapsPFT.jnl" $baseDir $run $year | bash

# Apply crop to images
vars=("bac" "cflx" "coc" "cru" "dia" "dpco2" "exp" "fer" "fix" "gel" "mes" "mix" "mld_ave" "mld_max" "no3" "o2" "pha" "pic" "po4" "ppint" "pro" "pte" "salin" "sil" "tchl" "temp")

for i in ${!vars[@]}; do
	FILE=${run}_${year}_${vars[$i]}.gif

	magick $FILE -crop 750x580+18+82 ${FILE::-4}.png
done

# Create difference plots in FERRET
echo "ferret -gif -script mapsDiff.jnl" $baseDir $run $year | bash

varsD=("dpco2" "no3" "o2" "po4" "salin" "sil" "tchl" "temp")

for i in ${!varsD[@]}; do
	FILE=${run}_${year}_difference_${varsD[$i]}.gif

	magick $FILE -crop 750x580+18+82 ${FILE::-4}.png
done

# Delete old .gif files
rm ${run}_${year}_*.gif

# Move files to save directory
mv *.png ${saveDir}
