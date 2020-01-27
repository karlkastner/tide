#!/bin/bash
# 2013-05-19 14:44:45
# Karl Kästner, Berlin
set -u +x

# generates TPXO input and runs TPXO

# TODO no magic file names
dir_tpxo=~/large/phd/gis/tide-tpxo/OTPS1_tpxo7
export PATH=$PATH:$dir_tpxo

if [ $# -gt 2 ]
then
	iname=$1;
	lat=$2;
	lon=$3;
	y0=$4
	ye=$5
else
	echo $@
	echo Usage: $0 filname lat lon
	exit -1;
#	iname="lat_lon_besar_mouth.csv";
fi

#dir_out=$(basename $iname)
dir_out=$iname
#${iname%.*}
echo $dir_out
#table-pontianak
mkdir -p $dir_out

script=$(readlink -f "$0")
dir=$(dirname "$script")
#ROOTFOLDER=~/phd/src/lib/tide/@Tidetable/

# generate table of input locations and dates
perl $dir/tidetable_generate_tpxo_input.pl $iname $dir_out $lat $lon $y0 $ye
#perl ~/phd/src/lib/tide/@Tidetable/tidetable_generate_tpxo_input.pl $iname $dir_out $lat $lon
# > $dir_out/lat_lon_time

if [ $? -eq 0 ]
then
	# run the tide prediction
	$dir_tpxo/predict_tide < $dir_out/level.inp 
	$dir_tpxo/predict_tide < $dir_out/current.inp
	$dir_tpxo/extract_HC   < $dir_out/hc.inp
	
	# extract tidal constituents
	tail -n +4 $dir_out/hc.out | sed -e 's/  *//' | sed -e 's/  */\n/g'  | sed 'N;s/\n/\t/' | tail -n +2 > $dir_out/hc_.out
	tail -n -2 $dir_out/hc.out | head -n 1 | sed -e 's/  *//' | sed -e 's/  */\n/g'  | sed 'N;s/\n/\t/' | cut -f 1  | tail -n +2 >> $dir_out/constituent_names.txt
	paste $dir/tidal-constituents-periods.csv $dir_out/hc_.out > $dir_out/tidal-constituents.dat
	paste $dir_out/constituent_names.txt $dir/tidal-constituents-periods.csv $dir_out/hc_.out > $dir_out/tidal-constituents-named.dat
	#paste ~/phd/src/lib/tide/@Tidetable/tidal-constituents-periods.dat $dir_out/hc_.out > $dir_out/tidal-constituents.dat
	# remove header
	tail -n +7 $dir_out/level.out > $dir_out/level_.out
	tail -n +7 $dir_out/current.out > $dir_out/current_.out
	# combine the level and current table
	paste $dir_out/level_.out $dir_out/current_.out > $dir_out/tide.out
	# separate date and time strings into separate columns
	sed 's/:/ /g' $dir_out/tide.out   | sed -s 's/\([0-9][0-9]\).\([0-9][0-9]\).\([0-9][0-9][0-9][0-9]\)/\1 \2 \3/g' > $dir_out/tide_.out
else
	echo "Error: Perl script failed";
	exit -2;	
fi

