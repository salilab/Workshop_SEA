#!/bin/bash

mkdir all_models.999

for i in `seq 0 1`; 
do
    for j in `seq 1 5`
    do
	rmf_file=`./Scripts/process_output.py -s 'rmf_file_full_path' -f ../kmeans_5_modeling$i/all_models.4/stat.0.out |\
	    grep -v ^# | awk -v lnum="$j" 'NR==lnum {print $2}'`

        rmf_frame=`./Scripts/process_output.py -s 'rmf_frame_index' -f ../kmeans_5_modeling${i}/all_models.4/stat.0.out | \
	    grep -v ^# | awk -v lnum="$j" 'NR==lnum {print $2}'`

	echo $rmf_file $rmf_frame

	rmf_slice $rmf_file  ./all_models.999/$i'_'$((j-1))'.rmf3' -f $rmf_frame

    done
done

