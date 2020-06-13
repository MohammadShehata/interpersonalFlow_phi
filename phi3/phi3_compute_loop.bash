#!/bin/bash

# Computes for every file in the given directory $1

source_dir=../tpms/
data_dir=4ch_diff_perSong_2t2_GP7Left/
out_dir=4ch_diff_perSong_2t2_GP7Left/

for file in ${source_dir}${data_dir}*; do # TODO: use input as directory location
	
	filename=$(basename $file)
	
	echo $file
	
	# Check current job list
	squeue -u aleung > job_list
	jobs=$(wc -l < job_list)
	while [ $jobs -ge 495 ]; do # Job limit is 500
		echo "too many jobs, sleeping"
		sleep 30s
		squeue -u aleung > job_list
		jobs=$(wc -l < job_list)
		echo "now there are $jobs jobs"
	done
	
	# Check if file exists, if not, submit job
	if [ ! -e ${out_dir}${filename} ]; then
		echo "submitting job $file"
		sbatch --job-name="${filename}" --output="logs/${filename}.out" --error="logs/${filename}.err" bash_loop_sbatch.bash $source_dir $data_dir $out_dir $filename
	fi
done
