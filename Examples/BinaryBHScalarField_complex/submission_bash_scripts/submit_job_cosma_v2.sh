#!/bin/bash
#
# This script should make a new directory in GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for COSMA6

work_dir=/cosma/home/dp174/dc-bamb1/GRChombo/Examples/BinaryBHScalarField_complex
cd $work_dir
data_directory=/cosma6/data/dp174/dc-bamb1/GRChombo_data/BinaryBHSF_complex

# specify the input params for each run I want to submit
# list for each is: mu, dt_mult, G, l, m

run0002=(params_10_orbits_gaussian.txt run0002_gaussian_kappa0.0125_mu0.34_L512_N128_G1_max_level8_boxsize16_Mcloud0.1M)
run0106=(params_6_orbits_G0.txt run0106_6_orbits_gaussian_kappa0.0125_mu0.373_L512_N128_G0_max_level8_boxsize16)


run_list=(
	run0002
	run0106
)

restart_hash="#"
restart_num="000000"

echo ${restart_hash}
echo ${restart_num}

for run in "${run_list[@]}"
do
  	cd $work_dir
        # extract parameters    
        val="$run[0]"; params_file="${!val}"
        val="$run[1]"; dir_name="${!val}"
	
	new_dir=${dir_name}
	   
        echo ${new_dir}
        new_dir_path=${data_directory}/${new_dir}
	#
	mkdir -p ${new_dir_path}
	# get latest checkpoint file
        cd ${new_dir_path}
	chk_file=BinaryBHSFChk_${restart_num}.3d.hdf5
	for chk in BinaryBHSFChk_*.hdf5
	do
	    chk_file=$chk
	done
	echo $chk_file
	#
	cd ${work_dir}
        cp slurm_files/slurm_submit_cosma7 ${new_dir_path}/slurm_submit

	cp ${params_file} ${new_dir_path}/params.txt
	cp BinaryBHLevel.cpp ${new_dir_path}/BinaryBHLevel.cpp.txt
	cp ScalarRotatingCloud.hpp ${new_dir_path}/ScalarRotatingCloud.hpp.txt
	cd ${new_dir_path}
        # add the input params to the input files
        sed -i "s|DATADIR|${new_dir_path}|" ${new_dir_path}/params.txt
        sed -i "s|DATASUBDIR|${new_dir}|" ${new_dir_path}/params.txt
        sed -i "s|JOBNAME|${run}BBH|" ${new_dir_path}/slurm_submit
	#
	mkdir -p outputs
        cd outputs
	pwd
        sbatch ../slurm_submit
        #
	cd ${work_dir}
done
