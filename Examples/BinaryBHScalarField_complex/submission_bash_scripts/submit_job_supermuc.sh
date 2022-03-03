#!/bin/bash
#
# This script should make a new directory in ~/GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for the KNL nodes

work_dir=/dss/dsshome1/04/di76bej/GRChombo/GRChombo/Examples/BinaryBHScalarField
cd $work_dir
data_directory=/hppfs/work/pn34tu/di76bej/GRChombo_data/BinaryBHScalarField

# specify the input params for each run I want to submit
# list for each is: mu, delay, dt, G

run0001=(1 0 0.125 0)
run0002=(0.08187607564 0 0.25 0)

params_file=params.txt

run_list=(
	run0001
)

params_file=params.txt
plot_interval=10
L=512
N1=64
box_size=16

for run in "${run_list[@]}"
do
  	cd $work_dir
        # extract parameters    
        val="$run[0]"; mu="${!val}"
        val="$run[1]"; delay="${!val}"
        val="$run[2]"; dt_mult="${!val}"
        val="$run[3]"; G="${!val}"

        # text_number=$(printf "%04d" ${run_number})
        new_dir=${run}_mu${mu}_delay${delay}_G${G}_ratio1
        #_L${L}_N$N1
        echo ${new_dir}
        new_dir_path=${data_directory}/${new_dir}
        #
	mkdir -p ${new_dir_path}
        cp slurm_submit_cosma ${new_dir_path}/slurm_submit
        cp ${params_file} ${new_dir_path}/params.txt

	cd ${new_dir_path}
        # add the input params to the input files
        sed -i "s|DATADIR|${new_dir_path}|" ${new_dir_path}/params.txt
        sed -i "s|DATASUBDIR|${new_dir}|" ${new_dir_path}/params.txt
        sed -i "s|JOBNAME|${run}KS|" ${new_dir_path}/slurm_submit
        sed -i "s|BOXLENGTH|${L}|" ${new_dir_path}/params.txt
        sed -i "s|BOXSIZE|${box_size}|" ${new_dir_path}/params.txt
        sed -i "s|CENTERX|$(($L/2))|" ${new_dir_path}/params.txt
        sed -i "s|CENTERY|$(($L/2))|" ${new_dir_path}/params.txt
        sed -i "s|MUVAL|${mu}|" ${new_dir_path}/params.txt
        sed -i "s|DELAYTIME|${delay}|" ${new_dir_path}/params.txt
        sed -i "s|DTMULT|${dt}|" ${new_dir_path}/params.txt
        sed -i "s|NBASIC|${N1}|" ${new_dir_path}/params.txt
	sed -i "s|NSPACE3|$(($N1/2))|" ${new_dir_path}/params.txt
        sed -i "s|PLOTINTERVAL|${plot_interval}|" ${new_dir_path}/params.txt
	#
	mkdir -p outputs
        cd outputs
        sbatch ../slurm_submit
        #
	cd ${work_dir}
done
