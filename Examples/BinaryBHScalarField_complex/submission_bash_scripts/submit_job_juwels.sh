#!/bin/bash
#
# This script should make a new directory in GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for COSMA6

work_dir=/p/home/jusers/bamber1/juwels/GRChombo/Examples/BinaryBHScalarField
cd $work_dir
data_directory=/p/project/pra116/bamber1/BinaryBHScalarField

# specify the input params for each run I want to submit
# list for each is: mu, delay, dt, G, BH mass ratio, l, m, Al

run0011=(1 0 0.0625 0 1 0 0 0)
run0012=(1 10000 0.0625 0 1 0 0 0)
run0013=(0.08187607564 0 0.25 0 1 0 0 0)
run0014=(1 0 0.0625 0 2 0 0 0)
run0015=(1 10000 0.0625 0 2 0 0 0)
run0016=(0.5 0 0.25 0 1 0 0 0)
run0016ICS=(0.5 0 0.25 0 1 0 0 0)
run0017=(0.5 10000 0.25 0 1 0 0 0)
run0018=(0.5 0 0.25 0.000001 1 0 0 0) # G = 10^{-6}
run0019=(1 0 0.0625 0 1 0 0 0) # resume from stationary BH distribution
run0020=(0.5 0 0.125 0 1 0 0 0) # resume from stationary BH distribution
run0021=(0.5 0 0.25 0.01 1 0 0 0)
run0022=(0.5 0 0.25 0 1 0 0 0) # this is a test, should be the same as 16
run0023=(0.5 0 0.25 0.0000000001 1 0 0 0) # G = 10^{-10}
run0023ICS=(0.5 0 0.25 0.0000000001 1 0 0 0) # G = 10^{-10}
run0024=(0.5 0 0.125 0.00000000000000000001 1 0 0 0) # G = 10^{-20}, restart as in run0020
run0025=(0.5 0 0.25 0 1 1 1 0)
run0026=(0.5 0 0.25 0.000000000000001 1 0 0 0) # G = 10^{-15}
run0027=(0.5 0 0.25 0.00000000000000000001 1 0 0 0) # G = 10^{-20}
run0028=(0.5 0 0.25 0.00000001 1 0 0 0) # G = 10^{-8}
run0029=(1 0 0.0625 0.0000000001 1 0 0 0) # G = 10^{-10} 
run0030=(0.5 0 0.25 0 1 2 2 0)
run0031=(0.5 0 0.25 0.0000000000000000000000001 1 0 0 0) # G = 10^{-25}
run0032=(0.5 0 0.25 0.000000000000000000000000000001 1 0 0 0) # G = 10^{-30} 
run0033=(0.5 0 0.25 0.000000001 1 0 0 0) # G = 10^{-9}
run0034=(0.5 0 0.25 0.000000000001 1 0 0 0) # G = 10^{-12}
run0035=(0.5 0 0.25 0.00000000000001 1 0 0 0) # G = 10^{-14}
run0035ICS=(0.5 0 0.25 0.00000000000001 1 0 0 0) # G = 10^{-14}
run0036=(0.5 0 0.25 0 1 1 -1 0)

run_Katy_test_v3=(1 0 0.0625 0.0000000001 1 0 0 0) # G = 10^{-10}

run_list=(
#    run0011
#    run0012
#    run0013
#    run0014
#    run0015
#    run0016ICS
#    run0017
#    run0018
#     run0019
#    run0020
#    run0025
    run0023ICS
#    run_Katy_test_v3
#    run0024
#    run0026
#    run0027
#    run0028
#    run0029
#    run0033
#   run0035ICS
#    run0036
)

plot_interval=1
L=512
N1=64
box_size=16
reflect_z=0

restart_hash=""
restart_num="000000"

echo ${restart_hash}
echo ${restart_num}

for run in "${run_list[@]}"
do
  	cd $work_dir
        # extract parameters    
        val="$run[0]"; mu="${!val}"
        val="$run[1]"; delay="${!val}"
        val="$run[2]"; dt_mult="${!val}"
        val="$run[3]"; G="${!val}"
        val="$run[4]"; ratio="${!val}"
	val="$run[5]"; l="${!val}"
	val="$run[6]"; m="${!val}"
	val="$run[7]"; Al="${!val}"
	
        # text_number=$(printf "%04d" ${run_number})
	if (( $l == 0 ))
	then
	   new_dir=${run}_mu${mu}_delay${delay}_G${G}_ratio${ratio}
	else
	   new_dir=${run}_mu${mu}_delay${delay}_G${G}_ratio${ratio}_l${l}_m${m}_Al${Al}
        fi
	   
        #_L${L}_N$N1
        echo ${new_dir}
        new_dir_path=${data_directory}/${new_dir}
	#
	mkdir -p ${new_dir_path}
	# get latest checkpoint file
        cd ${new_dir_path}
	chk_file=BinaryBHSFChk_${restart_num}.3d.hdf5
	#for chk in BinaryBHSFChk_*.hdf5
	#do
	#    chk_file=$chk
	#done
	echo $chk_file
	#
	cd ${work_dir}
        cp slurm_submit_juwels ${new_dir_path}/slurm_submit

	params_file=params_ratio${ratio}.txt
	#params_file=params_evolution_Katy_test_v3.txt
	
	cp ${params_file} ${new_dir_path}/params.txt
	cp BinaryBHLevel.cpp ${new_dir_path}/BinaryBHLevel.cpp.txt
	cd ${new_dir_path}
        # add the input params to the input files
        sed -i "s|DATADIR|${new_dir_path}|" ${new_dir_path}/params.txt
        sed -i "s|DATASUBDIR|${new_dir}|" ${new_dir_path}/params.txt
        sed -i "s|DTMULT|${dt_mult}|" ${new_dir_path}/params.txt
	sed -i "s|GVALUE|${G}|" ${new_dir_path}/params.txt
        sed -i "s|JOBNAME|${run}BBH|" ${new_dir_path}/slurm_submit
        sed -i "s|BOXLENGTH|${L}|" ${new_dir_path}/params.txt
        sed -i "s|BOXSIZE|${box_size}|" ${new_dir_path}/params.txt
        sed -i "s|CENTERX|$(($L/2))|" ${new_dir_path}/params.txt
        sed -i "s|CENTERY|$(($L/2))|" ${new_dir_path}/params.txt
        sed -i "s|MUVAL|${mu}|" ${new_dir_path}/params.txt
	sed -i "s|SCALARL|${l}|" ${new_dir_path}/params.txt
	sed -i "s|SCALARM|${m}|" ${new_dir_path}/params.txt
	sed -i "s|ALANGLE|${Al}|" ${new_dir_path}/params.txt
        sed -i "s|DELAYTIME|${delay}|" ${new_dir_path}/params.txt
        sed -i "s|NBASIC|${N1}|" ${new_dir_path}/params.txt
	sed -i "s|RESTARTHASH|${restart_hash}|" ${new_dir_path}/params.txt
	sed -i "s|RESTARTNUM|${restart_num}|" ${new_dir_path}/params.txt
	sed -i "s|CHKFILE|${chk_file}|" ${new_dir_path}/params.txt
	sed -i "s|PLOTINTERVAL|${plot_interval}|" ${new_dir_path}/params.txt
	#
	echo "reflect_z = ${reflect_z}"
	if (( $reflect_z==1 ))
	    then
		echo "reflect_z=1"
		sed -i "s|NSPACE3|$((${N1}/2))|" ${new_dir_path}/params.txt
		sed -i "s|CENTERZ|0|" ${new_dir_path}/params.txt
		sed -i "s|ZBOUND|2|" ${new_dir_path}/params.txt
	else
	        echo "reflect_z=0"
		sed -i "s|NSPACE3|${N1}|" ${new_dir_path}/params.txt
                sed -i "s|CENTERZ|$(($L/2))|" ${new_dir_path}/params.txt
                sed -i "s|ZBOUND|4|" ${new_dir_path}/params.txt
	fi
	#
	mkdir -p outputs
        cd outputs
        sbatch ../slurm_submit
        #
	cd ${work_dir}
done
