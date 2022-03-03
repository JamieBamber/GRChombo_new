#!/bin/bash
#
# This script should make a new directory in GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# this copy is for COSMA6

work_dir=/home/dc-bamb1/GRChombo/Examples/BinaryBHScalarField_complex
cd $work_dir
data_directory=/scratch/dp092/dc-bamb1/GRChombo_data/BinaryBHSF

# specify the input params for each run I want to submit
# list for each is: mu, dt_mult, G, l, m

run0001=(0.34 0.1 '1e-10' 0 0 homogeneous_noICS) # 10^{-10}
run0002=(0.34 0.1 '1e-12' 0 0 homogeneous_noICS) # 10^{-12}
run0003=(0.34 0.1 '1e-8' 0 0 homogeneous_noICS) # 10^{-8}

run0004=(0.34 0.1 '1e-10' 0 0 homogeneous_withICS) # 10^{-10}
run0005=(0.34 0.1 '1e-12' 0 0 homogeneous_withICS) # 10^{-12}
run0006=(0.34 0.1 '1e-8' 0 0 homogeneous_withICS) # 10^{-8}

run0021=(0.34 0.1 '1e-10' 0 0 gaussian_kappa0.0125_noICS) # 10^{-10}
run0022=(0.34 0.1 '1e-12' 0 0 gaussian_kappa0.0125_noICS) # 10^{-12}
run0023=(0.34 0.1 '1e-8' 0 0 gaussian_kappa0.0125_noICS) # 10^{-8}

run0024=(0.34 0.1 '1e-10' 0 0 gaussian_kappa0.0125_withICS) # 10^{-10}
run0025=(0.34 0.1 '1e-12' 0 0 gaussian_kappa0.0125_withICS) # 10^{-12}
run0026=(0.34 0.1 '1e-8' 0 0 gaussian_kappa0.0125_withICS) # 10^{-8}

run0051=(0.34 0.1 '1e-10' 0 0 Newtonian_run0031_wslope0.25_wradius50_n001500_ICS) # 10^{-10}
run0052=(0.34 0.1 '1e-12' 0 0 Newtonian_run0031_wslope0.25_wradius50_n001500_ICS) # 10^{-12}
run0053=(0.34 0.1 '1e-14' 0 0 Newtonian_run0031_wslope0.25_wradius50_n001500_ICS) # 10^{-8}

run0061=(0.34 0.1 '1e-10' 0 0 Newtonian_run0031_nowindow_n001500_ICS) # 10^{-10}
run0062=(0.34 0.1 '1e-12' 0 0 Newtonian_run0031_nowindow_n001500_ICS) # 10^{-12}
run0063=(0.34 0.1 '1e-14' 0 0 Newtonian_run0031_nowindow_n001500_ICS) # 10^{-8}

run0071=(0.34 0.1 '1e-10' 0 0 Newtonian_run0031_wslope0.25_wradius50_n001500_noICS) # 10^{-10}
run0072=(0.34 0.1 '1e-12' 0 0 Newtonian_run0031_wslope0.25_wradius50_n001500_noICS) # 10^{-12}
run0073=(0.34 0.1 '1e-14' 0 0 Newtonian_run0031_wslope0.25_wradius50_n001500_noICS) # 10^{-8}

run0036=(0.5 0.125 0 0 0)

run00Uli=(0.01 0.25 '1' 0 0 6orbit_parameters_no_SF)

run_list=(
	run00Uli
)

plot_interval=10
L=512
N1=64
box_size=8
reflect_z=1

delay=0
Al=0

restart_hash="#"
restart_num="000100"

echo ${restart_hash}
echo ${restart_num}

scalar_kappa=0.0
scalar_amplitude=0.0

for run in "${run_list[@]}"
do
  	cd $work_dir
        # extract parameters    
        val="$run[0]"; mu="${!val}"
        val="$run[1]"; dt_mult="${!val}"
        val="$run[2]"; G="${!val}"
	val="$run[3]"; l="${!val}"
	val="$run[4]"; m="${!val}"
	val="$run[5]"; suffix="${!val}"
	
	new_dir=${run}_mu${mu}_l${l}_m${m}_G${G}_${suffix}

	#_gaussian_kappa0.0125_noICS
	#gaussian_kappa0.0125_noICS
	   
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
        cp slurm_files/pbs_submit_dial3 ${new_dir_path}/slurm_submit

	#params_file=params_ratio${ratio}.txt
	params_file=params_ratio1_Uli.txt
	
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
	sed -i "s|SAMPLITUDE|${scalar_amplitude}|" ${new_dir_path}/params.txt
        sed -i "s|SKAPPA|${scalar_kappa}|" ${new_dir_path}/params.txt
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
