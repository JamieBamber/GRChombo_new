#!/bin/bash
#
# This script should make a new directory in ~/GRChombo_data for the data, and copy the slurm_submit file, params.txt file and output files to
# that directory, submit the job, then change back to the working directory

# load modules 
module load slurm_setup
module load intel mkl hdf5 gcc/7
module load valgrind/3.15.0-impi 

work_dir=/dss/dsshome1/04/di76bej/GRChombo/GRChombo/Examples/BinaryBHScalarField
cd $work_dir
#data_directory=/hppfs/work/pn34tu/di76bej/GRChombo_data/BinaryBHScalarField

run_number=1
suffix=valgrind

params_file=params_test.txt

new_dir=test_${suffix}
echo ${new_dir}
new_dir_path=${work_dir}/${new_dir}
#
mkdir -p ${new_dir_path}
cp ${params_file} ${new_dir_path}/params.txt

cd ${new_dir_path}
# add the location of the new directory to the input files
sed -i "s|DATADIR|${new_dir_path}|" ${new_dir_path}/params.txt
# 
#! Full path to application executable: 
application="/dss/dsshome1/04/di76bej/GRChombo/GRChombo/Examples/BinaryBHScalarField/Main_BinaryBH3d.Linux.64.mpiicpc.ifort.DEBUG.OPT.MPI.ex"

#! Run options for the application:
currentdir=$(pwd)
options="${currentdir}/params.txt"

#Run the program:
valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=valgrind-out.txt $application $options
#
cd ${work_dir}

