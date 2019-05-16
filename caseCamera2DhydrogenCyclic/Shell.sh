#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# hyperthreading off
#SBATCH --threads-per-core=1
 

# set max wallclock time
#SBATCH --time=6-0

# set name of job
#SBATCH --job-name=RDE

# set queue name
#SBATCH -p hydro

# run the application

export MPI_ROOT=/opt/software/intel/2017/compilers_and_libraries_2017.4.196/linux/mpi
source /opt/software/applied/hydro/OpenFOAM-6/etc/bashrc

mv ./system/setFieldsDict ./system/setFieldsDict-temp
mv ./system/setFieldsDict0 ./system/setFieldsDict
setFields
blockMesh
mv ./system/setFieldsDict ./system/setFieldsDict0
mv ./system/setFieldsDict-temp ./system/setFieldsDict
setFields
decomposePar -force
mpirun -np 32 RDE -parallel
reconstructPar
rm -r processor*