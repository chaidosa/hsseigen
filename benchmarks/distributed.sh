#!/bin/bash
#SBATCH -N 26
#SBATCH --ntasks-per-node=48
#SBATCH --time=01:00:00
#SBATCH --job-name=test_mpi
#SBATCH --error=job.tri.err.%j
#SBATCH --output=job.tri.out.%j
#SBATCH --partition=medium



spack load gcc@12.2.0%gcc@=11.2.0

HSS_HOME=/home/ext/extnikhil/abhishek/hsseigen

cd /home/ext/extnikhil/abhishek/hsseigen
BLAS_PATH=/home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/openblas-0.3.20-mooia7kq4d24w5ew7qlog6fst4j5o4zx
export LD_LIBRARY_PATH=$BLAS_PATH/lib:$LD_LIBRARY_PATH

cd /scratch/extnikhil

hostlist=`scontrol show hostnames $SLURM_JOB_NODELIST`
hosts=($hostlist)

touch hostfile_test_tri
true > hostfile_test_tri

for i in "${hosts[@]}"
do
echo $i >> hostfile_test_tri
done

cd /scratch/extnikhil

echo "65k"

Cutoff=(1 2 3 4 5 6 7 8 9 10)

for t in ${allThreads[@]}; do
    $HSS_HOME/Test /home/ext/extnikhil/abhishek/triad_65k.txt 66560 64 2 1 96 $t 
done
