#!/bin/bash
#SBATCH -N 26
#SBATCH --ntasks-per-node=48
#SBATCH --time=01:00:00
#SBATCH --job-name=test_mpi
#SBATCH --error=job.tri.err.%j
#SBATCH --output=job.tri.out.%j
#SBATCH --partition=medium



# spack load gcc@12.2.0%gcc@=11.2.0

HSS_HOME=/home/abhishek/hsseigen

cd $HSS_HOME
BLAS_PATH=/home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/openblas-0.3.20-mooia7kq4d24w5ew7qlog6fst4j5o4zx
export LD_LIBRARY_PATH=$BLAS_PATH/lib:$LD_LIBRARY_PATH

# cd /scratch/extnikhil

# hostlist=`scontrol show hostnames $SLURM_JOB_NODELIST`
# hosts=($hostlist)

# touch hostfile_test_tri
# true > hostfile_test_tri

# for i in "${hosts[@]}"
# do
# echo $i >> hostfile_test_tri
# done

# cd /scratch/extnikhil

echo "Tree height experiments:"
Cutoff=(5 6 7 8 9 10)
repetition=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
typeMat=(1 2) # 1 is tri 2 is banded 3 is ex4
mat2hss=2
band=1
blksize=64
for mat in "${typeMat[@]}"; do
	if [ $mat -eq 1 ]; then
		echo "Tridiagonal: "
		mat2hss=2
		band=1
		blksize=64
	fi		
	
	if [ $mat -eq 2 ]; then
		echo "Banded: "
		mat2hss=2
		band=5
		blksize=1024
	fi	
	
	if [ $mat -eq 3 ]; then
		echo "Ex4: "
		mat2hss=1
		band=1
		blksize=1024
	fi
	
	for t in ${Cutoff[@]}; do
    		echo "Cutoff $t" 
    		for rep in ${repetition[@]}; do
    			# $HSS_HOME/TestUntied /home/ext/extnikhil/abhishek/triad_65k.txt 66560 $blksize $mat2hss $band 96 $t 
				echo 5
    		done
	done
	echo " "
done

echo " "

echo "Comparision experiments:"

parallelType=(1 2 3) # 1 - Task ; 2 - Task untied ; 3 - Cilk version
prll=Task

typeMat=(1 2) # 1 is tri 2 is banded 3 is ex4
mat2hss=2
band=1
t=10

thrds=(8 16 32 64 96)

repetition=(1 2 3 4 5) 

for prlltype in "${parallelType[@]}"; do

	if [ $prlltype -eq 1 ]; then
		echo "Parallel Task tied: "
		prll=Task
	fi	

	if [ $prlltype -eq 2 ]; then
		prll=Untied
		echo "Parallel task untied: "
	fi	

	if [ $prlltype -eq 3 ]; then
		prll=Cilk
		echo "Parallel Cilk: "
	fi	

for mat in "${typeMat[@]}"; do

	if [ $mat -eq 1 ]; then
		echo "Tridiagonal: "
		mat2hss=2
		band=1
		blksize=64
	fi		
	
	if [ $mat -eq 2 ]; then
		echo "Banded: "
		mat2hss=2
		band=5
		blksize=1024
	fi	
	
	if [ $mat -eq 3 ]; then
		echo "Ex4: "
		mat2hss=1
		band=1
		blksize=1024
	fi	
for thrd in "${thrds[@]}"; do
	echo "N threads = $thrd"
	export CILK_NWORKERS=$thrd
    	for rep in ${repetition[@]}; do
    		$HSS_HOME/Test$prll /home/ext/extnikhil/abhishek/triad_65k.txt 66560 $blksize $mat2hss $band $thrd $t
    	done

	echo " "
done # thrd
	echo " "
done # mat
	echo " "
done # prlltype

echo " "
