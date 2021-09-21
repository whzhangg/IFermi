#!/bin/sh
#QSUB2 queue qS
#QSUB2 core 48
#QSUB2 mpi 48
#QSUB2 smp 1
#QSUB2 wtime 2:00:00
#PBS -N whzjob
cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh

module load comp1
module load intel-mkl/20.0.0
module load intel-mpi/2019.6
#module load fftw/2.1.5
#module load intelpython/3_2020.0

export PATH=$PATH:/home/whzhang/work/packages/qe-6.4.1/bin
export PATH=$PATH:/home/whzhang/work/packages/my_bin

mpirun pw.x -npool 12 -in scf.in > scf.out
mpirun pw.x -npool 12 -in nscf.in > nscf.out

