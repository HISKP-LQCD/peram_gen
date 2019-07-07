#!/bin/bash -x
#SBATCH --job-name=JOBNAME
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=bartosz_kostrzewa@fastmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem_bind=verbose
#SBATCH --time=06:00:00
#SBATCH --mem=666G
#SBATCH --gres=gpu:pascal:8
#SBATCH --partition=pascal

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/qbigwork2/bartek/libs/bleeding_edge/pascal/quda_develop-dynamic_clover/lib

p2p=P2P
gdr=GDR
rundir=RUNDIR
exe=EXEC
outfile=OUTFILE
infile=INFILE
export QUDA_RESOURCE_PATH=QUDA_RSC_PATH

if [ ! -d ${QUDA_RESOURCE_PATH} ]; then
  mkdir -p ${QUDA_RESOURCE_PATH}
fi

cd ${rundir}
date > ${outfile}
QUDA_RESOURCE_PATH=${QUDA_RESOURCE_PATH} OMP_NUM_THREADS=2 \
  QUDA_ENABLE_GDR=${gdr} QUDA_ENABLE_P2P=${p2p} QUDA_ENABLE_TUNING=1 \
  QUDA_ENABLE_DEVICE_MEMORY_POOL=0 \
  srun ${exe} -LapHsin ${infile} | tee -a ${outfile}

date >> ${outfile}

