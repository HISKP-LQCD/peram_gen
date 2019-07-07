#!/bin/bash

conf_start=200
conf_step=2
conf_end=220

for i in $( seq ${conf_start} ${conf_step} ${conf_end} ); do
  j=`printf %04d $i`
  wdir=cnfg${j}
  pushd ${wdir}
  jscr=quda.job.slurm.${j}.cmd
  sbatch ${jscr}
  popd
done
