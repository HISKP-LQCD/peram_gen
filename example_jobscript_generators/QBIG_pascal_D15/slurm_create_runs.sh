#!/bin/bash

conf_start=200
conf_step=2
conf_end=220
conf_skip=( 612 )

rv_start=0
rv_end=3
nb_rnd=$(( $rv_end - $rv_start + 1 ))

seeds=( 9204344 3900937 8146422 5730210 8383664 9030138 )

flavour="strange_1150"

gdr=0
p2p=3

RUNDIR="/hiskp4/bartek/peram_generation/nf211/D15.48/${flavour}\/cnfg"
EVDIR="/hiskp4/eigensystems/nf211/D15.48"
GCONFBASE="/hiskp4/gauges/nf211/D15.48/conf"
EXEC="/hadron/bartek/build/bleeding_edge/pascal/peram_gen_multigpu.tmLQCD.quda_work.quda_develop-dynamic_clover/peram_gen"
JOBNAME="D15.48_${flavour}"
QUDA_RSC_PATH="/qbigwork/bartek/quda_resources/pascal_9c0e0dc8e96d9beb8de56a0e58a406cb486ce300_gdr${gdr}_p2p${p2p}"

for i in $( seq ${conf_start} ${conf_step} ${conf_end} ); do
  skip=0
  for skip_id in ${conf_skip[@]}; do
    if [ $skip_id -eq $i ]; then
      skip=1
    fi
  done
  if [ $skip -eq 1 ]; then
    continue
  fi

  echo "creating config $i"
  j=`printf %04d $i`

  wdir=cnfg${j}

  mkdir -p ${wdir}/
  mkdir -p ${wdir}/outputs
  
  ifile=${wdir}/invert.input
  
  cp templates/quda.invert.input ${ifile}
  sed -i "s@NSTORE@${i}@g" ${ifile}
  sed -i "s@GCONFBASE@${GCONFBASE}@g" ${ifile}
  
  laphin=LapH_${j}.in

  jscr=${wdir}/quda.job.slurm.${j}.cmd
  outfile="outputs/run_${j}.out"
  
  cp templates/quda.job.slurm.cmd ${jscr}
  sed -i "s@GDR@${gdr}@g" ${jscr}
  sed -i "s@P2P@${p2p}@g" ${jscr}
  sed -i "s@RUNDIR@${RUNDIR}${j}/@g" ${jscr}
  sed -i "s@JOBNAME@${JOBNAME}_${j}@g" ${jscr}
  sed -i "s@INFILE@${laphin}@g" ${jscr}
  sed -i "s@OUTFILE@${outfile}@g" ${jscr}
  sed -i "s@EXEC@${EXEC}@g" ${jscr}
  sed -i "s@QUDA_RSC_PATH@${QUDA_RSC_PATH}@g" ${jscr}
  
  cp templates/quda.LapH.part1.in ${wdir}/${laphin}
  
  echo "nb_rnd = ${nb_rnd}" >> ${wdir}/${laphin}
  
  for rv in $( seq ${rv_start} ${rv_end} ); do
    seed=${seeds[${rv}]}
    echo "id ${rv} seed ${seed}" >> ${wdir}/${laphin}
  done

  cat templates/quda.LapH.part3.in >> ${wdir}/${laphin}
  
  sed -i "s@NSTORE@${i}@g" ${wdir}/${laphin} 
  sed -i "s@EVDIR@${EVDIR}@g" ${wdir}/${laphin}

done
