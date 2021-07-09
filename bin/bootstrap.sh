#!/bin/sh


NUMPROCS=1
LIMEX_TOL=0.001
PROBLEMID=footing2D_tri
SMOOTHER=vanka-ssc

function biot_exec {

CWD=$(pwd)

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua --limex-num-stages $2 --limex-tol $LIMEX_TOL --num-refs $3 --mg-smoother-type $SMOOTHER --problem-id $PROBLEMID
cd $CWD

}


NUMREFS=3
biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS
# biot_exec "vanka-nstages3-refs$NUMREFS" 3 $NUMREFS
# biot_exec "vanka-nstages4-refs$NUMREFS" 4 $NUMREFS

NUMREFS=4
biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS
# biot_exec "vanka-nstages3-refs$NUMREFS" 3 $NUMREFS
# biot_exec "vanka-nstages4-refs$NUMREFS" 4 $NUMREFS

NUMREFS=5
biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS