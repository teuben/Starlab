#!/bin/bash

# Make sure that STARLAB_INSTALL_PATH has been defined: export
#STARLAB_INSTALL_PATH="/home/starlab/manybody/starlab/usr" 
#Better do
# this automatically, but then the .cshrc.bask file has to be
# defined. Better to rewrite this script to a tcsh.

source ~/.bashrc.starlab

if [ "$LAMRANK" ]; then
  rank=$LAMRANK
else
  rank=$MPIRUN_RANK
fi

case $rank in
  0) out=out.$rank
     err=err.$rank
     ;;
  *) out=/dev/null
     err=/dev/null
     ;;
esac

$STARLAB_INSTALL_PATH/bin/kira "$@" 

#On the LISA cluster this script looks like:
# #!/bin/bash
# nodes=$1
# ppn=$2
# n=$3
# t=$4
# inputfile=$PWD/universe.$n
# (( rnodes = nodes ))
# #if [[ nodes -lt 2 ]] ; then
# #  (( rnodes++ ))
# #fi
# (( procs = nodes * ppn ))
# cat << eof > tmpjob
# #PBS -lnodes=$rnodes:ppn=$ppn:infiniband -lwalltime=4:00:00 -N jobkira-d10-$1-$2-$3-$4
# date
# echo $1 $2 $3 $4
# export GMON_OUT_PREFIX=kira
# 
# module load gnu-mpich-ib
# #module load mpicopy
# #mpicopy $inputfile
# 
# workdir=$PWD/kiraworkdirmpi-d10-$n-$nodes-$ppn-$t.\$(date +%F-%R-%S | sed s/://)
# exe=$PWD/run_kira.sh
# mkdir -p \$workdir || exit
# cd \$workdir || exit
# time mpiexec -n $procs \$exe -d 10 -t $t -R $inputfile
# eof
# qsub tmpjob
