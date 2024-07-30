#!/bin/sh
#BSUB -q scafellpikeSKL
#BSUB -o Logs/Simlogs/Simlog_%J_%I.test.log
#BSUB -e Logs/Errlogs/SimErrs_%J_%I.test.err
#BSUB -R "rusage[mem=10000] span[ptile=1]"
#BSUB -n 1
#BSUB -W 00:10
#BSUB -J "LibMeshAMSJob"


export WORK_DIREC="/lustre/scafellpike/local/HT04544/sht06/ssr82-sht06/SingularityCont"

module load use.dev
module load singularity/3.10.5

singularity exec --bind $WORK_DIREC/:/opt/scratch $WORK_DIREC/MOOSE_APP.simg bash 
