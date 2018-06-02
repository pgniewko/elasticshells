#! /bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -o fem_50Vb.out
#SBATCH -e fem_50Vb.err
#SBATCH -J fem_50Vb
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH -A m2282
#SBATCH -a 52-52

MY_PATH=`pwd`

cd $MY_PATH

#module unload PrgEnv-intel/6.0.3
#module unload PrgEnv-intel/5.2.56
#module load PrgEnv-gnu

module load PrgEnv-intel/6.0.3
module load gsl

N=50
NS=200
BSX=20.0
BSY=20.0
BSZ=20.0
BSXE=6.0
BSYE=6.0
BSZE=6.0
DX=-0.02
DY=-0.02
DZ=-0.02
ODIR=$SCRATCH/FEM/051717/output/
ODIR=$SCRATCH/FEM/072117/output/
IDIR=./input_$SLURM_ARRAY_TASK_ID/
BOXS=1
LOGS=1
RCUT=0.1
PBC="--pbc"
dP=0.25
ddP=0
DEPTH=8
DT=0.001
model=fem

export OMP_NUM_THREADS=1

PREFIX="Vconst_"$N"_"$model"_bend_P"$dP_"IDX_"$SLURM_ARRAY_TASK_ID"_pbc"

SEED=$RANDOM

cd $ODIR
if [ -f $PREFIX."box.xyz" ] && [ -s $PREFIX."box.xyz" ] ; then
    BSX=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print X}'`
    BSY=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Y}'`
    BSZ=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Z}'`
fi
cd -


#bf_software=/global/homes/p/pawelg/bin/biofilm/bin/biofilm
bf_software=/global/homes/p/pawelg/bin/cori_software/biofilm/bin/biofilm

time $bf_software -n $N --int=fire --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.5 --ecc 100 --ecw 200 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        $PBC --model $model --tt rnd \
        -d --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --no-bend --jam

