#! /bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -o rfem_50Vb.out
#SBATCH -e rfem_50Vb.err
#SBATCH -J rfem_50Vb
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH -A m2282
#SBATCH -a 51-51

#for sid in 1 2 3 4 5 6 7 8 9 10 11 12 15 16 18 19 20 21 22 23 24 25 27 28 30 31 32 34 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50; do
#    if [ $sid -eq $SLURM_ARRAY_TASK_ID  ]; then
#        echo  "NTH MORE TO DO"
#        exit 1
#    fi
#done

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

cd $IDIR
if [ ! -s schedule.config ]; then
    echo "Nth to do! Process $SLURM_ARRAY_TASK_ID terminates!"
    exit
fi
cd -

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


bf_software=/global/homes/p/pawelg/bin/cori_software/biofilm/bin/biofilm

time $bf_software -n $N --int=fire --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.5 --ecc 100 --ecw 200 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        $PBC --model $model --tt rnd \
        -d --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --no-bend --jam --restart

