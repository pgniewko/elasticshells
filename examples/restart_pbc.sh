#! /bin/bash -x

MY_PATH=`pwd`

cd $MY_PATH

N=10
NS=200
BSX=20.0
BSY=20.0
BSZ=20.0
BSXE=5.0
BSYE=5.0
BSZE=5.0
ODIR=$MY_PATH/output/
IDIR=$MY_PATH/input/
RCUT=0.20
PBC="--pbc"
dP=0.25
ddP=0
DT=0.001

PREFIX=$N"_shells"
SEED=1001

cd $IDIR
if [ ! -s schedule.config ]; then
    echo "Nth to do! Process $SLURM_ARRAY_TASK_ID terminates!"
    exit
fi
cd -

cd $ODIR
if [ -f $PREFIX."box.xyz" ] && [ -s $PREFIX."box.xyz" ] ; then
    BSX=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print X}'`
    BSY=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Y}'`
    BSZ=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Z}'`
fi
cd -

software=$MY_PATH/../bin/elasticshells

time $software -n $N --ns $NS --dt $DT --dp $dP --ddp=$ddP --th 0.1 --nu 0.5 --E-shell 100 --E-box 200 \
        --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --tt rnd -d --rv $RCUT --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --restart --jam $PBC 

