#! /bin/bash -x

MY_PATH=`pwd`

cd $MY_PATH

N=2
NS=200
BSX=6.5
BSY=6.5
BSZ=6.5
BSXE=1.2
BSYE=1.2
BSZE=1.2
DX=-0.02
DY=-0.02
DZ=-0.02
ODIR=$MY_PATH/output/
IDIR=$MY_PATH/input/
BOXS=1
LOGS=1
RCUT=0.20
PBC="--pbc"
PBC=
dP=0.25
ddP=0
DEPTH=6
DT=0.001

SEED=1000

cd $IDIR
if [ ! -s schedule.config ]; then
    echo "Nth to do! Process $SLURM_ARRAY_TASK_ID terminates!"
    exit
fi

cd -


PREFIX=$N"_shells"


cd $ODIR
if [ -f $PREFIX."box.xyz" ] && [ -s $PREFIX."box.xyz" ] ; then
    BSX=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print X}'`
    BSY=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Y}'`
    BSZ=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Z}'`
fi
cd -


bf_software=$MY_PATH/../bin/elasticshells


time $bf_software -n $N --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.5 --ecc 100 --ecw 200 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        $PBC --tt rnd \
        -d --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --bend --restart --jam

