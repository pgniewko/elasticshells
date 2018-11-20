#! /bin/bash -x

MY_PATH=`pwd`

cd $MY_PATH

N=10
NS=200
BSX=20.0
BSY=20.0
BSZ=20.0
BSXE=1.5
BSYE=1.5
BSZE=1.5
DX=-0.02
DY=-0.02
DZ=-0.02
ODIR=$MY_PATH/output/
IDIR=$MY_PATH/input/
BOXS=1
LOGS=1
RCUT=0.20
PBC="--pbc"
dP=0.25
ddP=0
DEPTH=6
DT=0.001
model=fem

PREFIX=$N"_shells"

SEED=$RANDOM

cd $ODIR
if [ -f $PREFIX."box.xyz" ] && [ -s $PREFIX."box.xyz" ] ; then
    BSX=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print X}'`
    BSY=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Y}'`
    BSZ=`tail -n 1 $PREFIX."box.xyz" | awk '{X=$1; Y=$2; Z=$3} END{print Z}'`
fi
cd -


bf_software=$MY_PATH/../bin/elasticshells


time $bf_software -n $N --int=fire --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.5 --ecc 100 --ecw 200 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        $PBC --model $model --tt rnd \
        -d --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --no-bend --jam

