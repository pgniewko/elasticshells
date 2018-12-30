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
ODIR=$MY_PATH/analyze/
IDIR=$MY_PATH/input/
RCUT=0.20
PBC="--pbc"
dP=0.25
ddP=0
DT=0.001

PREFIX=$N"_shells"
SEED=1001

software=$MY_PATH/../bin/elasticshells

time $software -n $N --ns $NS --dt $DT --dp $dP --ddp=$ddP --th 0.1 --nu 0.5 --E-shell 100 --E-box 200 \
        --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --tt rnd -d --rv $RCUT --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --analyze $PBC

