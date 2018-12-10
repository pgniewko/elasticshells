#! /bin/bash -x

MY_PATH=`pwd`

cd $MY_PATH

N=3
NS=200
BSX=7.0
BSY=7.0
BSZ=7.0
BSXE=2.5
BSYE=2.5
BSZE=2.5
DX=-0.02
DY=-0.02
DZ=-0.02
ODIR=$MY_PATH/analyze/
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

PREFIX=$N"_shells"
SEED=1000

bf_software=$MY_PATH/../bin/elasticshells

time $bf_software -n $N --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.5 --ecc 100 --ecw 200 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        $PBC --tt rnd \
        -d --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --analyze
