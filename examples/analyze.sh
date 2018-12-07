#! /bin/bash -x

cd $MY_PATH

N=3
NS=200
BSX=7.0
BSY=7.0
BSZ=7.0
BSXE=1.5
BSYE=1.5
BSZE=15.
DX=-0.0
DY=-0.0
DZ=-0.0
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

PREFIX=$N"_shells"
SEED=1000

bf_software=$MY_PATH/../bin/elasticshells

time $bf_software -n $N --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.5 --ecc 100 --ecw 200 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX --out-dir $ODIR --in-dir $IDIR  \
        $PBC --tt rnd \
        -d --ir 2.5 --ir2 2.5 --seed $SEED \
        --vol-f --bend --jam --analyze
