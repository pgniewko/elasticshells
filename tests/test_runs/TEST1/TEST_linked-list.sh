#! /bin/bash -x

N=8
NS=140
BSX=9.0
BSY=9.0
BSZ=9.0
BSXE=4.6
BSYE=4.6
BSZE=4.6
DX=-0.05
DY=-0.05
DZ=-0.05
ODIR=./output/
IDIR=./input/
BOXS=1
SS=2
LOGS=1
RCUT=0.4
PBC=
dP=10
ddP=0
DEPTH=5
DT=0.001
model=fem

export OMP_NUM_THREADS=1

cd $ODIR
RESULT=$?
if [ $RESULT -eq 0 ]; then
   rm *
   cd -
fi

simid=1
PREFIX1="LINKED_LIST_TEST"

time biofilm -n $N --int=cp --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.0 --ecc 1500 --ecw 2000 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --log-step $LOGS --save-step $SS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX1 --out-dir $ODIR --in-dir $IDIR  \
        $PBC --model $model \
        -d --ir 2.5 --ir2 2.5 --seed $simid --no-bend

PREFIX2="ON2_TEST"
time biofilm -n $N --int=cp --depth $DEPTH --ns $NS --dt $DT --dp $dP --ddp=$ddP  --th 0.1 --nu 0.0 --ecc 1500 --ecw 2000 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 0 --log-step $LOGS --save-step $SS --box-step $BOXS \
        --rv $RCUT --prefix $PREFIX2 --out-dir $ODIR --in-dir $IDIR  \
        $PBC --model $model \
        -d --ir 2.5 --ir2 2.5 --seed $simid --no-bend

cd $ODIR
echo "TRAJECTORY DIFFERENCES"
diff $PREFIX1".xyz"  $PREFIX2".xyz"

echo "OBSERVABLES DIFFERENCES"
diff $PREFIX1".out"  $PREFIX2".out"
