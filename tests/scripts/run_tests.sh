#!/bin/bash -x

N=25
NS=30000
BSX=9
BSY=9
BSZ=9
BSXE=4.5
BSYE=4.5
BSZE=4.5
DX=-0.01
DY=-0.01
DZ=-0.01
DIR=./data
BOXS=50
SS=50
VS=10
VR=3.0
RCUT=0.3


PROJECT_PATH=/home/pablo/Science/Biofilm
cd $PROJECT_PATH
make clean
make
cd -
cp ../../bin/biofilm ./bin/

PBC=

# VERLET-LIST
./bin/biofilm -n $N --int=hm --depth 5 --ns $NS --dt 0.001 --dp 10.0 -k 50 -a 100 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 1 --save-step $SS --box-step $BOXS \
        --rv $RCUT --vlist-step $VS -verlet-r $VR \
        -r $DIR/render_1b.py -t $DIR/traj_1b.xyz --surf $DIR/surf_1b.py -o $DIR/biofilm_1b.out \
        $PBC \
        -d

# LINKED-CELLS
./bin/biofilm -n $N --int=hm --depth 5 --ns $NS --dt 0.001 --dp 10.0 -k 50 -a 100 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --save-step $SS --box-step $BOXS \
        --rv $RCUT \
        -r $DIR/render_2b.py -t $DIR/traj_2b.xyz --surf $DIR/surf_2b.py -o $DIR/biofilm_2b.out \
        $PBC \
        -d

# NAIVE O(N^2)
./bin/biofilm -n $N --int=hm --depth 5 --ns $NS --dt 0.001 --dp 10.0 -k 50 -a 100 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 0 --save-step $SS --box-step $BOXS \
        --rv $RCUT \
        -r $DIR/render_0b.py -t $DIR/traj_0b.xyz --surf $DIR/surf_0b.py -o $DIR/biofilm_0b.out \
        $PBC \
        -d

PBC="--pbc"

# VERLET-LIST
./bin/biofilm -n $N --int=hm --depth 5 --ns $NS --dt 0.001 --dp 10.0 -k 50 -a 100 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 1 --save-step $SS --box-step $BOXS \
        --rv $RCUT --vlist-step $VS -verlet-r $VR \
        -r $DIR/render_1p.py -t $DIR/traj_1p.xyz --surf $DIR/surf_1p.py -o $DIR/biofilm_1p.out \
        $PBC \
        -d

# LINKED-CELLS
./bin/biofilm -n $N --int=hm --depth 5 --ns $NS --dt 0.001 --dp 10.0 -k 50 -a 100 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 2 --save-step $SS --box-step $BOXS \
        --rv $RCUT \
        -r $DIR/render_2p.py -t $DIR/traj_2p.xyz --surf $DIR/surf_2p.py -o $DIR/biofilm_2p.out \
        $PBC \
        -d

# NAIVE O(N^2)
./bin/biofilm -n $N --int=hm --depth 5 --ns $NS --dt 0.001 --dp 10.0 -k 50 -a 100 \
        --bsx $BSX --bsy $BSY --bsz $BSZ --bsdx $DX --bsdy $DY --bsdz $DZ --bsxe $BSXE --bsye $BSYE --bsze $BSZE \
        --nb 0 --save-step $SS --box-step $BOXS \
        --rv $RCUT \
        -r $DIR/render_0p.py -t $DIR/traj_0p.xyz --surf $DIR/surf_0p.py -o $DIR/biofilm_0p.out \
        $PBC \
        -d
