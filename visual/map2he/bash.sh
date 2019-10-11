#!/usr/bin/bash
imgdir="/home/alma/Documents/PhD/papers/STSC/rsc/mob_data/resized_imgs"
tdir="/home/alma/Documents/PhD/papers/STSC/rsc/mob_data/tmats"
wdir="/home/alma/Documents/PhD/papers/STSC/res/molb/try_1"
odir="/home/alma/Documents/PhD/papers/STSC/res/molb/try_1/he_overlay"


for ii in {9..12}; do
    ./map2he.py -i $imgdir/Rep${ii}_MOB.rz10.jpg -p $wdir/Rep${ii}_*/W.*.tsv -sf 0.1 -t $tdir/Rep${ii}_MOB_transformation.txt -o $odir -si
done
