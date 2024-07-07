#!/bin/bash

echo COMPILING PROGRAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meson setup builddir --buildtype=release
cd builddir
meson compile
cd ..
mkdir mackay_gold_clusters
echo GENERATING GOLD CLUSTERS WITH R_0=2.88428856046 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i in {1..15}
do
   ./builddir/milestones/mackay/mackay 1 2.88428856046 > mackay_gold_clusters/cluster_1.xyz
    echo "Cluster with $i layers has $(head -1 mackay_gold_clusters/cluster_$i.xyz) atoms"
done
