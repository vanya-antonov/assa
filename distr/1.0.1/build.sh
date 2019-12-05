#!/bin/bash
#
# Prepare (assemble) the distribution
# 

VER='1.0.1'
OUT="assa-${VER}"

mkdir -p  $OUT/src  $OUT/examples
cp -v  INSTALL.txt      $OUT
cp -v  ../data/*.fna    $OUT/examples
cp -rv ../../src        $OUT

#rm *.tar.gz
tar -czf $OUT.tar.gz  $OUT
zip  -r  $OUT.zip     $OUT
rm -r $OUT

