#!/bin/bash
#
# $Id$
#
# Prepare (assemble) the distribution
# 

VER='1.00'
OUT="assa-${VER}"

mkdir -p  $OUT/src  $OUT/examples
cp -v  INSTALL.txt                                  $OUT
cp -v  ../data/*.fna                                $OUT/examples
cp -v  ~/_my/Programming/Perl/scripts/SADB/assa.pl  $OUT/src
cp -rv ~/_my/Programming/Perl/lib/ASSA              $OUT/src

#rm *.tar.gz
tar -czf $OUT.tar.gz  $OUT
rm -r $OUT

