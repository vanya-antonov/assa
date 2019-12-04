
mkdir -p assa-0.0.2-source/ASSA
cp -v  ../INSTALL.txt ../seq1.fna ../seq2.fna       assa-0.0.2-source
cp -v  ~/_my/Programming/Perl/scripts/SADB/assa.pl  assa-0.0.2-source/ASSA
cp -rv ~/_my/Programming/Perl/lib/ASSA              assa-0.0.2-source/ASSA

zip -rv assa-0.0.2-source.zip assa-0.0.2-source
rm -rv assa-0.0.2-source

