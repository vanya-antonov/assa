
mkdir assa-0.0.1-source
cd assa-0.0.1-source
cp -v  ~/_my/Programming/Perl/scripts/SADB/README.txt .
mkdir ASSA
cp -v  ~/_my/Programming/Perl/scripts/SADB/assa.pl    ASSA
cp -rv ~/_my/Programming/Perl/lib/ASSA                ASSA

cd ..
zip -rv assa-0.0.1-source.zip assa-0.0.1-source
rm -rv assa-0.0.1-source

