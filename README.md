# Introduction
ASSA is a pipeline for large scale RNA-RNA interaction prediction that utilizes
the thermodymanics-based tool RNAup for the hybridization energy calcluation.

# Installation
Example of ASSA installation in Google Colab is available here: https://colab.research.google.com/drive/1Gem_i761us0XvfGLpFs9nZYSnaLLCUUV?usp=sharing

## Prerequisites
To use ASSA, you will need to install:
* Perl
* R: https://www.r-project.org/
* The 'Statistics::R' Perl module and (optionally) Parallel::ForkManager to enable parallel computing. On Ubuntu you can use `apt-get install libstatistics-r-perl` and `apt-get install libparallel-forkmanager-perl`, respectively
* The 'lastdb' and 'lastal' tools that are parts of the LAST package: http://last.cbrc.jp/
  * Can be installed from [bioconda](https://anaconda.org/bioconda/last): `conda install -c bioconda last`
* The 'RNAup' tool that is a part of the ViennaRNA package: https://www.tbi.univie.ac.at/RNA/
  * Can be installed from [bioconda](https://anaconda.org/bioconda/viennarna): `conda install -c bioconda viennarna`

## Download the [latest release](https://github.com/vanya-antonov/assa/releases) and uncompress it
Copy the content of the src/ directory to a folder from your $PATH, for example:
```bash
cp -rv  src/*  ~/bin
```

# Notes
 * The transcript sequences you provide as input to assa.pl must consist of ACGT symbols ONLY.
   No other letters (such as N) are allowed.

 * To test that ASSA is installed properly just type 'assa.pl' and you will see the usage message

 * Test run for the 7SL::TP53 interaction from [Abdelmohsen et al, 2014](https://www.ncbi.nlm.nih.gov/pubmed/25123665) using 4 CPUs:
```bash
assa.pl  --num_threads 4  examples/7SL.fna  examples/TP53.fna  >  output.txt
```

# Citation
If you use ASSA in your research, please cite:
> Antonov I, Marakhonov A, Zamkova M, Medvedeva Y.
> [ASSA: Fast identification of statistically significant interactions between long RNAs.](https://www.ncbi.nlm.nih.gov/pubmed/29375012)
> J Bioinform Comput Biol. 2018 Feb;16(1):1840001.
