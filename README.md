# SANS-KC: Phylogenomic Splits + *K*-mer Counting

This branch extends the main version of **SANS serif** with some basic pangenome *k*-mer counting capabilities.  
**Please note**: currently not all features of the master branch are supported and some parameters have different names.
Most notably, the `--norev` flag is now enabled by default and can be disabled using the new `--reverse` flag.
Please make sure to verify the correctness of the parameters when re-running your previous pipelines with this version.

**Symmetric Alignment-free phylogeNomic Splits**
* Reference-free
* Alignment-free
* Assembled genomes or reads as input
* Phylogenetic splits or tree as output

### Publications

* Rempel, A. and Wittler, R.: [SANS serif: alignment-free, whole-genome-based phylogenetic reconstruction](https://academic.oup.com/bioinformatics/article-pdf/37/24/4868/41726858/btab444.pdf).  
  Bioinformatics, 37(24), 4868-4870 (2021).
* Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](https://link.springer.com/content/pdf/10.1186/s13015-020-00164-3.pdf).  
  Algorithms for Molecular Biology, 15(1), 1-12 (2020).
* Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](http://drops.dagstuhl.de/opus/volltexte/2019/11032/pdf/LIPIcs-WABI-2019-2.pdf).  
  19th International Workshop on Algorithms in Bioinformatics (WABI). Schloss Dagstuhl-Leibniz-Zentrum für Informatik (2019).

## Table of Contents

* [Requirements](#requirements)
* [Compilation](#compilation)
* [Usage](#usage)
* [Examples](#examples)
* [Contact](#contact)
* [License](#license)

## Requirements

For the main program, there are no strict dependencies other than C++ version 14.  
However, there are some _optional_ features:
* To read in a **colored de Bruijn graph**, SANS uses the API of [Bifrost](https://github.com/pmelsted/bifrost).
* To convert the output into NEXUS format, the provided script requires Python 3.
* To visualize the splits, we recommend the tool [SplitsTree](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree).

## Compilation

```
cd <SANS directory>
make
```

By default, the installation creates:
* a binary (*SANS*)

You may want to make the binary (*SANS*) accessible via your *PATH* variable.

Optional: If Bifrost should be used, change the SANS makefile accordingly (easy to see how). Please note the installation instructions regarding the default maximum *k*-mer size of Bifrost in its README. If during the compilation, the Bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler. You may have to add `-I/usr/local/include` (with the corresponding folder) to the compiler flags in the makefile. We also recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).

In the *makefile*, two parameters are specified:
* *-DmaxK*: Maximum *k*-mer length that can be chosen when running SANS. Default: 32
* *-DmaxN*: Maximum number of input files for SANS. Default: 64

These values can simply be increased if necessary. To keep memory requirements small, do not choose these values unnecessarily large.

## Usage

```
SANS
```

displays the command line interface:
```
Usage: SANS [PARAMETERS]

  Input arguments:

    -i, --input   	 Input FASTA files: list of sequence files, one per line
    -j, --index   	 Input Index file: load a k-mer index, e.g. counts table
    -g, --graph   	 Input Graph file: load a Bifrost graph, filename prefix
    -s, --splits  	 Input Splits file: load an existing list of splits file

  Output arguments:

    -o, --output  	 Output TSV file: list of splits, sorted by weight desc.
    -n, --newick  	 Output Newick file: convert splits to a tree topology
    -c, --counts  	 Output K-mer file: list k-mer occurrence per input file
    -d, --diff    	 Print the difference between two index or splits files

  K-mer options:

    -k, --kmer    	 Length of k-mers (default: 31)
    -l, --gapped  	 Pattern of gapped k-mers (default: no gaps)
    -w, --window  	 Number of k-mers per minimizer window (default: 1)
    -x, --iupac   	 Extended IUPAC alphabet, resolve ambiguous bases
    -q, --qualify 	 Discard k-mers with lower coverage than a threshold
    -r, --reverse 	 Keep one repr. for reverse complement k-mers

  Filter options:

    -t, --top     	 Number of splits in the output list (default: all)
    -m, --mean    	 Mean weight function to handle asymmetric splits
    -f, --filter  	 Output a greedy maximum weight subset of splits

  Other settings:

    -p, --threads 	 Number of parallel threads (default: auto)
    -v, --verbose 	 Print information messages during execution
    -h, --help    	 Display an extended help page and quit
  
```

## Examples

1. **Determine splits from assemblies or read files**
   ```
   SANS -i list.txt -o sans.splits -k 31 -rv
   ```
   The 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt* (`-i list.txt`) are extracted. Splits are determined and written to *sans.splits* (`-o sans.splits`).

   To extract a tree (`-f strict`) in NEWICK format (`-n sans_greedytree.new`), use
   ```
   SANS -i list.txt -n sans_greedytree.new -k 31 -f strict -rv
   ```
   or filter from a set of splits (`-s sans.splits`)
   ```
   SANS -s sans.splits -n sans_greedytree.new -f strict -v
   ```

2. **Drosophila example data**
   ```
   # go to example directory
   cd <SANS directory>
   cd example_data/drosophila

   # download data
   ./download.sh

   # run SANS greedy tree
   cd fa
   SANS -i list.txt -o ../sans_greedytree.splits -n ../sans_greedytree.new -t 130 -f strict -rv
   cd ..

   # compare to reference
   ../../scripts/newick2sans.py Reference.new > Reference.splits
   ../../scripts/comp.py sans_greedytree.splits Reference.splits fa/list.txt
   ```

3. **Virus example data**
   ```
   # go to example directory
   cd <SANS directory>
   cd example_data/prasinoviruses

   # download data
   ./download.sh

   # run SANS
   cd fa
   SANS -i list.txt -o ../sans.splits -k 11 -t 130 -rv
   cd ..

   # compare to references
   ../../scripts/newick2sans.py Reference_Fig3.new > Reference_Fig3.splits
   ../../scripts/comp.py sans.splits Reference_Fig3.splits fa/list.txt
   ../../scripts/newick2sans.py Reference_Fig4.new > Reference_Fig4.splits
   ../../scripts/comp.py sans.splits Reference_Fig4.splits fa/list.txt
   ```

## Contact

For any question, feedback, or problem, please feel free to file an issue on this Git repository or write an email, and we will get back to you as soon as possible:
[sans-service@cebitec.uni-bielefeld.de](mailto:sans-service@cebitec.uni-bielefeld.de).  
SANS is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/).
We would appreciate if you participate in the evaluation of SANS by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans).

## License

* SANS is licensed under the [GNU general public license](https://gitlab.ub.uni-bielefeld.de/gi/sans/blob/master/LICENSE).
* The Bifrost library is licensed under the [BSD-2 license](https://github.com/pmelsted/bifrost/blob/master/LICENSE).
* The sparse-map library is licensed under the [MIT license](https://github.com/Tessil/sparse-map/blob/master/LICENSE).
* The concurrent-queue library is licensed under [BSD license](https://github.com/cameron314/concurrentqueue/blob/master/LICENSE.md).
<img src="https://piwik.cebitec.uni-bielefeld.de/matomo.php?idsite=12&rec=1&action_name=VisitGitLab&url=https://gitlab.ub.uni-bielefeld.de/gi/sans" style="border:0;" alt="" />