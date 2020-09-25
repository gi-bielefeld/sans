# SANS *serif*

**Symmetric Alignment-free phylogeNomic Splits**  
***--- Space and time Efficient Re-Implementation including Filters***

* Reference-free
* Alignment-free
* Assembled genomes or reads as input
* Phylogenetic splits as output

### Publications

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs.](https://pub.uni-bielefeld.de/download/2942421/2942423/s13015-020-00164-3.wittler.pdf)
Algorithms for Molecular Biology. 15: 4 (2020).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs.](http://drops.dagstuhl.de/opus/volltexte/2019/11032/pdf/LIPIcs-WABI-2019-2.pdf)
In: Huber, K. and Gusfield, D. (eds.) Proceedings of WABI 2019. LIPIcs. 143, Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, Dagstuhl, Germany (2019).

## Table of Contents

* [Requirements](https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/serif#requirements)
* [Compilation](https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/serif#compilation)
* [Usage](https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/serif#usage)
* [FAQ](https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/serif#faq)
* [Contact](https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/serif#contact)
* [License](https://gitlab.ub.uni-bielefeld.de/gi/sans/tree/serif#license)

## Requirements

For the main program, there are no strict dependencies. However, there are some optional features:
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

If Bifrost should be used, change the SANS makefile accordingly (easy to see how).
Please note the installation instructions regarding the default maximum *k*-mer size of Bifrost in its README.

If during the compilation, the Bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler.
You may have to add `-I/usr/local/include` (with the corresponding folder) to the compiler flags in the makefile.

## Usage:

```
SANS
```

displays the command line interface:
```
Usage: SANS [PARAMETERS]

  Required arguments:

    -i, --input   	 Input file: list of sequence files, one per line

    -g, --graph   	 Graph file: load a Bifrost graph, file name prefix
                  	 (at least --input or --graph must be provided, or both)

    -s, --splits  	 Splits file: load an existing list of splits file
                  	 (allows to filter -t/-f, other arguments are ignored)

                  	 (either --input and/or --graph, or --splits must be provided)
                  	 
    -o, --output  	 Output file: list of splits, sorted by weight desc.

    -N, --newick  	 Output newick file
                  	 (only applicable in combination with -f strict or -f n-tree)

                  	 (at least --output or --newick must be provided, or both)

    Optional arguments:

    -k, --kmer    	 Length of k-mers (default: 31)

    -t, --top     	 Number of splits (default: all)

    -m, --mean    	 Mean weight function to handle asymmetric splits
                  	 options: arith: arithmetic mean
                  	          geom:  geometric mean (default)
                  	          geom2: geometric mean with pseudo-counts

    -f, --filter  	 Output a greedy maximum weight subset
                  	 options: strict: compatible to a tree
                  	          weakly: weakly compatible network
                  	          n-tree: compatible to a union of n trees
                  	                  (where n is an arbitrary number)

    -x, --iupac   	 Extended IUPAC alphabet, resolve ambiguous bases
                  	 Specify a number to limit the k-mers per position
                  	 between 1 (no ambiguity) and 4^k (allows NNN...N)

    -n, --norev   	 Do not consider reverse complement k-mers

    -v, --verbose 	 Print information messages during execution

    -h, --help    	 Display this help page and quit
```

### Examples

1. **Determine splits from assemblies or read files**
   ```
   SANS  -k 31 -i list.txt -o sans.splits
   ```
   The 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt* (`-i list.txt`) are extracted. Splits are determined and written to *sans.splits* (`-o sans.splits`).

   To extract a tree (`-f strict`) in NEWICK format (`-N sans_greedytree.new`), use 
   ```
   SANS -i list.txt -k 31 -f strict -N sans_greedytree.new 
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
   SANS -i list.txt -f strict -o ../sans_greedytree.splits -N sans_greedytree.new -t 130 -v
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
   SANS -i list.txt -o ../sans.splits -k 11 -t 130 -v
   cd ..

   # compare to references
   ../../scripts/newick2sans.py Reference_Fig3.new > Reference_Fig3.splits
   ../../scripts/comp.py sans.splits Reference_Fig3.splits fa/list.txt
   ../../scripts/newick2sans.py Reference_Fig4.new > Reference_Fig4.splits
   ../../scripts/comp.py sans.splits Reference_Fig4.splits fa/list.txt
   ```

## FAQ

We recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).

## Contact

For any question, feedback, or problem, please feel free to file an issue on this Git repository and we will get back to you as soon as possible.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)
* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)
* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)
* The kseq library is copyrighted by Heng Li and released under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)
* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)
* Bifrost is BSD-2 licensed (https://github.com/pmelsted/bifrost)
* SANS is under GNU general public license (https://gitlab.ub.uni-bielefeld.de/gi/sans/blob/serif/LICENSE)
