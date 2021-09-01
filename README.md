# SANS *serif*

**Symmetric Alignment-free phylogeNomic Splits**  
***--- Space and time Efficient Re-Implementation including Filters***

* Reference-free
* Alignment-free
* Assembled genomes or reads as input
* Phylogenetic splits or tree as output
* **NEW:** Coding sequences / amino acid sequences as input (see --code and --amino), implemented by Marco Sohn

### Publications

Rempel, A., Wittler, R.: [SANS serif: alignment-free, whole-genome based phylogenetic reconstruction](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab444/6300510). Bioinformatics. (2021).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](https://pub.uni-bielefeld.de/download/2942421/2942423/s13015-020-00164-3.wittler.pdf).
Algorithms for Molecular Biology. 15: 4 (2020).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](http://drops.dagstuhl.de/opus/volltexte/2019/11032/pdf/LIPIcs-WABI-2019-2.pdf).
In: Huber, K. and Gusfield, D. (eds.) Proceedings of WABI 2019. LIPIcs. 143, Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, Dagstuhl, Germany (2019).

## Table of Contents

* [Requirements](#requirements)
* [Compilation](#compilation)
* [Usage](#usage)
* [Details on filter option](#details-on-filter-option)
* [Examples](#examples)
* [Performance evaluation on predicted open reading frames](#performance-evaluation-on-predicted-open-reading-frames)
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

## Usage

```
SANS
```

displays the command line interface:
```
Usage: SANS [PARAMETERS]

  Input arguments:

    -i, --input   	 Input file: file of files format
                  	 Either: list of files, one genome per line (space-separated for multifile genomes)
                  	 Or: kmtricks input format (see https://github.com/tlemane/kmtricks)

    -g, --graph   	 Graph file: load a Bifrost graph, file name prefix
                  	 (requires compiler flag -DuseBF, please edit makefile)

    -s, --splits  	 Splits file: load an existing list of splits file
                  	 (allows to filter -t/-f, other arguments are ignored)

    (either --input and/or --graph, or --splits must be provided)

  Output arguments:

    -o, --output  	 Output TSV file: list of splits, sorted by weight desc.

    -N, --newick  	 Output Newick file
                  	 (only applicable in combination with -f strict or n-tree)

    (at least --output or --newick must be provided, or both)

  Optional arguments:

    -k, --kmer    	 Length of k-mers (default: 31, or 10 for --amino and --code)

    -t, --top     	 Number of splits in the output list (default: all).
                  	 Use -t <integer>n to limit relative to number of input files, or
                  	 use -t <integer> to limit by absolute value.

    -m, --mean    	 Mean weight function to handle asymmetric splits
                  	 options: arith: arithmetic mean
                  	          geom:  geometric mean
                  	          geom2: geometric mean with pseudo-counts (default)

    -f, --filter  	 Output (-o, -N) is a greedy maximum weight subset (see README)
                  	 options: strict: compatible to a tree
                  	          weakly: weakly compatible network
                  	          n-tree: compatible to a union of n trees
                  	                  (where n is an arbitrary number, e.g. 2-tree)

    -x, --iupac   	 Extended IUPAC alphabet, resolve ambiguous bases or amino acids
                  	 Specify a number to limit the k-mers per position between 
                  	 1 (no ambiguity) and 4^k respectively 22^k (allows NNN...N)
                  	 Without --iupac respective k-mers are ignored

    -n, --norev   	 Do not consider reverse complement k-mers

    -a, --amino   	 Consider amino acids: --input provides amino acid sequences
                  	 Implies --norev and a default k of 10

    -c, --code    	 Translate DNA: --input provides coding sequences
                  	 Implies --norev and a default k of 10
                  	 optional: ID of the genetic code to be used
                  	 Default: 1 (The Standard Code)
                  	 Use 11 for Bacterial, Archaeal, and Plant Plastid Code
                  	 (See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.)

    -v, --verbose 	 Print information messages during execution

    -h, --help    	 Display this help page and quit
  
```

## Details on filter option

The sorted list of splits is greedily filtered, i.e., splits are iterated from strongest to weakest and a split is kept if and only if the filter criterion is met.

* `strict`: a split is kept if it is compatible to all previously filtered splits, i.e., the resulting set of splits is equivalent to a tree.
* `weakly`: a split is kept if it is weakly compatible to all previously filtered splits (see publication for definition of "weak compatibility").
* `n-tree`: several sets of compatible splits (=trees) are maintained. A split is added to the first, second, ... *n*-th set if possible (compatible).

<!--**Cleanliness:** The filtered set of splits is compared to the originally given (`--splits`) or computed (`--input`, `--graph`) set of splits to obtain a measure of how many incompatible splits have been filtered out. Consider a list of splits *S* that has been filtered to the sublist *F*, both sorted in descending order. To make the measure robust against low weighting splits and the choice of paramter `--top`, we truncate *S* to "just contain F": let *S'* be the shortest prefix of *S* that contains *F*. Then the **cleanliness** is the ratio of the length of *F* to the length of *S'*. The **weighted cleanliness** is the ratio of the corresponding sum of weights of splits in *F* and *S'*, resp. 
-->


## Examples

1. **Determine splits from assemblies or read files**
   ```
   SANS -i list.txt -o sans.splits -k 31
   ```
   The 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt* (`-i list.txt`) are extracted. Splits are determined and written to *sans.splits* (`-o sans.splits`).

   To extract a tree (`-f strict`) in NEWICK format (`-N sans_greedytree.new`), use
   ```
   SANS -i list.txt -k 31 -f strict -N sans_greedytree.new
   ```
   or filter from a set of splits (`-s sans.splits`)
   ```
   SANS -s sans.splits -f strict -N sans_greedytree.new
   ```

2. **Drosophila example data**
   ```
   # go to example directory
   cd <SANS directory>
   cd example_data/drosophila

   # download data: whole genome and coding sequences
   ./download_WG.sh
   ./download_CDS.sh

   # run SANS greedy tree
   ../../SANS -i WG/list.txt -o sans_greedytree_WG.splits -f strict -N sans_greedytree_WG.new -v
   ../../SANS -i CDS/list.txt -o sans_greedytree_CDS.splits -f strict -N sans_greedytree_CDS.new -v -c
   
   # compare to reference
   ../../scripts/newick2sans.py Reference.new > Reference.splits
   ../../scripts/comp.py sans_greedytree_WG.splits Reference.splits WG/list.txt > sans_greedytree_WG.comp
   ../../scripts/comp.py sans_greedytree_CDS.splits Reference.splits CDS/list.txt > sans_greedytree_CDS.comp
   ```

3. **Virus example data**
   ```
   # go to example directory
   cd <SANS directory>
   cd example_data/prasinoviruses

   # download data
   ./download.sh

   # run SANS
   ../../SANS -i fa/list.txt -o sans.splits -k 11 -t 130 -v
   
   # compare to references
   ../../scripts/newick2sans.py Reference_Fig3.new > Reference_Fig3.splits
   ../../scripts/comp.py sans.splits Reference_Fig3.splits fa/list.txt > fig3.comp
   ../../scripts/newick2sans.py Reference_Fig4.new > Reference_Fig4.splits
   ../../scripts/comp.py sans.splits Reference_Fig4.splits fa/list.txt > fig4.comp
   ```

## Performance evaluation on predicted open reading frames
SANS-serif can predict phylogenies based on amino acid sequences as input. Processing coding sequences is faster than processing whole genome data. Experiments show that the reconstruction quality is comparable. 

If you want to use selected marker genes, the number of genes should be as high as possible to provide sufficient sequence information for extracting phylogenetic signals.

If the genomes at hand are not annotated, you can use a tool to predict open reading frames. The following experiment shows that the reconstruction quality does not suffer from such a simple pre-processing or even improves, while saving total running time, especially because genomes can easily be pre-processed in parallel. 

The following tools have been used with the parameters shown in the table.

| Tool | Parameters| Reference |
|:--|:--|:--|
| SANS | -k (see below) -m geom2 -t 10n -filter strict [-a (if preprocessed)]| |
| Getorf (EMBOSS)| -find (see below) -t 11 | Gary Williams, 2000 |
| ORFfinder | -n true -g 11 | NCBI |
| Prodigal | -q | V2.6.3,Doug Hyatt, 2016 |

For estimating the reconstruction accuracy, the (weighted) F1 score has been determined as follows.

| Measure ||
|--|--|
|F1-score  | harmonic mean of precision and recall |
|precision | (number of called splits that are also in reference) <br> / (total number of called splits) |
|recall    | (number of reference splits that are also called) <br> / (total number of reference splits) |
|weighted precision | (total weight of called splits that are also in reference) <br> / (total weight of all called splits) |
|weighted recall | (total weight of reference splits that are also called) <br> / (total weight of all reference splits) |

Further information on the datasets can be found in the initial publication of SANS (Wittler, 2019), see above.

### *Salmonella enterica* Para C
220 genomes, k=31

| Preprocessing | Running time <br> pre-processing | Running time <br>SANS | Running time both <br>(parall. pre-proc.,<br> factor 16) | F1-Score <br> (weighted) |
|:--|--:|--:|--:|--:|
| none (whole genome) | -- | **610s** | -- | **0.878 <br> (0.999)** |
| Getorf (-find 0) | 247s | 459s | 706s <br> (474s) | 0.881 <br> (0.999) |
| Getorf (-find 1) | 199s | 339s | 538s <br> (351s) | 0.868 <br> (0.999) |
| ORFfinder | 5195s | 166s | 5361s <br> (491s) | 0.858 <br> (0.997) |
| Prodigal | 6292s | 134s | 6326s <br> (521s) |0.853 <br> (0.998)  |

### *Salmonella enterica* subspecies enterica
2964 genomes, k=21

| Preprocessing | Running time <br> pre-processing | Running time <br>SANS | Running time both <br>(parall. pre-proc.,<br> factor 16) | F1-Score <br> (weighted) |
|:--|--:|--:|--:|--:|
| none (whole genome) | -- | **190min** | -- | **0.587 <br> (0.792)** |
| Getorf (-find 0) | 55min | 220min | 274min <br> (223min) | 0.624 <br> (0.807) |
| Getorf (-find 1) | 48min | 165min | 213min <br> (168min) | 0.620 <br> (0.799) |
| ORFfinder | 1621min | 78min | 1698min <br> (179min) | 0.594 <br> (0.766) |
| Prodigal | 1322min | 50min | 1372min <br> (133min) |0.587 <br> (0.792)  |

## Contact

For any question, feedback, or problem, please feel free to file an issue on this Git repository or write an email and we will get back to you as soon as possible.

[sans-service@cebitec.uni-bielefeld.de](mailto:sans-service@cebitec.uni-bielefeld.de)

SANS is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/). We would appriciate if you would participate in the evaluation of SANS by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans).


## License

* The sparse-map library is licensed under the [MIT license](https://github.com/Tessil/sparse-map/blob/master/LICENSE).
* The Bifrost library is licensed under the [BSD-2 license](https://github.com/pmelsted/bifrost/blob/master/LICENSE).
* SANS is licensed under the [GNU general public license](https://gitlab.ub.uni-bielefeld.de/gi/sans/blob/master/LICENSE).

<img src="https://piwik.cebitec.uni-bielefeld.de/matomo.php?idsite=12&rec=1&action_name=VisitGitLab&url=https://gitlab.ub.uni-bielefeld.de/gi/sans" style="border:0;" alt="" />



