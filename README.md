# SANS

### Symmetric Alignment-free phylogeNomic Splits

* Reference-free
* Alignment-free
* Assembled genomes or reads as input
* Phylogenetic splits as output

### Publication

Preprint: Wittler, R.: Alignment- and reference-free phylogenomics with colored de-Bruijn graphs. [arXiv:1905.04165](https://arxiv.org/abs/1905.04165). (2019). 


## Table of Contents

* [Requirements](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#requirements)
* [Compilation](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#compilation)
* [Usage](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#usage)
* [FAQ](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#faq)
* [Contact](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#contact)
* [License](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans#license)

## Requirements

From the given genomes, a **colored de-Bruijn graph** is built to efficiently extract common subsequences. To this end, SANS uses the API of [Bifrost](https://github.com/pmelsted/bifrost). Apart from the requirements of Bifrost (c++ and cmake), there are no further strict dependencies.

* To convert the output into NEXUS format, the provided script requires Phython 3.

* To filter (greedy tree/compatible, greedy weakly compatible) or visualize the splits, the tool [SplitsTree](http://www.splitstree.org/) can be used.


## Compilation

```
cd <sans_directory>
make
```

By default, the installation creates:
* a binary (*SANS*)

You may want to make the binary (*SANS*) accessible via your (*PATH*) variable.

Please note the installation instructions regsarding the default maximum *k*-mer size of Bifrost in its README.
If your Bifrost libraries have been compiled for 64 bit, change the SANS makefile accordingly (easy to see how).


## Usage:

```
SANS
```

displays the command line interface:
```
Usage: SANS [PARAMETERS]

[PARAMETERS]:

  > Mandatory with required argument:

  -s, --input-seq-files   Input sequence files (FASTA/FASTQ possibly gzipped)
                          Input files can be provided as a list in a TXT file (one file per line)
                          K-mers with exactly 1 occurrence in the input files will be discarded
  -r, --input-ref-files   Input reference files (FASTA/FASTQ possibly gzipped and GFA)
                          Input files can be provided as a list in a TXT file (one file per line)
                          All k-mers of the input reference files are used
  -o, --output-file       name of output file

  > Optional with required argument:

  -t, --threads           Number of threads (default: 1)
  -T, --top               Output the top T splits sorted by weight descending (default: all)
  -k, --kmer-length       Length of k-mers (default: 31)
  -m, --min-length        Length of minimizers (auto-adjusted by default: see verbose output)
  -b, --bloom-bits        Number of Bloom filter bits per k-mer with 1+ occurrences (default: 14)
  -B, --bloom-bits2       Number of Bloom filter bits per k-mer with 2+ occurrences (default: 14)
  -l, --load-mbbf         Input Blocked Bloom Filter file, skips filter step (default: no input)
  -w, --write-mbbf        Output Blocked Bloom Filter file (default: no output)
  -u, --chunk-size        Read chunk size per thread (default: 64)

  > Optional with no argument:

  -i, --clip-tips         Clip tips shorter than k k-mers in length
  -d, --del-isolated      Delete isolated contigs shorter than k k-mers in length
  -y, --keep-mercy        Keep low coverage k-mers connecting tips
  -v, --verbose           Print information messages during execution


...
```

### Examples

1. **Determine splits from assemblies**
   ```
   SANS -t 4 -k 31 -o splits.txt -r list.txt
   ```
   The colored de-Bruijn graph is built with Bifrost using 4 threads (`-t 4`) from the 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt*  (`-r list.txt`). Splits are determined and written to *splits.txt* (`-o splits.txt`)

2. **Determine splits from read files**
   ```
   SANS -t 4 -k 31 -o splits.txt -s list.txt
   ```
   The colored de-Bruijn graph is built with Bifrost using 4 threads (`-t 4`) from the 31-mers (`-k 31`) of those fasta or fastq files listed in *list.txt*  (`-s list.txt`). By using parameter `-s`, all files are filtered: k-mers occurring exactly once in a file are discarded from the construction.  Splits are determined and written to *splits.txt* (`-o splits.txt`).

3. **Drosophila example data**
   ```
   # go to example directory
   cd <SANS dir>
   cd example_data/drosophila
   
   # download data
   ./download.sh
   
   # run SANS
   cd fa
   SANS -r list.txt -o ../sans.splits -T 130 -t 4 -v
   cd ..
   
   # greedy tree
   ../../scripts/sans2nexus.py sans.splits fa/list.txt > sans.nexus
   # open sans.nexus in Splitstree (splitstree.org)
   # -> Data -> Greedily Make Compatible
   # -> File -> Save As "sans_greedytree.nexus" (or export the tree in newick format)
   # results is prepared in folder
   ../../scripts/nexus2sans.py sans_greedytree.nexus > sans_greedytree.splits
      
   #compare to reference
   ../../scripts/newick2sans.py Reference.new > Reference.splits
   ../../scripts/comp.py sans_greedytree.splits Reference.splits fa/list.txt
   ```

4. **Virus example data**
   ```
   # go to example directory
   cd <SANS dir>
   cd example_data/prasinoviruses
      
   # download data
   ./download.sh
   
   # run SANS
   cd fa
   SANS -r list.txt -o ../sans.splits -T 130 -t 4 -v -k11
   cd ..

   #compare to references
   ../../scripts/newick2sans.py Reference_Fig3.new > Reference_Fig3.splits
   ../../scripts/comp.py sans.splits Reference_Fig3.splits fa/list.txt
   ../../scripts/newick2sans.py Reference_Fig4.new > Reference_Fig4.splits
   ../../scripts/comp.py sans.splits Reference_Fig4.splits fa/list.txt
   ```


## FAQ

We recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).


## Contact

For any question, feedback or problem, please feel free to file an issue on this Git repository and we will get back to you as soon as possible.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

* The kseq library is copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)

* Bifrost is BSD-2 licensed (https://github.com/pmelsted/bifrost)

* SANS is BSD-2 licensed [LICENSE](https://gitlab.ub.uni-bielefeld.de/roland.wittler/sans/blob/master/LICENSE)

