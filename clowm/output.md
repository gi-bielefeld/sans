# Output Files

<img src="https://raw.githubusercontent.com/gi-bielefeld/sans/master/example_data/prasinoviruses/weakly.splits.nexus.png" style="border:0;" alt="Example network" align="right" width="25%"/>

* `sans_splitnetwork.pdf` shows the phylogeny, generated with [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).
* `sans_splitnetwork.nexus` can be opened in [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html) to explore the phylogeny interactively.
* `sans_splitnetwork.tsv` is a tab separated file. Each line corresponds to one split. The first column specifies the split weight. The remaining columns specify a set of genomes that is split from the others. Splits are odered by weight.


* `sans_tree.pdf` shows the phylogenetic tree, generated with [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).
* `sans_tree.newick` contains the phylogenetic tree in newick format.
* `sans_tree.tsv` analog to the tsv file above but containing only the tree splits.


* `sans.log` shows the logging output of the SANS run including the actual parameter settings.


Depending on selected advanced parameters, further output files can be generated. See "Usage" for details.
