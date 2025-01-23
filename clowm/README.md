<img src="https://raw.githubusercontent.com/gi-bielefeld/sans/master/example_data/drosophila/WG_weakly_groups.png" style="border:0;" alt="Example network" align="right"/>

# SANS

**Symmetric Alignment-free phylogeNomic Splits**  


SANS is a whole-genome based, alignment- and reference-free approach
for reconstructing phylogenies. It does not rely on a pairwise comparison of the genomes.
In a pangenomic approach, evolutionary relationships are determined based on the
similarity of the whole sequences. Sequence segments (*k*-mers) shared by a subset
of genomes are interpreted as a phylogenetic split indicating the closeness of these
genomes and their separation from the other genomes. The resulting splits can be
visualized as a phylogenetic tree or network.


* Reference-free
* Alignment-free
* Input: assembled genomes / reads, or coding sequences / amino acid sequences
* Output: phylogenetic splits or tree


## Dos and Don'ts

* The genomes should not be too diverged. SANS works well on species level.
* Be careful with outliers and outgroups (for the reason above).
* The sequences should not be too short. Provide whole-genome data or as many coding sequences as possible.
* Be careful with viruses (for the reasons above).
* Have a look at the network (weakly compatible or 2-tree). It does not make much sense to extract a tree, if the split network is a hairball.
* Reconstructed phylogenies are unrooted, even though a Newick file suggests a root.
* In case of problems, contact us (see below).


## Local installation

SANS can be also easily installed locally from our [git repository](https://github.com/gi-bielefeld/sans). 


## Contact

For any question, feedback, or problem, please feel free to file an issue on this Git repository or write an email and we will get back to you as soon as possible.

[sans-service@cebitec.uni-bielefeld.de](mailto:sans-service@cebitec.uni-bielefeld.de)

SANS is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/). We would appreciate if you would participate in the evaluation of SANS by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=sans).



## Publications

Rempel, A., Wittler, R.: [SANS serif: alignment-free, whole-genome based phylogenetic reconstruction](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab444/6300510). Bioinformatics. (2021).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](https://pub.uni-bielefeld.de/download/2942421/2942423/s13015-020-00164-3.wittler.pdf).
Algorithms for Molecular Biology. 15: 4 (2020).

Wittler, R.: [Alignment- and reference-free phylogenomics with colored de Bruijn graphs](http://drops.dagstuhl.de/opus/volltexte/2019/11032/pdf/LIPIcs-WABI-2019-2.pdf).
In: Huber, K. and Gusfield, D. (eds.) Proceedings of WABI 2019. LIPIcs. 143, Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik, Dagstuhl, Germany (2019).


## License

* The sparse-map library is licensed under the [MIT license](https://github.com/Tessil/sparse-map/blob/master/LICENSE).
* The Bifrost library is licensed under the [BSD-2 license](https://github.com/pmelsted/bifrost/blob/master/LICENSE).
* SANS uses gzstream, licensed under the [LGPL license](https://github.com/gi-bielefeld/sans/blob/master/src/gz/COPYING.LIB).
* SANS is licensed under the [GNU general public license](https://github.com/gi-bielefeld/sans/blob/master/LICENSE).

