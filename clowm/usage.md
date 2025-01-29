## Parameter Views

---

To ease the application of SANS for unexperienced users, the CloWM version of SANS provides a slightly different parameter handling compared to a local installation from our [git repository](https://github.com/gi-bielefeld/sans).

|  |||
|:--|:--|:--|
| <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/simple.png" style="border:0;" alt="simple" height="30"/> | | Select the input and output folders and run SANS with default parameters. |
| <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/advanced.png" style="border:0;" alt="advanced" height="30"/> | | If your input are **reads** or **coding sequences**, or if you want to **beautify the output**, switch the parameter view to "advanced". |
| <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/expert.png" style="border:0;" alt="expert" height="30"/> | | This parameter view provides further options, such as bootstrapping. |


</br>

## <i class="fa fa-file-text" aria-hidden="true"></i> Input / Output   <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/simple.png" style="border:0;" alt="simple" align="right"/>

---

To compute a phylogeny with SANS, the only requirement are the input sequences.

#### Input File Format

Input sequences are read from Fasta or Fastq files.
* Each file can be gzipped.
* Each file can contain multiple sequence entries, e.g., contigs.

If your input are reads or coding sequences, switch the parameter view to "advanced".


#### Upload

Use the menu "Files", "My Data Buckets" and
* upload all files into one folder individually, or
* upload one zip or tar.gz file containing all input files (without folder structure).

You can also transfer data using an S3 management software such as [S3cmd](https://s3tools.org/s3cmd) or provided by [AWS](https://aws.amazon.com/cli/) or [MinIO](https://min.io/docs/minio/linux/reference/minio-mc.html). For example, for using the minIO client, look up the S3 endpoint, access key, and secret key under "Files, S3 Bucket Keys", create an alias with `mc alias set sans-clowm <endpoint> <access_key> <secret_key> --api "s3v4" --path "auto"`, and then upload with `mc cp <files> sans-clowm/<bucket_name>`.


#### Output Files

<img src="https://raw.githubusercontent.com/gi-bielefeld/sans/master/example_data/prasinoviruses/weakly.splits.nexus.png" style="border:0;" alt="Example network" align="right" width="25%"/>

* `sans_splitnetwork.pdf` shows the phylogeny, generated with [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).
* `sans_splitnetwork.nexus` can be opened in [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html) to explore the phylogeny interactively.
* `sans_splitnetwork.tsv` is a tab separated file. Each line corresponds to one split. The first column specifies the split weight. The remaining columns specify a set of genomes that is split from the others. Splits are ordered by weight.

* `sans_tree.pdf` shows the phylogenetic tree, generated with [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).
* `sans_tree.newick` contains the phylogenetic tree in newick format.
* `sans_tree.tsv` analog to the tsv file above but containing only the tree splits.

* `sans.log` shows the logging output of the SANS run including the actual parameter settings.

#### Download

Use the menu "Files", "My Data Buckets" to access the output files or an S3 command line tool (see "Upload").

</br>

## <i class="fa fa-pencil" aria-hidden="true"></i> Advanced Input Arguments   <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/advanced.png" style="border:0;" alt="advanced" align="right"/>

---

#### Abundance threshold
When analyzing read data, a common preprocessing step is to filter out low coverage
*k*-mers that typically arise from sequencing errors. SANS includes the
option `--qualify` to perform such a filtering step, allowing raw read data to be analyzed without
the need to run another tool first. A minimum coverage threshold can be specified,
i.e., `--qualify 2` filters out all *k*-mers that occur less than 2 times per genome.

#### File-of-files
By default, 
* each genome is identified by the name of its input file,
* one genome corresponds to one file, and
* the same abundance filter (see "qualify") is used for each genome.

By providing a `file-of-files` in the following format, you can 
* assign custom identifiers to the genomes,
* specify multiple input files for each genome, and/or
* specify individual abundance thresholds for each genome (overwriting the global qualify-value).

The format has been introduced by the developers of [kmtricks](https://github.com/tlemane/kmtricks/wiki/Input-data) and is specified as follows.

  ```
  <Identifier> : <File1> ; ... ; <FileN> ! <MinAbundance>
  ...
  ```

Example:

  ```
  genome_A : reads_a_forward.fa ; reads_a_reverse.fa ! 2
  genome_B : genome_b_chr_1.fa ; genome_b_chr_2.fa ! 1
  ...
  ```

#### Coding sequences as input
Even though SANS is originally developed to process whole genome DNA data, it also offers to process protein sequences, either translated (using parameter `--amino`), or
untranslated employing automatic translation (using parameter `--translate`). Reverse complement *k*-mers are not considered and the default *k*-mer length is 10. By default, the standard genetic code will be used for translation. See "Genetic code" for further options. 



</br>

## <i class="fa fa-wrench" aria-hidden="true"></i> Advanced Parameters   <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/advanced.png" style="border:0;" alt="advanced" align="right"/>

---


#### *k*-mer length
The default *k*-mer length is 31 for DNA sequences, and 10 for amino acid sequences. 
You may want to try different values. On shorter or rather heterogeneous sequences, use a smaller *k*, e.g., `-k 15`. If you want to use a *k*-mer length larger than 63, please use a local installation of [SANS](https://github.com/gi-bielefeld/sans).


#### Labeled output
<img src="https://raw.githubusercontent.com/gi-bielefeld/sans/master/example_data/drosophila/WG_weakly_groups.png" style="border:0;" alt="Example network" align="right" width="33%"/>
To depict the phylogeny on a higher level, taxa can be assigned to groups. Each group is then represented by a color and individual text labels of taxa are replaced by colored circles accordingly as shown in the figure on the right.
 
Use option `--label` to provide a mapping of genome identifiers to group names. The file needs to be tab-separated with genome identifiers in the first column and group names in the second. Not all genomes need to be mapped. Group names can be arbitrary strings. 

Colors are selected automatically. Optionally, you can use `--label_colors` to specify custom color assignments to groups using an additional tab-separated file with group names in the first and colors (rgb values, e.g. 90 0 255) in the second column.


#### Number of splits
For large data sets, the list of splits can become very long. We recommend to restrict the output for *n* genomes as input to the 10*n* strongest splits in the output using `--top 10n`. Use an integer (without *n*) to limit by an absolute value.


#### Tree or network
By default, the CloWM version of SANS first generates a phylogenetic split network (weakly compatible splits), and then also filters the splits to obtain a tree (strictly compatible splits). The tree filtering step can be turned off by setting `--tree` to false.


#### PDF output
By default, the CloWM version of SANS generates a PDF of the phylogeny using [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html). This can be turned off by setting `--pdf` to false.

#### SVG output
Generates an SVG file of the phylogeny using [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).

</br>


## <i class="fa fa-graduation-cap" aria-hidden="true"></i> Expert Parameters   <img src="https://raw.githubusercontent.com/gi-bielefeld/sans/clowm-integration/clowm/expert.png" style="border:0;" alt="expert" align="right"/>

---


#### Filter criteria
The sorted list of splits is greedily filtered, i.e., splits are iterated from strongest to weakest and a split is kept if and only if the `--filter` criterion is met.
* strict: a split is kept if it is compatible to all previously filtered splits, i.e., the resulting set of splits is equivalent to a tree
* weakly: a split is kept if it is weakly compatible to all previously filtered splits
* planar: a split is kept if the resulting set of splits can be displayed in the plain without any edges crossing (a.k.a. circular compatible, outer-labeled planar)
* *n*-tree: multiple sets of compatible splits (=trees) are maintained. A split is added to the first, second, ... *n*-th set if possible (compatible).
* default: By default, the CloWM version of SANS first generates a phylogenetic split network (weakly compatible splits), and then also filters the splits to obtain a tree (strictly compatible splits).
* none: apply no filter (not recommended)

Any choice other than "default" deactivates the two-step procedure (see "Tree or network").


#### Bootstrapping
To assess the robustness of reconstructed splits with respect to noise in the input data, bootstrap replicates can be constructed by randomly varying the observed *k*-mer content. To compare the originally determined splits or tree to, e.g., 1000 bootstrap replicates, use `--bootstrap 1000`. Bootstrap support values will be integrated into the nexus/newick output file and can, e.g., be visualized in [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html). Further, an additional output file `sans_*.bootstrap` containing the bootstrap support values will be created. This option requires to select a filter criterion other than "none" and "default", see "Filter criteria".

The bootstrap support values can also be used for filtering splits:
* To filter out low support splits, e.g., those appearing in less than 75% of the bootstrap replicates, use `--support 0.75`. (Note that choosing filter "tree" in combination with a threshold of 50% or higher implies compatibility of the filtered splits, i.e., a tree.)
* To filter all splits (greedily according to their support value) to obtain a consensus split network or tree, select a filter for `--consensus`. You could, e.g., use support values from tree replicates (`--filter strict`) to filter for a consensus tree (`--consensus strict`) or for a split network (`--consensus weakly`).


#### IUPAC characters
By default, *k*-mers spanning ambiguous characters like 'N' are ignored/skipped. Instead, using `--iupac <threshold>` all variants of such *k*-mers can be considered, replacing the ambiguous character by all the characters it represents. For instance, for DNA input, one occurrence of 'N' results in four variants. If a *k*-mer spans multiple ambiguous characters, the number of resulting combinations grows quickly. Thus, the number of considered *k*-mers is limited. If resolving the next read ambiguous character would exceed `<threshold>`, the *k*-mer is skipped instead.

#### Canonical *k*-mers
By default, each *k*-mer is compared to its reverse complement and the lexicografically smaller is chosen as a canonical representative.
Set `--norev` to true to not consider reverse complement *k*-mers.

#### Genetic code
Use `--code` to select the ID of the genetic code to be used for translation, if `--translate` is used.
Default: 1 (standard genetic code). Use 11 for Bacterial, Archaeal, and Plant Plastid Code. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details.

#### Mean function
If *a* *k*-mers are present in genome set *A* but not *B*, and *b* *k*-mers are present in genome set *B* but not *A*, both counts are combined by a mean function to obtain a final weight for the split *{A,B}*. Option `--mean` offers:
* arith: arithmetic mean, *(a+b)/2*
* geom:  geometric mean, *sqrt(a b)*
* geom2: geometric mean with pseudo-counts, *sqrt((a+1)(b+1))* (default)

#### Output core k-mers
Use `--core` to output all core *k*-mers, i.e., *k*-mers appearing in all input genomes, in a Fasta file `sans_core.fasta`.

#### Ignore certain k-mers
Use `--blacklist` to provide a Fasta or Fastq file containing *k*-mers to be ignored when reading input files.

