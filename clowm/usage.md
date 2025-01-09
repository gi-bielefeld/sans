## Input

To ease the application of SANS for unexperiences users, the CloWM vesion of SANS provides a simplified parameter handling compared to a local installation from our [git repository](https://github.com/gi-bielefeld/sans).

### Format

Input sequences are read from Fasta or Fastq files.
* Each file can be gzipped.
* Each file can contain multiple sequence enties, e.g., contigs.


### Upload

Use the menu "Files", "My Data Buckets" and
* upload all files into one folder individually, or
* join all files in one zip or tar.gz file (no folder structure).

You can also transfer data using an S3 management software such as provided by [AWS](https://aws.amazon.com/cli/) or [minIO](https://min.io/docs/minio/linux/reference/minio-mc.html).

## Output

### Files

* `sans_splitnetwork.pdf` shows the phylogeny, generated with [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).
* `sans_splitnetwork.nexus` can be opened in [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html) to explore the phylogeny interactively.
* `sans_splitnetwork.tsv` is a tab separated file. Each line corresponds to one split. The first column specifies the split weight. The remaining columns specify a set of genomes that is split from the others. Splits are odered by weight.

* `sans_tree.pdf` shows the phylogenetic tree, generated with [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html).
* `sans_tree.newick` contains the phylogenetic tree in newick format.
* `sans_tree.tsv` analog to the tsv file above but containing only the tree splits.

* `sans.log` shows the logging output of the SANS run including the actual parameter settings.

### Download

Use the menu "Files", "My Data Buckets" to acces the ouput files or an S3 command line tool (see "Upload").

## Optional input options

### Abundance threshold
When analyzing read data, a common preprocessing step is to filter out low coverage
*k*-mers that typically arise from sequencing errors. SANS includes the
option `--qualify` to perform such a filtering step, allowing raw read data to be analyzed without
the need to run another tool first. A minimum coverage threshold can be specified,
i.e., `--qualify 2` filters out all *k*-mers that occur less than 2 times per genome.

### File-of-files
By default, 
* each genome is identified by the name of its input file,
* one genome corresponds to one file, and
* the same abundance filter (see "qualify") is used for each genome.

By providing a `file-of-files` in the following format, you can 
* assign custom identifyers to the genomes,
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

### Coding sequences as input
Even though SANS is originally developed to process whole genome DNA data, it also offers to process protein sequences, either translated (using parameter `--amino`), or
untranslated employing automatic translation (using parameter `--code`). Reverse complement *k*-mers are not considered and the default *k*-mer length is 10. By default, the standard genetic code will be used for translation. (A local installation of [SANS](https://github.com/gi-bielefeld/sans) allows to specify different code tables.) 


## Further optional options


### *k*-mer length
You may want to try different values for the *k*-mer length. On shorter or rather heterogeneous sequences, use a smaller *k*, e.g., `-k 15`.

### Labeled output
To depict the phylogeny on a higher level, taxa can be assigned to groups. Each group is then represented by a color and individual text labels of taxa are replaced by colored circles accordingly. An example is shown on the "Description" tab.
 
Use option `--label` to provide a mapping of genome identifiers to group names. The file needs to be tab-separated with genome identifyers in the first column and group names in the second. Not all genomes need to be mapped. Group names can be arbitrary strings. Colors are selected automatically. (A local installation of [SANS](https://github.com/gi-bielefeld/sans) allows to specify custom color assignments to groups.) 

### Number of splits
For large data sets, the list of splits can become very long. We recommend to restrict the output for *n* genomes as input to the 10*n* strongest splits in the output using `--top 10n`. Use an integer (without *n*) to limit by an absolute value.


### Tree or network
By default, the CloWM version of SANS first generates a phylogenetic split network (weakly compatible splits), and then also filters the splits to obtain a tree (strictly compatible spltis). The tree filtering step can be turned off by setting `--tree` to false.


### PDF output
By default, the CloWM version of SANS generates a PDF of the phylogeny using [SplitsTree 4](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html). This can be turned off by setting `--pdf` to false.
