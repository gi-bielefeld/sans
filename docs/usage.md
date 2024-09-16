
[comment]: <> (Add the new SANS logo here)
<br>
<img src="resource/logo.png" alt="Description of the image">

***

<br>
The <b><i> SANS ambages </i></b> software provides an easy way of phylogenomic reconstruction, including a pipeline for immediate output rendering via 
the tool <i>Splits Tree </i>. 
The software comes with an extensive set of parameters and features, allowing for applicability in a wide range of scenarios.
For a simplistic run however, it is sufficient to provide the the input 
and output location of the target data. 
<br>
<br>

- A minimal example for such a scenario is illustrated in the   [`Quickstart Guide`](#QickstartChap) . 

<br>

- Additional parameters are explained in detail in the
[Parameter Collection](#ParamChap).

<br>

- Convinced by SANS? An installation guide for your system is provided in the
['Installation Guide'](#InstalChap).

<br>

***

<br>
<br>

 # [`Quickstart Guide`](#QickstartChap)

- <b>  [`Uploading Data`](#UploadSec) </b>

- <b> Preparing the Workflow </b>

- <b> Executing the Workflow</b>

- <b> Accessing the Results </b>


<br>
<br>


##  [Uploading Data](#UploadSec)


<br>
<br>



## Preparing the Job



## Executing 



# Accessing the Results

***

**Input files**

Specify your input by `-i <list>` where `<list>` is either a file-of-files or in kmtricks format. Each file can be in fasta, multiple fasta or fastq format.
- **File-of-files:**
  ```
  genome_a.fa
  genome_b.fa
  ...
  ```
  Files can be in subfolders and/or compressed:
  ```
  dataset_1/genome_a.fa.gz
  dataset_1/genome_b.fa.gz
  ...
  ```
  One genome can also be composed of several files (the first one will be used as identifier in the output):
  ```
  reads_a_forward.fa reads_a_reverse.fa
  genome_b_chr_1.fa genome_b_chr_2.fa
  ...
  ```
- **kmtricks format:**
  In this format, you can specify individual identifiers and, optionally, abundance thresholds (see "read data as input"):
  ```
  genome_A : reads_a_forward.fa ; reads_a_reverse.fa ! 2
  genome_B : genome_b_chr_1.fa ; genome_b_chr_2.fa ! 1
  ...
  ```

**Input paramters**

- genomes/assemblies as input: just use `-i <list>`
- read data as input: to filter out *k*-mers of low abundance, either use `-q 2` (or higher thresholds) to specify a global threshold for all input files, or use the kmtricks file-of-files format to specify (individual) thresholds.
- mix of assemblies and read data as input: use the kmtricks file-of-files format to specify individual thresholds.
- coding sequences as input: add `-a` if input is provided as translated sequences, or add `-c` if translation is required. See usage information (`SANS --help`) for further details.


**Output**
- The output file, specified by `-o <split-file>`, lists all splits line by line, sorted by their weight in a tab-separated format where the first column is the split weight and the remaining entries are the identifiers of genomes split from the others.
- For large data sets, the list of splits can become very long. We recommend to restrict the output for *n* genomes as input to the *10n* strongest splits in the output using `-t 10n`.
- We recommend to filter the splits using `-f <filter>`. Then, the sorted list of splits is greedily filtered, i.e., splits are iterated from strongest to weakest and a split is kept if and only if the filter criterion is met.

***
## [Parameter Collection](#ParamChap)
***


***
## Installing SANS on your sytem
***