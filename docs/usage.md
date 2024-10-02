
# SANS *ambages*: Usage

[comment]: <> (Add the new SANS logo here:: PID=1302)

![SANS LOGO](https://gitlab.ub.uni-bielefeld.de/gi/sans/-/blob/clowm-integration-extensions/docs/resource/logo.png)

![SANS LOGO](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Flogo.png/raw?ref=clowm-integration-extensions)
<br>

The ***SANS ambages*** software provides an easy way of phylogenomic reconstruction, including a pipeline for immediate output rendering via 
the tool <i>Splits Tree </i>. 
The software comes with an extensive set of parameters and features, allowing for applicability to a wide range of scenarios.
In the following, we provide a set of explanations and examples of how to configure and stage a SANS workflow, using the CloWM environment. 

[comment]: <> (F.K. Added-3 break-spacing between sections ) 
<br>
<br>
<br>

# Contents
***
<br>

[**Quickstart Guide**](#quickstart-guide)

- [Preparing the Data](#preparing-the-data)
- Preparing the Workflow
- Accessing Results

### Parameter Tweeking

- Adjusting k-mer Length
- Trees and Networks
- Graph Generation
- Amino Acid Data
- Filtering Singletons

### Convinced by SANS? An Installation Guide 

- Downloading SANS
- Installing

<br><br>

## Quickstart Guide

<br>

Great! You've made it to the SANS quickstart guide!

In this section we would like to give a minimal example of the SANS *ambages* workflow to help inexperienced users get started without requireing any prior knowledge of the software or the CloWM environment. 

The three steps below will show you how to upload data from your local system, run a minimalist workflow and download the computed phylogeny.

<br>

### Preparing the Data

![QS_Preparing_Data: 1_Goto_Data](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2FQuickstart%2F1_Preparing_Data%2F1_Goto_Data.png/raw?ref=clowm-integration-extensions)

Starting out, look for the to the `Files` tab at the top left of your browser and select `My Data Buckets`.

<br>
<br>

![QS_Preparing_Data: 1_Goto_Data](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2FQuickstart%2F1_Preparing_Data%2F_Goto_Data.png/raw?ref=clowm-integration-extensions)

## Preparing the Job





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