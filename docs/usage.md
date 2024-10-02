
[comment]: <> (Author: Fabian Kolesch)

![SANS LOGO](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Flogo.png/raw?ref=clowm-integration-extensions)
<br>

The ***SANS ambages*** software provides an easy way of phylogenomic reconstruction, including a pipeline for immediate output rendering via 
the tool <i>Splits Tree </i>. 
The software comes with an extensive set of parameters and features, allowing for applicability to a wide range of scenarios.
In the following, we provide a set of explanations and examples of how to configure and stage a SANS workflow, using the CloWM environment. 

[comment]: <> (Format: 3 break between sections, 2 br between paragraphs) 
<br>
<br>
<br>

# Contents
***

<br>

[**Quickstart Guide**](#quickstart-guide)

- [Preparing the Data](#preparing-the-data)
- [Preparing the Workflow](#preparing-the-workflow)
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

<br>
<br>

# Quickstart Guide

***

<br>

Great! You've made it to the SANS quickstart guide!

In this section we would like to give a minimal example of the SANS *ambages* workflow to help inexperienced users get started without requireing any prior knowledge of the software or the CloWM environment. 

The three steps below will show you how to access data from your local system, run a minimalist workflow and download the computed phylogeny.

<br>

## Preparing the Data
Starting out, we need to upload the input data, to make it accessible in the CloWM environment.
To achieve this we will take a look at how to [create buckets](#1-creating-buckets), [add folders](#2-creating-folders) to them, and [upload local data](#uploading-local-data). 

### 1 Creating Buckets  

<br>

![QS_Preparing_Data: 1_1](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_1.png/raw?ref=clowm-integration-extensions)

Look for the to the `Files` tab at the top left of your browser and select `My Data Buckets`.

<br>
<br>

![QS_Preparing_Data: 1_2](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_2.png/raw?ref=clowm-integration-extensions)

In this tab, you can inspect and mange your data buckets. For the first execution of this workflow, we want to create a new empty bucket. To do so, click an on the `+` sign at the top left of your browser.

<br>
<br>

![QS_Preparing_Data: 1_3](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_3.png/raw?ref=clowm-integration-extensions)

The occuring pop-up window allows to select a name for the new bucket. Please notice that this name has to be unique to ensure that the bucket can be identified afterwards. 
The entry below allows to add an description of our new containers content.
Once we are satisfied, clicking the `save` button will create the new bucket.

<br>
<br>

![QS_Preparing_Data: 1_4](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_4.png/raw?ref=clowm-integration-extensions)

The newly created bucket is automatically selected in the bucket overview.
At this point, it is possible to start uploading our data into the new bucket. However, to sepparate the workflows input data and produced output, it is usefull to create sepparate folders first.

<br>
<br>

### 2 Creating Folders

<br>

To create a new folder, click on the `+Folder` button on the top right of your browser.

<br>
<br>

![QS_Preparing_Data: 1_5](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_5.png/raw?ref=clowm-integration-extensions)

The pop-up allows to set a name for the folder.
As this folder will hold our input data, it is named
*Quickstart_Data*, here. Again, clicking the `save` button will create the folder.

<br>
<br>

![QS_Preparing_Data: 1_6](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_6.png/raw?ref=clowm-integration-extensions)

Repeating this process, we create a folder for the output of the workflow. It is named *Quickstart_Output* here.

<br>
<br>

![QS_Preparing_Data: 1_7](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_7.png/raw?ref=clowm-integration-extensions)

Both new folders should should appear in the bucket overview. With these folders set up, we can now move on to uploading the target files.

<br>
<br>

### 3 Uploading Local Data

<br>

To upload local files into the newly created *Quickstart_Data* folder, we need to select it, by clicking on its entry in the folder list.

![QS_Preparing_Data: 1_8](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_8.png/raw?ref=clowm-integration-extensions)

To check whether we are at the right position, we can take a look at the center top of the browser, where our current position in the file tree is shown. It shows the selected *sans-tutorial* bucket 
and the *Quickstart_Data* folder that we are currently inspecting.

 Ensured that we are at the right folder, we can start uploading, by clicking on the `UploadFile` button at the top right of the browser.

<br>
<br>

![QS_Preparing_Data: 1_9](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_9.png/raw?ref=clowm-integration-extensions)

Inside of the pop-up window, we click on the browse option to open the file selector. 

<br>
<br>

![QS_Preparing_Data: 1_10](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_10.png/raw?ref=clowm-integration-extensions)

By navigating to the first of our input files and clicking the `open` button, we can start the upload of the first input file. In this example, the file uses the *fasta* format. More details on other supported formats is provided in the parameter section.

<br>
<br>

![QS_Preparing_Data: 1_11](https://gitlab.ub.uni-bielefeld.de/api/v4/projects/1302/repository/files/docs%2Fresource%2Fquickstart%2F1_11.png/raw?ref=clowm-integration-extensions)



By repeating this process for the other input files, we can add all of our target data to the folder, completing the data preparation.

<br>
<br>
<br>

## Preparing the Workflow

Important sentence!!



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