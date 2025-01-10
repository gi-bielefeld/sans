fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.input + "/" + it },type : "file")

params.label = "$projectDir/clowm/NO_FILE"
params.file_of_files = "$projectDir/clowm/NO_FILE2"

process sans {
  container "ghcr.io/gi-bielefeld/sans:latest"
  publishDir params.outdir, mode: 'symlink'
  debug false
  label 'highmemMedium'

  input:
    path inputFiles
    file label
    file fof

  output:
    path 'sans_splitnetwork.pdf', optional: true
    path 'sans_splitnetwork.nexus', optional: true
    path 'sans_splitnetwork.tsv', optional: true
    path 'sans_tree.pdf', optional: true
    path 'sans_tree.newick', optional: true
    path 'sans_tree.tsv', optional: true
    path 'sans_tree.tsv.bootstrap', optional: true
    path 'sans_splitnetwork.tsv.bootstrap', optional: true
    path 'sans.log'

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  
  if [ ${params.pdf ? "1" : "0"} -eq 1 ]; then
    /usr/bin/Xvfb &
  fi
  
  if [ ${params.bootstrapping} != null ]; then
    if [ ${params.filter} == "none" ] || [ ${params.filter} == "default" ]; then
      echo "ERROR: For bootstrapping, you have to choose a filter criterioin using --filter." > sans.log;
      exit 1
    fi
  fi  
  
  SANS \
  ${ params.filter == 'strict' ? "--output sans_tree.tsv" : '--output sans_splitnetwork.tsv' } \
  ${ params.consensus == 'strict' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus =! null && params.consensus != 'strict' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == null && params.filter == 'strict' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == null && params.filter != 'strict' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } \
  ${ fof.name != 'NO_FILE2' ? "--input $fof" : '--input genomeList.txt' } \
  ${ params.amino ? "--amino" : "" } \
  ${ params.translate ? "--code ${ params.code }" : "" } \
  --kmer ${ params.kmer } \
  ${ label.name != 'NO_FILE' && label_colors == 'NO_FILE' ? "--label $label" : '' } \
  ${ label.name != 'NO_FILE' && label_colors != 'NO_FILE' ? "--label $label $label_colors" : '' } \
  ${ params.top != null ? "--top ${ params.top }" : "" } \
  ${ params.pdf ? "--pdf sans_splitnetwork.pdf" : "" } \
  ${ params.filter == 'strict' ? "--filter strict" : '' } \
  ${ params.filter == 'weakly' ? "--filter weakly" : '' } \
  ${ params.filter == '2-tree' ? "--filter 2-tree" : '' } \
  ${ params.filter == '3-tree' ? "--filter 3-tree" : '' } \
  ${ params.filter == 'default' ? "--filter weakly" : '' } \
  ${ params.bootstrapping != null ? "--bootstrapping ${ params.bootstrapping } ${ params.support }" : "" } \
  ${ params.consensus != null && params.consensus != "none" ? "--consensus ${ params.consensus }" : "" } \
  ${ params.iupac != 0 ? "--iupac ${ params.iupac }" : "" } \
  ${ params.norev ? "--norev" : "" } \
  ${ params.mean != "geom2" ? "--mean ${ params.mean }" : "" } \
  --verbose \
  --threads ${ task.cpus } \
  2>&1 | grep -v "(genome" | grep -v "%)" > sans.log
  
  
  if [ ${params.filter} == "default" ] && [ ${params.tree} ? "1" : "0" -eq 1 ]; then
    SANS \
    --splits sans_splitnetwork.tsv \
    --output sans_tree.tsv \
    --newick sans_tree.newick \
    ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } \
    ${ fof.name != 'NO_FILE2' ? "--input $fof" : '--input genomeList.txt' } \
    ${ params.amino ? "--amino" : "" } \
    ${ params.translate ? "--code ${ params.code }" : "" } \
    --kmer ${ params.kmer } \
    ${ label.name != 'NO_FILE' && label_colors == 'NO_FILE' ? "--label $label" : '' } \
    ${ label.name != 'NO_FILE' && label_colors != 'NO_FILE' ? "--label $label $label_colors" : '' } \
    ${ params.top != null ? "--top ${ params.top }" : "" } \
    ${ params.pdf ? "--pdf sans_tree.pdf" : "" } \
    --filter strict
    ${ params.iupac != 0 ? "--iupac ${ params.iupac }" : "" } \
    ${ params.norev ? "--norev" : "" } \
    ${ params.mean != "geom2" ? "--mean ${ params.mean }" : "" } \
    --verbose \
    --threads ${ task.cpus } \
    2>&1 | grep -v "(genome" | grep -v "%)" >> sans.log
  fi

  rm -f gegnomeList.txt
  """
}

process unzip {

  container 'python:3.12'
  input:
  path zipgenomes
  output:
  path 'output/*'
  
  script:
  """
  #!/usr/local/bin/python
  
  import zipfile
  with zipfile.ZipFile("${zipgenomes }", 'r') as zip_ref:
    zip_ref.extractall("output")
  """
}

process untargz {

  container 'python:3.12'
  input:
  path zipgenomes
  output:
  path 'output/*'
  
  script:
  """
  #!/usr/local/bin/python

  import tarfile 
  file = tarfile.open("${zipgenomes }") 
  file.extractall("output") 
  file.close()
  """
}


workflow {
  opt_label = file(params.label, checkIfExists:true)
  opt_fof = file(params.file_of_files, checkIfExists:true)
  if (params.input.endsWith(".zip")) {
    unzip(params.input)
    sans(unzip.output,opt_label,opt_fof)
  } else if (params.input.endsWith(".tar.gz")) {
    untargz(params.input)
    sans(untargz.output,opt_label,opt_fof)
  } else {
    sans(inputChannel.collect(),opt_label,opt_fof)
  }
}



