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
    path 'sans_splitnetwork.pdf'
    path 'sans_splitnetwork.nexus'
    path 'sans_splitnetwork.tsv'
    path 'sans_tree.pdf', optional: true
    path 'sans_tree.newick', optional: true
    path 'sans_tree.tsv', optional: true
    path 'sans.log'

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  if [ ${params.pdf ? "1" : "0"} -eq 1 ]; then
    /usr/bin/Xvfb &
  fi
  SANS ${ fof.name != 'NO_FILE2' ? "--input $fof" : '--input genomeList.txt' } -o sans_splitnetwork.tsv --verbose  --kmer ${ params.kmer } ${ params.top != null ? "--top ${ params.top }" : "" } --filter weakly ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } --threads ${ task.cpus } ${ params.pdf ? "--pdf sans_splitnetwork.pdf" : "" } ${ params.amino ? "--amino" : "" } ${ params.code ? "--code" : "" } ${ label.name != 'NO_FILE' ? "--label $label" : '' } --nexus sans_splitnetwork.nexus 2>>&1 | grep -v "(genome" | grep -v "%)" > sans.log
  if [ ${params.tree ? "1" : "0"} -eq 1 ]; then
    SANS ${ fof.name != 'NO_FILE2' ? "--input $fof" : '--input genomeList.txt' } -s sans_splitnetwork.tsv -o sans_tree.tsv --verbose  --kmer ${ params.kmer }  --filter strict --threads ${ task.cpus } ${ params.pdf ? "--pdf sans_tree.pdf" : "" } ${ label.name != 'NO_FILE' ? "--label $label" : '' } --newick sans_tree.newick 2>>&1 | grep -v "(genome" | grep -v "%)" >> sans.log
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



