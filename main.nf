fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.input + "/" + it },type : "file")

params.label = "$projectDir/clowm/NO_FILE"
params.file_of_files = "$projectDir/clowm/NO_FILE"

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
    path 'splits.tsv'
    path 'genomeList.txt'
    path 'splits.pdf', optional: true

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  if [ ${params.pdf ? "1" : "0"} -eq 1 ]; then
    /usr/bin/Xvfb &
  fi
  SANS ${ fof.name != 'NO_FILE' ? "--input $fof" : '--input genomeList.txt' } -o splits.tsv ${ params.verbose ? "--verbose" : "" } --mean ${ params.mean } --kmer ${ params.kmer } ${ params.top != null ? "--top ${ params.top }" : "" } ${ params.filter != null ? "--filter ${ params.filter }" : "" } ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } --threads ${ task.cpus } ${ params.pdf ? "--pdf splits.pdf" : "" } ${ params.amino ? "--amino" : "" } ${ params.code ? "--code" : "" } ${ label.name != 'NO_FILE' ? "--label $label" : '' }
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
  with ipfile.ZipFile("${zipgenomes }", 'r') as zip_ref:
    zip_ref.extractall("output")
  """
}

workflow {
  opt_label = file(params.label, checkIfExists:true)
  opt_fof = file(params.file_of_files, checkIfExists:true)
  if (params.input.endsWith(".zip")) {
    unzip(params.input)
    sans(unzip.output,opt_label,opt_fof)
  } else {
    sans(inputChannel.collect(),opt_label,opt_fof)
  }
}



