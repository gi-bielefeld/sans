fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.inputdir + "/" + it },type : "file")

params.label = "$projectDir/clowm/NO_FILE"


process sans {
  container "ghcr.io/gi-bielefeld/sans:latest"
  publishDir params.outdir, mode: 'symlink'
  debug false
  label 'highmemMedium'

  input:
    path inputFiles
    file label
  output:
    tuple path("splits.*"), path("genomeList.txt")

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  if [ ${params.pdf ? "1" : "0"} -eq 1 ]; then
    /usr/bin/Xvfb &
  fi
  #if [ ${params.label} ]; then
  #  touch ${params.label}
  #fi
  SANS -i genomeList.txt -o splits.tsv ${ params.verbose ? "--verbose" : "" } --mean ${ params.mean } --kmer ${ params.kmer } ${ params.top != null ? "--top ${ params.top }" : "" } ${ params.filter != null ? "--filter ${ params.filter }" : "" } ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } --threads ${ task.cpus } ${ params.pdf ? "--pdf splits.pdf" : "" } ${ params.amino ? "--amino" : "" } ${ params.code ? "--code" : "" } ${ label.name != 'NO_FILE' ? "--label $opt" : '' }
  """
}

workflow {
  opt_label = file(params.label, checkIfExists:true)
  sans(inputChannel.collect(),opt_label)
}
