fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.inputdir + "/" + it },type : "file")

process sans {
  container "ghcr.io/gi-bielefeld/sans:latest"
  publishDir params.outdir, mode: 'symlink'
  debug false
  label 'highmemMedium'

  input:
    path inputFiles
  output:
    tuple path("splits.*"), path("genomeList.txt")

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  if [ ${params.pdf ? "1" : "0"} -eq 1 ]; then
    /usr/bin/Xvfb &
  fi
  SANS -i genomeList.txt -o splits.tsv ${ params.verbose ? "-v" : "" } --mean ${ params.mean } --kmer ${ params.kmer } ${ params.top != null ? "--top ${ params.top }" : "" } ${ params.filter != null ? "--filter ${ params.filter }" : "" } ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } --threads ${ task.cpus } ${ params.pdf ? "--pdf splits.pdf" : "" }
  """
}

workflow {
  sans(inputChannel.collect())
}
