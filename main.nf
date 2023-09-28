fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.inputdir + "/" + it },type : "file")

process sans {
  container "gitlab.ub.uni-bielefeld.de:4567/gi/sans:sans-latest"
  publishDir params.outdir, mode: 'symlink'
  debug false

  input:
    path inputFiles
  output:
    tuple path("splits.tsv"), path("genomeList.txt")

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  SANS -i genomeList.txt -o splits.tsv ${ params.verbose ? "-v" : "" } --mean ${ params.mean } --kmer ${ params.kmer } ${ params.top != null ? "--top ${ params.top }" : "" } ${ params.filter != null ? "--filter ${ params.filter }" : "" } ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" }
  """
}

process splitstree {
  container "gitlab.ub.uni-bielefeld.de:4567/gi/sans:splitstree-latest"
  publishDir params.outdir, mode: 'symlink'
  debug false

  input:
    path inputFiles
  output:
    path "splits.tsv.nexus*"
  
  script:
  """
  /usr/bin/Xvfb &
  sans2pdf.py ${inputFiles[0]} ${inputFiles[1]}
  """
}

workflow {
  sans(inputChannel.collect())
  if ( params.visualization ) {
    splitstree(sans.out)
  }
}
