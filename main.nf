fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.inputdir + it },type : "file")

process sans {
  container "gitlab.ub.uni-bielefeld.de:4567/gi/sans:latest"
  //container "localhost/sans:latest"
  publishDir params.outdir, mode: 'symlink'
  shell '/bin/sh', '-euo', 'pipefail'

  input:
    path inputFiles
  output:
    path "splits.tsv"

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > tempfile.txt
  SANS -i tempfile.txt -o splits.tsv ${ params.verbose ? "-v" : "" } --mean ${ params.mean } --kmer ${ params.kmer } ${ params.top != null ? "--top ${ params.top }" : "" } ${ params.filter != null ? "--filter ${ params.filter }" : "" } ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" }
  """
}

workflow {
   inputChannel | collect | sans
}