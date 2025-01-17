fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.input + "/" + it },type : "file")

params.label = "$projectDir/clowm/NO_FILE"
params.label_colors = "$projectDir/clowm/NO_FILE2"
params.file_of_files = "$projectDir/clowm/NO_FILE3"
params.blacklist = "$projectDir/clowm/NO_FILE4"

process sans {
  container "ghcr.io/gi-bielefeld/sans:latest"
  publishDir params.outdir, mode: 'symlink'
  debug false
  label 'highmemMedium'

  input:
    path inputFiles
    file label
    file label_colors
    file fof
    file blacklist

  output:
    path 'sans_splitnetwork.pdf', optional: true
    path 'sans_splitnetwork.nexus', optional: true
    path 'sans_splitnetwork.tsv', optional: true
    path 'sans_tree.pdf', optional: true
    path 'sans_tree.newick', optional: true
    path 'sans_tree.tsv', optional: true
    path 'sans_tree.tsv.bootstrap', optional: true
    path 'sans_splitnetwork.tsv.bootstrap', optional: true
    path 'sans_core.fasta', optional: true
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


  
  SANS_PARAMS=\"\
  ${ params.consensus == "none" && params.filter == 'strict' ? "--output sans_tree.tsv" : '' } \
  ${ params.consensus == "none" && params.filter != 'strict' ? "--output sans_splitnetwork.tsv" : '' } \
  ${ params.consensus == "strict" ? "--output sans_tree.tsv" : '' } \
  ${ params.consensus != "none" && params.consensus != "strict"? "--output sans_splitnetwork.tsv" : '' } \
  ${ params.consensus == "none" && params.filter == 'strict' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "none" && params.filter == '2-tree' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "none" && params.filter == '3-tree' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "tree" ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "2-tree" ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "3-tree" ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "none" && params.filter == 'weakly' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "none" && params.filter == '2-tree' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "none" && params.filter == '3-tree' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "weakly" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "2-tree" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "3-tree" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } \
  ${ fof.name != 'NO_FILE3' ? "--input $fof" : '--input genomeList.txt' } \
  ${ params.amino ? "--amino" : "" } \
  ${ params.translate ? "--code ${ params.code }" : "" } \
  --kmer ${ params.kmer } \
  ${ label.name != 'NO_FILE' && label_colors.name == 'NO_FILE2' ? "--label $label" : '' } \
  ${ label.name != 'NO_FILE' && label_colors.name != 'NO_FILE2' ? "--label $label $label_colors" : '' } \
  ${ params.top != null ? "--top ${ params.top }" : "" } \
  ${ params.pdf ? "--pdf sans_splitnetwork.pdf" : "" } \
  ${ params.filter == 'strict' ? "--filter strict" : '' } \
  ${ params.filter == 'weakly' ? "--filter weakly" : '' } \
  ${ params.filter == '2-tree' ? "--filter 2-tree" : '' } \
  ${ params.filter == '3-tree' ? "--filter 3-tree" : '' } \
  ${ params.filter == 'default' ? "--filter weakly" : '' } \
  ${ params.bootstrapping != null ? "--bootstrapping ${ params.bootstrapping } ${ params.support }" : "" } \
  ${ params.consensus != "none" ? "--consensus ${ params.consensus }" : "" } \
  ${ params.iupac != 0 ? "--iupac ${ params.iupac }" : "" } \
  ${ params.norev ? "--norev" : "" } \
  ${ params.mean != "geom2" ? "--mean ${ params.mean }" : "" } \
  ${ params.core ? "--core sans_core.fasta" : "" } \
  ${ blacklist.name != 'NO_FILE4' ? "--blacklist $blacklist" : "" } \
  --verbose \
  --threads ${ task.cpus }\"

  echo SANS \$SANS_PARAMS > sans.log

  SANS \$SANS_PARAMS 2>&1 | grep -v \"Fontconfig error\" | awk -F \"\r\" '{print \$NF}' >> sans.log
  
  if [ ${params.filter} == "default" ] && [ ${params.tree} ? "1" : "0" -eq 1 ]; then

    echo \"\" >> sans.log

    SANS_PARAMS=\"--splits sans_splitnetwork.tsv \
    --output sans_tree.tsv \
    --newick sans_tree.newick \
    ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } \
    ${ fof.name != 'NO_FILE3' ? "--input $fof" : '--input genomeList.txt' } \
    ${ params.amino ? "--amino" : "" } \
    ${ params.translate ? "--code ${ params.code }" : "" } \
    --kmer ${ params.kmer } \
    ${ label.name != 'NO_FILE' && label_colors.name == 'NO_FILE2' ? "--label $label" : '' } \
    ${ label.name != 'NO_FILE' && label_colors.name != 'NO_FILE2' ? "--label $label $label_colors" : '' } \
    ${ params.top != null ? "--top ${ params.top }" : "" } \
    ${ params.pdf ? "--pdf sans_tree.pdf" : "" } \
    --filter strict
    ${ params.iupac != 0 ? "--iupac ${ params.iupac }" : "" } \
    ${ params.norev ? "--norev" : "" } \
    ${ params.mean != "geom2" ? "--mean ${ params.mean }" : "" } \
    ${ params.core ? "--core sans_core.fasta" : "" } \
    ${ blacklist.name != 'NO_FILE4' ? "--blacklist $blacklist" : "" } \
    --verbose \
    --threads ${ task.cpus }\"
    
    echo SANS \$SANS_PARAMS >> sans.log
 
    SANS \$SANS_PARAMS 2>&1 | grep -v \"Fontconfig error\" | awk -F \"\r\" '{print \$NF}' >> sans.log

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
  opt_label_colors = file(params.label_colors, checkIfExists:true)
  opt_fof = file(params.file_of_files, checkIfExists:true)
  opt_blacklist = file(params.blacklist, checkIfExists:true)
  if (params.input.endsWith(".zip")) {
    unzip(params.input)
    sans(unzip.output,opt_label,opt_label_colors,opt_fof,opt_blacklist)
  } else if (params.input.endsWith(".tar.gz")) {
    untargz(params.input)
    sans(untargz.output,opt_label,opt_label_colors,opt_fof,opt_blacklist)
  } else {
    sans(inputChannel.collect(),opt_label,opt_label_colors,opt_fof,opt_blacklist)
  }
}



