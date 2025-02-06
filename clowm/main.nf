fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.faa", "*.faa.gz", "*.fna", "*.fna.gz", "*.mpfa", "*.mpfa.gz", "*.ffn", "*.ffn.gz", "*.faa.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

inputChannel=Channel.fromPath(fileEndingList.collect { params.input + "/" + it },type : "file")

params.label = "$projectDir/NO_FILE"
params.label_colors = "$projectDir/NO_FILE2"
params.file_of_files = "$projectDir/NO_FILE3"
params.blacklist = "$projectDir/NO_FILE4"

process sans {
  container "ghcr.io/gi-bielefeld/sans:v1.0.2"
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
    path 'sans_splitnetwork.svg', optional: true
    path 'sans_splitnetwork.nexus', optional: true
    path 'sans_splitnetwork.tsv', optional: true
    path 'sans_tree.pdf', optional: true
    path 'sans_tree.svg', optional: true
    path 'sans_tree.newick', optional: true
    path 'sans_tree.tsv', optional: true
    path 'sans_tree.tsv.bootstrap', optional: true
    path 'sans_splitnetwork.tsv.bootstrap', optional: true
    path 'sans_core.fasta', optional: true
    path 'sans.log', optional: true
    path 'sans.err', optional: true
    

  script:
  """
  touch $inputFiles
  echo "$inputFiles" | tr " " "\n" > genomeList.txt
  
  if [ ${params.pdf ? "1" : "0"} -eq 1 ] || [ ${params.svg ? "1" : "0"} -eq 1 ]  || [ label.name != 'NO_FILE' ]; then
    /usr/bin/Xvfb &
  fi

  if [ \$( cat genomeList.txt | wc -l) -lt 2 ]; then
    echo \"ERROR: The input is either empty or fewer than two files.\" > sans.err;
    exit 0
  fi

 
  if [ ${params.bootstrapping} != null ]; then
    if [ ${params.filter} == "none" ] || [ ${params.filter} == "default" ]; then
      echo \"ERROR: For bootstrapping, you have to choose a filter criterion using --filter.\" > sans.err;
      exit 0
    fi
  fi  
  if [ ${params.consensus} != "none" ] && [ ${params.bootstrapping} == null ]; then
    echo \"ERROR: Filter on bootstrap values (--consensus) can only be chosen in combination with bootstrapping (--boostrapping).\" > sans.err;
    exit 0
  fi
  if [ ${params.support} != null ] && [ ${params.support} != "0" ] && [ ${params.bootstrapping} == null ]; then
    echo \"ERROR: Bootstrap support filter (--support) can only be chosen in combination with bootstrapping (--boostrapping).\" > sans.err;
    exit 0
  fi

  if [ ${params.pdf ? "1" : "0"} -eq 1 ] && [ ${params.svg ? "1" : "0"} -eq 0 ]; then
    echo "NOTE: In some cases, the PDF might not be properly readable. Thus, the CloWM version generates the PDF output by first generating an SVG file that is then converted to PDF." >> sans.log
  fi

  
  SANS_PARAMS=\"\
  ${ params.consensus == "none" && params.filter == 'strict' ? "--output sans_tree.tsv" : '' } \
  ${ params.consensus == "none" && params.filter != 'strict' ? "--output sans_splitnetwork.tsv" : '' } \
  ${ params.consensus == "strict" ? "--output sans_tree.tsv" : '' } \
  ${ params.consensus != "none" && params.consensus != "strict" ? "--output sans_splitnetwork.tsv" : '' } \
  ${ params.consensus == "none" && params.filter == 'strict' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "none" && params.filter == '2-tree' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "none" && params.filter == '3-tree' ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "tree" ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "2-tree" ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "3-tree" ? "--newick sans_tree.newick" : '' } \
  ${ params.consensus == "none" && params.filter == 'weakly' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "none" && params.filter == 'planar' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "none" && params.filter == '2-tree' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "none" && params.filter == '3-tree' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "none" && params.filter == 'default' ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "weakly" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "planar" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "2-tree" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.consensus == "3-tree" ? "--nexus sans_splitnetwork.nexus" : '' } \
  ${ params.qualify != null ? "--qualify ${ params.qualify }" : "" } \
  ${ fof.name != 'NO_FILE3' ? "--input $fof" : '--input genomeList.txt' } \
  ${ params.amino ? "--amino" : "" } \
  ${ params.translate ? "--code ${ params.code }" : "" } \
  ${ params.kmer != null ? "--kmer ${ params.kmer }" : "" } \
  ${ label.name != 'NO_FILE' && label_colors.name == 'NO_FILE2' ? "--label $label" : '' } \
  ${ label.name != 'NO_FILE' && label_colors.name != 'NO_FILE2' ? "--label $label $label_colors" : '' } \
  ${ params.top != null ? "--top ${ params.top }" : "" } \
  ${ params.pdf  && params.consensus == "none" && params.filter == 'strict' ? "--svg sans_tree.svg" : "" } \
  ${ params.pdf  && params.consensus == "none" && params.filter != 'strict' ? "--svg sans_splitnetwork.svg" : "" } \
  ${ params.pdf  && params.consensus == "tree" ? "--svg sans_tree.svg" : "" } \
  ${ params.pdf  && params.consensus != "none" && params.consensus != "strict" ? "--svg sans_splitnetwork.svg" : "" } \
  ${ params.svg  && params.consensus == "none" && params.filter == 'strict' ? "--svg sans_tree.svg" : "" } \
  ${ params.svg  && params.consensus == "none" && params.filter != 'strict' ? "--svg sans_splitnetwork.svg" : "" } \
  ${ params.svg  && params.consensus == "tree" ? "--svg sans_tree.svg" : "" } \
  ${ params.svg  && params.consensus != "none" && params.consensus != "strict" ? "--svg sans_splitnetwork.svg" : "" } \
  ${ params.filter == 'strict' ? "--filter strict" : '' } \
  ${ params.filter == 'weakly' ? "--filter weakly" : '' } \
  ${ params.filter == 'planar' ? "--filter weakly" : '' } \
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

  echo SANS \$SANS_PARAMS >> sans.log

  SANS-autoN.sh \$SANS_PARAMS 2>&1 | grep -v \"Fontconfig error\" | awk -F \"\r\" '{print \$NF}' >> sans.log
  
  if [ ${params.filter} == "default" ] && [ ${params.tree ? "1" : "0"} -eq 1 ]; then

    echo \"\" >> sans.log

    SANS_PARAMS=\"--splits sans_splitnetwork.tsv \
    --output sans_tree.tsv \
    --newick sans_tree.newick \
    ${ fof.name != 'NO_FILE3' ? "--input $fof" : '--input genomeList.txt' } \
    ${ label.name != 'NO_FILE' && label_colors.name == 'NO_FILE2' ? "--label $label" : '' } \
    ${ label.name != 'NO_FILE' && label_colors.name != 'NO_FILE2' ? "--label $label $label_colors" : '' } \
    ${ params.top != null ? "--top ${ params.top }" : "" } \
    ${ params.pdf ? "--svg sans_tree.svg" : "" } \
    ${ params.svg ? "--svg sans_tree.svg" : "" } \
    --filter strict
    ${ params.mean != "geom2" ? "--mean ${ params.mean }" : "" } \
    --verbose \
    --threads ${ task.cpus }\"
    
    echo SANS \$SANS_PARAMS >> sans.log
 
    SANS-autoN.sh \$SANS_PARAMS 2>&1 | grep -v \"Fontconfig error\" | awk -F \"\r\" '{print \$NF}' >> sans.log

  fi

  if [ ${params.pdf ? "1" : "0"} -eq 1 ] && [ ${params.svg ? "1" : "0"} -eq 0 ]; then
    if [ -f "sans_tree.svg" ]; then cairosvg sans_tree.svg -o sans_tree.pdf; rm sans_tree.svg; fi
    if [ -f "sans_splitnetwork.svg" ]; then cairosvg sans_splitnetwork.svg -o sans_splitnetwork.pdf; rm sans_splitnetwork.svg; fi
  fi

  rm -f genomeList.txt
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
    sans(inputChannel.ifEmpty(file("$projectDir/NO_FILE5")).collect(),opt_label,opt_label_colors,opt_fof,opt_blacklist)
  }
}



