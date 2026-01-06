// --------------------------------------------------
// SANS process
// --------------------------------------------------
process SANS {
    label 'highmemMedium'
    container "ghcr.io/gi-bielefeld/sans:v1.0.8"
    debug false

    input:
    path inputFiles
    path label
    path label_colors
    path fof
    path blacklist

    output:
    path 'sans_*', optional: true, emit: results
    path 'sans.*', optional: true, emit: logs

    script:
    """
    touch ${inputFiles}
    echo "${inputFiles}" | tr " " "\\n" > genomeList.txt

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
    ${ params.consensus == "strict" ? "--newick sans_tree.newick" : '' } \
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
    ${ params.pdf  && params.consensus == "strict" ? "--svg sans_tree.svg" : "" } \
    ${ params.pdf  && params.consensus != "none" && params.consensus != "strict" ? "--svg sans_splitnetwork.svg" : "" } \
    ${ params.svg  && params.consensus == "none" && params.filter == 'strict' ? "--svg sans_tree.svg" : "" } \
    ${ params.svg  && params.consensus == "none" && params.filter != 'strict' ? "--svg sans_splitnetwork.svg" : "" } \
    ${ params.svg  && params.consensus == "strict" ? "--svg sans_tree.svg" : "" } \
    ${ params.svg  && params.consensus != "none" && params.consensus != "strict" ? "--svg sans_splitnetwork.svg" : "" } \
    ${ params.filter == 'strict' ? "--filter strict" : '' } \
    ${ params.filter == 'weakly' ? "--filter weakly" : '' } \
    ${ params.filter == 'planar' ? "--filter planar" : '' } \
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
    --threads ${ task.cpus } \
    --stats sans.stats\"

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

// --------------------------------------------------
// ZIP extractor
// --------------------------------------------------
process UNZIP {
    container 'python:3.14-slim'

    input:
    path zipgenomes

    output:
    path 'output/*', emit: unzipped

    script:
    """
    #!/usr/bin/env python

    import zipfile
    with zipfile.ZipFile("${zipgenomes}") as z:
        z.extractall("output")
    """
}

// --------------------------------------------------
// TAR extractor
// --------------------------------------------------
process UNTAR {
    container 'python:3.14-slim'

    input:
    path targenomes

    output:
    path 'output/*', emit: untarred

    script:
    """
    #!/usr/bin/env python

    import tarfile
    with tarfile.open("${targenomes}") as t:
        t.extractall("output")
    """
}

// --------------------------------------------------
// Workflow definition
// --------------------------------------------------
workflow {
    main:
    opt_label        = params.label == null ? projectDir.resolve("NO_FILE") : file(params.label, checkIfExists:true)
    opt_label_colors = params.label_colors == null ? projectDir.resolve("NO_FILE2") : file(params.label_colors, checkIfExists:true)
    opt_fof          = params.file_of_files == null ? projectDir.resolve("NO_FILE3") : file(params.file_of_files, checkIfExists:true)
    opt_blacklist    = params.blacklist == null ? projectDir.resolve("NO_FILE4") : file(params.blacklist, checkIfExists:true)

    def fileEndingList = [
        "*.fa","*.fa.gz","*.fasta","*.fasta.gz",
        "*.faa","*.faa.gz","*.fna","*.fna.gz",
        "*.mpfa","*.mpfa.gz","*.ffn","*.ffn.gz",
        "*faa.gz", "*.fastq","*.fastq.gz","*.fq","*.fq.gz"
    ]

    if (params.input.endsWith(".zip")) {
        UNZIP(file(params.input, checkIfExists: true))
        SANS(UNZIP.out.unzipped, opt_label, opt_label_colors, opt_fof, opt_blacklist)
    }
    else if (params.input.endsWith(".tar.gz")) {
        UNTAR(file(params.input, checkIfExists: true))
        SANS(UNTAR.out.untarred, opt_label, opt_label_colors, opt_fof, opt_blacklist)
    }
    else {
        inputChannel = Channel
        .fromPath(fileEndingList.collect { "${params.input}/${it}" }, type: "file")
        .ifEmpty(projectDir.resolve("NO_FILE5"))
        .collect()
        SANS(inputChannel, opt_label, opt_label_colors, opt_fof, opt_blacklist)
    }

    publish:
    sans_results = SANS.out.results
    sans_logs = SANS.out.logs
}

output {
    sans_results {
        path '.'
    }
    sans_logs {
        path '.'
    }
}
