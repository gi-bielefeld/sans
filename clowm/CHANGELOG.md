# SANS Changelog



## v1.0.4 [2025-02-11]

### Fixed

- last codon in coding sequence translation got lost

### Added

- example data set


## v1.0.3 [2025-02-07]

### Fixed

- consensus planar was missing in parameter definitions
- some filtering combinations with tree/strict were not handled properly

## v1.0.2 [2025-02-0]

### Fixed

- parameter translate did not work (due to typo)

## v1.0.1 [2025-02-06]

### Added

- support further fasta file endings

## v1.0.0 [2025-01-30]

### Added

- made public

### Fixed

- too many parameters for second filter step in default mode
- check if input is empty


## v0.9.4 [2025-01-29]

### Added

- new filter option: planar

## v0.9.3 [2025-01-24]

### Added

- optional SVG output

### Fixed

- PDF generated via SVG to obtain a PDF that is MacOS compatible

### Changed

- no SANS version check anymore

## v0.9.2 [2025-01-22]

### Fixed

- paramter descriptions, validation and error handling


## v0.9.1 [2025-01-20]

### Added

- further parameter defaults in nextflow.config

### Changed

- moved main.nf into clowm subfolder

### Fixed

- graphics in usage.md



## v0.9 [2025-01-20]

Second release of SANS ambage on CloWM staging.

### Added

- SANS parameters: file-of-files, amino, translate, code, label, colors, nexus, newick, filter, bootstrapping, consensus, iupac, noref, core, blacklist, threads
- log ouput
- automated re-compiling to adjust maximum number of input files
- input files in a folder or as zip or tar.gz
- simple default mode: network und tree

### Changed

- pdf generation: Splitstree called from SANS
- documentation
