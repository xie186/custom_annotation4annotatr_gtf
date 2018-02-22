extract-transcript-regions
================================

For bioinformatics on genome annotation sets. Ensembl-derived GTF into transcript regions (i.e. exons, introns, UTRs and CDS).

This script is modified from https://github.com/stephenfloor/extract-transcript-regions to output files used for generating custome annotation information for annotatr (bioconductor.org/packages/annotatr).

### Input

Ensembl GTF. For this to work properly the knownGenes or input GTF *must* define the CDS either using thickStart/thickEnd (knownGenes) or with CDS/start_codon/stop_codon directives (GTF). 
  
### Output
  Region files in .bed format 
  
### Usage
```
 python ../../../../scripts/custom_annotation4annotatr_gtf/extract_transcript_regions.py  -h
 ----------------------------------
| Extract Regions from annotations |
 ----------------------------------


usage: extract_transcript_regions.py [-h] -i INPUT -o OUTPUT [-s SPECIES]

Create transcript regions (5' UTR/CDS/3'UTR etc) from knownGenes or a GTF

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input filename
  -o OUTPUT, --output OUTPUT
                        output basename
  -s SPECIES, --species SPECIES
                        species name

```
