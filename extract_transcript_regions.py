#!/usr/bin/env python

# Stephen N. Floor
# 7 October 2014
# floor@berkeley.edu


### Modified by Shaojun Xie
### 2/21/2018


import sys, os, argparse
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import GTF

from collections import defaultdict 

from Transcript import *

print " ----------------------------------"
print "| Extract Regions from annotations |"
print "|  snf        Fall 2014            |"
print " ----------------------------------\n\n"


# ------ ARGUMENT PARSING ----------

parser = argparse.ArgumentParser(description="Create transcript regions (5' UTR/CDS/3'UTR etc) from knownGenes or a GTF") 

parser.add_argument("-i", "--input", help="input filename", required=True)
parser.add_argument("-o", "--output", help="output basename", required=True) 
parser.add_argument("-s", "--species", help="species name", default="mm10")

args = parser.parse_args() 
tx_id = args.species

#if ( (not (args.ucsc or args.gtf)) or (args.ucsc and args.gtf)):
#    sys.exit("FATAL: must set one but not both of --ucsc and --gtf") 

# output filenames:
utr5FName = args.output + "_5utr.bed"
utr5StartFName = args.output + "_5utr_start.bed"
cdsFName = args.output + "_cds.bed"
utr3FName = args.output + "_3utr.bed"
exonFName = args.output + "_exons.bed"
intronFName = args.output + "_introns.bed"
codingExonFName = args.output + "_codingexons.bed"
codingIntronFName = args.output + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA
noncodingExonFName = args.output + "_noncodingexons.bed"
noncodingIntronFName = args.output + "_noncodingintrons.bed"

promoterFName = args.output + "_promoter.bed"
p15kFName = args.output + "_p15k.bed"

#keep track of where we are
genesRead = 0


#terminate if output files exist

if os.path.exists(utr5FName) or os.path.exists(utr5StartFName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) or os.path.exists(promoterFName) or os.path.exists(p15kFName) \
        or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
    sys.exit("ERROR: output basename %s files already exist" % args.output)

#process the file

with open(utr5FName, "w") as utr5File, open(utr5StartFName, "w") as utr5StartFile, open(cdsFName, "w") as cdsFile, \
        open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
        open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
        open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile, \
        open(promoterFName, "w") as promoterFile, open(p15kFName, "w") as p15kFile:

    def writeOutput(gene):
        if(gene.coding):
            for entry in gene.bedFormat(region="5utr", tx_id=tx_id):
                utr5File.write(entry + "\n")
            for entry in gene.bedFormat(region="5utr_start", tx_id=tx_id):
                utr5StartFile.write(entry + "\n")
            for entry in gene.bedFormat(region="cds", tx_id=tx_id):
                cdsFile.write(entry + "\n")
            for entry in gene.bedFormat(region="3utr", tx_id=tx_id):
                utr3File.write(entry + "\n")

            for entry in gene.bedFormat(region="exons", tx_id=tx_id):
                exonFile.write(entry + "\n")
                codingExonFile.write(entry + "\n")

            for entry in gene.bedFormat(region="introns", tx_id=tx_id):
                intronFile.write(entry + "\n")
                codingIntronFile.write(entry + "\n")

        else: # noncoding transcripts just have exons and introns
            for entry in gene.bedFormat(region="exons", tx_id=tx_id):
                exonFile.write(entry + "\n")
                noncodingExonFile.write(entry + "\n")

            for entry in gene.bedFormat(region="introns", tx_id=tx_id):
                intronFile.write(entry + "\n")
                noncodingIntronFile.write(entry + "\n")

        for entry in gene.bedFormat(region="promoter", tx_id=tx_id):
            promoterFile.write(entry + "\n")
        for entry in gene.bedFormat(region="p1-5k", tx_id=tx_id):
            p15kFile.write(entry + "\n")


    txDict = defaultdict(list) 

    print "Building GTF dictionary..." 

    # the issue here is that lines for various transcripts may be interleaved, so can either create lots of objects, or a giant dict. opted for giant dict. 
    for line in GTF.lines(args.input): 
        # only want to read in lines corresponding to these features
        if line["feature"] in ["exon", "CDS", "start_codon", "stop_codon"]:
            txDict[line["transcript_id"]].append(line)
            genesRead += 1

            if (not genesRead % 25000):
                print "\tProcessed %d lines..." %  genesRead

    print "Dictionary built." 

    print "Writing transcript properties."
    genesRead = 0
    
    # now create a Transcript object for each transcript and output it 

    for key in txDict: 

        #print key

        tx = createGTFTranscript(txDict[key])

        #print tx 
        writeOutput(tx)
        genesRead += 1
        
        if (not genesRead % 2500):
            print "\tProcessed %d entries..." %  genesRead


print "Processed %d entries." %  genesRead


