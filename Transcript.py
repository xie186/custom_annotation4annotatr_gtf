# general class containing info about a transcribed region.  Can come from UCSC Knowngenes (BED) or Ensembl GTF files currently 

# Stephen N. Floor
# Fall 2014 

class Transcript: 
    def __init__(self):
        #properties defined in UCSC knowngenes 
        self.name = ''
        self.symbol = '' 
        self.chrom = ''
        self.strand = '' 
        self.txStart = 0
        self.txEnd = 0 
        self.cdsStart = 0
        self.cdsEnd = 0
        self.exonCt = 0
        self.exonStarts = []
        self.exonEnds = []
        self.exonLengths = []
        
        #meta properties to be computed during construction.  these are lists of BED first four field tuples with the exception of Len terms which are the length of the total region for the gene 
        self.utr5 = []
        self.utr5Len = 0
        self.utr5start = []
        self.utr5startLen = 0
        self.cds = []
        self.cdsLen = 0 
        self.utr3 = []
        self.utr3Len = 0
        self.exons = []
        self.exonsLen = 0
        self.introns = []
        self.intronsLen = 0
        self.promoter = []  ## promoter 1kb
        self.p15k = []  ### 1-5k upstream of promoter
        self.d1k = []  ## downstream 1k  
        self.coding = False


    def __str__(self):  #currently roughly knownGenes format with a second line containing metadata 
        return "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s" % (self.name, self.chrom, self.strand, self.txStart, self.txEnd, self.cdsStart, self.cdsEnd, self.exonCt, self.exonStarts, self.exonEnds, self.utr5, self.utr5Len, self.cds, self.cdsLen, self.utr3, self.utr3Len, self.exons, self.exonsLen, self.introns, self.intronsLen, self.coding)
    
#BED format output is goal.  Fields are optional after featureEnd 
# chrom    featureStart   featureEnd   nameOfLine   score(0-1000)   strand   thickStart  thickEnd  itemRGBtuple  blockCount  blockSizes   blockStarts 

#this function returns a list of BED-formatted strings for the feature passed as region with multiple entries per region possible, one for each primitive (exon/intron) 
    def bedFormat(self, region="exons", tx_id = "mm10"):
        if (not self.coding and (region == "5utr" or region == "cds" or region == "3utr")):
            print "Transcript.py bedFormat error: noncoding transcripts do not have 5utr/cds/3utr"
            return []

        returnVal = []

        if (region == "5utr"):
            for chunk in self.utr5:
                returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (chunk[0], chunk[1], chunk[2], self.name, self.symbol, tx_id))

        elif (region == "5utr_start"):
            for chunk in self.utr5start:
                returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (chunk[0], chunk[1], chunk[2], self.name, self.symbol, tx_id))

        elif (region == "cds"):
            for chunk in self.cds:
                returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (chunk[0], chunk[1], chunk[2], self.name, self.symbol, tx_id))

        elif (region == "3utr"):
            for chunk in self.utr3:
                returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (chunk[0], chunk[1], chunk[2], self.name, self.symbol, tx_id))

        elif (region == "exons"):
            for chunk in self.exons:
                returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (chunk[0], chunk[1], chunk[2], self.name, self.symbol, tx_id))

        elif (region == "introns"):
            for chunk in self.introns:
                returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (chunk[0], chunk[1], chunk[2], self.name, self.symbol, tx_id))
        elif (region == "promoter"):
            #print("%s\t%d\t%d\t%s\t%s\t%s" % (self.chrom, self.promotStart, self.promotEnd, self.name, self.symbol, tx_id))
            returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (self.chrom, self.promotStart, self.promotEnd, self.name, self.symbol, tx_id))
        elif (region == "p1-5k"):
            returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (self.chrom, self.p15kStart, self.p15kEnd, self.name, self.symbol, tx_id))
        elif (region == "d1k"):
            returnVal.append("%s\t%d\t%d\t%s\t%s\t%s" % (self.chrom, self.d1kStart, self.d1kEnd, self.name, self.symbol, tx_id))

        else:
            print "Transcript.py bedFormat error: currently only regions 5utr/cds/3utr/exons/introns are supported"
            

        return returnVal

    def computeMetadata(self): 
    # -- begin computing metadata -- 

    # -- note: chose clarity of code and conditionals here over most efficient computation (i.e. some clauses may be redundant)

        if (self.strand == "+"): 
            ### Promoter 
            self.promotStart = self.txStart -1000 if self.txStart -1000 > 0 else 0
            self.promotEnd = self.txStart -1 if self.txStart -1 > 0 else 0
            ### 1-5K 
            self.p15kStart = self.txStart -5000 if  self.txStart -5000 > 0 else 0
            self.p15kEnd = self.txStart -1000 if  self.txStart -1000 >0 else 0
            ### d1k
            self.d1kStart = self.txEnd + 1
            self.d1kEnd = self.txEnd + 1000
        #print ("DBUG - exonCt %d i %d exonEnds[i] %d cdsStart %d exonStarts[i] %d cdsEnd %d") % \
            #    (self.exonCt, i, self.exonEnds[i], self.cdsStart, self.exonStarts[i], self.cdsEnd)
            for i in range (self.exonCt): 
                if (self.cdsStart != self.cdsEnd): # if this is a coding transcript
                    self.coding = True
                # -- first compute 5'utr, CDS, 3'utr regions --
                #case 1 - exon spans 5' UTR/CDS/3' UTR
                    if (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] > self.cdsEnd):
                        self.utr5.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr5Len += self.cdsStart - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name)) # for now just append the 5' utr exons to the utr5start 
                        self.utr5startLen += self.cdsStart - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.cdsStart
                        self.utr3.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.cdsEnd
                #case 2 - exon spans 5' UTR/CDS junction
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] >= self.cdsStart):
                        self.utr5.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr5Len += self.cdsStart - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name)) 
                        self.utr5startLen += self.cdsStart  - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i]- self.cdsStart
                #case 3 - exon spans CDS/3'UTR junction 
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonStarts[i] <= self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.exonStarts[i]
                        self.utr3.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.cdsEnd
                #case 4 - exon is 5' UTR only 
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] < self.cdsStart): 
                        self.utr5.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name)) 
                        self.utr5startLen += self.exonEnds[i] - self.exonStarts[i]
                #case 5 - exon is CDS only
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonEnds[i] <= self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i] - self.exonStarts[i]
                #case 6 - exon is 3' UTR only 
                    elif (self.exonStarts[i] > self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.utr3.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.exonStarts[i]
                    else: 
                        print "Thar be dragons - Transcript computeMetadata + stranded gene region parsing" 


            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
                self.exons.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                self.exonsLen += self.exonEnds[i] - self.exonStarts[i]
            
            #print "DBUG2: i %d self.exonCt-1 %d self.exonEnds %s self.exonStarts %s" % (i, self.exonCt-1, self.exonEnds, self.exonStarts)
        
                if (i < self.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                    self.introns.append((self.chrom, self.exonEnds[i], self.exonStarts[i+1], self.name))
                    self.intronsLen += self.exonStarts[i+1] - self.exonEnds[i] 


        elif (self.strand == "-"):
            ### Promoter
            self.promotStart = self.txEnd
            self.promotEnd = self.txEnd + 1000
            ### 1-5K
            self.p15kStart = self.txEnd + 1000
            self.p15kEnd = self.txEnd + 5000
            ### downstread 1k
            self.d1kStart = self.txStart - 1000 if  self.txStart - 1000 > 0 else 0
            self.d1kEnd = self.txStart - 1 if  self.txStart - 1 >0 else 0
     #uc001ach.2	    chr1    -	    910578  917473  911551  916546  5	    910578,911878,914260,916516,917444,	    911649,912004,916037,916553,917473,	    Q5SV97  uc001ach.2
            #	name		chrom	strand	txStart txEnd	cdsStart self.cdsEnd exonCt	exonStarts		exonEnds		proteinID  alignID 
            # for the minus strand everything is the same except the order of encountering regions is reversed
            # i.e. 3' UTR -> CDS -> 5' UTR 
            
            for i in range (self.exonCt): 
            #print ("DBUG - exonCt %d i %d self.exonEnds[i] %d self.cdsStart %d exonStarts[i] %d self.cdsEnd %d") % \
                #    (self.exonCt, i, self.exonEnds[i], self.cdsStart, self.exonStarts[i], self.cdsEnd)
                
                if (self.cdsStart != self.cdsEnd):
                    self.coding = True
                    #if self.symbol == "Mrpl15":
                    #    print self.cdsStart + "\t" + self.cdsEnd
                # -- first compute 5'utr, CDS, 3'utr regions --
                # -- this is the same as for + sense except 5' UTR and 3' UTR are swapped throughout
                #case 1 - exon spans 3' UTR/CDS/5' UTR
                    if (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] > self.cdsEnd):
                        self.utr3.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr3Len += self.cdsStart - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.cdsStart
                        self.utr5.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.cdsEnd
                        self.utr5start.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5startLen += self.exonEnds[i] - (self.cdsEnd)
                #case 2 - exon spans 3' UTR/CDS junction
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] >= self.cdsStart):
                        self.utr3.append((self.chrom, self.exonStarts[i], self.cdsStart, self.name))
                        self.utr3Len += self.cdsStart - self.exonStarts[i]
                        self.cds.append((self.chrom, self.cdsStart, self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i]- self.cdsStart
                #case 3 - exon spans CDS/5'UTR junction 
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonStarts[i] <= self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.cdsEnd, self.name))
                        self.cdsLen += self.cdsEnd - self.exonStarts[i]
                        self.utr5.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.cdsEnd
                        self.utr5start.append((self.chrom, self.cdsEnd, self.exonEnds[i], self.name))
                        self.utr5startLen += self.exonEnds[i] - (self.cdsEnd)
                #case 4 - exon is 3' UTR only 
                    elif (self.exonStarts[i] < self.cdsStart and self.exonEnds[i] < self.cdsStart): 
                        self.utr3.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr3Len += self.exonEnds[i] - self.exonStarts[i]
                #case 5 - exon is CDS only
                    elif (self.exonStarts[i] >= self.cdsStart and self.exonEnds[i] <= self.cdsEnd):
                        self.cds.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.cdsLen += self.exonEnds[i] - self.exonStarts[i]
                #case 6 - exon is 5' UTR only 
                    elif (self.exonStarts[i] > self.cdsEnd and self.exonEnds[i] > self.cdsEnd):
                        self.utr5.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                        self.utr5Len += self.exonEnds[i] - self.exonStarts[i]
                        self.utr5start.append((self.chrom, self.exonStarts[i] , self.exonEnds[i], self.name))
                        self.utr5startLen += self.exonEnds[i] - self.exonStarts[i]
                    else: 
                        print "Thar be dragons - Transcript computeMetadata - stranded gene region parsing" 
                    
            #else: 
            #    print "- strand noncoding transcript"
                

            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
                self.exons.append((self.chrom, self.exonStarts[i], self.exonEnds[i], self.name))
                self.exonsLen += self.exonEnds[i] - self.exonStarts[i]
            
                if (i < self.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                    self.introns.append((self.chrom, self.exonEnds[i], self.exonStarts[i+1], self.name))
                    self.intronsLen += self.exonStarts[i+1] - self.exonEnds[i] 
                
        else:
            print "Thar be dragons - Transcript computeMetadata strand does not match + or -"
        

# input to createGTFTranscript below must be a list of dictionaries for each line of the input GTF file 
# these are created inside knowngenes_to_transcript_regions.py 

# example input: 

#[{'gene_name': 'DDX11L1', 'seqname': '1', 'end': '12227', 'start': '11869', 'frame': None, 'transcript_source': 'havana', 'feature': 'exon', 'exon_number': '1', 'exon_id': 'ENSE00002234944', 'tss_id': 'TSS15145', 'source': 'processed_transcript', 'gene_source': 'ensembl_havana', 'score': None, 'gene_biotype': 'pseudogene', 'gene_id': 'ENSG00000223972', 'transcript_id': 'ENST00000456328', 'transcript_name': 'DDX11L1-002', 'strand': '+'}, {'seqname': '1', 'end': '14409', 'start': '11869', 'frame': None, 'transcript_source': 'havana', 'feature': 'transcript', 'gene_id': 'ENSG00000223972', 'tss_id': 'TSS15145', 'source': 'processed_transcript', 'gene_source': 'ensembl_havana', 'score': None, 'gene_biotype': 'pseudogene', 'gene_name': 'DDX11L1', 'transcript_id': 'ENST00000456328', 'transcript_name': 'DDX11L1-002', 'strand': '+'}]

# keys for each dict:
#  gene_name
#  seqname
#  start
#  end
#  frame
#  transcript_source
#  feature
#  exon_number
#  exon_id
#  tss_id
#  source
#  gene_source
#  score
#  gene_biotype
#  gene_id
#  transcript_id
#  transcript_name
#  strand

def createGTFTranscript(gtfLines):
    foo = Transcript()

    
    # these properties (better be) all identical for each entry in the list of dicts 
    
    first = gtfLines[0] 

    foo.name = first["transcript_id"]
    foo.symbol = first["gene_name"] if "gene_name" in first else first["transcript_id"]
    foo.chrom = first["seqname"]
    foo.strand = first["strand"]

    # now process all lines for this transcript ID 

    for dict in gtfLines: 
        
        # ensembl GTFs have special lines where feature = "transcript" and feature = "CDS" that define the transcript and CDS start/ends, respectively 
    
        # GTF files are closed intervals while BED are right-open-left-closed, so --- 
        #   need to subtract one from all start coordinates? seems counterintuitive maybe the input genome.fa is zero based? 

        if (dict["feature"] == "exon"):
            ### Get the start and end position of transcript
            if (foo.txStart == 0 or int(dict["start"]) < foo.txStart):
                foo.txStart = int(dict["start"]) - 1
            if (foo.txEnd == 0 or int(dict["end"]) > foo.txEnd):
                foo.txEnd = int(dict["end"])
            
            ### 
            foo.exonCt += 1 
            foo.exonStarts.append(int(dict["start"]) - 1)
            foo.exonEnds.append(int(dict["end"]))


        if (dict["feature"] == "CDS"):
            #print dict
            if (foo.cdsStart == 0 or int(dict["start"]) < foo.cdsStart):
                foo.cdsStart = int(dict["start"]) - 1
            if (foo.cdsEnd== 0 or int(dict["end"]) > foo.cdsEnd):
                foo.cdsEnd = int(dict["end"])
            #if foo.symbol == "Mrpl15":
            #    print str(foo.cdsStart) +"\t"+str(foo.cdsEnd) 
            
        ##if (dict["feature"] == "exon"):  comment out by xie186
        #    foo.exonCt += 1
        #
        #    foo.exonStarts.append(int(dict["start"]) - 1)
        #    foo.exonEnds.append(int(dict["end"]))

    foo.exonStarts = sorted(foo.exonStarts)
    foo.exonEnds = sorted(foo.exonEnds) 

    foo.computeMetadata() 
    return foo 

