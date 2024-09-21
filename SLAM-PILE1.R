library('Rsamtools')

## ScanBamParam constructor
sbp <- ScanBamParam()

## PileupParam constructor

phred2ASCIIOffset(10, "Illumina 1.8+")

p_param <- PileupParam(max_depth=250, 
                       min_base_quality=13, 
                       min_mapq=0, 
                       min_nucleotide_depth=1, 
                       min_minor_allele_depth=0, 
                       distinguish_strands=TRUE, 
                       distinguish_nucleotides=TRUE, 
                       ignore_query_Ns=TRUE, 
                       include_deletions=TRUE, 
                       include_insertions=FALSE, 
                       left_bins=NULL, 
                       query_bins=NULL)

## Designate destination folder
destination = '/folder_with_bam_and_bai_files/'
## Designate name of .bam and .bai files.
file_name = 'file_name'

bamFile = paste(destination, file_name, '.bam', sep='')
baiFile = paste(bamFile, '.bai', sep='')
pileFile = paste(bamFile, '.pileup.csv', sep='')
print(bamFile)
print(baiFile)
print(pileFile)
## run and save fileup
res <- pileup(bamFile, index=baiFile, scanBamParam=sbp, pileupParam=p_param)
write.csv(res, file = pileFile)
print('Done.')