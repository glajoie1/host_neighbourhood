# Analyzing sequence data for
# "Host neighbourhood shapes bacterial community assembly and specialization on tree species across a latitudinal gradient"
# Geneviève Lajoie 2019-2020

# For the first part: processing with Dada2 (adjusted from Ben Callahan's tutorial)

# Load librarires

library(dada2)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(msa)
library(phangorn)
library(ips)


# Check if all files can be read
path<-'C:/Genevieve/SeqData2018/idemp-demultiplexed'
list.files(path)

# Match lists of forward and reverse files
fnFs <- sort(list.files(path, pattern="_R2_trim.fastq", full.names=TRUE)) # F is R2
fnRs <- sort(list.files(path, pattern="_R1_trim.fastq", full.names=TRUE)) # R is R1

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

####################
# Quality profiles #
####################

# Sample
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

###################
# Filter and trim #
###################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,170),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs=TRUE) # On Windows set multithread=FALSE

head(out2)
# Only 1 sample (N3) did not pass the filter at all.
# remove N3 from filtFs
# remove N3 from filtRs

filtFs<-filtFs[-which(filtFs=='C:/Genevieve/SeqData2018/idemp-demultiplexed/filtered/N3_F_filt.fastq.gz')]
filtRs<-filtRs[-which(filtRs=='C:/Genevieve/SeqData2018/idemp-demultiplexed/filtered/N3_R_filt.fastq.gz')]

# saved as dada2env2.RData
# If need be:
# load("~/Genevieve/SeqData2018/idemp-demultiplexed/dada2env2.RData")

#########################
# Learn the error rates #
#########################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize error rates
plotErrors(errF, nominalQ=TRUE)

##################
# De-replication #
##################

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[-which(sample.names=='N3')]
names(derepRs) <- sample.names[-which(sample.names=='N3')]

####################
# Sample inference #
####################

dadaFs <- dada(derepFs, err=errF, multithread=T, pool="pseudo")
dadaRs <- dada(derepRs, err=errR, multithread=T, pool="pseudo")

# Check up results
dadaFs[[1]]
dadaRs[[1]]

######################
# Merge paired reads #
######################

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE) # From the output, most sequences seem to have merged

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

############################
# Construct sequence table #
############################

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) # From 210 to 368

# Remove non-target sequences
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(293,322)]

###################
# Remove chimeras #
###################

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Proportion
sum(seqtab.nochim)/sum(seqtab)

###############
# Track reads #
###############

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

###################
# Assign taxonomy #
###################

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("C:/Genevieve/SeqData2018/idemp-demultiplexed/tax/SILVA_SSU_r132_March2018.RData")
# this dataset downloaded from http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData
ids <- IdTaxa(dna, trainingSet, strand="both", processors=1, verbose=FALSE) # Processors=NULL crashes R - Only works with processors=1 / Prend entre 8-12h ? runner
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks 
rownames(taxid) <- getSequences(seqtab.nochim)

taxa<-taxid

# Examining output
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


############
# Phyloseq #
############

theme_set(theme_bw())

# Import metadata
metadata<-read.csv('C:/Genevieve/SeqData2018/metadata_plotspp.csv', sep=',')
rownames(metadata)<-metadata$Sample_ID

# Construct phyloseq object
ps<-phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(metadata), tax_table(taxa))

# Save object as dadaphyl.RData
# If need be:
# load('~/Genevieve/SeqData2018/idemp-demultiplexed/dadaphyl.RData')

#####################
# Phylogenetic tree #
#####################

# Prepare data for alignment
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
seq.fa<-DNAStringSet(seqs)
writeXStringSet(seq.fa, '/mnt/storage/genevieve/seq.fa') # Need to transfer the file to phyllo server where Qiime is installed

# Perform alignment, masking and tree-building within the console
# Script is in file /mnt/storage/genevieve/Field2017/postdada2.txt

# For checking: import alignment
    TS3alignment <- read.dna(file="/mnt/storage/genevieve/seq_aligned.fasta",format="fasta",as.matrix=TRUE)
    # Run diagnostics
    checkAlignment(TS3alignment)
    
# Lost 2 seqs

# Import phylogenetic tree
tree <- read_tree("/mnt/storage/genevieve/FastTree2/FT_tree")

# Add to the phyloseq object
ps <- merge_phyloseq(ps, tree)

#############
# ASV names #
#############

# Rename ASVs while conserving sequences
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

############################
# Load ecological datasets #
############################

##########
# Rarefy #
##########

# Examine distribution of number of sequences
max(sample_sums(ps)) # Min:15 # Max: 72419
sort(sample_sums(ps))

# Rarefaction curve
# rarecurve(otu_table(ps.epi), step=100, cex=0.5, label=F)

# Rarefy (in total, jumping from 323 samples to 298 samples)
ps.epi7K<-rarefy_even_depth(ps, sample.size=7000, replace=F, rngseed=711)

# Remove the samples from the plots where only one species remains after the rarefaction
# VER-HI-P3, MSH-HI-P2: AS20 et B31
ps.epi7K <- subset_samples(ps.epi7K, Sample_ID != c('AS20'))
ps.epi7K <- subset_samples(ps.epi7K, Sample_ID != c('B31'))

# Remove these plots from other objects
site.UTM<-site.UTM[-which(site.UTM$Site.name%in%c('VER-HI-P3','MSH-HI-P2')),]
comp<-comp[-which(rownames(comp)%in%c('VER-HI-P3','MSH-HI-P2')),]
metadata$Site.name<-paste(metadata$Site,metadata$Type,metadata$Plot, sep='-')
metadata<-metadata[-which(metadata$Site.name%in%c('VER-HI-P3','MSH-HI-P2')),]



##########
# Export #
##########

# Can now export ps object or save within the RData workspace
save.image(file="/mnt/storage/genevieve/phyloseq_obj.RData")

