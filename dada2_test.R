### Test du package DADA2
### Date de création : 30.01.2020
### 
### Objectif : tester le package dada2 avec le jeu de données fourni par les développeurs puis avec les jeux de données de metabarcode
### Packages requis : devtools, dada2

rm(list - ls())

##Install packages
if (require(devtools))
{
    install.packages("devtools")
}

devtools::install_github("benjjneb/dada2", ref="v1.14") # change the ref argument to get other versions

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ShortRead")

library(ShortRead)

library(dada2)

                                                    ##############################################################
                                                    ##                                                          ##
                                                    ##                       Dataset test                       ##
                                                    ##                                                          ##
                                                    ##############################################################
path <- "C:/Users/quent/Desktop/essai DADA2/2e essai/echantillon 2017 meta37/raw(trim sur gapA)/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
file_names_Fwd_seq <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
file_names_Rev_seq <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- sapply(strsplit(basename(file_names_Fwd_seq), "_"), `[`, 1)

trimRs <- file.path(path, "trimmed", paste0(sample_names, "_R_trim.fastq.gz"))

file_to_trim <- removePrimers(fn = file_names_Rev_seq[i], fout = trimRs, primer.fwd = "TCRTACCARGAAACCAGTT", max.mismatch = 5)

# Check read quality along the sequence
x11() ; plotQualityProfile(file_names_Fwd_seq)

x11() ; plotQualityProfile(trimRs)


####################################################
# Filter sequences
####################################################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "trimmed", "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "trimmed", "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

#Filtering
out <- filterAndTrim(file_names_Fwd_seq, filtFs, trimRs, filtRs, truncLen=c(280,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


#test <- readFastq(paste(path, "/trimmed", sep=""), "_R_trim.fastq.gz")
#sread(test)

####################################################
# Learn errors
####################################################

errF <- learnErrors(filtFs, multithread=TRUE)
x11() ; plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)
x11() ; plotErrors(errR, nominalQ=TRUE)


####################################################
# Sample Inference
####################################################

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


####################################################
# Merge paired reads
####################################################

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers)


####################################################
# Construct sequence table
####################################################

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

####################################################
# Remove chimera
####################################################

seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_nochim)

sum(seqtab_nochim)/sum(seqtab)


####################################################
# Track reads through the pipeline
####################################################

getN <- function(x) sum(getUniques(x))

track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)

ratio <- function(x,y) x/y

track_prctg <- cbind(ratio(out[,2], out[,1]), ratio(track[,3], out[,1]), ratio(track[,4], out[,1]), ratio(track[,5], out[,1]), ratio(track[,6], out[,1]))
colnames(track_prctg) <- c("filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_prctg) <- sample_names
head(track_prctg)


####################################################
# Assign taxonomy
####################################################

taxa <- assignTaxonomy(seqtab_nochim, "C:/Users/quent/Desktop/essai DADA2/2e essai/gapA2.fa.gz", taxLevels = c("Genus", "species", "strain"), minBoot = 80, outputBootstraps = TRUE,multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
row.names(taxa.print) <- NULL
head(taxa.print)


####################################################
# Evaluate accuracy
####################################################

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

####################################################
# Handoff to phyloseq
####################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq"); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(seqtab_nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

install.packages("seqinr")
require(seqinr)

id = character(dim(seqtab_nochim)[2])
for (i in 1:length(id)) {
    id[i] <- paste(taxa[[1]][i,1], "(", taxa[[2]][i,1], ") ", taxa[[1]][i,2], "(", taxa[[2]][i,2], ") ",
                   taxa[[1]][i,3], "(", taxa[[2]][i,3], ")","-", seqtab_nochim[i], sep = "")
}

for (i in 1:length(sample_names)) {
    write.fasta(as.list(colnames(seqtab_nochim)[which(seqtab_nochim != 0)]), names = id, 
                file.out = paste(path, "/", sample_names[i], "-ASV-sequences.fa" , sep = ""), open = "w")
}
