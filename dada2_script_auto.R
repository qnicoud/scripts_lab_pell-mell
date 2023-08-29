### Applpication of DADA2 package to mutliple datasets
### Created  09.12.2020
### 
### Objective: automatized samples analysis
### Required packages: devtools, dada2, seqinr

# options(encoding = "UTF-8")

# rm(list = ls())
# 
# ##Install packages
# if (!require(devtools))
# {
#   install.packages("devtools")
# }
# z
# devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ShortRead")
# 
# # a lancer a chaque debut de session
# #library(ShortRead)
# 
# library(dada2)

##############################################################
##                                                          ##
##                       Dataset test                       ##
##                                                          ##
##############################################################

# path <- "C:/Users/jeremy.cigna/Desktop/metabarcod gapA/2020/metabarcod/19_67/seq_sept2020"

get_files_list <- function(path) {
    path <- gsub("\\\\", "/", path)
    
    # Forward and reverse fastq filenames have format: SAMPLENAME_trim_map_r1.fastq and SAMPLENAME_trim_map_r2.fastq
    file_names <- sort(list.files(path = path, pattern=".fastq", full.names = TRUE))
    
    
    # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq simplifie 
    sample_names <- sapply(strsplit(basename(file_names), "_"), `[`, 1)
    
    directionnality <- rep("F", times = length(file_names))
    directionnality[grep("r2", file_names)] <- "R"
    
    files_table <- data.frame(sampleNames = sample_names, fileNames = file_names, dir = directionnality)
    
    return(files_table)
}

files_transfo <- function(files) {
    new_files <- list()
    
    for (i in 1:dim(files)[1]) {
        new_files[[i]] <- files[i,]
    }
    return(new_files)
}

plotQuality <- function(path_to_file) {
    dir.create("quality_check", showWarnings = FALSE)

    pdf(file = paste("./quality_check/", path_to_file[1], "_", path_to_file[3], ".pdf", sep = ""))
    plotQualityProfile(path_to_file)
    dev.off()
}


####################################################
# Filter sequences
####################################################

filtering_set_up <- function(files) {
    
    # Place filtered files in filtered/ subdirectory
    filt <- file.path("./filtered", paste0(files$sampleNames, "_", files$dir, ".fastq.gz"))
    names(filt) <- files$sampleNames
    
    # Generate a frame that contains location of different files with filtered and raw sequences
    filtering_list <- data.frame(filesF = files$fileNames[which(files$dir %in% "F")],
                                 filtF = filt[which(files$dir %in% "F")],
                                 filesR = files$fileNames[which(files$dir %in% "R")],
                                 filtR = filt[which(files$dir %in% "R")])
    
    return(filtering_list)
}


filtering <- function(filtering_list, trunc_location, nMax = 0, eeMax = c(2,2), qTrunc = 2, compress = TRUE) {

    trunc <- openxlsx::read.xlsx(trunc_location)

    # Filtering
    out <- list()
    
    if (sum(class(eeMax) == "numeric") && dim(filtering_list)[1] != 1)
        eeMax <- matrix(data = rep(eeMax, times = dim(filtering_list)[1]), ncol = 2, byrow = TRUE, dimnames(list(names(filtering_list), c("F", "R"))))
    
    for ( i in 1:dim(filtering_list)[1]) {
        out[i] <- filterAndTrim(filtering_list$filesF[i], filtering_list$filtF[i], 
                                filtering_list$filesR[i], filtering_list$filtF[i], 
                                truncLen = trunc[i,],
                                maxN = nMax, maxEE = eeMax, truncQ = qTrunc, rm.phix = TRUE,
                                compress = compress, multithread = FALSE) # On Windows set multithread=FALSE
    }
    
    return(out)
} 

 

####################################################
# Learn errors
####################################################

findErrors <- function(filt) {
    dir.create("error_rates", showWarnings = FALSE)
    
    errors <- list()
    
    for ( i in c("F", "R")) {
        if ( i == "F")
            ind <- 2
        else
            ind <- 4
        
        err <- learnErrors(filt[ind], multithread=TRUE)
        
        pdf(file = paste("./error_rates/", i, ".pdf", sep = ""))
        plotErrors(err, nominalQ=TRUE)
        dev.off()
        
        errors[i] <- err
    }
    
    return(errors)
}

####################################################
# Sample Inference
####################################################

getN <- function(x) sum(getUniques(x))

ratio <- function(x,y) x/y

goDada <- function(filt, errors, chim_method = "consensus") {
    
    # Sample Inference
    dadaFs <- dada(filt$filtF, err = errors$"F", multithread = TRUE)
    dadaRs <- dada(filt$filtR, err = errors$"R", multithread = TRUE)
    
    # Merge paired reads
    mergers <- mergePairs(dadaFs, filt$filtF, dadaRs, filt$filtR, verbose = TRUE)
    
    # Construct sequence table
    seqtab <- makeSequenceTable(mergers)
    
    # Remove chimera
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = chim_method, multithread = TRUE, verbose = TRUE)
    
    # Track reads through the pipeline
    track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab_nochim))
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample_names
    
    track_prctg <- cbind(ratio(out[,2], out[,1]), ratio(track[,3], out[,1]), ratio(track[,4], out[,1]), ratio(track[,5], out[,1]), ratio(track[,6], out[,1]))
    colnames(track_prctg) <- c("filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track_prctg) <- sample_names
    
    res <- list(dadaFs, dadaRs, mergers, seqtab, seqtab_nochim, track, track_prctg)
    
    return(res)
}



####################################################
# Assign taxonomy
####################################################

taxa <- assignTaxonomy(seqtab_nochim, "C:/Users/jeremy.cigna/Desktop/metabarcod gapA/2017/essai DADA2/base de données gapA- v2.fa", taxLevels = c("Genus", "species", "strain"), minBoot = 80, outputBootstraps = TRUE,multithread=TRUE)







# ####################################################
# # Evaluate accuracy
# ####################################################
#
# unqs_mock <- seqtab_nochim["Mock",]
# unqs_mock <- sort(unqs_mock[unqs_mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
# cat("DADA2 inferred", length(unqs_mock), "sample sequences present in the Mock community.\n")
#
# mock_ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
# match_ref <- sum(sapply(names(unqs_mock), function(x) any(grepl(x, mock_ref))))
# cat("Of those,", sum(match_ref), "were exact matches to the expected reference sequences.\n")

####################################################
# Handoff to phyloseq
####################################################
# 
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

wd <- "C:/Users/quent/Desktop/dada2/"

setwd(wd)

files <- get_files_list("./")

new_files <- files_transfo(files)

lapply(new_files, plotQuality)

filtering_table <- filtering_set_up(files)

filt <- filtering(filtering_table, "truncation.xlsx", )

errors <- findErrors(filt)

output <- goDada(filt, errors, chim_method = "consensus")