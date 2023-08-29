### Applpication of DADA2 package to mutliple datasets
### Created  09.12.2020
### 
### Objective: automatized samples analysis
### Required packages: devtools, dada2, seqinr

# options(encoding = "UTF-8")

rm(list = ls())

##Install packages
if (!require(devtools))
{
  install.packages("devtools")
}

#Bugué
# devtools::install_github("benjjneb/dada2", ref="v1.18") # change the ref argument to get other versions

# Fonctionne mais nécéssite de télécharger manuellement le package
install.packages("C:/Users/jeremy.cigna/Desktop/dada2-1.18/dada2-1.18",
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")

# a lancer a chaque debut de session
#library(ShortRead)

if (!require(openxlsx))
{
  install.packages("openxlsx")
}

if (!require(seqinr))
{
  install.packages("seqinr")
}

library(openxlsx)

library(seqinr)

library(dada2)

# Ne pas oublier de charger les fonctions en les exécutant (ne donne pas de résultats, enregistre seulement les fonctions dans l'environment)

##############################################################
##                                                          ##
##                       Dataset test                       ##
##                                                          ##
##############################################################

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
  print(plotQualityProfile(path_to_file[2]))
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
  
  truncation_table <- matrix(data = 0, ncol = 2, nrow = dim(filtering_list)[1], dimnames = list(row.names(filtering_list), c("F", "R")))
  openxlsx::write.xlsx(truncation_table, file = "truncation.xlsx", row.names = TRUE)
  
  return(filtering_list)
}


filtering <- function(filtering_list, trunc_location, nMax = 0, eeMax = c(2,2), qTrunc = 2, compress = TRUE) {
  
  trunc <- openxlsx::read.xlsx(trunc_location, rowNames = TRUE)
  
  # Filtering
  out <- list()
  
  for ( i in 1:dim(filtering_list)[1]) {
    out[[i]] <- filterAndTrim(filtering_list$filesF[i], filtering_list$filtF[i], 
                            filtering_list$filesR[i], filtering_list$filtR[i], 
                            truncLen = trunc[i,],
                            maxN = nMax, maxEE = eeMax, truncQ = qTrunc, rm.phix = TRUE,
                            compress = compress, multithread = FALSE) # On Windows set multithread=FALSE
  }
  
  return(out)
} 



####################################################
# Learn errors
####################################################

findErrors <- function() {
  dir.create("error_rates", showWarnings = FALSE)
  
  errors <- list()
  
  for ( i in c("F", "R")) {
    # if ( i == "F")
    #   ind <- 2
    # else
    #   ind <- 4
    filtered <- list.files("./filtered", pattern = paste(i, ".fastq.gz", sep = ""), full.names = TRUE)
    sample.names <- sapply(strsplit(basename(filtered), "_"), `[`, 1)
    names(filtered) <- sample.names
    
    err <- learnErrors(filtered, multithread=TRUE)
    
    pdf(file = paste("./error_rates/", i, ".pdf", sep = ""))
    print(plotErrors(err, nominalQ=TRUE))
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

goDada <- function(out, errors, chim_method = "consensus") {
  
  files <- list()
  
  for ( i in c("F", "R")) {
    filt_files <- list.files("./filtered", pattern = paste(i, ".fastq.gz", sep = ""), full.names = TRUE)
    sample.names <- sapply(strsplit(basename(filt_files), "_"), `[`, 1)
    names(filt_files) <- sample.names
    
    files[[i]] <- filt_files
    
    # Sample Inference 
        assign(paste("dada", i, sep = ""), dada(filt_files, err = errors[[i]], multithread = TRUE))
  }
  
  # Merge paired reads
  mergers <- mergePairs(dadaF, files[[1]], dadaR, files[[2]], verbose = TRUE)
  
  # Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  
  # Remove chimera
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = chim_method, multithread = TRUE, verbose = TRUE)
  
  # Track reads through the pipeline
  out <- do.call("rbind", out)
  track <- cbind(out, lapply(dadaF, getN), lapply(dadaR, getN), lapply(mergers, getN), rowSums(seqtab_nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names

  track_prctg <- cbind(ratio(out[,2], out[,1]), ratio(do.call("rbind", track[,3]), out[,1]), ratio(do.call("rbind", track[,4]), out[,1]), ratio(do.call("rbind", track[,5]), out[,1]), ratio(do.call("rbind", track[,6]), out[,1]))
  colnames(track_prctg) <- c("filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track_prctg) <- sample.names

  res <- list(ddFwd = dadaF, ddRev = dadaR, merg = mergers, seqTab = seqtab, chimRm = seqtab_nochim, trackRd = track, trackRdPerc = track_prctg)
  # Add write report (xlsx)
  
  return(res)
}

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

wd <- "C:/Users/jeremy.cigna/Desktop/metabarcod gapA/2020/metabarcod/all_sample_sept"

setwd(wd)

files <- get_files_list("./")

new_files <- files_transfo(files)

lapply(new_files, plotQuality)

filtering_table <- filtering_set_up(files)

filt <- filtering(filtering_table, "truncation.xlsx")

errors <- findErrors()

output <- goDada(filt, errors, chim_method = "consensus")

# Enregistrer les données de l'environnement dans un fichier de données:
# save.image(file = "all_sept.RData")
# Charger le fichier de données:
# load(file = "up_to_output.RData")

taxa <- assignTaxonomy(output$chimRm, "C:/Users/jeremy.cigna/Desktop/metabarcod gapA/2017/essai DADA2/base de données gapA- v3.fa", taxLevels = c("Genus", "species", "strain"), minBoot = 80, outputBootstraps = TRUE,multithread=TRUE)

# Sauvegarder le fichier fasta final
for ( i in 1:dim(output$chimRm)[1] ) {
  cond_name <- dimnames(output$chimRm)[[1]][i]
  sub_out <- output$chimRm[i,]
  sub_out <- sub_out[which(sub_out != 0)]
  sub_out <- sort(sub_out, decreasing = TRUE)
  
  id <- character(length(sub_out))
  for (j in 1:length(sub_out)) {
    asv <- names(sub_out)[j]
    id[j] <- paste(taxa[[1]][asv,1], "(", taxa[[2]][asv,1], ") ", taxa[[1]][asv,2], "(", taxa[[2]][asv,2], ") ",
                   taxa[[1]][asv,3], "(", taxa[[2]][asv,3], ")","-", sub_out[asv], sep = "")
  }
  
  dir.create("./output", showWarnings = FALSE)
  write.fasta(as.list(names(sub_out)), names = id, 
                file.out = paste("./output/", cond_name, "-ASV-sequences.fa" , sep = ""),
                open = "w")
}



############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
## A TESTER

####################################################
# Handoff to phyloseq
####################################################
# 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq"); packageVersion("phyloseq")
require(phyloseq)
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(output$chimRm)

ps <- phyloseq(otu_table(output$chimRm, taxa_are_rows=FALSE),
               tax_table(taxa$tax))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, measures=c("Shannon", "Simpson"))

# # Transform data to proportions as appropriate for Bray-Curtis distances
# ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
# ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
# 
# plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Genus")
plot_bar(ps.top20, fill="species")
