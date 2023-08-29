rm(list = ls())

load("ortholog_objects")

#Load Packages
{
    if (! require(openxlsx)) {
        install.packages("openxlsx")
        library(openxlsx)
    }
    
    if (! require(xlsx)) {
        install.packages("xlsx")
        library(xlsx)
    }
    
    #Package required for removing duplicates in tables.
    if( ! require(dplyr)) {
        install.packages("dplyr")
        library(dplyr)
    }
    
    if( ! require(ggplot2)) {
        install.packages("ggplot2")
        library(ggplot2)
    }
    
    if( ! require(ggpmisc)) {
        install.packages("ggpmisc")
        library(ggpmisc)
    }
    
    if( ! require(ggrepel)) {
        install.packages("ggrepel")
        library(ggrepel)
    }
    
    # Package required to format tables in the long-table format
    if( ! require(tidyr)) {
        install.packages("tidyr")
        library(tidyr)
    }
     
    #Packages required to make multipannel figures but cowplot won't be installed...
    if( ! require(gridExtra)) {
        install.packages("gridExtra")
        library(gridExtra)
    }
    
    # if( ! require(cowplot)) {
    #     devtools::install_github("wilkelab/cowplot")
    #     library(cowplot)
    # }
    
    # library
    if( ! require(pheatmap)) {
        install.packages("pheatmap")
        library(pheatmap)
    }
    
    if( ! require(corrplot)) {
        install.packages("corrplot")
        library(corrplot)
    }
    
    if( ! require(ggplot2)) {
        install.packages("ggplot2")
        library(ggplot2)
    }
    
    if( ! require(kohonen)) {
        install.packages("kohonen")
        library(kohonen)
    }
    
    if( ! require(reshape2)) {
        install.packages("reshape2")
        library(reshape2)
    }
    
    if( ! require(NbClust)) {
        install.packages("NbClust")
        library(NbClust)
    }
    
    if( ! require(edgeR)) {
        install.packages("edgeR")
        library(edgeR)
    }
}

# At the end of the session
save(list = ls(), file = "ortholog_objects")

#####
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |              Get data                |      /
#                                               \      |                                      |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\
#####

## Gather transcriptomic data
ors <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/summary_ORS_USDA_transcriptome.xlsx", sheet = 1)
usda <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/summary_ORS_USDA_transcriptome.xlsx", sheet = 2)

## Get ortholog list (phyloprofiles)
orth <- read.csv("C:/Users/quent/Desktop/compa_phyloprofile0050_CIDsansFiltre.csv", header = TRUE, sep =";", stringsAsFactors = FALSE)
orth <- orth[orth$X == 3 , c(1,6)]

corres <- openxlsx::read.xlsx("C:/Users/quent/Desktop/ORS285_v1-v2.xlsx")

orth <- orth[which(orth$ï..Label %in% corres$Previous.label),]
a <- corres[which(corres$Previous.label %in% orth$ï..Label),]

for (i in 1:dim(orth)[1]) {
    if (sum(orth$ï..Label[i] %in% a$Previous.label)) {
        orth$ï..Label[i] <- a$Current.label[which(a$Previous.label == orth$ï..Label[i])]
    } else {
        orth$ï..Label[i] <- ""
    }
}

#####
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |             Format data              |      /
#                                               \      |                                      |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\
#####


# Clean the ortholog table by removing lines with "No Hit" or empty cells
orth <- orth[orth$Bradyrhizobium.japonicum.USDA.110 != "",]
orth <- orth[orth$Bradyrhizobium.japonicum.USDA.110 != "No Hit",]

# If a list of genes is given as orthologs, remove all but the first (with higher identity and e-value)
orth_cln <- data.frame(Label = orth$ï..Label, Bjap = gsub("\\,.*","",orth$Bradyrhizobium.japonicum.USDA.110))
colnames(orth_cln) <- c("Label", "Bjap")

# Test, restrain transcriptomic tables only to the genes that are in the ortholog table
ors_orth <- ors[which(ors$Gene.accession %in% orth_cln$Label), c(1,16)]
usda_orth <- usda[which(usda$Label %in% orth_cln$Bjap), c(1,12)]

# Problem, tables don't have the same dimensions
dim(ors_orth)
dim(usda_orth)

# Test, shorten the ortholog table so only the genes that are only in each transcriptomic table remain 
orth_shrt <- orth_cln[which(orth_cln$Bjap %in% usda_orth$Label),]
orth_shrt <- orth_shrt[which(orth_shrt$Label %in% ors_orth$Gene.accession),]

# Test, restrain transcriptomic tables only to the genes that are in the short ortholog table
ors_orth_shrt <- ors[which(ors$Gene.accession %in% orth_shrt$Label),]
usda_orth_shrt <- usda[which(usda$Label %in% orth_shrt$Bjap),]

# THis time it is still not fine --> presence of duplicates
dim(ors_orth_shrt)
dim(usda_orth_shrt)


# Removed not-unique genes from ortholog list
orth_shrt_unq <- (orth_cln %>% distinct(Bjap, .keep_all = TRUE))
dim(orth_shrt_unq)

# Applied this to the data
ors_orth_shrt_unq <- ors[which(ors$Gene.accession %in% orth_shrt_unq$Label),]
usda_orth_shrt_unq <- usda[which(usda$Label %in% orth_shrt_unq$Bjap),]

# Still not good
dim(ors_orth_shrt_unq)
dim(usda_orth_shrt_unq)

# Seems that there are differences in gene contents in each dataset.
ors_spec <- orth_shrt_unq$Label[which(usda_orth_shrt_unq$Label %in% orth_shrt_unq$Bjap)]
usda_spec <- orth_shrt_unq$Bjap[which(ors_orth_shrt_unq$Gene.accession %in% orth_shrt_unq$Label)]

#Just to have an idea
length(ors_spec)
length(usda_spec)

#Once I got the list I remove the genes lacking in the datasets from ortholog list
orth_spec <- orth_shrt_unq[which(orth_shrt_unq$Label %in% ors_spec) ,]
orth_spec <- orth_spec[which(orth_spec$Bjap %in% usda_spec) ,]
dim(orth_spec)

#And also get them out of data 
ors_fin <- ors_orth_shrt_unq[which(ors_orth_shrt_unq$Gene.accession %in% ors_spec),]
usda_fin <- usda_orth_shrt_unq[which(usda_orth_shrt_unq$Label %in% usda_spec),]

usda_fin <- usda_fin[!is.na(usda_fin$AA_vs_YM_LFC),]

orth_spec <- orth_spec[which(orth_spec$Bjap %in% usda_fin$Label),]
orth_spec <- orth_spec[which(orth_spec$Label %in% ors_fin$Gene.accession),]

ors_fin <- ors_fin[which(ors_fin$Gene.accession %in% orth_spec$Label),]
usda_fin <- usda_fin[which(usda_fin$Label %in% orth_spec$Bjap),]

#Should be OK
dim(ors_fin)
dim(usda_fin)
dim(orth_spec)



#=======================================#
#  Order the data
#=======================================#

# Perepare data for the final table (need to get them in the right order)
usda_id <- usda_fin$Label[order(usda_fin$Label)][rank(orth_spec$Bjap)]
usda_lfc <- usda_fin$AA_vs_YM_LFC[order(usda_fin$Label)][rank(orth_spec$Bjap)][!is.na(usda_id)]
usda_id <- usda_id[!is.na(usda_id)]

ors_id <- ors_fin$Gene.accession[order(ors_fin$Gene.accession)][rank(orth_spec$Label)]
ors_lfc <- ors_fin$AA_vs_YM.LFC[order(ors_fin$Gene.accession)][rank(orth_spec$Label)][!is.na(ors_id)]
ors_id <- ors_id[!is.na(ors_id)]


# Final data.frame
fin_fin <- data.frame(ors_id = ors_id, usda_id = usda_id, 
                      ors_lfc = ors_lfc, usda_lfc = usda_lfc)


#=======================================#
#  Filtering on FDR
#=======================================#

sub_fdr_ors <- ors$Gene.accession %in% fin_fin$ors_id
sub_fdr_usda <- usda$Label %in% fin_fin$usda_id

sub_fdr_ors_filtered <- ors$AA_vs_YM.FDR[which(sub_fdr_ors)]#[order(ors$Gene.accession[which(sub_fdr_ors)])][rank(fin_fin$ors_id)]
sub_fdr_usda_filtered <- usda$AA_vs_YM_FDR[which(sub_fdr_usda)]#[order(usda$Label[which(sub_fdr_usda)])][rank(fin_fin$usda_id)]

fdr_thresh <- 0.01
fdr_trim_ors <- ors$Gene.accession[which(sub_fdr_ors)][which(sub_fdr_ors_filtered < fdr_thresh)]
fdr_trim_usda <- usda$Label[which(sub_fdr_usda)][which(sub_fdr_usda_filtered < fdr_thresh)]

sub_ok <- as.logical(fin_fin$ors_id %in% fdr_trim_ors * fin_fin$usda_id %in% fdr_trim_usda)

fin_fdr_filt <- fin_fin[sub_ok,]



#####
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |            Visualize data            |      /
#                                               \      |                                      |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\
#####

disp <- ggplot(data = fin_fdr_filt, aes(x=ors_lfc, y=usda_lfc, color=a)) + 
    geom_point(size = 1) +
    scale_color_manual(values = c("black", "red", "orange", "#ff6f00", "grey", "#248f25", "#728772", "#a67777", "#f5d902")) + 
    stat_ellipse() +
    theme(legend.position = "bottom")
disp

inter <- 0.285*7


lfc_thresh <- 1.58

fin_fdr_lfc_filt <- fin_fdr_filt[which((fin_fdr_filt$ors_lfc > lfc_thresh | fin_fdr_filt$ors_lfc < -lfc_thresh) & (fin_fdr_filt$usda_lfc > lfc_thresh | fin_fdr_filt$usda_lfc < -lfc_thresh)),]
png(filename = "~/disp_fin2.png", width = 1036, height = 808)
disp <- ggplot(data = fin_fdr_lfc_filt, aes(x=ors_lfc, y=usda_lfc)) + 
    geom_smooth(data = fin_fdr_lfc_filt,
                method = "lm", formula = y ~ x, linetype="dashed",
                color="darkred", fill="blue") +
    geom_point(size = 2) +
    labs(x = "ORS285 LFC - AA vs YM", y = "USDA110 LFC - AA vs YM") + 
    stat_ellipse(data = fin_fdr_lfc_filt[intersect(which(fin_fdr_lfc_filt$ors_lfc > 7), which(fin_fdr_lfc_filt$usda_lfc > 7)),], type = "t", colour = "red") +
    theme(legend.position = "bottom", 
          axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
          plot.background = element_rect(fill = "#f5f5f5")) #+ 
    # geom_label_repel(data = fin_fdr_lfc_filt[intersect(which(fin_fdr_lfc_filt$ors_lfc > 7), which(fin_fdr_lfc_filt$usda_lfc > 7)),],
    #                  aes(label = ors_id, size = 0.5),
    #                  box.padding = unit(0.35, "lines"),
    #                  point.padding = unit(0.3, "lines"))
    # stat_poly_eq(data = fin_fdr_lfc_filt,
    #              formula = y ~ x,
    #              eq.with.lhs = "italic(hat(y))~`=`~",
    #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    #              parse = TRUE,
    #              label.x = 0.025,
    #              label.y = 0.995) +
    # geom_abline(slope = 0.6, intercept = inter, color = "darkred", size = 1) +
    # geom_abline(slope = 0.6, intercept = -inter, color = "darkred", size = 1) +

    
disp
dev.off()


### WARNNG: this uses the table generated at the very end of this script !
table <- cbind(table, a)
disp <- ggplot(data = table, aes(x=ors_lfc, y=usda_lfc, color=a)) + 
    geom_point(size = 1) +
    scale_color_manual(values = c("black", "red", "orange", "#ff6f00", "grey", "#248f25", "#728772", "#a67777", "#f5d902")) + 
    theme(legend.position = "None") + 
    stat_ellipse()
disp


# disp + geom_label_repel(data = subset(fin_fin, a == "u_u" | a == "u_l"),
#                         aes(label = ors_id, size = 0.5),
#                         box.padding = unit(0.35, "lines"),
#                         point.padding = unit(0.3, "lines"))



##===========================================##
##   Quick ggplot2 heatmap #### USELESS !!!  ##
##===========================================##

#Convert from short-table to long-table
fin_fin2 <- fin_fdr_filt[which(fin_fdr_filt$ors_lfc < -lfc_thresh),]
fin_fin2 <- rbind(fin_fin2, fin_fdr_filt[which(fin_fdr_filt$ors_lfc > lfc_thresh),])
fin_fin2 <- rbind(fin_fin2, fin_fdr_filt[which(fin_fdr_filt$usda_lfc < -lfc_thresh),])
fin_fin2 <- rbind(fin_fin2, fin_fdr_filt[which(fin_fdr_filt$usda_lfc > lfc_thresh),])

fin_long <- gather(data = fin_fin, key = org, value = lfc, -ors_id, -usda_id, -a, -gene_type)
fin_long2 <- gather(data = fin_fin2, key = org, value = lfc, -ors_id, -usda_id, -a, -gene_type)




heat_map_orthologs <- ggplot(data = fin_long2, aes(x=org, y=usda_id, fill=lfc)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("blue", "black", "yellow"), limits = c(-5, 5))
heat_map_orthologs


#####
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |      Get subsets for hm w/clust      |      /
#                                               \      |                                      |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\
#####


#============================#
#Full phyloprofile data
#============================#
hm_clust_name <- "fin_fin"
fin_hm_clust <- (fin_fin %>% distinct())
#------------------------------------------DONE

#============================#
#FDR filtering
#============================#
hm_clust_name <- "fin_fdr_filt"
fin_hm_clust <- (fin_fdr_filt %>% distinct())
#------------------------------------------DONE


## Required for followings
sub_fdr_ors <- ors$Gene.accession %in% fin_fin$ors_id
sub_fdr_usda <- usda$Label %in% fin_fin$usda_id
sub_fdr_ors_filtered <- ors$AA_vs_YM.FDR[which(sub_fdr_ors)]#[order(ors$Gene.accession[which(sub_fdr_ors)])][rank(fin_fin$ors_id)]
sub_fdr_usda_filtered <- usda$AA_vs_YM_FDR[which(sub_fdr_usda)]#[order(usda$Label[which(sub_fdr_usda)])][rank(fin_fin$usda_id)]
#===================================#
#FDR filtering ors not-DEG usda not-DEG
#===================================#
fdr_trim_ors <- ors$Gene.accession[which(sub_fdr_ors)][which(sub_fdr_ors_filtered > fdr_thresh)]
fdr_trim_usda <- usda$Label[which(sub_fdr_usda)][which(sub_fdr_usda_filtered > fdr_thresh)]

sub_ok <- as.logical(fin_fin$ors_id %in% fdr_trim_ors * fin_fin$usda_id %in% fdr_trim_usda)

fin_fdr_filt_not <- fin_fin[sub_ok,]

hm_clust_name <- "fin_fdr_filt_not"
fin_hm_clust <- (fin_fdr_filt_not %>% distinct())
#------------------------------------------DONE

#===================================#
#FDR filtering ors DEG usda not-DEG
#===================================#
fdr_trim_ors <- ors$Gene.accession[which(sub_fdr_ors)][which(sub_fdr_ors_filtered < fdr_thresh)]
fdr_trim_usda <- usda$Label[which(sub_fdr_usda)][which(sub_fdr_usda_filtered > fdr_thresh)]

sub_ok <- as.logical(fin_fin$ors_id %in% fdr_trim_ors * fin_fin$usda_id %in% fdr_trim_usda)

fin_fdr_filt_ors <- fin_fin[sub_ok,]

hm_clust_name <- "fin_fdr_filt_ors"
fin_hm_clust <- (fin_fdr_filt_ors %>% distinct())
#------------------------------------------DONE

#===================================#
#FDR filtering ors not-DEG usda DEG 
#===================================#
fdr_trim_ors <- ors$Gene.accession[which(sub_fdr_ors)][which(sub_fdr_ors_filtered > fdr_thresh)]
fdr_trim_usda <- usda$Label[which(sub_fdr_usda)][which(sub_fdr_usda_filtered < fdr_thresh)]

sub_ok <- as.logical(fin_fin$ors_id %in% fdr_trim_ors * fin_fin$usda_id %in% fdr_trim_usda)

fin_fdr_filt_usda <- fin_fin[sub_ok,]
hm_clust_name <- "fin_fdr_filt_usda"
fin_hm_clust <- (fin_fdr_filt_usda %>% distinct())
#------------------------------------------DONE


#==============================#
#FDR and LFC filtering for all
#==============================#
hm_clust_name <- "fin_fdr_lfc_filt"
fin_hm_clust <- (fin_fdr_lfc_filt %>% distinct())
#------------------------------------------DONE


# Now useless
{
    # ## Extracting subsets
    # sbset <- list(out_up = fin_fin[which(fin_fin$usda_lfc > (fin_fin$ors_lfc * 0.6 + inter)),],
    #               out_dwn = fin_fin[which(fin_fin$usda_lfc < (fin_fin$ors_lfc * 0.6 - inter)),])
    # 
    # gene_type <- character(length = dim(fin_fin)[1])
    # 
    # gene_type[which(fin_fin$ors_id %in% sbset[[1]]$ors_id)] <- "out_up"
    # gene_type[which(fin_fin$ors_id %in% sbset[[2]]$ors_id)] <- "out_dwn"
    # gene_type[which(gene_type == "")] <- "in"
    # 
    # fin_fin$gene_type <- as.factor(gene_type)
    # 
    # ## Extracting other subsets
    # {
    # sbset <- list( u_u = fin_fin[which(fin_fin$ors_lfc > 2.5 & fin_fin$usda_lfc > 2.5),],
    #                u_z = fin_fin[which((fin_fin$ors_lfc > -2.5 & fin_fin$ors_lfc < 2.5) & fin_fin$usda_lfc > 3),],
    #                u_l = fin_fin[which(fin_fin$ors_lfc < -2.5 & fin_fin$usda_lfc > 2.5),],
    #                z_u = fin_fin[which((fin_fin$ors_lfc > 2.5 & (fin_fin$usda_lfc < 2.5 & fin_fin$usda_lfc > -2.5))),],
    #                z_l = fin_fin[which((fin_fin$ors_lfc < -2.5 & (fin_fin$usda_lfc < 2.5 & fin_fin$usda_lfc > -2.5))),],
    #                l_u = fin_fin[which(fin_fin$ors_lfc > 2.5 & fin_fin$usda_lfc < -2.5),],
    #                l_z = fin_fin[which((fin_fin$ors_lfc > -2.5 & fin_fin$ors_lfc < 2.5) & fin_fin$usda_lfc < -2.5),],
    #                l_l = fin_fin[which(fin_fin$ors_lfc < -2.5 & fin_fin$usda_lfc < -2.5),])
    # 
    # a <- character(dim(fin_fin)[1])
    # levels(a) <- c("u_u", "u_z", "u_l", "z_u", "z_l", "l_u", "l_z", "l_l")
    # 
    # a[which(fin_fin$ors_id %in% sbset$u_u$ors_id)] <- "u_u"
    # a[which(fin_fin$ors_id %in% sbset$u_z$ors_id)] <- "u_z"
    # a[which(fin_fin$ors_id %in% sbset$u_l$ors_id)] <- "u_l"
    # a[which(fin_fin$ors_id %in% sbset$z_u$ors_id)] <- "z_u"
    # a[which(fin_fin$ors_id %in% sbset$z_l$ors_id)] <- "z_l"
    # a[which(fin_fin$ors_id %in% sbset$l_u$ors_id)] <- "l_u"
    # a[which(fin_fin$ors_id %in% sbset$l_z$ors_id)] <- "l_z"
    # a[which(fin_fin$ors_id %in% sbset$l_l$ors_id)] <- "l_l"
    # a[which(a == "")] <- "z"
    # 
    # 
    # fin_fin <- cbind(fin_fin, a)
    # }
}


#####
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |              Heat map                |      /
#                                               \      |            w/ clustering             |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\
#####


fin_clust <- fin_hm_clust[,3:4]
#fin_clust[,2] <- -fin_clust[,2]
row.names(fin_clust) <- make.names(as.character(paste(fin_hm_clust$ors_id, fin_hm_clust$usda_id, sep = "_-_")), unique =TRUE)

##Temp
n <- as.character(paste(fin_hm_clust$ors_id, fin_hm_clust$usda_id, table$ors_annot[which(table$ors_id %in% fin_hm_clust$ors_id)][order(table$ors_id[which(table$ors_id %in% fin_hm_clust$ors_id)])][rank(fin_hm_clust$ors_id)], sep = "_-_"))
n <- substr(n, 1, 50)
row.names(fin_clust) <- make.names(n, unique =TRUE)

fin_clust <- fin_clust[!is.na(fin_clust$usda_lfc),]


#============================#
# functions
#============================#

## fonction qui créée un fichier par cluster : extrait les cluster de "count_data"
split_cluster_table = function(nb_cluster, som_data, count_data){
    for (cluster in 1:nb_cluster){
        print(cluster)
        index_cluster = which(som_data$unit.classif == cluster)
        table = data.frame(ors = numeric(), usda = numeric())
        for (i in 1:length(index_cluster)){
            index = index_cluster[i]
            table[i,] <- count_data[index,] 
        }
        one_cluster = assign(paste("cluster", cluster , sep = "."), table, envir=parent.frame())
    }
}

# creation of cluster profiles
plot_cluster = function(nb_cluster, som_data, count_data){
    for (cluster in 1:nb_cluster){
        print(cluster)
        index_cluster = which(som_data$unit.classif == cluster)
        table = data.frame(ors = numeric(), usda = numeric())
        for (i in 1:length(index_cluster)){
            index = index_cluster[i]
            table[i,] <- count_data[index,] 
        }
        
        table_cluster = as.data.frame(table)
        dlong = melt(data.frame(gene=row.names(table_cluster), table_cluster))
        plot_cls = plot_cls = ggplot(dlong, aes(x=variable , y=value, color=gene, group=gene)) +
            geom_line(show.legend = FALSE, alpha=0.4)
        name_plot = paste("plot", cluster, sep = "_")
        assign(name_plot, plot_cls, envir=parent.frame())
    }
}

#===================================#
# estimation nombre de cluster
#===================================#

#Now useless
{
    # # mclust BIC
    # par(mfrow=c(1,1))
    # library(mclust)
    # BIC = mclustBIC(mean_DE_cpm_log2_CR, G=c(1:50))
    # plot(BIC)
    # 
    # #summary(BIC)
    # #mod1 = Mclust(relative_DE_cpm, x = BIC)
    # #summary(mod1, parameters = TRUE)
    # #plot(mod1, what = "classification")
    # #table(class, mod1$classification)
    # 
    # 
    # ## nbClust several indexes
    # #cindex ou dindex
    # nb_cpm = NbClust(mean_DE_cpm_log2_CR, diss=NULL, distance = "euclidean", 
    #                  min.nc=2, max.nc=50, method = "kmeans", 
    #                  index = "dindex", alphaBeale = 0.1)
    # hist(nb_cpm$Best.nc, breaks = max(na.omit(nb_cpm$Best.nc)))
    
    ## tuto elbow (A privilÃ©gier)
}

wss <- (nrow(fin_clust)-1)*sum(apply(fin_clust,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(fin_clust,
                                     iter.max = 100, 
                                     centers=i)$withinss)

par(mfrow=c(2,2))
plot(1:50, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#==================================================================================================#
# clustering
#==================================================================================================#
nb_clust <- 7

set.seed(7)
som_CR = som(as.matrix(fin_clust[,1:2]),  grid= somgrid(1,nb_clust, "hexagonal"), rlen = 100)
table(som_CR$unit.classif)
par(mfrow=c(2,2))
plot(som_CR, type = "counts")
plot(som_CR, type = "changes")
plot(som_CR, type = "codes")
plot(som_CR, type = "dist.neighbours")


# rm(cluster.1, cluster.2, cluster.3, cluster.4, cluster.5)
# rm(som_CR)

#===========================================================#
# creation of cluster file (all and splited) 
#===========================================================#
# split cluster
split_cluster_table(nb_cluster = nb_clust, som_data = som_CR, count_data = fin_clust)

# clusters tables
total_cluster_id = data.frame(ors = numeric(), usda = numeric(), cluster = numeric())

nb_cluster = nb_clust
list_object = ls()
cluster_list = grep(paste("cluster.*[0-", nb_clust, "]$", sep = ""), list_object, value = TRUE)
cluster_list
# creation of file containing all DE gene and the correspondant cluster number
for (i in 1:length(cluster_list)){
    cluster_name = get(cluster_list[i])
    cluster_name = cbind(cluster_name, rep(cluster_list[i], nrow(cluster_name)))
    colnames(cluster_name) = c("ors","usda","cluster")
    total_cluster_id = rbind(total_cluster_id, cluster_name)
    assign(paste(i, "cluster_plot", sep ="_"), cluster_name)
    write.csv(cluster_name, file = cluster_list[i])
}
write.csv(total_cluster_id, file = "cluster_CR_all_4")
# rm(cluster.1, cluster.2, cluster.3, cluster.4,cluster.5,cluster.6,cluster.7)

#===========================================================#
# visualization
#===========================================================#

#Now useless
{
    # ## heatmap with cluster details
    # my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 399)
    # colors = c(seq(-5,-3.501,length=100),seq(-3.5,0.5,length=100),seq(0.501,5,length=100),seq(5.001,10,length=100))
    # #attention a la colonne cluster a preciser
    # cluster_label = as.data.frame(total_cluster_id[,3])
    # row.names(cluster_label) = row.names(total_cluster_id)
    # 
    # pheatmap(total_cluster_id[,1:2], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors, 
    #          color = my_palette, show_rownames = FALSE, border_color = NA)
    # 
    # ## heatmap without cluster details
    # pheatmap(total_cluster_id[,1:2], cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors, 
    #          color = my_palette, show_rownames = TRUE, fontsize_row = 1, border_color = NA)
    # 
    
}

##========================================================##
# 4 paper - heatmap without cluster details and gene names #
##========================================================##

# my_palette <- colorRampPalette(c("violet", "blue", "black", "yellow", "orange"))(n = 399)
# colors = c(seq(-30,-5.501,length=100),seq(-5.5,0.5,length=100),seq(0.501,5,length=100),seq(5.001,30,length=100))

my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 399)
colors = c(seq(-10,-2.501,length=100),seq(-2.5,0.5,length=100),seq(0.501,2.5,length=100),seq(2.501,10,length=100))
cluster_label = as.data.frame(total_cluster_id[,3])
row.names(cluster_label) = row.names(total_cluster_id)

dev.off()
bmp(file = paste("C:/Users/quent/Desktop/ortholog_heatmap_", hm_clust_name, ".bmp", sep = ""), width = 400, height = 1200, units = "px")
pheatmap(total_cluster_id[,1:2], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors, 
         color = my_palette, show_rownames = FALSE, border_color = NA)
dev.off()
pdf(file = paste("C:/Users/quent/Desktop/ortholog_heatmap_", hm_clust_name, "_with_names.pdf", sep = ""), width = 10, height = 30)
pheatmap(total_cluster_id[,1:2], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors, 
         color = my_palette, show_rownames = TRUE, fontsize = 2, border_color = NA)
dev.off()


# Now useless
{
    # # plot profiles
    # plot_cluster(nb_cluster = 7, som_CR, fin_clust)
    # 
    # grid.arrange(plot_1, plot_2, plot_3, plot_4, plot_5, plot_6)
        
    # #Specific heatmaps
    # total_cluster_id <- total_cluster_id[order(row.names(total_cluster_id)),][rank(row.names(fin_clust)),]
    # 
    # my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 399)
    # colors = c(seq(-15,-7.501,length=100),seq(-7.5,0.5,length=100),seq(0.501,7,length=100),seq(7.001,15,length=100))
    # 
    # for (test in levels(a)) {
    #     nb <- dim(total_cluster_id[which(fin_fin3$a == test),1:2])[1]
    #     pdf(file = paste("C:/Users/quent/Desktop/ortholog_heatmap_subset_", test, ".pdf", sep = ""), width = 7, height = 1+0.14*nb)
    #     pheatmap(total_cluster_id[which(fin_fin3$a == test),1:2], cluster_rows = TRUE, cluster_cols = FALSE, breaks = colors, 
    #              color = my_palette, show_rownames = TRUE, fontsize_row = 5, border_color = NA)
    #     dev.off()
    # }
    
    
    # #Pour refaire les heatmap et cluster à partir de cluster_CR_all
    # #===========================================================#
    # # download cluster from file
    # #===========================================================#
    # cluster_table = read.table("cluster_CR_all", row.names = 1, header = TRUE, sep = ",")
    # 
    # ## heatmap with cluster details
    # my_palette <- colorRampPalette(c("green", "black", "red", "yellow"))(n = 399)
    # colors = c(seq(-2.5,-1,length=100),seq(-0.9999,0.5,length=100),seq(0.501,1.5,length=100),seq(1.501,2.5,length=100))
    # #attention a la colonne cluster a preciser
    # cluster_label = as.data.frame(cluster_table[,15])
    # row.names(cluster_label) = row.names(cluster_table)
    # pheatmap(cluster_table[,1:14], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors, 
    #          color = my_palette, show_rownames = FALSE, border_color = NA)
    # 
    # # profil
    # dlong = melt(data.frame(gene=row.names(`1_cluster_plot`), `1_cluster_plot`))
    # plot  = ggplot(dlong, aes(x=variable , y=value, color="black", group=gene)) +
    #     geom_line(show.legend = FALSE, alpha=0.4)
    # plot + stat_summary(aes(group = dlong$cluster), 
    #                     fun.y=mean, geom="line", colour="black")
}

#####
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |             Export data              |      /
#                                               \      |                                      |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\
#####



##=====================================##
##             Whole Table             ##
##=====================================##

ors_mean_YM <- ors$YM.mean.reads[which(ors$Gene.accession %in% fin_fin$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fin$ors_id)])][rank(fin_fin$ors_id)]
ors_mean_AA <- ors$AA.mean.reads[which(ors$Gene.accession %in% fin_fin$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fin$ors_id)])][rank(fin_fin$ors_id)]
ors_annot <- ors$Product[which(ors$Gene.accession %in% fin_fin$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fin$ors_id)])][rank(fin_fin$ors_id)]
ors_gene_name <- ors$Gene.name[which(ors$Gene.accession %in% fin_fin$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fin$ors_id)])][rank(fin_fin$ors_id)]

ors_fdr <- ors$AA_vs_YM.FDR[which(ors$Gene.accession %in% fin_fin$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fin$ors_id)])][rank(fin_fin$ors_id)]

ors_lfc_calc <- log(x = ors_mean_AA / ors_mean_YM, base = 2)


usda_mean_YM <- usda$YM.reads[which(usda$Label %in% fin_fin$usda_id)][order(usda$Label[which(usda$Label %in% fin_fin$usda_id)])][rank(fin_fin$usda_id)]
usda_mean_AA <- usda$AA.reads[which(usda$Label %in% fin_fin$usda_id)][order(usda$Label[which(usda$Label %in% fin_fin$usda_id)])][rank(fin_fin$usda_id)]
usda_annot <- usda$EugenePP_annotation[which(usda$Label %in% fin_fin$usda_id)][order(usda$Label[which(usda$Label %in% fin_fin$usda_id)])][rank(fin_fin$usda_id)]
usda_gene_name <- usda$gene_name[which(usda$Label %in% fin_fin$usda_id)][order(usda$Label[which(usda$Label %in% fin_fin$usda_id)])][rank(fin_fin$usda_id)]

usda_fdr <- usda$AA_vs_YM_FDR[which(usda$Label %in% fin_fin$usda_id)][order(usda$Label[which(usda$Label %in% fin_fin$usda_id)])][rank(fin_fin$usda_id)]

usda_lfc_calc <- log(x = usda_mean_AA / usda_mean_YM, base = 2)


table <- data.frame(ors_id = fin_fin$ors_id,
                    ors_gene = ors_gene_name,
                    usda_id = fin_fin$usda_id,
                    usda_gene = usda_gene_name,
                    ors_annot = ors_annot,
                    usda_annot = usda_annot,
                    ors_mean_YM = ors_mean_YM,
                    ors_mean_AA = ors_mean_AA,
                    #ors_lfc_QN = ors_lfc_calc,
                    ors_lfc = fin_fin$ors_lfc,
                    ors_fdr = ors_fdr,
                    usda_mean_YM = usda_mean_YM,
                    usda_mean_AA = usda_mean_AA,
                    #usda_lfc_QN = usda_lfc_calc,
                    usda_lfc = fin_fin$usda_lfc,
                    usda_fdr = usda_fdr)

#[order(fin_fin$usda_lfc)][rank(fin_fin$usda_id)]


openxlsx::write.xlsx(table, file= "C:/Users/quent/Desktop/flo_Orth_data_orthologs_ORS_USDA_YM_vs_AA.xlsx")


#Now useless
{
    # ##=====================================##
    # ##            Genes subsets            ##
    # ##=====================================##
    # 
    # for (cond in levels(fin_fdr_filt$gene_type)) {
    #     curr <- fin_fdr_filt[which(fin_fdr_filt$gene_type == cond),]
    #     
    #     ors_annot_s <- ors$Product[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)])][rank(fin_fdr_filt$ors_id)]
    #     ors_annot_s <- ors_annot_s[which(fin_fdr_filt$gene_type == cond)]
    #     
    #     usda_annot_s <- usda$EugenePP_annotation[which(usda$Label %in% fin_fdr_filt$usda_id)][order(usda$Label[which(usda$Label %in% fin_fdr_filt$usda_id)])][rank(fin_fdr_filt$usda_id)]
    #     usda_annot_s <- usda_annot_s[which(fin_fdr_filt$gene_type == cond)]
    #     
    #     curr <- cbind(curr, ors_annot_s, usda_annot_s)
    #     
    #     xlsx::write.xlsx(curr, file = "C:/Users/quent/Desktop/ortholog_genes_subsets.xlsx", sheetName = cond, append = TRUE)
    # }
}

if(!require(xlsx)) {
    instal.packages("xlsx")
    library(xlsx)
}


ors_mean_YM <- ors$YM.mean.reads[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)])][rank(fin_fdr_filt$ors_id)]
ors_mean_AA <- ors$AA.mean.reads[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)])][rank(fin_fdr_filt$ors_id)]
ors_annot <- ors$Product[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)])][rank(fin_fdr_filt$ors_id)]
ors_gene_name <- ors$Gene.name[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)])][rank(fin_fdr_filt$ors_id)]

ors_fdr <- ors$AA_vs_YM.FDR[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)][order(ors$Gene.accession[which(ors$Gene.accession %in% fin_fdr_filt$ors_id)])][rank(fin_fdr_filt$ors_id)]

usda_mean_YM <- usda$YM.reads[which(usda$Label %in% fin_fdr_filt$usda_id)][order(usda$Label[which(usda$Label %in% fin_fdr_filt$usda_id)])][rank(fin_fdr_filt$usda_id)]
usda_mean_AA <- usda$AA.reads[which(usda$Label %in% fin_fdr_filt$usda_id)][order(usda$Label[which(usda$Label %in% fin_fdr_filt$usda_id)])][rank(fin_fdr_filt$usda_id)]
usda_annot <- usda$EugenePP_annotation[which(usda$Label %in% fin_fdr_filt$usda_id)][order(usda$Label[which(usda$Label %in% fin_fdr_filt$usda_id)])][rank(fin_fdr_filt$usda_id)]
usda_gene_name <- usda$gene_name[which(usda$Label %in% fin_fdr_filt$usda_id)][order(usda$Label[which(usda$Label %in% fin_fdr_filt$usda_id)])][rank(fin_fdr_filt$usda_id)]

usda_fdr <- usda$AA_vs_YM_FDR[which(usda$Label %in% fin_fdr_filt$usda_id)][order(usda$Label[which(usda$Label %in% fin_fdr_filt$usda_id)])][rank(fin_fdr_filt$usda_id)]


xlsx::write.xlsx(fin_fdr_filt, file = "C:/Users/quent/Desktop/DEG-notDEG_subsets.xlsx", sheetName = "all_DEG")
xlsx::write.xlsx(fin_fdr_filt_not, file = "C:/Users/quent/Desktop/DEG-notDEG_subsets.xlsx", sheetName = "no_DEG", append = TRUE)
xlsx::write.xlsx(fin_fdr_filt_ors, file = "C:/Users/quent/Desktop/DEG-notDEG_subsets.xlsx", sheetName = "ors_DEG", append = TRUE)
xlsx::write.xlsx(fin_fdr_filt_usda, file = "C:/Users/quent/Desktop/DEG-notDEG_subsets.xlsx", sheetName = "usda_DEG", append = TRUE)
