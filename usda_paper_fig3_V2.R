rm(list = ls())

#Load Packages
{
    if (! require(openxlsx)) {
        install.packages("openxlsx")
        library(openxlsx)
    }
    
    # if (! require(xlsx)) {
    #     install.packages("xlsx")
    #     library(xlsx)
    # }
    # 
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
        if (!requireNamespace("BiocManager", quietly = TRUE))
           install.packages("BiocManager")
      
        BiocManager::install("edgeR")
        library(edgeR)
    }
}

# Get and Format data =============================================================================================================================
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |         Get and Format data          |      /
#                                               \      |                                      |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\



# data <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/USDA_mSystems_SuppTables_1704_ qn.xlsx", sheet = "supp table 1")
# colnames(data) <- make.names(data[1,]) 
# data <- data[2:dim(data)[1],]
# # data <- data[which(rowSums(is.na(data[12:14])) == 0),]
# # data <- data[which(rowSums(is.na(data[23:25])) == 0),]
# row.names(data) <- make.names(data$Label)
# 
# 
# data_rna <- data[,c(3,4,8,12:14)]
# data_rna <- data_rna[which(rowSums(is.na(data_rna[4:6])) == 0),]
# data_rna_data <- as.data.frame(lapply(data_rna[4:6], as.numeric))
# data_rna_rpkm <- (data_rna_data / (colSums(data_rna_data) / 1000000))  / as.numeric(data_rna$Length)
# data_rna[4:6] <- data_rna_rpkm
# 
# 
# data_sc <- data[,c(3,4,8,23:25)]
# data_sc <- data_sc[which(rowSums(is.na(data_sc[4:6])) == 0),]
# data_sc_data <- as.data.frame(lapply(data_sc[4:6], as.numeric))
# data_sc_rpkm <- (data_sc_data / (colSums(data_sc_data) / 1000))  / as.numeric(data_sc$Length)
# data_sc[4:6] <- data_sc_rpkm
# 
# 
# data_int <- intersect(data_sc$Label, data$Label)
# data_int <- intersect(data_rna$Label, data_int)
# 
# data_fin <- data[data_int, 3:4]
# data_fin <- cbind(data_fin, data_sc[,4:6], data_rna[data_int,4:6])
# 
# dap_deg_nod_raw <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/USDA_mSystems_SuppTables_1704_ qn.xlsx", sheet = "DEG_DAP_nod")
# colnames(dap_deg_nod_raw) <- make.names(dap_deg_nod_raw[1,]) 
# dap_deg_nod_raw <- dap_deg_nod_raw[2:dim(dap_deg_nod_raw)[1],]
# names_later <- paste(dap_deg_nod_raw$Label, dap_deg_nod_raw$gene_name, sep = "_")
# dap_deg_nod_raw <- dap_deg_nod_raw$Label
# dap_deg_nod <- data_fin[dap_deg_nod_raw, 3:8]
# row.names(dap_deg_nod) <- make.names(names_later)
# 
# aa_gm_raw <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/USDA_mSystems_SuppTables_1704_ qn.xlsx", sheet = "AA>GM")
# colnames(aa_gm_raw) <- make.names(aa_gm_raw[1,]) 
# aa_gm_raw <- aa_gm_raw[2:dim(aa_gm_raw)[1],]
# names_later <- paste(aa_gm_raw$Label, aa_gm_raw$gene_name, sep = "_")
# aa_gm_raw <- aa_gm_raw$Label
# aa_gm <- data_fin[aa_gm_raw, 3:8]
# row.names(aa_gm) <- make.names(names_later)
# 
# gm_aa_raw <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/USDA_mSystems_SuppTables_1704_ qn.xlsx", sheet = "GM>AA")
# colnames(gm_aa_raw) <- make.names(gm_aa_raw[1,])
# gm_aa_raw <- gm_aa_raw[2:dim(gm_aa_raw)[1],]
# names_later <- paste(gm_aa_raw$Label, gm_aa_raw$gene_name, sep = "_")
# gm_aa_raw <- gm_aa_raw$Label
# gm_aa <- data_fin[gm_aa_raw, 3:8]
# row.names(gm_aa) <- make.names(names_later)
# 
# wf <- list(dap_deg_nod = dap_deg_nod, aa_gm = aa_gm, gm_aa = gm_aa)

dap_deg_nod_raw <- openxlsx::read.xlsx(xlsxFile = "C:/Users/quent/Desktop/Global_t+p_table_with_all_genes_2020.xlsx", sheet = "DAP DEG Bacteroids vs YM")
dap_deg_nod <- cbind(dap_deg_nod_raw[,c(3,4)], as.data.frame(lapply(dap_deg_nod_raw[,c(15,16,17,42,43,44)], as.numeric), stringsAsFactors = FALSE))


aa_gm_raw <- openxlsx::read.xlsx(xlsxFile = "D:/Work/temp_desktop/Global_t+p_table_with_all_genes_2020.xlsx", sheet = "DAP DEG AA VS GM")
aa_gm <- cbind(aa_gm_raw[,c(3,4)], as.data.frame(lapply(aa_gm_raw[,c(15,16,17,42,43,44)], as.numeric), stringsAsFactors = FALSE))

aa_gm_sup <- aa_gm[which(aa_gm_raw$AA_vs_GM_LFC > 0),]
prot_sup <- (aa_gm_sup$AA_SC+1) / (aa_gm_sup$GM_SC+1)
aa_gm_sup_trim <- aa_gm_sup[which(prot_sup >1),]

aa_gm_inf <- aa_gm[which(aa_gm_raw$AA_vs_GM_LFC < 0),]
prot_inf <- (aa_gm_inf$GM_SC+1) / (aa_gm_inf$AA_SC+1)
aa_gm_inf_trim <- aa_gm_inf[which(prot_inf >1),]

all_res <- list(dap_deg_nod = dap_deg_nod, aa_gm = aa_gm_sup_trim, gm_aa = aa_gm_inf_trim)
wf <- lapply(all_res, format_USDA_datasets)
  

# Heat map w/ clustering =============================================================================================================================
#                                                         _______________________________________
#                                               _________|                                      |_______
#                                               \       |              Heat map                |      /
#                                               \      |            w/ clustering             |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\


lapply(wf, nb_cluster_estimation)

nb_clust <- 4

som_cr <- lapply(wf, clustering, nb_clust)

total_clust_id <- list()
for (i in 1:length(som_cr)) {
  total_clust_id[[i]] <- creation_cluster(som_cr[[i]], wf[[i]],  nb_clust)
  
  wf_temp <- wf[[i]][order(row.names(wf[[i]])),][rank(row.names(total_clust_id[[i]])),]
  total_clust_id[[i]] <- cbind(total_clust_id[[i]][,1:3], wf_temp[,4:6], cluster = total_clust_id[[i]][,4])
}
names(total_clust_id) <- make.names(names(som_cr))


for (i in 1:length(total_clust_id)) {
  draw_hm_clust(total_clust_id[[i]], names(total_clust_id[i]))
}

# for (i in 1:length(total_clust_id)) {
#   draw_hm_clust_v2(total_clust_id[[i]], names(total_clust_id[i]))
# }



all_res$dap_deg_nod[all_res$dap_deg_nod$Label == "bll6927",]
all_res$dap_deg_nod[243,]
wf$dap_deg_nod[243,]
wf$dap_deg_nod["bll6927_hypE_2",]
x_fin["blr6979_groEL_6",]
x_fin["bll6927_hypE_2",]




# Heat map w/ clustering - Functions =============================================================================================================================
#                                                         _______________________________________
#                                               _________|              Heat map                |_______
#                                               \       |           w/ clustering              |      /
#                                               \      |              Functions               |     /
#                                              /      |______________________________________|     \
#                                            /__________)                                (_________\


format_USDA_datasets <- function(x) {
  
  if (sum( x[,3:5] == 0))
    x[,3:5] <- x[,3:5] + 1
  
  if (sum( x[,6:8] == 0))
    x[,6:8] <- x[,6:8] + 1
  
  x_sub <- x[,c(3,4,5,6,7,8)]
  
  # x_sub <- x_sub[which(rowSums(is.na(x_sub)) == 0),]  
  
  x_fin_prot <- (x_sub[4:6])
  x_fin_prot <- log2(x_fin_prot)
  # x_fin_prot <- scale(x = t(x_fin_prot), center = TRUE, scale = FALSE)
  # x_fin_prot <- t(x_fin_prot)
  
  x_fin_rna <- (x_sub[,1:3])
  x_fin_rna <- log2(x_fin_rna)
  # x_fin_rna <- scale(x = t(x_fin_rna), center =, TRUE, scale = FALSE)
  # x_fin_rna <- t(x_fin_rna)

  x_fin <- cbind(x_fin_prot, x_fin_rna) 
  
  row.names(x_fin) <- make.names(paste(x$Label, x$gene_name, sep = "_"), unique = TRUE) 
  
  x_fin <- x_fin[!is.na(x_fin[,1]),]
  
  return(as.data.frame(x_fin))
}


# estimation nombre de cluster
nb_cluster_estimation <- function(clust_ready) {
    wss <- (nrow(clust_ready)-1)*sum(apply(clust_ready,2,var))
    for (i in 2:50) wss[i] <- sum(kmeans(clust_ready,
                                         iter.max = 100, centers=i)$withinss)
    x11()
    par(mfrow=c(2,2))
    plot(1:50, wss, type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares")
} 

## fonction qui créée un fichier par cluster : extrait les cluster de "count_data"
split_cluster_table = function(nb_cluster, som_data, count_data){
    for (cluster in 1:nb_cluster){
        print(cluster)
        index_cluster = which(som_data$unit.classif == cluster)
        table = data.frame(AA_SC = numeric(), GM_SC = numeric(), YM_SC = numeric())
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
        table = data.frame(AA_SC = numeric(), GM_SC = numeric(), YM_SC = numeric())
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

# clustering

clustering <- function(data, nb_clust) {
    set.seed(7)
    som_CR <- som(as.matrix(data[,1:3]),  grid= somgrid(1,nb_clust, "hexagonal"), rlen = 100)
    table(som_CR$unit.classif)
    
    x11()
    par(mfrow=c(2,2))
    plot(som_CR, type = "counts")
    plot(som_CR, type = "changes")
    plot(som_CR, type = "codes")
    plot(som_CR, type = "dist.neighbours")
    
    return(som_CR)
    
}

# creation of cluster file (all and splited) 
creation_cluster <- function(som_CR, data, nb_clust) {
    data <- data[,1:3]
    # split cluster
    split_cluster_table(nb_cluster = nb_clust, som_data = som_CR, count_data = data)
    
    # clusters tables
    total_cluster_id = data.frame(AA_SC = numeric(), GM_SC = numeric(), YM_SC = numeric(), AA.reads = numeric(), GM.reads = numeric(), YM.reads = numeric(), cluster = numeric())
    
    list_object = ls()
    cluster_list = grep(paste("cluster.*[0-", nb_clust, "]$", sep = ""), list_object, value = TRUE)

    # creation of file containing all DE gene and the correspondant cluster number
    for (i in 1:length(cluster_list)){
        cluster_name = get(cluster_list[i])
        cluster_name = cbind(cluster_name, rep(cluster_list[i], nrow(cluster_name)))
        colnames(cluster_name) = c("AA_SC", "GM_SC", "YM_SC", "cluster")
        total_cluster_id = rbind(total_cluster_id, cluster_name)
        assign(paste(i, "cluster_plot", sep ="_"), cluster_name)
        write.csv(cluster_name, file = cluster_list[i])
    }
    write.csv(total_cluster_id, file = "cluster_CR_all_4")
    
    return(total_cluster_id)
}

# heatmap without cluster details and gene names 
draw_hm_clust <- function(total_cluster_id, name) {
    my_palette_prot <- colorRampPalette(c("black", "yellow"))(n = 399)
    colors_prot = c(seq(0,1.66,length=100),seq(1.6601,3.33,length=100), seq(3.330001,5,length=100))
    
    my_palette_rna <- colorRampPalette(c("black", "yellow"))(n = 399)
    colors_rna = c(seq(0,7.5,length=100),seq(7.50001,10,length=100), seq(10.0001,15,length=100))
    
    cluster_label = as.data.frame(total_cluster_id[,7])
    row.names(cluster_label) = row.names(total_cluster_id)
    
    while (!is.null(dev.list()))  dev.off()
    bmp(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_prot.bmp", sep = ""), width = 400, height = 1200, units = "px")
    pheatmap(total_cluster_id[,1:3], annotation_row = NULL, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_prot, 
             color = my_palette_prot, show_rownames = FALSE, show_colnames = FALSE, border_color = "black", legend =FALSE, annotation_legend = FALSE)
    dev.off()
    bmp(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_rna.bmp", sep = ""), width = 400, height = 1200, units = "px")
    pheatmap(total_cluster_id[,4:6], annotation_row = NULL, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_rna, 
             color = my_palette_rna, show_rownames = FALSE, show_colnames = FALSE, border_color = "black", legend =FALSE, annotation_legend = FALSE)
    dev.off()
    
    
    pdf(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_prot_with_names.pdf", sep = ""), width = 10, height = 30)
    pheatmap(total_cluster_id[,1:3], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_prot, 
             color = my_palette_prot, show_rownames = TRUE, fontsize = 2, border_color = NA)
    dev.off()
    pdf(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_rna_with_names.pdf", sep = ""), width = 10, height = 30)
    pheatmap(total_cluster_id[,4:6]/2, annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_rna, 
             color = my_palette_rna, show_rownames = TRUE, fontsize = 2, border_color = NA)
    dev.off()
} 
# # heatmap without cluster details and gene names 
# draw_hm_clust <- function(total_cluster_id, name) {
#   my_palette_prot <- colorRampPalette(c("black", "orange"))(n = 399)
#   colors_prot = c(seq(0,1.66,length=100),seq(1.6601,3.33,length=100), seq(3.330001,5,length=100))
#   
#   my_palette_rna <- colorRampPalette(c("black", "yellow"))(n = 399)
#   colors_rna = c(seq(0,4,length=100),seq(4.0001,6,length=100), seq(6.0001,7.5,length=100))
#   
#   cluster_label = as.data.frame(total_cluster_id[,7])
#   row.names(cluster_label) = row.names(total_cluster_id)
#   
#   while (!is.null(dev.list()))  dev.off()
#   bmp(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_prot.bmp", sep = ""), width = 400, height = 1200, units = "px")
#   pheatmap(total_cluster_id[,1:3], annotation_row = NULL, cluster_rows = TRUE, clustering_method = "mcquitty", cluster_cols = FALSE, breaks = colors_prot, 
#            color = my_palette_prot, show_rownames = FALSE, show_colnames = FALSE, border_color = "black", legend =FALSE, annotation_legend = FALSE)
#   dev.off()
#   bmp(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_rna.bmp", sep = ""), width = 400, height = 1200, units = "px")
#   pheatmap(total_cluster_id[,4:6]/2, annotation_row = NULL, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_rna, 
#            color = my_palette_rna, show_rownames = FALSE, show_colnames = FALSE, border_color = "black", legend =FALSE, annotation_legend = FALSE)
#   dev.off()
#   
#   
#   pdf(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_prot_with_names.pdf", sep = ""), width = 10, height = 30)
#   pheatmap(total_cluster_id[,1:3], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_prot, 
#            color = my_palette_prot, show_rownames = TRUE, fontsize = 2, border_color = NA)
#   dev.off()
#   pdf(file = paste("C:/Users/quent/Desktop/fig_3_heatmap_", name, "_rna_with_names.pdf", sep = ""), width = 10, height = 30)
#   pheatmap(total_cluster_id[,4:6]/2, annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_rna, 
#            color = my_palette_rna, show_rownames = TRUE, fontsize = 2, border_color = NA)
#   dev.off()
# } 
# 
# draw_hm_clust_v2 <- function(total_cluster_id, name) {
#   my_palette <- colorRampPalette(c("black", "yellow"))(n = 399)
#   colors_rna = c(seq(0.00,0.25,length=100), seq(0.2501,0.5,length=100), seq(0.5001,1,length=100))
#   colors_prot = c(seq(-1,-0.501,length=100),seq(-0.5,0,length=100),seq(0.01,0.5,length=100), seq(0.5001,1,length=100))
#   cluster_label = as.data.frame(total_cluster_id[,7])
#   row.names(cluster_label) = row.names(total_cluster_id)
#   
#   while (!is.null(dev.list()))  dev.off()
#   bmp(file = paste("C:/Users/quent/Desktop/ortholog_heatmap_", name, ".bmp", sep = ""), width = 400, height = 1200, units = "px")
#   pheatmap(total_cluster_id[,1:6], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_rna, 
#            color = my_palette, show_rownames = FALSE, border_color = NA)
#   dev.off()
#   pdf(file = paste("C:/Users/quent/Desktop/ortholog_heatmap_", name, "_with_names.pdf", sep = ""), width = 10, height = 30)
#   pheatmap(total_cluster_id[,1:6], annotation_row = cluster_label, cluster_rows = FALSE, cluster_cols = FALSE, breaks = colors_rna, 
#            color = my_palette, show_rownames = TRUE, fontsize = 2, border_color = NA)
#   dev.off()
# } 

