rm(list = ls())

if (!require(xlsx)){
    install.packages("xlsx")
    library(xlsx)
}

if (!require(tidyr)){
    install.packages("tidyr")
    library(tidyr)
}

table <- read.csv(file = "C:/Users/quent/Desktop/elArtist_data.csv", sep = ";")
table <- table[c(1:4,6:8,10,12)]

bnm_bnm_ess_genes <- table[which(table[,6] %in% c(2,3)),c(1,8,9)]
write.xlsx(bnm_bnm_ess_genes, file = "c:/Users/quent/Desktop/ess_bnm.xlsx")

cog_annot <- read.xlsx(file = "C:/Users/quent/Desktop/BRAD285.2-CDSs-COG.xlsx", sheetIndex = 1)
cog_classes <- data.frame(id = sapply(cog_annot[,1], function(x) {return(strsplit(as.character(x), split="|", fixed = TRUE)[[1]][1])}), 
                          cog = cog_annot$S)

####################################################################################################################################################
# Functions

genEssTable <- function(table, data, multiple) {
    
    if (multiple < 1 || multiple > 3) stop("'multiple' must take the value 1, 2 or 3")
    
    if (multiple == 1) {
        ess_table = matrix(data = FALSE, nrow = dim(table)[1], ncol = length(data))
        
        for ( i in 1:dim(table)[1] ) {
            
            ess_table[i,] = idSpeEss(table[i,data])
        }
    } else if (multiple == 2) {
        comb <- combn(data, 2)
        res = matrix(data = FALSE, nrow = dim(table)[1], ncol = dim(comb)[2])
        
        for ( i in 1:dim(table)[1] ) {
            
            res[i,] = idBiEss(table[i,data])
        }
        col_names <- c(paste(colnames(table)[comb[1,1]], colnames(table)[comb[2,1]], sep = "_"), 
                       paste(colnames(table)[comb[1,2]], colnames(table)[comb[2,2]], sep = "_"),
                       paste(colnames(table)[comb[1,3]], colnames(table)[comb[2,3]], sep = "_"))
        colnames(res) <- col_names
        ess_table = list(comb = combn(data, 2), results = res)
    } else {
        ess_table <- logical(dim(table)[1])
        
        for ( i in 1:dim(table)[1] ) {
            
            ess_table[i] = idAllEss(table[i,data])
        }
    }
    
    return(ess_table)
}


idSpeEss <- function(reg) {
    
    reg_state = logical(3)
    
    for (i in 1:3) {
        
        switch(EXPR = i,
               "1" = {other = c(2,3)},
               "2" = {other = c(1,3)},
               "3" = {other = c(1,2)})
        
        
        if ((reg[1,i] == 2 || reg[1,i] == 3) && 
            (reg[1,other[1]] != 2 && reg[1,other[2]] != 2) && 
            (reg[1,other[1]] != 3 && reg[1,other[2]] != 3)) {
            reg_state[i] = 1
        } else {
            reg_state[i] = 0
        }
        
    }
    return(reg_state)
}


idBiEss <- function(reg) {
    
    reg_state = logical(dim(combn(1:length(reg), 2))[2])
    
    for (i in 1:dim(combn(1:length(reg), 2))[2]) {
        
        comb = combn(1:length(reg), 2)
        
        if ((reg[1,comb[1,i]] == 2 || reg[1,comb[1,i]] == 3) && 
            (reg[1,comb[2,i]] == 2 || reg[1,comb[2,i]] == 3) && 
            (reg[1,-comb[,i]] != 2 && reg[1,-comb[,i]] != 3)) {
            reg_state[i] = 1
        } else {
            reg_state[i] = 0
        }
        
    }
    return(reg_state)
}


idAllEss <- function(reg) {
    
    reg_state <- logical(1)
    
    if ((reg[1,1] == 2 || reg[1,1] == 3) && (reg[1,2] == 2 || reg[1,2] == 3) && (reg[1,3] == 2 || reg[1,3] == 3)) {
        reg_state <- TRUE
    } else {
        reg_state <- FALSE
    }
    
    return(reg_state)
}


assignCog <- function(id, cog_table) {
    if (!(id %in% cog_table$id)) {
        return(NA)
    } else {
        return(as.character(cog_table$cog[which(cog_table$id %in% id)]))
    }
}

####################################################################################################################################################
# YM selection media

a <- genEssTable(table, 2:4, 1)

mix_ym <- data.frame(reg = table[which(as.logical(a[,1])),1], gene = table[which(as.logical(a[,1])),8], annot = table[which(as.logical(a[,1])),9] )
ym_ym <- data.frame(reg = table[which(as.logical(a[,2])),1], gene = table[which(as.logical(a[,2])),8], annot = table[which(as.logical(a[,2])),9] )
bnm_ym <- data.frame(reg = table[which(as.logical(a[,3])),1], gene = table[which(as.logical(a[,3])),8], annot = table[which(as.logical(a[,3])),9] )


all_a <- list(mix_ym,ym_ym,bnm_ym)
len_a <- lapply(all_a, dim)

len_ini_a = list(length(which(table[,2] %in% c(2,3))), length(which(table[,3] %in% c(2,3))), length(which(table[,4] %in% c(2,3))))
len_ini_a_shrt = list(length(which(table[which(table$Annotation != ""),2] %in% c(2,3))), 
                      length(which(table[which(table$Annotation != ""),3] %in% c(2,3))), 
                      length(which(table[which(table$Annotation != ""),4] %in% c(2,3))))


write.xlsx(mix_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "mix_ym")
write.xlsx(ym_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "ym_ym", append = TRUE)
write.xlsx(bnm_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "bnm_ym", append = TRUE)


mix_ym_short <- mix_ym[which(mix_ym$annot != ""),]
mix_ym_short <- cbind(mix_ym_short, cog = as.character(rep(x = "a", length.out = dim(mix_ym_short)[1])))
mix_ym_short$cog <- sapply(mix_ym_short$reg, assignCog, cog = cog_classes)

ym_ym_short <- ym_ym[which(ym_ym$annot != ""),]
ym_ym_short <- cbind(ym_ym_short, cog = rep(x = "a", length.out = dim(ym_ym_short)[1]))
ym_ym_short$cog <- sapply(ym_ym_short$reg, assignCog, cog = cog_classes)

bnm_ym_short <- bnm_ym[which(bnm_ym$annot != ""),]
bnm_ym_short <- cbind(bnm_ym_short, cog = rep(x = "a", length.out = dim(bnm_ym_short)[1]))
bnm_ym_short$cog <- sapply(bnm_ym_short$reg, assignCog, cog = cog_classes)


all_a_shrt <- list(mix_ym_short,ym_ym_short,bnm_ym_short)
len_a_shrt <- lapply(all_a_shrt, dim)

write.xlsx(mix_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "mix_ym")
write.xlsx(ym_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "ym_ym", append = TRUE)
write.xlsx(bnm_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "bnm_ym", append = TRUE)


a_bis <- genEssTable(table, 2:4, 2)
a_bis <- a_bis[[2]]

mix_ym__ym_ym <- data.frame(reg = table[which(as.logical(a_bis[,1])),1], 
                            gene = table[which(as.logical(a_bis[,1])),8], 
                            annot = table[which(as.logical(a_bis[,1])),9] )

mix_ym__bnm_ym <- data.frame(reg = table[which(as.logical(a_bis[,2])),1], 
                             gene = table[which(as.logical(a_bis[,2])),8], 
                             annot = table[which(as.logical(a_bis[,2])),9] )

ym_ym__bnm_ym <- data.frame(reg = table[which(as.logical(a_bis[,3])),1], 
                           gene = table[which(as.logical(a_bis[,3])),8], 
                           annot = table[which(as.logical(a_bis[,3])),9] )

all_a_bis <- list(mix_ym__ym_ym,mix_ym__bnm_ym,ym_ym__bnm_ym)
len_a_bis <- lapply(all_a_bis, dim)

write.xlsx(mix_ym__ym_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "mix_ym__ym_ym", append = TRUE)
write.xlsx(mix_ym__bnm_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "mix_ym__bnm_ym", append = TRUE)
write.xlsx(ym_ym__bnm_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "ym_ym_bnm_ym", append = TRUE)


mix_ym__ym_ym_short <- mix_ym__ym_ym[which(mix_ym__ym_ym$annot != ""),]
mix_ym__ym_ym_short <- cbind(mix_ym__ym_ym_short, cog = rep(x = "a", length.out = dim(mix_ym__ym_ym_short)[1]))
mix_ym__ym_ym_short$cog <- sapply(mix_ym__ym_ym_short$reg, assignCog, cog = cog_classes)

mix_ym__bnm_ym_short <- mix_ym__bnm_ym[which(mix_ym__bnm_ym$annot != ""),]
mix_ym__bnm_ym_short <- cbind(mix_ym__bnm_ym_short, cog = rep(x = "a", length.out = dim(mix_ym__bnm_ym_short)[1]))
mix_ym__bnm_ym_short$cog <- sapply(mix_ym__bnm_ym_short$reg, assignCog, cog = cog_classes)

ym_ym__bnm_ym_short <- ym_ym__bnm_ym[which(ym_ym__bnm_ym$annot != ""),]
ym_ym__bnm_ym_short <- cbind(ym_ym__bnm_ym_short, cog = rep(x = "a", length.out = dim(ym_ym__bnm_ym_short)[1]))
ym_ym__bnm_ym_short$cog <- sapply(ym_ym__bnm_ym_short$reg, assignCog, cog = cog_classes)


all_a_bis_shrt <- list(mix_ym__ym_ym_short,mix_ym__bnm_ym_short,ym_ym__bnm_ym_short)
len_a_bis_shrt <- lapply(all_a_bis_shrt, dim)


write.xlsx(mix_ym__ym_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "mix_ym__ym_ym", append = TRUE)
write.xlsx(mix_ym__bnm_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "mix_ym__bnm_ym", append = TRUE)
write.xlsx(ym_ym__bnm_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "ym_ym__bnm_ym", append = TRUE)



a_all <- genEssTable(table, 2:4, 3)
a_all_shrt <- a_all[which(table$Annotation != "")]


common_ess_all_ym <- data.frame(reg = table[which(a_all),1], gene = table[which(a_all),8], annot = table[which(a_all),9] )

write.xlsx(common_ess_all_ym, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "ym_all_cond_ess_common", append = TRUE)

common_ess_all_ym_short <- common_ess_all_ym[which(common_ess_all_ym$annot != ""),]

write.xlsx(common_ess_all_ym_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "ym_all_cond_ess_common", append = TRUE)

( sum(a_all) + len_a[[1]][1] + len_a_bis[[1]][1] + len_a_bis[[2]][1] ) / len_ini_a[[1]][1]

## Summary table

summary_a <- data.frame(totMix = len_ini_a [[1]][1], totYM = len_ini_a [[2]][1], totBNM = len_ini_a [[3]][1], common_ess = sum(a_all), 
                        mix_ym = len_a_bis[[1]][1], mix_bnm = len_a_bis[[2]][1], ym_bnm = len_a_bis[[3]][1], 
                        mix = len_a[[1]][1], ym = len_a[[2]][1], bnm = len_a[[3]][1])

summary_a_shrt <- data.frame(totMix = len_ini_a_shrt [[1]][1], totYM = len_ini_a_shrt [[2]][1], totBNM = len_ini_a_shrt [[3]][1], common_ess = sum(a_all_shrt), 
                             mix_ym = len_a_bis_shrt [[1]][1], mix_bnm = len_a_bis_shrt [[2]][1], ym_bnm = len_a_bis_shrt [[3]][1], 
                             mix = len_a_shrt[[1]][1], ym = len_a_shrt[[2]][1], bnm = len_a_shrt[[3]][1])

####################################################################################################################################################
# BNM selection media

b <- genEssTable(table, 5:7, TRUE)

mix_bnm <- data.frame(reg = table[which(as.logical(b[,1])),1], gene = table[which(as.logical(b[,1])),8], annot = table[which(as.logical(b[,1])),9] )
ym_bnm <- data.frame(reg = table[which(as.logical(b[,2])),1], gene = table[which(as.logical(b[,2])),8], annot = table[which(as.logical(b[,2])),9] )
bnm_bnm <- data.frame(reg = table[which(as.logical(b[,3])),1], gene = table[which(as.logical(b[,3])),8], annot = table[which(as.logical(b[,3])),9] )

all_b <- list(mix_bnm,ym_bnm,bnm_bnm)
len_b <- lapply(all_b, dim)

len_ini_b = list(length(which(table[,5] %in% c(2,3))), length(which(table[,6] %in% c(2,3))), length(which(table[,7] %in% c(2,3))))
len_ini_b_shrt = list(length(which(table[which(table$Annotation != ""),5] %in% c(2,3))), length(which(table[which(table$Annotation != ""),6] %in% c(2,3))), length(which(table[which(table$Annotation != ""),7] %in% c(2,3))))


write.xlsx(mix_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "mix_bnm", append = TRUE)
write.xlsx(ym_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "ym_bnm", append = TRUE)
write.xlsx(bnm_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "bnm_bnm", append = TRUE)


mix_bnm_short <- mix_bnm[which(mix_bnm$annot != ""),]
mix_bnm_short <- cbind(mix_bnm_short, cog = rep(x = "a", length.out = dim(mix_bnm_short)[1]))
mix_bnm_short$cog <- sapply(mix_bnm_short$reg, assignCog, cog = cog_classes)

ym_bnm_short <- ym_bnm[which(ym_bnm$annot != ""),]
ym_bnm_short <- cbind(ym_bnm_short, cog = rep(x = "a", length.out = dim(ym_bnm_short)[1]))
ym_bnm_short$cog <- sapply(ym_bnm_short$reg, assignCog, cog = cog_classes)

bnm_bnm_short <- bnm_bnm[which(bnm_bnm$annot != ""),]
bnm_bnm_short <- cbind(bnm_bnm_short, cog = rep(x = "a", length.out = dim(bnm_bnm_short)[1]))
bnm_bnm_short$cog <- sapply(bnm_bnm_short$reg, assignCog, cog = cog_classes)


all_b_shrt <- list(mix_bnm_short,ym_bnm_short,bnm_bnm_short)
len_b_shrt <- lapply(all_b_shrt, dim)

write.xlsx(mix_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "mix_bnm", append = TRUE)
write.xlsx(ym_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "ym_bnm", append = TRUE)
write.xlsx(bnm_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "bnm_bnm", append = TRUE)



b_bis <- genEssTable(table, 5:7, 2)
b_bis <- b_bis[[2]]

mix_bnm__ym_bnm <- data.frame(reg = table[which(as.logical(b_bis[,1])),1], gene = table[which(as.logical(b_bis[,1])),8], annot = table[which(as.logical(b_bis[,1])),9] )
mix_bnm__bnm_bnm <- data.frame(reg = table[which(as.logical(b_bis[,2])),1], gene = table[which(as.logical(b_bis[,2])),8], annot = table[which(as.logical(b_bis[,2])),9] )
ym_bnm__bnm_bnm <- data.frame(reg = table[which(as.logical(b_bis[,3])),1], gene = table[which(as.logical(b_bis[,3])),8], annot = table[which(as.logical(b_bis[,3])),9] )

all_b_bis <- list(mix_bnm__ym_bnm,mix_bnm__bnm_bnm,ym_bnm__bnm_bnm)
len_b_bis <- lapply(all_b_bis, dim)

write.xlsx(mix_bnm__ym_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "mix_bnm__ym_bnm", append = TRUE)
write.xlsx(mix_bnm__bnm_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "mix_bnm__bnm_bnm", append = TRUE)
write.xlsx(ym_bnm__bnm_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "ym_bnm__bnm_bnm", append = TRUE)

mix_bnm__ym_bnm_short <- mix_bnm__ym_bnm[which(mix_bnm__ym_bnm$annot != ""),]
mix_bnm__ym_bnm_short <- cbind(mix_bnm__ym_bnm_short, cog = rep(x = "a", length.out = dim(mix_bnm__ym_bnm_short)[1]))
mix_bnm__ym_bnm_short$cog <- sapply(mix_bnm__ym_bnm_short$reg, assignCog, cog = cog_classes)

mix_bnm__bnm_bnm_short <- mix_bnm__bnm_bnm[which(mix_bnm__bnm_bnm$annot != ""),]
mix_bnm__bnm_bnm_short <- cbind(mix_bnm__bnm_bnm_short, cog = rep(x = "a", length.out = dim(mix_bnm__bnm_bnm_short)[1]))
mix_bnm__bnm_bnm_short$cog <- sapply(mix_bnm__bnm_bnm_short$reg, assignCog, cog = cog_classes)

ym_bnm__bnm_bnm_short <- ym_bnm__bnm_bnm[which(ym_bnm__bnm_bnm$annot != ""),]
ym_bnm__bnm_bnm_short <- cbind(ym_bnm__bnm_bnm_short, cog = rep(x = "a", length.out = dim(ym_bnm__bnm_bnm_short)[1]))
ym_bnm__bnm_bnm_short$cog <- sapply(ym_bnm__bnm_bnm_short$reg, assignCog, cog = cog_classes)


all_b_bis_shrt <- list(mix_bnm__ym_bnm_short,mix_bnm__bnm_bnm_short,ym_bnm__bnm_bnm_short)
len_b_bis_shrt <- lapply(all_b_bis_shrt, dim)

write.xlsx(mix_bnm__ym_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "mix_bnm__ym_bnm", append = TRUE)
write.xlsx(mix_bnm__bnm_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "mix_bnm__bnm_bnm", append = TRUE)
write.xlsx(ym_bnm__bnm_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "ym_bnm__bnm_bnm", append = TRUE)


b_all <- genEssTable(table, 5:7, 3)
b_all_shrt <- b_all[which(table$Annotation != "")]

common_ess_all_bnm <- data.frame(reg = table[which(b_all),1], gene = table[which(b_all),8], annot = table[which(b_all),9] )

write.xlsx(common_ess_all_bnm, file = "C:/Users/quent/Desktop/res_ess.xlsx", sheetName = "bnm_all_cond_ess_common", append = TRUE)

common_ess_all_bnm_short <- common_ess_all_bnm[which(common_ess_all_bnm$annot != ""),]

write.xlsx(common_ess_all_bnm_short, file = "C:/Users/quent/Desktop/res_ess_short.xlsx", sheetName = "bnm_all_cond_ess_common", append = TRUE)

( sum(b_all) + len_b[[1]][1] + len_b_bis[[1]][1] + len_b_bis[[2]][1] ) / len_ini_b[[1]][1]

summary_b <- data.frame(totMix = len_ini_b [[1]][1], totYM = len_ini_b [[2]][1], totBNM = len_ini_b [[3]][1], common_ess = sum(b_all), 
                        mix_ym = len_b_bis[[1]][1], mix_bnm = len_b_bis[[2]][1], ym_bnm = len_b_bis[[3]][1], 
                        mix = len_b[[1]][1], ym = len_b[[2]][1], bnm = len_b[[3]][1])

summary_b_shrt <- data.frame(totMix = len_ini_b_shrt [[1]][1], totYM = len_ini_b_shrt [[2]][1], totBNM = len_ini_b_shrt [[3]][1], 
                             common_ess = sum(b_all_shrt), mix_ym = len_b_bis_shrt[[1]][1], mix_bnm = len_b_bis_shrt[[2]][1], ym_bnm = len_b_bis_shrt[[3]][1], 
                             mix = len_b_shrt[[1]][1], ym = len_b_shrt[[2]][1], bnm = len_b_shrt[[3]][1])

####################################################################################################################################################

summary_tot <- rbind(summary_a, summary_b)

rownames(summary_tot) <- c("selection ym", "selection bnm")

write.xlsx(summary_tot, file = "C:/Users/quent/Desktop/summary_ess.xlsx")


summary_tot_shrt <- rbind(summary_a_shrt, summary_b_shrt)

rownames(summary_tot_shrt) <- c("selection ym", "selection bnm")

write.xlsx(summary_tot_shrt, file = "C:/Users/quent/Desktop/summary_ess.xlsx")


####################################################################################################################################################

require(ggplot2)

summary_tot <- cbind(select = rownames(summary_tot), summary_tot)
new_summary <- as.data.frame(gather(data = summary_tot, key = intrsct, value = cnt, totMix:bnm, factor_key = TRUE))

cnt_ess_tot <- ggplot(data = new_summary, aes(x = intrsct, y = cnt, fill = select)) +
    geom_bar(stat= "identity", position = position_dodge()) + 
    scale_x_discrete(name = "Libraries", limits = c("totMix", "totYM", "totBNM", "", "common_ess")) + 
    scale_y_continuous(name = "Number of essential genes") +
    scale_fill_brewer(palette ="Dark2") +
    geom_text(aes(label=cnt, y = 1), position = position_dodge(width = 0.9), vjust=-0.5, size=3) +
    theme(legend.position = "bottom")
cnt_ess_tot

cnt_ess_int <- ggplot(data = new_summary, aes(x = intrsct, y = cnt, fill = select)) +
    geom_bar(stat= "identity", position = position_dodge()) + 
    scale_x_discrete(name = "Intersection                                                Libraries", 
                     limits = c("mix_ym", "mix_bnm", "ym_bnm", "", "mix", "ym", "bnm")) + 
    scale_y_continuous(name = "Number of essential genes", limits = c(0,150)) +
    scale_fill_brewer(palette ="Dark2") +
    geom_text(aes(label=cnt, y = 1), position = position_dodge(width = 0.9), vjust=-0.5, size=3) +
    theme(legend.position = "bottom", panel.background = element_rect("grey92"))
cnt_ess_int



summary_tot_shrt <- cbind(select = rownames(summary_tot_shrt), summary_tot_shrt)
new_summary_shrt <- as.data.frame(gather(data = summary_tot_shrt, key = intrsct, value = cnt, totMix:bnm, factor_key = TRUE))

cnt_ess_tot_shrt <- ggplot(data = new_summary_shrt, aes(x = intrsct, y = cnt, fill = select)) +
    geom_bar(stat= "identity", position = position_dodge()) + 
    scale_x_discrete(name = "Libraries", limits = c("totMix", "totYM", "totBNM", "", "common_ess")) + 
    scale_y_continuous(name = "Number of essential genes") +
    scale_fill_brewer(palette ="Dark2") +
    geom_text(aes(label=cnt, y = 1), position = position_dodge(width = 0.9), vjust=-0.5, size=3) +
    theme(legend.position = "bottom")
cnt_ess_tot_shrt

cnt_ess_int_shrt <- ggplot(data = new_summary_shrt, aes(x = intrsct, y = cnt, fill = select)) +
    geom_bar(stat= "identity", position = position_dodge()) + 
    scale_x_discrete(name = "Intersection                                                Libraries", 
                     limits = c("mix_ym", "mix_bnm", "ym_bnm", "", "mix", "ym", "bnm")) + 
    scale_y_continuous(name = "Number of essential genes", limits = c(0,100)) +
    scale_fill_brewer(palette ="Dark2") +
    geom_text(aes(label=cnt, y = 1), position = position_dodge(width = 0.9), vjust=-0.5, size=3) +
    theme(legend.position = "bottom", panel.background = element_rect("grey92"))
cnt_ess_int_shrt

####################################################################################################################################################
## Pi chart generation for COG classes

lvl <-levels(cog_classes$cog)

a <- as.data.frame(table(ym_bnm__bnm_bnm_short$cog))

bp<- ggplot(a, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(width = 1, stat = "identity", color ="black") +
    coord_polar("y", start = 0) +
    geom_text(aes()) + ## Add coordinates based on central position of slices.
    theme(axis.text = element_blank())
bp

####################################################################################################################################################
## Correlation essentiality
#options(java.parameters = "-Xmx2048m")  
if(!require(openxlsx)) {
    install.packages("openxlsx")
    library(openxlsx)
}

corr <- read.xlsx("C:/Users/quent/Desktop/elArtist_data.xlsx", sheet = 3)
corr_head <- corr[1:5,]
corr <- corr[5:dim(corr)[1],]
corrr <- corr[,2:dim(corr)[2]]
corrr <- as.factor(corrr)

anova(lm(corr[,2] ~corr[,3]))

table(corr$"3")

test <- cor.test(as.numeric(corr$"16"), as.numeric(corr$"17"), method = "pearson")
