rm(list= ls())

source("C:/Users/quent/Documents/GitHub/SAARA/rScript.R")

require(openxlsx)
require(ggplot2)
require(tidytext)
require(multcompView)


mass_data <- list()

mass_data$Ms <- read.xlsx("C:/Users/quent/Desktop/masses_Ms.xlsx")
names(mass_data$Ms) <- c("condition_name", "id", "empty", "full")
mass_data$Mt <- read.xlsx("C:/Users/quent/Desktop/masses_Mt.xlsx")
names(mass_data$Mt) <- c("condition_name", "id", "empty", "full")


calc <- list()
for (i in names(mass_data)){
    calc[[i]] <- cbind(mass_data[[i]], nmolC2H4_H_plant = mass_data[[i]]$full - mass_data[[i]]$empty)
    
    calc[[i]]$nmolC2H4_H_plant[which(calc[[i]]$nmolC2H4_H_plant < 0 | calc[[i]]$nmolC2H4_H_plant > 2)] <- NA
    
    boxplot(calc[[i]]$nmolC2H4_H_plant ~ calc[[i]]$cond)
}

a <- check_normality(calc)
b <- check_var_h(calc)
test <- check_means(calc, a, b)
result <- calc

ref <- list(
    Mt = matrix(c("NI", "WT", "bacA", "rpoH1", "lpsB", "lpxXL", "yjeA", "yjeE", "yjeF", 
                  "Controls", "Controls", "Controls", "Enveloppe functions", "Enveloppe functions", "Enveloppe functions", "yej transporter", "yej transporter", "yej transporter",
                  1, 2, 3, 1, 2, 3, 1, 2, 3), 
                ncol = 3),
    Ms =  matrix(c("NI", "WT", "bacA", "rpoH1", "lpsB", "lpxXL", "yejA", "yejE", "yejF", 
                   "Controls", "Controls", "Controls", "Enveloppe functions", "Enveloppe functions", "Enveloppe functions", "yej transporter", "yej transporter", "yej transporter",
                   1, 2, 3, 1, 2, 3, 1, 2, 3), 
                 ncol = 3))

fun <- function(result_ref) {
    x <- c()
    y <- c()
    
    for (i in 1:dim(result_ref$result)[1]) {
        x[i] <- result_ref$ref[which(result_ref$ref == result_ref$result$condition_name[i]),2]
        y[i] <- result_ref$ref[which(result_ref$ref == result_ref$result$condition_name[i]),3]
        
    }
    result_ref$result <- cbind(result_ref$result, class = x)
    result_ref$result <- cbind(result_ref$result, ord = as.numeric(y))

    return(result_ref)
}

resultTest <- lapply(list( Mt = list(result = result$Mt, ref = ref$Mt), Ms = list(result = result$Ms, ref = ref$Ms)), fun)
result$Mt <- resultTest$Mt$result
result$Ms <- resultTest$Ms$result


get_high_values <- function(x, lim) {
    y <- as.data.frame(table(x[which(x$nmolC2H4_H_plant > lim),"condition_name"]))
    if (dim(y)[1] != 0) {
        colnames(y) <- c("condition_name", "nmolC2H4_H_plant")
        y$nmolC2H4_H_plant <- gsub( "[0-9]", 0.5, y$nmolC2H4_H_plant)
    }
    else
        y <- NULL
    
    return(y)
}

high <- lapply(result, get_high_values, lim_max)
high_ok <- lapply(list( Mt = list(result = high$Mt, ref = ref$Mt), Ms = list(result = high$Ms, ref = ref$Ms)), fun)
high_ok <- lapply(list( Ms = list(result = high$Ms, ref = ref$Ms)),  fun)

trim_data_graphs <- function(data, lim) {
    data <- data[!is.na(data$nmolC2H4_H_plant),]
    data <- data[!(data$nmolC2H4_H_plant >= lim),]
    
    return(data)
}

result <- lapply(result, trim_data_graphs, lim = lim_max)


# Temp ####

get_stat_lab <- function(res) 
{
    names(res$stat$P.adjusted) <- gsub(" ", "", res$stat$comparisons)
    assigned_class <- multcompLetters(res$stat$P.adjusted)["Letters"]
    condition_name <- names(assigned_class[["Letters"]])
    
    boxplot.df <- ddply(res$res, .(condition_name), function(x) max(fivenum(x$nmolC2H4_H_plant)+0.2*(x$nmolC2H4_H_plant), na.rm = TRUE))
    #boxplot.df <- boxplot.df[, c("condition_name", "V1")]
    plot.levels <- data.frame(condition_name, labels = assigned_class[['Letters']],
                              stringsAsFactors = FALSE)
    
    labels.df <- merge(plot.levels, boxplot.df, by= "condition_name", sort = FALSE)
    
    return(labels.df)
}

lim_max <- 0.5

labels_stat <- lapply(list(Ms = list( stat = test$Ms, res = result$Ms), Mt = list(stat = test$Mt, res = result$Mt)), get_stat_lab)
####
label_ok <- lapply(list( Mt = list(result = labels_stat$Mt, ref = ref$Mt), Ms = list(result = labels_stat$Ms, ref = ref$Ms)), fun)
label_ok <- list(Ms = label_ok$Ms$result, Mt = label_ok$Mt$result)

plots <- list( Mt =
                   ggplot(result$Mt, aes(x = reorder_within(condition_name, ord, class), y = nmolC2H4_H_plant, fill = condition_name)) +
                   geom_boxplot(size = 1.5, outlier.colour="red", outlier.shape=8, outlier.size=4) +
                   facet_grid(~class, scales = "free") + 
                   scale_x_reordered() +
                   scale_fill_manual(values = c("black", "#4CB4BE", "#93aa00", "grey", "#619cff", "white", "#FFBA00", "#FF9223", "#FF5123")) +
                   labs(x = "", y = "Plant dry mass (mg)") +
                   theme(legend.position = "none",
                         strip.text = element_text(size = 24, face="bold"),
                         axis.title.y = element_text(size = 22),
                         axis.text.y = element_text(size = 20)) + 
                   geom_text(data = label_ok$Mt, aes(x = reorder_within(condition_name, ord, class), y = V1, label = labels),size=8)
               , Ms = 
                   ggplot(result$Ms, aes(x = reorder_within(condition_name, ord, class), y = nmolC2H4_H_plant, fill = condition_name)) + 
                   geom_boxplot(size = 1.5, outlier.colour="red", outlier.shape=8, outlier.size=4) +
                   facet_grid(~class, scales = "free") + 
                   scale_x_reordered() +
                   scale_fill_manual(values = c("black", "#4CB4BE", "#93aa00", "grey", "#619cff", "white", "#FFBA00", "#FF9223", "#FF5123")) +
                   labs(x = "", y = "Plant dry mass (mg)") +
                   theme(legend.position = "none",
                         strip.text = element_text(size = 24, face="bold"),
                         axis.title.y = element_text(size = 22),
                         axis.text.y = element_text(size = 20)) + 
                   geom_text(data = label_ok$Ms, aes(x = reorder_within(condition_name, ord, class), y = V1, label = labels),size=8)
)

plots

save_bmp(plots, awidth = 1100, aheight = 500)
