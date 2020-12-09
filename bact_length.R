rm(list = ls())

bact_length <- list()
bact_length$Area <- openxlsx::read.xlsx("C:/Users/quent/Desktop/usda/All_Stats_r.xlsx", "Area")
bact_length$Length <- openxlsx::read.xlsx("C:/Users/quent/Desktop/usda/All_Stats_r.xlsx", "Length")
bact_length$Width <- openxlsx::read.xlsx("C:/Users/quent/Desktop/usda/All_Stats_r.xlsx", "Width")

value <- lapply(bact_length, stack)
levels(value$Area$ind)

dunn.test::dunn.test(x = value[[3]]$values, g = as.factor(value[[3]]$ind), alpha = 0.001, method = "by")

# linear model for anova
anova_model <- lm( unlist(value[[2]]['values']) ~ as.factor(unlist(value[[2]]['ind'])))

# Graphical Verification of residues normality
plot(anova_model)
car::qqPlot(anova_model, simulate = TRUE, id = "y", id.n = 2, main = "Q-Q plot with confidence enveloppe")
dev.off()

# Verify residues normality
res <- anova_model$residuals
resNorm <- numeric()
isResNorm <- logical()
condNames <- levels(value[[1]]$ind)

for (j in 1:length(condNames))
{
    resNorm[j] <- shapiro.test(res[which(value[[1]]$ind == condNames[j])])$p.value
    if (resNorm[j] >= normalityThreshold)
    {
        isResNormal[j] <- TRUE          
    }
}

# AnOVa
car::Anova(anova_model, type = 2)

tukey <- TukeyHSD(aov(value[[2]]$values ~ value[[2]]$ind, data = value[[2]] ), conf.level = 0.95)
print(plot(tukey))


require(stats)
# kolgorow smirnoff
ks.test()
# anderson daring



leng_wid <- data.frame(ind = as.character(values$Length$ind), width = values$Width$values, length = values$Length$values)

require(ggplot2)

ggplot(value$Area, aes(x = ind , y = values, fill = ind)) +
    geom_boxplot(size = 1.2) + 
    scale_fill_manual(values = c("#F8AC9C", "#F8766D", "#90ee90", "#00BA38", "#BF9000")) +
    theme(legend.position = "none")
ggplot(value$Area, aes(x = ind , y = values, fill = ind)) +
    geom_violin(trim = F, size = 1, width = 1.2) +
    geom_boxplot(width = 0.05, size = 1, outlier.shape = NA) +
    scale_fill_manual(values = c("#F8AC9C", "#F8766D", "#90ee90", "#00BA38", "#BF9000"))


ggplot(value$Length, aes(x = ind , y = values, fill = ind)) +
    geom_boxplot(size = 1.2) + 
    scale_fill_manual(values = c("#F8AC9C", "#F8766D", "#90ee90", "#00BA38", "#BF9000")) +
    theme(legend.position = "none")
ggplot(value$Length, aes(x = ind , y = values, fill = ind)) +
    geom_violin(trim = F, size = 1, width = 1.2) +
    geom_boxplot(width = 0.05, size = 1, outlier.shape = NA) +
    scale_fill_manual(values = c("#F8AC9C", "#F8766D", "#90ee90", "#00BA38", "#BF9000"))

ggplot(value$Width, aes(x = ind , y = values, fill = ind)) +
    geom_boxplot(size = 1.2) + 
    scale_fill_manual(values = c("#F8AC9C", "#F8766D", "#90ee90", "#00BA38", "#BF9000"))
ggplot(value$Width, aes(x = ind , y = values, fill = ind)) +
    geom_violin(trim = F, size = 1, width = 1.2) +
    geom_boxplot(width = 0.05, size = 1, outlier.shape = NA) +
    scale_fill_manual(values = c("#F8AC9C", "#F8766D", "#90ee90", "#00BA38", "#BF9000"))

ggplot(leng_wid, aes(x = Width, y = length, fill = ind)) +
    geom_point()
    