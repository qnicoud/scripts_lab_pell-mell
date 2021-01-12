rm(list = ls())

tab <- read.csv(file = "C:/Users/quent/Desktop/NCR_data_R.csv", header = TRUE, sep=";", stringsAsFactors = TRUE)

new_tab <- pivot_longer(data = tab, cols = c(rep1, rep2, rep3), names_to = "rep")

for (i in levels(new_tab$treatment)) {
        print(dunn.test::dunn.test(x = as.numeric(new_tab$value[which(new_tab$treatment == i)]), g= new_tab$strain[which(new_tab$treatment == i)]))
}
