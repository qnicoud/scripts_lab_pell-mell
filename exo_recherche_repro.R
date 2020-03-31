## Exercices prérequis formation INRAE "Recherche reproductible"
## Quentin Nicoud
## Date de création: 26.02.2019


## Exo 1
data("ChickWeight")

time_8 <- ChickWeight[ChickWeight$Time == 8,]
time_8$Diet <- as.factor(time_8$Diet)

model <- lm(formula = weight ~ Diet, data = time_8)
    coef(model)
    fitted(model)
    residuals(model)
    predict(model, interval="confidence")
    
anova(model)


## Exo 2
install.packages("nycflights13")
require(nycflights13)
data("flights")
carriers <- table(flights$carrier)
sub_carriers <- carriers[carriers > 10000]
sub_flights <- subset(flights, flights$carrier %in% names(sub_carriers))

arrival_delays <- sub_carriers
for ( i in 1:length(sub_carriers)) {
    arrival_delays[i] <- sum(flights$arr_delay[which(flights$carrier == names(sub_carriers[i]))], na.rm = TRUE)/sub_carriers[i]
}

require(tibble)
arrival_delays_per_month <- tibble(carrier="a", month=0, avDelayPerMonth=0)
cnt <- 0
for ( i in 1:length(sub_carriers) ) {
    for ( j in 1:12 ) {
        cnt <- cnt + 1
        arrival_delays_per_month[cnt,1] <- names(sub_carriers[i])
        arrival_delays_per_month[cnt,2] <- j
        a <- flights[which(flights$carrier == names(sub_carriers[i])),]
        arrival_delays_per_month[cnt,3] <- sum(a$arr_delay[which(a$month == j)], na.rm = TRUE)/(sum(a$month[which(a$month == j)], na.rm = TRUE)/j)
    }
}

require(ggplot2)

p <- ggplot(data = arrival_delays_per_month, mapping = aes(x = month, y = avDelayPerMonth, color = carrier)) + geom_line() + theme_gray()


## Exo 3
syracuse <- function(a) {
    res <- numeric(0)
    res[1] <- a
    for (i in 1:2000) {
        if (res[i] == 1) {
            break
        } else if ((res[i] %% 2) == 0) {
            res[i+1] <- res[i]/2
        } else if ((res[i] %% 2) == 1) {
            res[i+1] <- res[i] * 3 + 1
        } 
    }
    return(res)
}

syr_res_list <- lapply(1:1000, syracuse)
sum(unlist(lapply(syr_res_list, min)))
