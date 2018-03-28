# R code to fit bud depths model

# load packages
library(FREE)

# load data
source("./code/load-data.R")

# prepare data set
data_set <- list(y = data_hist,
                 x = factor(data_info$TREATMENT),
                 z = cbind(as.integer(data_info$YEAR),
                           as.integer(as.factor(data_info$SPP))))

# fit model
mod <- FREEfit(y ~ x, groups = data_set$z, data = data_set,
               n.iters = 1000, n.burnin = 500, n.chains = 2)

# cross validate fitted model
# mod_cv <- FREEfitCV(y ~ x, groups = data$z, data = data_set,
#                     n.cv = 10, n.iters = 10, n.burnin = 5,
#                     n.chains = 1)

# summarise fitted model

