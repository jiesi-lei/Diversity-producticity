# Fit logistic growth curve
# https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html#input-data
library(growthcurver)

# read data #
normal <- read.csv("normal.csv", header = TRUE, stringsAsFactors = FALSE)
salt <- read.csv("salt.csv", header = TRUE, stringsAsFactors = FALSE)
hungry <- read.csv("hungry.csv", header = TRUE, stringsAsFactors = FALSE)

normal.trim <- normal[1:21,]

gc_fit <- SummarizeGrowth(normal.$time,normal$X0101A1.4)
gc_fit
plot(gc_fit)
