# code to load and clean buds data

# load raw data
raw_data <- read.csv("./data/raw/compiled buddepths_for_Jian.csv",
                     stringsAsFactors = FALSE)

# identify columns that contain depth measurements
data_cols <- c(6:ncol(raw_data))

# convert depth measurements to a single vector
depth_long <- unlist(raw_data[, data_cols])

# add treatment info to depth measurements
data_long <- data.frame(do.call("rbind", lapply(seq_len(length(data_cols)),
                                                function(x) raw_data[, 1:5])),
                        depths_mm = depth_long)
data_long <- data_long[-which(is.na(data_long$depths_mm)), ]

# prepare binned data

## CLEAN UP plant + spp columns to give an ID and SPCODE

hist_fun <- function(x) hist(x, plot = FALSE)[c("counts", "mids")]
data_hist <- tapply(data_long$depths_mm,
                    list(data_long$plant, data_long$spp,
                         data_long$yr, data_long$treat),
                    hist_fun)

