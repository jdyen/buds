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

# sort data by species
data_long <- data_long[order(data_long$sp), ]

# create species-specific IDs for individual plants
data_long$plant_code <- as.integer(as.factor(paste0(data_long$spp, data_long$plant, data_long$yr)))

# prepare binned data
breaks_set <- seq(-25, 95, by = 1)
hist_fun <- function(x) hist(x,
                             breaks = breaks_set,
                             plot = FALSE)[c("counts", "mids")]
data_hist <- matrix(NA,
                    nrow = length(unique(data_long$plant_code)),
                    ncol = (length(breaks_set) - 1))
data_info <- matrix(NA, nrow = nrow(data_hist), ncol = 4)
for (i in seq_along(unique(data_long$plant_code))) {

  data_subset <- which(data_long$plant_code == unique(data_long$plant_code)[i])
  hist_tmp <- hist_fun(data_long$depths_mm[data_subset])
  data_hist[i, ] <- hist_tmp$counts
  data_info[i, ] <- c(unique(data_long$plant_code[data_subset]),
                      unique(data_long$treat[data_subset]),
                      unique(data_long$yr[data_subset]) + 1,
                      unique(data_long$spp[data_subset]))
  
}
data_mids <- hist_tmp$mids
data_info <- as.data.frame(data_info,
                           stringsAsFactors = FALSE)
colnames(data_info) <- c("CODE", "TREATMENT", "YEAR", "SPP")
data_info$CODE <- as.integer(data_info$CODE)
data_info$YEAR <- as.integer(data_info$YEAR)
