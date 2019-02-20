# code to load and clean buds data

# load raw data
raw_data <- read.csv("./data/compiled-bud-depths.csv", stringsAsFactors = FALSE)

# load p_resprout data
resprout <- read.csv("data/survival-data.csv", stringsAsFactors = FALSE)

# remove c2 treatment
raw_data <- raw_data[raw_data$treat %in% c("c", "b"), ]

# load trait data
trait_data <- read.csv("./data/buds-traits-data.csv")

# identify columns that contain depth measurements
data_cols <- c(6:ncol(raw_data))

# convert depth measurements to a single vector
depth_long <- unlist(raw_data[, data_cols])

# add treatment info to depth measurements
data_long <- data.frame(do.call(rbind, lapply(seq_len(length(data_cols)),
                                              function(x) raw_data[, 1:5])),
                        depths_mm = depth_long)
data_long <- data_long[-which(is.na(data_long$depths_mm)), ]

# sort data by species
data_long <- data_long[order(data_long$sp), ]

# create species-specific IDs for individual plants
data_long$plant_code <- as.integer(as.factor(paste0(data_long$spp, data_long$plant, data_long$yr)))

# prepare binned data
breaks_set <- seq(-15, 95, length = 21)
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

# fill NAs in trait data with mean values
# for (i in which(apply(trait_data, 2, function(x) any(is.na(x)))))
#   trait_data[which(is.na(trait_data[, i])), i] <- mean(trait_data[, i], na.rm = TRUE)

# alternative: remove any species with missing trait data
trait_data <- trait_data[!apply(trait_data, 1, function(x) any(is.na(x))), ]

# standardise trait data
trait_data$MAXHTstd <- scale(trait_data$MAXHT)
trait_data$STEMCATSstd <- scale(trait_data$STEMCATS)
trait_data$SLAstd <- scale(trait_data$SLA)
trait_data$LWCstd <- scale(trait_data$LWC)
trait_data$BAstd <- scale(trait_data$BA)
trait_data$CAstd <- scale(trait_data$CA)

# add trait data to data_info
trait_rows <- match(data_info$SPP, trait_data$CODE)
data_info$GFORM <- trait_data$GFORM[trait_rows]
data_info$MAXHTstd <- trait_data$MAXHTstd[trait_rows]
data_info$STEMCATSstd <- trait_data$STEMCATSstd[trait_rows]
data_info$SLAstd <- trait_data$SLAstd[trait_rows]
data_info$LWCstd <- trait_data$LWCstd[trait_rows]
data_info$BAstd <- trait_data$BAstd[trait_rows]
data_info$CAstd <- trait_data$CAstd[trait_rows]

# add trait data to resprout data
resprout <- data.frame(resprout,
                       trait_data[match(resprout$spp, trait_data$CODE), ])

# remove NAs from data_hist, data_info and resprout if they exist
if (any(is.na(data_info))) {
  data_hist <- data_hist[!apply(data_info, 1, function(x) any(is.na(x))), ]
  data_info <- data_info[!apply(data_info, 1, function(x) any(is.na(x))), ]
}
if (any(is.na(resprout)))
  resprout <- resprout[!apply(resprout, 1, function(x) any(is.na(x))), ]
