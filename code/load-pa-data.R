# R code for loading and preparing all GB data for beta diversity analyses

# load survey data
cent.bird.survey <- read.csv("./data/cent-birds-survey-data.csv", header = TRUE)
west.bird.survey <- read.csv("./data/west-birds-survey-data.csv", header = TRUE)
cent.bfs.survey <- read.csv("./data/cent-bfs-survey-data.csv", header = TRUE)
west.bfs.survey <- read.csv("./data/west-bfs-survey-data.csv", header = TRUE)

# list of spp to remove
exclusions <- read.csv("./data/spp_exclusions.csv")
exclusions$Common.name <- gsub(pattern = "Black-crowned Night-heron",
                               replacement = "Black-crowned Night-Heron",
                               exclusions$Common.name)
exclusions$Common.name <- gsub(pattern = "Northern Saw-whet Owl",
                               replacement = "Northern Saw-Whet Owl",
                               exclusions$Common.name)
exclusions$Common.name <- gsub(pattern = "Sagebrush Sparrow",
                               replacement = "Sage Sparrow",
                               exclusions$Common.name)
exclusions_list <- exclusions$Common.name[which(exclusions$Poorly.sampled.by.point.counts == "Y")]

# rename spp with misformatted commas
cent.bird.survey$sppname <- gsub(pattern = "\xd5", "\\'", cent.bird.survey$sppname)
west.bird.survey$sppname <- gsub(pattern = "\xd5", "\\'", west.bird.survey$sppname)
west.bird.survey$sppname <- gsub(pattern = "Grossbeak", "Grosbeak'", west.bird.survey$sppname)

# exclude those species
cent.bird.survey <- cent.bird.survey[which(is.na(match(cent.bird.survey$sppname, exclusions_list))), ]
west.bird.survey <- west.bird.survey[which(is.na(match(west.bird.survey$sppname, exclusions_list))), ]

# load covariate data
cent.bird.cov.pts <- read.csv("./data/cent-birds-cov-pt.csv", header = TRUE)
cent.bird.cov.can <- read.csv("./data/cent-birds-cov-can.csv", header = TRUE)
west.bird.cov.pts <- read.csv("./data/west-birds-cov-pt.csv", header = TRUE)
west.bird.cov.can <- read.csv("./data/west-birds-cov-can.csv", header = TRUE)
cent.bfs.cov.pts <- read.csv("./data/cent-bfs-cov-tran.csv", header = TRUE)
cent.bfs.cov.can <- read.csv("./data/cent-bfs-cov-can.csv", header = TRUE)
west.bfs.cov.pts <- read.csv("./data/west-bfs-cov-tran.csv", header = TRUE)
west.bfs.cov.can <- read.csv("./data/west-bfs-cov-can.csv", header = TRUE)

# load guild data
bird.guilds <- read.csv("./data/bird-guilds.csv", row.names = 1)
bfs.guilds <- read.csv("./data/bfs-guilds.csv", row.names = 1)

# compile western and central bird surveys into one file
bird.survs <- rbind(cent.bird.survey, west.bird.survey)

# clean up inconsistent bird species names
bird.survs$sppname <- gsub("'", "", bird.survs$sppname)
bird.survs$sppname <- tolower(bird.survs$sppname)
bird.survs$sppname <- gsub("-", "", bird.survs$sppname)
bird.survs$sppname <- gsub(" ", "", bird.survs$sppname)
bird.survs$sppname <- gsub("lewisswoodpecker", "lewiswoodpecker", bird.survs$sppname)

# tidy spp names from guild data
bird.guilds$sppname <- gsub("'", "", bird.guilds$Bird.species)
bird.guilds$sppname <- gsub("\xd5", "", bird.guilds$sppname)
bird.guilds$sppname <- tolower(bird.guilds$sppname)
bird.guilds$sppname <- gsub("-", "", bird.guilds$sppname)
bird.guilds$sppname <- gsub(" ", "", bird.guilds$sppname)
bird.guilds$sppname <- gsub("lewisswoodpecker", "lewiswoodpecker", bird.guilds$sppname)
bird.guilds$sppname <- gsub("sagebrushsparrow", "sagesparrow", bird.guilds$sppname)

# exclude spp from bird guilds
bird.guilds <- bird.guilds[-which(is.na(match(bird.guilds$sppname, bird.survs$sppname))), ]

# exlucde spp repeats in bird guilds
bird.guilds <- bird.guilds[-c(8, 13, 15, 22, 27, 53, 87), ]

# add sppcodes for butterflies
bfs.splist <- rbind(cent.bfs.survey[, 1:10], west.bfs.survey[, 1:10])
bfs.guilds$sppcode <- bfs.splist$sppcode[match(bfs.guilds$Butterfly.species, bfs.splist$sppname)]

# compile annual survey data from detection lists for birds
bird.annual.survs <- tapply(rep(1, nrow(bird.survs)),
                            list(bird.survs$sppname, bird.survs$point, bird.survs$year),
                          sum,
                            na.rm = TRUE)
all.bird.surv.combos <- unique(paste(bird.survs$year, bird.survs$point, sep = "-"))
bird.annual.survs <- ifelse(is.na(bird.annual.survs), 0, bird.annual.survs)
bird.annual.can.survs <- tapply(rep(1, nrow(bird.survs)),
                                list(bird.survs$sppname, bird.survs$canyon, bird.survs$year),
                                sum,
                                na.rm = TRUE)
all.bird.can.combos <- unique(paste(bird.survs$year, bird.survs$canyon, sep = "-"))
bird.annual.can.survs <- ifelse(is.na(bird.annual.can.survs), 0, bird.annual.can.survs)

# convert non-surveys to NAs rather than zeros
bird.year.list <- dimnames(bird.annual.survs)[[3]]
for (i in seq(along = bird.year.list)) {
  bird.site.list <- sapply(strsplit(all.bird.surv.combos[grep(bird.year.list[i], all.bird.surv.combos)], "-"), function(x) x[2]) 
  bird.annual.survs[, which(is.na(match(colnames(bird.annual.survs), bird.site.list))), i] <- rep(NA, nrow(bird.annual.survs))
}
bird.can.year.list <- dimnames(bird.annual.can.survs)[[3]]
for (i in seq(along = bird.can.year.list)) {
  bird.can.list <- sapply(strsplit(all.bird.can.combos[grep(bird.can.year.list[i], all.bird.can.combos)], "-"), function(x) x[2]) 
  bird.annual.can.survs[, which(is.na(match(colnames(bird.annual.can.survs), bird.can.list))), i] <- rep(NA, nrow(bird.annual.can.survs))
}

# compile annual survey data from detection lists for butterflies
bfs.survs <- rbind(cent.bfs.survey[, 1:11], west.bfs.survey[, 1:11])
bfs.annual.survs <- tapply(rep(1, nrow(bfs.survs)),
                           list(bfs.survs$sppcode, bfs.survs$segment, bfs.survs$year),
                           sum,
                           na.rm = TRUE)
all.bfs.surv.combos <- unique(paste(bfs.survs$year, bfs.survs$segment, sep = "-"))
bfs.annual.survs <- ifelse(is.na(bfs.annual.survs), 0, bfs.annual.survs)
bfs.annual.can.survs <- tapply(rep(1, nrow(bfs.survs)),
                                list(bfs.survs$sppcode, bfs.survs$canyon, bfs.survs$year),
                                sum,
                                na.rm = TRUE)
all.bfs.can.combos <- unique(paste(bfs.survs$year, bfs.survs$canyon, sep = "-"))
bfs.annual.can.survs <- ifelse(is.na(bfs.annual.can.survs), 0, bfs.annual.can.survs)

# convert non-surveys to NAs rather than zeros
bfs.year.list <- dimnames(bfs.annual.survs)[[3]]
for (i in seq(along = bfs.year.list)) {
  bfs.site.list <- sapply(strsplit(all.bfs.surv.combos[grep(bfs.year.list[i], all.bfs.surv.combos)], "-"), function(x) x[2]) 
  bfs.annual.survs[, which(is.na(match(colnames(bfs.annual.survs), bfs.site.list))), i] <- rep(NA, nrow(bfs.annual.survs))
}
bfs.can.year.list <- dimnames(bfs.annual.can.survs)[[3]]
for (i in seq(along = bfs.can.year.list)) {
  bfs.can.list <- sapply(strsplit(all.bfs.can.combos[grep(bfs.can.year.list[i], all.bfs.can.combos)], "-"), function(x) x[2]) 
  bfs.annual.can.survs[, which(is.na(match(colnames(bfs.annual.can.survs), bfs.can.list))), i] <- rep(NA, nrow(bfs.annual.can.survs))
}

# convert arrays to matrices for all bird point-level survey data
bird.ann.pt.data <- matrix(NA,
                           nrow = (dim(bird.annual.survs)[2] * dim(bird.annual.survs)[3]),
                           ncol = dim(bird.annual.survs)[1])
bird.ann.pt.info <- data.frame(year = rep(dimnames(bird.annual.survs)[[3]], each = dim(bird.annual.survs)[2]),
                               point = rep(dimnames(bird.annual.survs)[[2]], times = dim(bird.annual.survs)[3]))
for (i in 1:dim(bird.annual.survs)[3]) {
  bird.ann.pt.data[(((i - 1) * dim(bird.annual.survs)[2]) + 1):(i * dim(bird.annual.survs)[2]), ] <- t(bird.annual.survs[, , i])
}
colnames(bird.ann.pt.data) <- dimnames(bird.annual.survs)[[1]]

# convert arrays to matrices for bird canyon-level data
bird.ann.can.data <- matrix(NA,
                            nrow = (dim(bird.annual.can.survs)[2] * dim(bird.annual.can.survs)[3]),
                            ncol = dim(bird.annual.can.survs)[1])
bird.ann.can.info <- data.frame(year = rep(dimnames(bird.annual.can.survs)[[3]], each = dim(bird.annual.can.survs)[2]),
                                canyon = rep(dimnames(bird.annual.can.survs)[[2]], times = dim(bird.annual.can.survs)[3]))
for (i in 1:dim(bird.annual.can.survs)[3]) {
  bird.ann.can.data[(((i - 1) * dim(bird.annual.can.survs)[2]) + 1):(i * dim(bird.annual.can.survs)[2]), ] <- t(bird.annual.can.survs[, , i])
}
colnames(bird.ann.can.data) <- dimnames(bird.annual.can.survs)[[1]]

# remove bird non-surveys from data set
bird.na.survs <- which(apply(bird.ann.pt.data, 1, function(x) all(is.na(x))))
bird.ann.pt.data <- bird.ann.pt.data[-bird.na.survs, ]
bird.ann.pt.info <- bird.ann.pt.info[-bird.na.survs, ]
bird.ann.pt.data <- ifelse(bird.ann.pt.data > 0, 1, 0)
bird.na.can.survs <- which(apply(bird.ann.can.data, 1, function(x) all(is.na(x))))
bird.ann.can.data <- bird.ann.can.data[-bird.na.can.survs, ]
bird.ann.can.info <- bird.ann.can.info[-bird.na.can.survs, ]
bird.ann.can.data <- ifelse(bird.ann.can.data > 0, 1, 0)

# convert arrays to matrices for all butterfly transect-level survey data
bfs.ann.tran.data <- matrix(NA,
                            nrow = (dim(bfs.annual.survs)[2] * dim(bfs.annual.survs)[3]),
                            ncol = dim(bfs.annual.survs)[1])
bfs.ann.tran.info <- data.frame(year = rep(dimnames(bfs.annual.survs)[[3]], each = dim(bfs.annual.survs)[2]),
                                transect = rep(dimnames(bfs.annual.survs)[[2]], times = dim(bfs.annual.survs)[3]))
for (i in 1:dim(bfs.annual.survs)[3]) {
  bfs.ann.tran.data[(((i - 1) * dim(bfs.annual.survs)[2]) + 1):(i * dim(bfs.annual.survs)[2]), ] <- t(bfs.annual.survs[, , i])
}
colnames(bfs.ann.tran.data) <- dimnames(bfs.annual.survs)[[1]]

# convert arrays to matrices for butterfly canyon-level data
bfs.ann.can.data <- matrix(NA,
                           nrow = (dim(bfs.annual.can.survs)[2] * dim(bfs.annual.can.survs)[3]),
                           ncol = dim(bfs.annual.can.survs)[1])
bfs.ann.can.info <- data.frame(year = rep(dimnames(bfs.annual.can.survs)[[3]], each = dim(bfs.annual.can.survs)[2]),
                               canyon = rep(dimnames(bfs.annual.can.survs)[[2]], times = dim(bfs.annual.can.survs)[3]))
for (i in 1:dim(bfs.annual.can.survs)[3]) {
  bfs.ann.can.data[(((i - 1) * dim(bfs.annual.can.survs)[2]) + 1):(i * dim(bfs.annual.can.survs)[2]), ] <- t(bfs.annual.can.survs[, , i])
}
colnames(bfs.ann.can.data) <- dimnames(bfs.annual.can.survs)[[1]]

# remove butterfly non-surveys from data set
bfs.na.survs <- which(apply(bfs.ann.tran.data, 1, function(x) all(is.na(x))))
bfs.ann.tran.data <- bfs.ann.tran.data[-bfs.na.survs, ]
bfs.ann.tran.info <- bfs.ann.tran.info[-bfs.na.survs, ]
bfs.ann.tran.data <- ifelse(bfs.ann.tran.data > 0, 1, 0)
bfs.na.can.survs <- which(apply(bfs.ann.can.data, 1, function(x) all(is.na(x))))
bfs.ann.can.data <- bfs.ann.can.data[-bfs.na.can.survs, ]
bfs.ann.can.info <- bfs.ann.can.info[-bfs.na.can.survs, ]
bfs.ann.can.data <- ifelse(bfs.ann.can.data > 0, 1, 0)

# compile environmental data for birds at points with appropriate years
bird.cov.pts <- rbind(cent.bird.cov.pts[, 1:10], west.bird.cov.pts[, 1:10])
bird.cov.all <- cbind(bird.ann.pt.info, bird.cov.pts[match(bird.ann.pt.info$point, bird.cov.pts$point), ])
bird.years <- unique(bird.cov.all$year)
bird.cov.all$precip <- rep(NA, nrow(bird.cov.all))
bird.cov.all$maxtemp <- rep(NA, nrow(bird.cov.all))
bird.cov.all$mean_maxtemp <- rep(NA, nrow(bird.cov.all))
bird.cov.all$mintemp <- rep(NA, nrow(bird.cov.all))
bird.cov.all$mean_mintemp <- rep(NA, nrow(bird.cov.all))
for (i in seq(along = bird.years)) {
  year.set <- which(bird.cov.all$year == bird.years[i])
  cent.subset.rows <- match(cent.bird.cov.pts$point,
                            bird.cov.all$point[year.set])
  west.subset.rows <- match(west.bird.cov.pts$point,
                            bird.cov.all$point[year.set])
  if (all(is.na(west.subset.rows))) {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bird.cov.all$precip[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                    grep(paste0(bird.years[i], "_total_precip"),
                                                                         colnames(cent.bird.cov.pts))]
    bird.cov.all$maxtemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                     grep(paste0(bird.years[i], "_maxtemp"),
                                                                          colnames(cent.bird.cov.pts))]
    bird.cov.all$mintemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                     grep(paste0(bird.years[i], "_mintemp"),
                                                                          colnames(cent.bird.cov.pts))]
    bird.cov.all$mean_maxtemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                          grep(paste0(bird.years[i], "_mean_maxtemp"),
                                                                               colnames(cent.bird.cov.pts))]
    bird.cov.all$mean_mintemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                          grep(paste0(bird.years[i], "_mean_mintemp"),
                                                                               colnames(cent.bird.cov.pts))]
  } else {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bird.cov.all$precip[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                    grep(paste0(bird.years[i], "_total_precip"),
                                                                         colnames(cent.bird.cov.pts))]
    bird.cov.all$maxtemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                     grep(paste0(bird.years[i], "_maxtemp"),
                                                                          colnames(cent.bird.cov.pts))]
    bird.cov.all$mintemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                     grep(paste0(bird.years[i], "_mintemp"),
                                                                          colnames(cent.bird.cov.pts))]
    bird.cov.all$mean_maxtemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                          grep(paste0(bird.years[i], "_mean_maxtemp"),
                                                                               colnames(cent.bird.cov.pts))]
    bird.cov.all$mean_mintemp[year.set][cov.row.set] <- cent.bird.cov.pts[which(!is.na(cent.subset.rows)),
                                                                          grep(paste0(bird.years[i], "_mean_mintemp"),
                                                                               colnames(cent.bird.cov.pts))]
    west.cov.row.set <- west.subset.rows[which(!is.na(west.subset.rows))]
    bird.cov.all$precip[year.set][west.cov.row.set] <- west.bird.cov.pts[which(!is.na(west.subset.rows)),
                                                                         grep(paste0(bird.years[i], "_total_precip"),
                                                                              colnames(west.bird.cov.pts))]
    bird.cov.all$maxtemp[year.set][west.cov.row.set] <- west.bird.cov.pts[which(!is.na(west.subset.rows)),
                                                                          grep(paste0(bird.years[i], "_maxtemp"),
                                                                               colnames(west.bird.cov.pts))]
    bird.cov.all$mintemp[year.set][west.cov.row.set] <- west.bird.cov.pts[which(!is.na(west.subset.rows)),
                                                                          grep(paste0(bird.years[i], "_mintemp"),
                                                                               colnames(west.bird.cov.pts))]
    bird.cov.all$mean_maxtemp[year.set][west.cov.row.set] <- west.bird.cov.pts[which(!is.na(west.subset.rows)),
                                                                               grep(paste0(bird.years[i], "_mean_maxtemp"),
                                                                                    colnames(west.bird.cov.pts))]
    bird.cov.all$mean_mintemp[year.set][west.cov.row.set] <- west.bird.cov.pts[which(!is.na(west.subset.rows)),
                                                                               grep(paste0(bird.years[i], "_mean_mintemp"),
                                                                                    colnames(west.bird.cov.pts))]
  }
}

# compile environmental data for birds in canyons with appropriate years
bird.cov.can <- rbind(cent.bird.cov.can[, 1:14], west.bird.cov.can[, 1:14])
bird.cov.can.all <- cbind(bird.ann.can.info, bird.cov.pts[match(bird.ann.can.info$canyon, bird.cov.can$canyon), ])
bird.years.can <- unique(bird.cov.can.all$year)
bird.cov.can.all$precip <- rep(NA, nrow(bird.cov.can.all))
bird.cov.can.all$maxtemp <- rep(NA, nrow(bird.cov.can.all))
bird.cov.can.all$mean_maxtemp <- rep(NA, nrow(bird.cov.can.all))
bird.cov.can.all$mintemp <- rep(NA, nrow(bird.cov.can.all))
bird.cov.can.all$mean_mintemp <- rep(NA, nrow(bird.cov.can.all))
for (i in seq(along = bird.years.can)) {
  year.set <- which(bird.cov.can.all$year == bird.years[i])
  cent.subset.rows <- match(cent.bird.cov.can$canyon,
                            bird.cov.can.all$canyon[year.set])
  west.subset.rows <- match(west.bird.cov.can$canyon,
                            bird.cov.can.all$canyon[year.set])
  if (all(is.na(west.subset.rows))) {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bird.cov.can.all$precip[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                        grep(paste0(bird.years[i], "_total_precip"),
                                                                             colnames(cent.bird.cov.can))]
    bird.cov.can.all$maxtemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                         grep(paste0(bird.years[i], "_maxtemp"),
                                                                              colnames(cent.bird.cov.can))]
    bird.cov.can.all$mintemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                         grep(paste0(bird.years[i], "_mintemp"),
                                                                              colnames(cent.bird.cov.can))]
    bird.cov.can.all$mean_maxtemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                              grep(paste0(bird.years[i], "_mean_maxtemp"),
                                                                                   colnames(cent.bird.cov.can))]
    bird.cov.can.all$mean_mintemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                              grep(paste0(bird.years[i], "_mean_mintemp"),
                                                                                   colnames(cent.bird.cov.can))]
  } else {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bird.cov.can.all$precip[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                        grep(paste0(bird.years[i], "_total_precip"),
                                                                             colnames(cent.bird.cov.can))]
    bird.cov.can.all$maxtemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                         grep(paste0(bird.years[i], "_maxtemp"),
                                                                              colnames(cent.bird.cov.can))]
    bird.cov.can.all$mintemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                         grep(paste0(bird.years[i], "_mintemp"),
                                                                              colnames(cent.bird.cov.can))]
    bird.cov.can.all$mean_maxtemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                              grep(paste0(bird.years[i], "_mean_maxtemp"),
                                                                                   colnames(cent.bird.cov.can))]
    bird.cov.can.all$mean_mintemp[year.set][cov.row.set] <- cent.bird.cov.can[which(!is.na(cent.subset.rows)),
                                                                              grep(paste0(bird.years[i], "_mean_mintemp"),
                                                                                   colnames(cent.bird.cov.can))]
    west.cov.row.set <- west.subset.rows[which(!is.na(west.subset.rows))]
    bird.cov.can.all$precip[year.set][west.cov.row.set] <- west.bird.cov.can[which(!is.na(west.subset.rows)),
                                                                             grep(paste0(bird.years[i], "_mean_total_precip"),
                                                                                  colnames(west.bird.cov.can))]
    bird.cov.can.all$maxtemp[year.set][west.cov.row.set] <- west.bird.cov.can[which(!is.na(west.subset.rows)),
                                                                              grep(paste0(bird.years[i], "_maxtemp"),
                                                                                   colnames(west.bird.cov.can))]
    bird.cov.can.all$mintemp[year.set][west.cov.row.set] <- west.bird.cov.can[which(!is.na(west.subset.rows)),
                                                                              grep(paste0(bird.years[i], "_mintemp"),
                                                                                   colnames(west.bird.cov.can))]
    bird.cov.can.all$mean_maxtemp[year.set][west.cov.row.set] <- west.bird.cov.can[which(!is.na(west.subset.rows)),
                                                                                   grep(paste0(bird.years[i], "_mean_maxtemp"),
                                                                                        colnames(west.bird.cov.can))]
    bird.cov.can.all$mean_mintemp[year.set][west.cov.row.set] <- west.bird.cov.can[which(!is.na(west.subset.rows)),
                                                                                   grep(paste0(bird.years[i], "_mean_mintemp"),
                                                                                        colnames(west.bird.cov.can))]
  }
}

# compile environmental data for butterflies at points with appropriate years
colnames(west.bfs.cov.pts)[3] <- "range"
bfs.cov.tran <- rbind(cent.bfs.cov.pts[, 1:6], west.bfs.cov.pts[, 1:6])
bfs.cov.all <- cbind(bfs.ann.tran.info, bfs.cov.tran[match(bfs.ann.tran.info$transect, bfs.cov.tran$transect), ])
bfs.years <- unique(bfs.cov.all$year)
bfs.cov.all$precip <- rep(NA, nrow(bfs.cov.all))
bfs.cov.all$maxtemp <- rep(NA, nrow(bfs.cov.all))
bfs.cov.all$nectar <- rep(NA, nrow(bfs.cov.all))
bfs.cov.all$mintemp <- rep(NA, nrow(bfs.cov.all))
bfs.cov.all$mud <- rep(NA, nrow(bfs.cov.all))
for (i in seq(along = bfs.years)) {
  year.set <- which(bfs.cov.all$year == bfs.years[i])
  cent.subset.rows <- match(cent.bfs.cov.pts$transect,
                            bfs.cov.all$transect[year.set])
  west.subset.rows <- match(west.bfs.cov.pts$transect,
                            bfs.cov.all$transect[year.set])
  if (all(is.na(west.subset.rows))) {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bfs.cov.all$precip[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                  grep(paste0("precip_", bfs.years[i]),
                                                                       colnames(cent.bfs.cov.pts))]
    bfs.cov.all$maxtemp[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                   grep(paste0("maxtemp_", bfs.years[i]),
                                                                        colnames(cent.bfs.cov.pts))]
    bfs.cov.all$mintemp[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                   grep(paste0("mintemp_", bfs.years[i]),
                                                                        colnames(cent.bfs.cov.pts))]
    if (length(grep(paste0("nectar_", bfs.years[i]), colnames(cent.bfs.cov.pts))) > 0) {
      bfs.cov.all$nectar[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                    grep(paste0("nectar_", bfs.years[i]),
                                                                         colnames(cent.bfs.cov.pts))]
    }
    if (length(grep(paste0("mud_", bfs.years[i]), colnames(cent.bfs.cov.pts))) > 0) {
      bfs.cov.all$mud[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                 grep(paste0("mud_", bfs.years[i]),
                                                                      colnames(cent.bfs.cov.pts))]
    }
  } else {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bfs.cov.all$precip[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                  grep(paste0("precip_", bfs.years[i]),
                                                                       colnames(cent.bfs.cov.pts))]
    bfs.cov.all$maxtemp[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                   grep(paste0("maxtemp_", bfs.years[i]),
                                                                        colnames(cent.bfs.cov.pts))]
    bfs.cov.all$mintemp[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                   grep(paste0("mintemp_", bfs.years[i]),
                                                                        colnames(cent.bfs.cov.pts))]
    if (length(grep(paste0("nectar_", bfs.years[i]), colnames(cent.bfs.cov.pts))) > 0) {
      bfs.cov.all$nectar[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                    grep(paste0("nectar_", bfs.years[i]),
                                                                         colnames(cent.bfs.cov.pts))]
    }
    if (length(grep(paste0("mud_", bfs.years[i]), colnames(cent.bfs.cov.pts))) > 0) {
      bfs.cov.all$mud[year.set][cov.row.set] <- cent.bfs.cov.pts[which(!is.na(cent.subset.rows)),
                                                                 grep(paste0("mud_", bfs.years[i]),
                                                                      colnames(cent.bfs.cov.pts))]
    }
    west.cov.row.set <- west.subset.rows[which(!is.na(west.subset.rows))]
    if (length(grep(paste0("precip_", bfs.years[i]), colnames(west.bfs.cov.pts))) > 0) {
      bfs.cov.all$precip[year.set][west.cov.row.set] <- west.bfs.cov.pts[which(!is.na(west.subset.rows)),
                                                                       grep(paste0("precip_", bfs.years[i]),
                                                                            colnames(west.bfs.cov.pts))]
    }
    if (length(grep(paste0("maxtemp_", bfs.years[i]), colnames(west.bfs.cov.pts))) > 0) {
      bfs.cov.all$maxtemp[year.set][west.cov.row.set] <- west.bfs.cov.pts[which(!is.na(west.subset.rows)),
                                                                        grep(paste0("maxtemp_", bfs.years[i]),
                                                                             colnames(west.bfs.cov.pts))]
    }
    if (length(grep(paste0("mintemp_", bfs.years[i]), colnames(west.bfs.cov.pts))) > 0) {
      bfs.cov.all$mintemp[year.set][west.cov.row.set] <- west.bfs.cov.pts[which(!is.na(west.subset.rows)),
                                                                         grep(paste0("mintemp_", bfs.years[i]),
                                                                              colnames(west.bfs.cov.pts))]
    }
    if (length(grep(paste0("nectar.", bfs.years[i]), colnames(west.bfs.cov.pts))) > 0) {
      bfs.cov.all$nectar[year.set][west.cov.row.set] <- west.bfs.cov.pts[which(!is.na(west.subset.rows)),
                                                                    grep(paste0("nectar.", bfs.years[i]),
                                                                         colnames(west.bfs.cov.pts))]
    }
    if (length(grep(paste0("mud.", bfs.years[i]), colnames(west.bfs.cov.pts))) > 0) {
      bfs.cov.all$mud[year.set][west.cov.row.set] <- west.bfs.cov.pts[which(!is.na(west.subset.rows)),
                                                                 grep(paste0("mud.", bfs.years[i]),
                                                                      colnames(west.bfs.cov.pts))]
    }
  }
}

# compile environmental data for butterflies in canyons with appropriate years
colnames(west.bfs.cov.can)[2] <- "range"
colnames(cent.bfs.cov.can)[2] <- "range"
cent.bfs.cov.can$canRange <- paste0(cent.bfs.cov.can$canyon, cent.bfs.cov.can$range)
west.bfs.cov.can$canRange <- paste0(west.bfs.cov.can$canyon, west.bfs.cov.can$range)
bfs.cov.can <- rbind(cent.bfs.cov.can[, c(1:11, ncol(cent.bfs.cov.can))],
                     west.bfs.cov.can[, c(1:11, ncol(west.bfs.cov.can))])
bfs.ann.can.info$canyon <- gsub(" ", "", bfs.ann.can.info$canyon)
bfs.cov.can.all <- cbind(bfs.ann.can.info, bfs.cov.can[match(bfs.ann.can.info$canyon, bfs.cov.can$canyon), ])
bfs.years <- unique(bfs.cov.can.all$year)
bfs.cov.can.all$mean_precip <- rep(NA, nrow(bfs.cov.can.all))
bfs.cov.can.all$mean_maxtemp <- rep(NA, nrow(bfs.cov.can.all))
bfs.cov.can.all$mean_mintemp <- rep(NA, nrow(bfs.cov.can.all))
for (i in seq(along = bfs.years)) {
  year.set <- which(bfs.cov.can.all$year == bfs.years[i])
  cent.subset.rows <- match(cent.bfs.cov.can$canRange,
                            bfs.cov.can.all$canRange[year.set])
  west.subset.rows <- match(west.bfs.cov.can$canRange,
                            bfs.cov.can.all$canRange[year.set])
  if (all(is.na(west.subset.rows))) {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    bfs.cov.can.all$mean_precip[year.set][cov.row.set] <- cent.bfs.cov.can[which(!is.na(cent.subset.rows)),
                                                                  grep(paste0("mean_precip_", bfs.years[i]),
                                                                       colnames(cent.bfs.cov.can))]
    bfs.cov.can.all$mean_maxtemp[year.set][cov.row.set] <- cent.bfs.cov.can[which(!is.na(cent.subset.rows)),
                                                                   grep(paste0("mean_maxtemp_", bfs.years[i]),
                                                                        colnames(cent.bfs.cov.can))]
    bfs.cov.can.all$mean_mintemp[year.set][cov.row.set] <- cent.bfs.cov.can[which(!is.na(cent.subset.rows)),
                                                                   grep(paste0("mean_mintemp_", bfs.years[i]),
                                                                        colnames(cent.bfs.cov.can))]
  } else {
    cov.row.set <- cent.subset.rows[which(!is.na(cent.subset.rows))]
    if (length(grep(paste0("mean_precip_", bfs.years[i]), colnames(cent.bfs.cov.can))) > 0) {
      bfs.cov.can.all$mean_precip[year.set][cov.row.set] <- cent.bfs.cov.can[which(!is.na(cent.subset.rows)),
                                                                             grep(paste0("mean_precip_", bfs.years[i]),
                                                                                  colnames(cent.bfs.cov.can))]
    }
    if (length(grep(paste0("mean_maxtemp_", bfs.years[i]), colnames(cent.bfs.cov.can))) > 0) {
      bfs.cov.can.all$mean_maxtemp[year.set][cov.row.set] <- cent.bfs.cov.can[which(!is.na(cent.subset.rows)),
                                                                            grep(paste0("mean_maxtemp_", bfs.years[i]),
                                                                                 colnames(cent.bfs.cov.can))]
    }
    if (length(grep(paste0("mean_mintemp_", bfs.years[i]), colnames(cent.bfs.cov.can))) > 0) {
      bfs.cov.can.all$mean_mintemp[year.set][cov.row.set] <- cent.bfs.cov.can[which(!is.na(cent.subset.rows)),
                                                                            grep(paste0("mean_mintemp_", bfs.years[i]),
                                                                                 colnames(cent.bfs.cov.can))]
    }
    west.cov.row.set <- west.subset.rows[which(!is.na(west.subset.rows))]
    bfs.cov.can.all$mean_precip[year.set][west.cov.row.set] <- west.bfs.cov.can[which(!is.na(west.subset.rows)),
                                                                           grep(paste0("mean_precip_", bfs.years[i]),
                                                                                colnames(west.bfs.cov.can))]
    bfs.cov.can.all$mean_maxtemp[year.set][west.cov.row.set] <- west.bfs.cov.can[which(!is.na(west.subset.rows)),
                                                                            grep(paste0("mean_maxtemp_", bfs.years[i]),
                                                                                 colnames(west.bfs.cov.can))]
    bfs.cov.can.all$mean_mintemp[year.set][west.cov.row.set] <- west.bfs.cov.can[which(!is.na(west.subset.rows)),
                                                                            grep(paste0("mean_mintemp_", bfs.years[i]),
                                                                                 colnames(west.bfs.cov.can))]
  }
}

# clear workspace
rm("all.bfs.can.combos", "all.bfs.surv.combos", "all.bird.can.combos", "all.bird.surv.combos",
   "bfs.annual.can.survs", "bfs.annual.survs", "bfs.splist",
   "bfs.can.list", "bfs.can.year.list", "bfs.cov.can", "bfs.cov.tran", "bfs.na.can.survs",
   "bfs.na.survs", "bfs.site.list", "bfs.survs", "bfs.year.list", "bfs.years",
   "bird.annual.can.survs", "bird.annual.survs",
   "bird.can.list", "bird.can.year.list", "bird.cov.can", "bird.cov.pts",
   "bird.na.can.survs", "bird.na.survs", "bird.site.list", "bird.survs",
   "bird.year.list", "bird.years", "bird.years.can", "cent.bfs.cov.can",
   "cent.bfs.cov.pts", "cent.bfs.survey", "cent.bird.cov.can", "cent.bird.cov.pts",
   "cent.bird.survey", "cent.subset.rows", "cov.row.set",
   "exclusions", "exclusions_list", "i",
   "west.bfs.cov.can", "west.bfs.cov.pts", "west.bfs.survey",
   "west.bird.cov.can", "west.bird.cov.pts", "west.bird.survey", "west.cov.row.set",
   "west.subset.rows", "year.set")