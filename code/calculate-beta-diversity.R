# R code to calculate temporal and spatial components of beta diversity in birds and butterflies
#   from the Great Basin.
# Calculations return total beta diversity (Sørensen index), beta diversity due to
#   turnover (Simpson"s index), and beta diversity due to nestedness (Baselga"s index)

# load coordinates for all points, transects, and canyons
cent.pt.coords <- read.csv("./data/cent-pt-coords.csv")
west.pt.coords <- read.csv("./data/west-pt-coords.csv")
cent.can.coords <- read.csv("./data/cent-bird-can-coords.csv")[, c(1:4)]
west.can.coords <- read.csv("./data/west-bird-can-coords.csv")[, c(1:4)]
cent.bfs.coords <- read.csv("./data/cent-bfs-tran-coords.csv")[, c(1, 15:16, 18:19)]
west.bfs.coords <- read.csv("./data/west-bfs-tran-coords.csv")[, c(1, 6:9)]

# calculate temporal beta diversity values at canyon scale for birds and butterflies
beta.all.can <- vector("list", length = 2)
beta.null.can <- vector("list", length = 2)

# birds
beta.can.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bird.ann.can.info$canyon)))
beta.can.null.tmp <- vector("list", length = length(unique(bird.ann.can.info$canyon)))
for (i in 1:length(unique(bird.ann.can.info$canyon))) {
  data.sub <- bird.ann.can.data[which(bird.ann.can.info$canyon == unique(bird.ann.can.info$canyon)[i]), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bird.ann.can.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.can.tmp[i, ] <- beta.set$out.beta
    beta.can.null.tmp[[i]] <- beta.set$out.null
  } else {
    beta.can.tmp[i, ] <- rep(NA, ncol(beta.can.tmp))
    beta.can.null.tmp[[i]] <- NA
  }
}
rownames(beta.can.tmp) <- unique(bird.ann.can.info$canyon)
colnames(beta.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.all.can[[1]] <- beta.can.tmp
beta.null.can[[1]] <- beta.can.null.tmp

# butterflies
beta.can.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bfs.ann.can.info$canyon)))
beta.can.null.tmp <- vector("list", length = length(unique(bfs.ann.can.info$canyon)))
for (i in 1:length(unique(bfs.ann.can.info$canyon))) {
  data.sub <- bfs.ann.can.data[which(bfs.ann.can.info$canyon == unique(bfs.ann.can.info$canyon)[i]), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bfs.ann.can.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.can.tmp[i, ] <- beta.set$out.beta
    beta.can.null.tmp[[i]] <- beta.set$out.null
  } else {
    beta.can.tmp[i, ] <- rep(NA, ncol(beta.can.tmp))
    beta.can.null.tmp[[i]] <- NA
  }
}
rownames(beta.can.tmp) <- unique(bfs.ann.can.info$canyon)
colnames(beta.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.all.can[[2]] <- beta.can.tmp
beta.null.can[[2]] <- beta.can.null.tmp

# calculate temporal beta diversity values at point/transect scale for birds and butterflies
beta.all.pt <- vector("list", length = 2)
beta.null.pt <- vector("list", length = 2)

# birds
beta.pt.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bird.ann.pt.info$point)))
beta.pt.null.tmp <- vector("list", length = length(unique(bird.ann.pt.info$point)))
for (i in 1:length(unique(bird.ann.pt.info$point))) {
  data.sub <- bird.ann.pt.data[which(bird.ann.pt.info$point == unique(bird.ann.pt.info$point)[i]), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bird.ann.pt.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.pt.tmp[i, ] <- beta.set$out.beta
    beta.pt.null.tmp[[i]] <- beta.set$out.null
  } else {
    beta.pt.tmp[i, ] <- rep(NA, ncol(beta.pt.tmp))
    beta.pt.null.tmp[[i]] <- NA
  }
}
rownames(beta.pt.tmp) <- unique(bird.ann.pt.info$point)
colnames(beta.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.all.pt[[1]] <- beta.pt.tmp
beta.null.pt[[1]] <- beta.pt.null.tmp

# butterflies
beta.pt.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bfs.ann.tran.info$transect)))
beta.pt.null.tmp <- vector("list", length = length(unique(bfs.ann.tran.info$transect)))
for (i in 1:length(unique(bfs.ann.tran.info$transect))) {
  data.sub <- bfs.ann.tran.data[which(bfs.ann.tran.info$transect == unique(bfs.ann.tran.info$transect)[i]), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bfs.ann.tran.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.pt.tmp[i, ] <- beta.set$out.beta
    beta.pt.null.tmp[[i]] <- beta.set$out.null
  } else {
    beta.pt.tmp[i, ] <- rep(NA, ncol(beta.pt.tmp))
    beta.pt.null.tmp[[i]] <- NA
  }
}
rownames(beta.pt.tmp) <- unique(bfs.ann.tran.info$transect)
colnames(beta.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.all.pt[[2]] <- beta.pt.tmp
beta.null.pt[[2]] <- beta.pt.null.tmp

# calculate spatial beta diversity values at canyon scale for birds and butterflies
beta.spat.can <- vector("list", length = 2)
beta.spat.null.can <- vector("list", length = 2)

# birds
vals.list <- expand.grid(unique(bird.cov.can.all$mtnrange),
                         unique(bird.cov.can.all$year))
beta.spat.can.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
beta.spat.null.can.tmp <- vector("list", length = nrow(vals.list))
for (i in 1:nrow(vals.list)) {
  data.sub <- bird.ann.can.data[which((bird.cov.can.all$mtnrange == vals.list[i, 1]) &
                                        (bird.cov.can.all$year == vals.list[i, 2])), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bird.ann.can.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.spat.can.tmp[i, ] <- beta.set$out.beta
    beta.spat.null.can.tmp[[i]] <- beta.set$out.null
  } else {
    beta.spat.can.tmp[i, ] <- rep(NA, ncol(beta.spat.can.tmp))
    beta.spat.null.can.tmp[[i]] <- NA
  }
}
rownames(beta.spat.can.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
colnames(beta.spat.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.spat.can[[1]] <- beta.spat.can.tmp
beta.spat.null.can[[1]] <- beta.spat.null.can.tmp

# butterflies
vals.list <- expand.grid(unique(bfs.cov.can.all$range),
                         unique(bfs.cov.can.all$year))
beta.spat.can.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
beta.spat.null.can.tmp <- vector("list", length = nrow(vals.list))
for (i in 1:nrow(vals.list)) {
  data.sub <- bfs.ann.can.data[which((bfs.cov.can.all$range == vals.list[i, 1]) &
                                        (bfs.cov.can.all$year == vals.list[i, 2])), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bfs.ann.can.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.spat.can.tmp[i, ] <- beta.set$out.beta
    beta.spat.null.can.tmp[[i]] <- beta.set$out.null
  } else {
    beta.spat.can.tmp[i, ] <- rep(NA, ncol(beta.spat.can.tmp))
    beta.spat.null.can.tmp[[i]] <- NA
  }
}
rownames(beta.spat.can.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
colnames(beta.spat.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.spat.can[[2]] <- beta.spat.can.tmp
beta.spat.null.can[[2]] <- beta.spat.null.can.tmp

# calculate spatial beta diversity values at point/transect scale for birds and butterflies
beta.spat.pt <- vector("list", length = 2)
beta.spat.null.pt <- vector("list", length = 2)

# birds
vals.list <- expand.grid(unique(bird.cov.all$canyon),
                         unique(bird.cov.all$year))
beta.spat.pt.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
beta.spat.null.pt.tmp <- vector("list", length = nrow(vals.list))
for (i in 1:nrow(vals.list)) {
  data.sub <- bird.ann.pt.data[which((bird.cov.all$canyon == vals.list[i, 1]) &
                                       (bird.cov.all$year == vals.list[i, 2])), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bird.ann.pt.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.spat.pt.tmp[i, ] <- beta.set$out.beta
    beta.spat.null.pt.tmp[[i]] <- beta.set$out.null
  } else {
    beta.spat.pt.tmp[i, ] <- rep(NA, ncol(beta.spat.pt.tmp))
    beta.spat.null.pt.tmp[[i]] <- NA
  }
}
rownames(beta.spat.pt.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
colnames(beta.spat.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.spat.pt[[1]] <- beta.spat.pt.tmp
beta.spat.null.pt[[1]] <- beta.spat.null.pt.tmp

# butterflies
vals.list <- expand.grid(unique(bfs.cov.all$canyon),
                         unique(bfs.cov.all$year))
beta.spat.pt.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
beta.spat.null.pt.tmp <- vector("list", length = nrow(vals.list))
for (i in 1:nrow(vals.list)) {
  data.sub <- bfs.ann.tran.data[which((bfs.cov.all$canyon == vals.list[i, 1]) &
                                        (bfs.cov.all$year == vals.list[i, 2])), ]
  data.sub <- ifelse(data.sub > 0, 1, 0)
  if (length(data.sub) > ncol(bfs.ann.tran.data)) {
    beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
    beta.spat.pt.tmp[i, ] <- beta.set$out.beta
    beta.spat.null.pt.tmp[[i]] <- beta.set$out.null
  } else {
    beta.spat.pt.tmp[i, ] <- rep(NA, ncol(beta.spat.pt.tmp))
    beta.spat.null.pt.tmp[[i]] <- NA
  }
}
rownames(beta.spat.pt.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
colnames(beta.spat.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
beta.spat.pt[[2]] <- beta.spat.pt.tmp
beta.spat.null.pt[[2]] <- beta.spat.null.pt.tmp

rm("beta.can.null.tmp", "beta.can.tmp", "beta.pt.null.tmp",
   "beta.pt.tmp", "beta.set", "beta.spat.can.tmp", "beta.spat.null.can.tmp",
   "beta.spat.null.pt.tmp", "beta.spat.pt.tmp", "data.sub",
   "i", "vals.list")


# calculate beta diversity for specific bird and butterfly guilds

# calculate temporal beta diversity values at canyon scale for bird nesting guilds
beta.nest.can <- vector("list", length = length(unique(bird.guilds$nesting)))
beta.nest.null.can <- vector("list", length = length(unique(bird.guilds$nesting)))
for (i in seq(along = unique(bird.guilds$nesting))) {
  bird.sub <- bird.ann.can.data[, match(bird.guilds$sppname[which(bird.guilds$nesting ==
                                                                    unique(bird.guilds$nesting)[i])],
                                        colnames(bird.ann.can.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.can.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bird.ann.can.info$canyon)))
  beta.can.null.tmp <- vector("list", length = length(unique(bird.ann.can.info$canyon)))
  for (j in 1:length(unique(bird.ann.can.info$canyon))) {
    data.sub <- bird.sub[which(bird.ann.can.info$canyon == unique(bird.ann.can.info$canyon)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.can.tmp[j, ] <- beta.set$out.beta
      beta.can.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.can.tmp[j, ] <- rep(NA, ncol(beta.can.tmp))
      beta.can.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.can.tmp) <- unique(bird.ann.can.info$canyon)
  colnames(beta.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.nest.can[[i]] <- beta.can.tmp
  beta.nest.null.can[[i]] <- beta.can.null.tmp
}
names(beta.nest.can) <- unique(bird.guilds$nesting)
names(beta.nest.null.can) <- unique(bird.guilds$nesting)

# calculate temporal beta diversity values at point scale for bird nesting guilds
beta.nest.pt <- vector("list", length = length(unique(bird.guilds$nesting)))
beta.nest.null.pt <- vector("list", length = length(unique(bird.guilds$nesting)))
for (i in seq(along = unique(bird.guilds$nesting))) {
  bird.sub <- bird.ann.pt.data[, match(bird.guilds$sppname[which(bird.guilds$nesting ==
                                                                   unique(bird.guilds$nesting)[i])],
                                       colnames(bird.ann.pt.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.pt.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bird.ann.pt.info$point)))
  beta.pt.null.tmp <- vector("list", length = length(unique(bird.ann.pt.info$point)))
  for (j in 1:length(unique(bird.ann.pt.info$point))) {
    data.sub <- bird.sub[which(bird.ann.pt.info$point == unique(bird.ann.pt.info$point)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.pt.tmp[j, ] <- beta.set$out.beta
      beta.pt.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.pt.tmp[j, ] <- rep(NA, ncol(beta.pt.tmp))
      beta.pt.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.pt.tmp) <- unique(bird.ann.pt.info$point)
  colnames(beta.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.nest.pt[[i]] <- beta.pt.tmp
  beta.nest.null.pt[[i]] <- beta.pt.null.tmp
}
names(beta.nest.pt) <- unique(bird.guilds$nesting)
names(beta.nest.null.pt) <- unique(bird.guilds$nesting)

# calculate temporal beta diversity values at canyon scale for bird riparian guilds
beta.ripa.can <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
beta.ripa.null.can <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
for (i in seq(along = unique(bird.guilds$riparian.dependence))) {
  bird.sub <- bird.ann.can.data[, match(bird.guilds$sppname[which(bird.guilds$riparian.dependence ==
                                                                    unique(bird.guilds$riparian.dependence)[i])],
                                        colnames(bird.ann.can.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.can.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bird.ann.can.info$canyon)))
  beta.can.null.tmp <- vector("list", length = length(unique(bird.ann.can.info$canyon)))
  for (j in 1:length(unique(bird.ann.can.info$canyon))) {
    data.sub <- bird.sub[which(bird.ann.can.info$canyon == unique(bird.ann.can.info$canyon)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.can.tmp[j, ] <- beta.set$out.beta
      beta.can.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.can.tmp[j, ] <- rep(NA, ncol(beta.can.tmp))
      beta.can.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.can.tmp) <- unique(bird.ann.can.info$canyon)
  colnames(beta.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.ripa.can[[i]] <- beta.can.tmp
  beta.ripa.null.can[[i]] <- beta.can.null.tmp
}
names(beta.ripa.can) <- unique(bird.guilds$riparian.dependence)
names(beta.ripa.null.can) <- unique(bird.guilds$riparian.dependence)

# calculate temporal beta diversity values at point scale for bird riparian guilds
beta.ripa.pt <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
beta.ripa.null.pt <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
for (i in seq(along = unique(bird.guilds$riparian.dependence))) {
  bird.sub <- bird.ann.pt.data[, match(bird.guilds$sppname[which(bird.guilds$riparian.dependence ==
                                                                   unique(bird.guilds$riparian.dependence)[i])],
                                       colnames(bird.ann.pt.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.pt.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bird.ann.pt.info$point)))
  beta.pt.null.tmp <- vector("list", length = length(unique(bird.ann.pt.info$point)))
  for (j in 1:length(unique(bird.ann.pt.info$point))) {
    data.sub <- bird.sub[which(bird.ann.pt.info$point == unique(bird.ann.pt.info$point)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.pt.tmp[j, ] <- beta.set$out.beta
      beta.pt.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.pt.tmp[j, ] <- rep(NA, ncol(beta.pt.tmp))
      beta.pt.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.pt.tmp) <- unique(bird.ann.pt.info$point)
  colnames(beta.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.ripa.pt[[i]] <- beta.pt.tmp
  beta.ripa.null.pt[[i]] <- beta.pt.null.tmp
}
names(beta.ripa.pt) <- unique(bird.guilds$riparian.dependence)
names(beta.ripa.null.pt) <- unique(bird.guilds$riparian.dependence)

# calculate spatial beta diversity values at canyon scale for bird nesting guilds
beta.spat.nest.can <- vector("list", length = length(unique(bird.guilds$nesting)))
beta.spat.nest.null.can <- vector("list", length = length(unique(bird.guilds$nesting)))
vals.list <- expand.grid(unique(bird.cov.can.all$mtnrange),
                         unique(bird.cov.can.all$year))
for (i in seq(along = unique(bird.guilds$nesting))) {
  bird.sub <- bird.ann.can.data[, match(bird.guilds$sppname[which(bird.guilds$nesting ==
                                                                    unique(bird.guilds$nesting)[i])],
                                        colnames(bird.ann.can.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.spat.can.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.can.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bird.sub[which((bird.cov.can.all$mtnrange == vals.list[j, 1]) &
                                 (bird.cov.can.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.can.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.can.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.can.tmp[j, ] <- rep(NA, ncol(beta.spat.can.tmp))
      beta.spat.null.can.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.can.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.nest.can[[i]] <- beta.spat.can.tmp
  beta.spat.nest.null.can[[i]] <- beta.spat.null.can.tmp
}

# calculate spatial beta diversity values at point scale for bird nesting guilds
beta.spat.nest.pt <- vector("list", length = length(unique(bird.guilds$nesting)))
beta.spat.nest.null.pt <- vector("list", length = length(unique(bird.guilds$nesting)))
vals.list <- expand.grid(unique(bird.cov.all$canyon),
                         unique(bird.cov.all$year))
for (i in seq(along = unique(bird.guilds$nesting))) {
  bird.sub <- bird.ann.pt.data[, match(bird.guilds$sppname[which(bird.guilds$nesting ==
                                                                   unique(bird.guilds$nesting)[i])],
                                       colnames(bird.ann.pt.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.spat.pt.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.pt.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bird.sub[which((bird.cov.all$canyon == vals.list[j, 1]) &
                                 (bird.cov.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.pt.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.pt.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.pt.tmp[j, ] <- rep(NA, ncol(beta.spat.pt.tmp))
      beta.spat.null.pt.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.pt.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.nest.pt[[i]] <- beta.spat.pt.tmp
  beta.spat.nest.null.pt[[i]] <- beta.spat.null.pt.tmp
}

# calculate spatial beta diversity values at canyon scale for bird riparian guilds
beta.spat.ripa.can <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
beta.spat.ripa.null.can <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
vals.list <- expand.grid(unique(bird.cov.can.all$mtnrange),
                         unique(bird.cov.can.all$year))
for (i in seq(along = unique(bird.guilds$riparian.dependence))) {
  bird.sub <- bird.ann.can.data[, match(bird.guilds$sppname[which(bird.guilds$riparian.dependence ==
                                                                    unique(bird.guilds$riparian.dependence)[i])],
                                        colnames(bird.ann.can.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.spat.can.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.can.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bird.sub[which((bird.cov.can.all$mtnrange == vals.list[j, 1]) &
                                 (bird.cov.can.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.can.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.can.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.can.tmp[j, ] <- rep(NA, ncol(beta.spat.can.tmp))
      beta.spat.null.can.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.can.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.ripa.can[[i]] <- beta.spat.can.tmp
  beta.spat.ripa.null.can[[i]] <- beta.spat.null.can.tmp
}

# calculate spatial beta diversity values at point scale for bird riparian guilds
beta.spat.ripa.pt <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
beta.spat.ripa.null.pt <- vector("list", length = length(unique(bird.guilds$riparian.dependence)))
vals.list <- expand.grid(unique(bird.cov.all$canyon),
                         unique(bird.cov.all$year))
for (i in seq(along = unique(bird.guilds$riparian.dependence))) {
  bird.sub <- bird.ann.pt.data[, match(bird.guilds$sppname[which(bird.guilds$riparian.dependence ==
                                                                   unique(bird.guilds$riparian.dependence)[i])],
                                       colnames(bird.ann.pt.data))]
  if (!is.matrix(bird.sub)) {
    bird.sub <- matrix(bird.sub, ncol = 1)
  }
  beta.spat.pt.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.pt.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bird.sub[which((bird.cov.all$canyon == vals.list[j, 1]) &
                                 (bird.cov.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    if (length(data.sub) > ncol(bird.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.pt.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.pt.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.pt.tmp[j, ] <- rep(NA, ncol(beta.spat.pt.tmp))
      beta.spat.null.pt.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.pt.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.ripa.pt[[i]] <- beta.spat.pt.tmp
  beta.spat.ripa.null.pt[[i]] <- beta.spat.null.pt.tmp
}

# calculate temporal beta diversity values at canyon scale for butterfly overwintering guilds
guild.list <- c("pupae", "adults", "larvae", "eggs")
beta.owint.can <- vector("list", length = length(guild.list))
beta.owint.null.can <- vector("list", length = length(guild.list))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.can.data[, match(bfs.guilds$sppcode[which(bfs.guilds$overwintering.stage ==
                                                                 guild.list[i])],
                                      colnames(bfs.ann.can.data))]
  beta.can.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bfs.ann.can.info$canyon)))
  beta.can.null.tmp <- vector("list", length = length(unique(bfs.ann.can.info$canyon)))
  for (j in 1:length(unique(bfs.ann.can.info$canyon))) {
    data.sub <- bfs.sub[which(bfs.ann.can.info$canyon == unique(bfs.ann.can.info$canyon)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.can.tmp[j, ] <- beta.set$out.beta
      beta.can.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.can.tmp[j, ] <- rep(NA, ncol(beta.can.tmp))
      beta.can.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.can.tmp) <- unique(bfs.ann.can.info$canyon)
  colnames(beta.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.owint.can[[i]] <- beta.can.tmp
  beta.owint.null.can[[i]] <- beta.can.null.tmp
}
names(beta.owint.can) <- guild.list
names(beta.owint.null.can) <- guild.list

# calculate temporal beta diversity values at transect scale for butterfly overwintering guilds
guild.list <- c("pupae", "adults", "larvae", "eggs")
beta.owint.tran <- vector("list", length = length(guild.list))
beta.owint.null.tran <- vector("list", length = length(guild.list))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.tran.data[, match(bfs.guilds$sppcode[which(bfs.guilds$overwintering.stage ==
                                                                  guild.list[i])],
                                       colnames(bfs.ann.tran.data))]
  beta.tran.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bfs.ann.tran.info$transect)))
  beta.tran.null.tmp <- vector("list", length = length(unique(bfs.ann.tran.info$transect)))
  for (j in 1:length(unique(bfs.ann.tran.info$transect))) {
    data.sub <- bfs.sub[which(bfs.ann.tran.info$transect == unique(bfs.ann.tran.info$transect)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.tran.tmp[j, ] <- beta.set$out.beta
      beta.tran.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.tran.tmp[j, ] <- rep(NA, ncol(beta.tran.tmp))
      beta.tran.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.tran.tmp) <- unique(bfs.ann.tran.info$transect)
  colnames(beta.tran.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.owint.tran[[i]] <- beta.tran.tmp
  beta.owint.null.tran[[i]] <- beta.tran.null.tmp
}
names(beta.owint.tran) <- guild.list
names(beta.owint.null.tran) <- guild.list

# calculate temporal beta diversity values at canyon scale for butterfly vagility guilds
guild.list <- unique(bfs.guilds$vagility)
beta.vgty.can <- vector("list", length = length(guild.list))
beta.vgty.null.can <- vector("list", length = length(guild.list))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.can.data[, match(bfs.guilds$sppcode[which(bfs.guilds$vagility ==
                                                                 guild.list[i])],
                                      colnames(bfs.ann.can.data))]
  beta.can.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bfs.ann.can.info$canyon)))
  beta.can.null.tmp <- vector("list", length = length(unique(bfs.ann.can.info$canyon)))
  for (j in 1:length(unique(bfs.ann.can.info$canyon))) {
    data.sub <- bfs.sub[which(bfs.ann.can.info$canyon == unique(bfs.ann.can.info$canyon)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.can.tmp[j, ] <- beta.set$out.beta
      beta.can.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.can.tmp[j, ] <- rep(NA, ncol(beta.can.tmp))
      beta.can.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.can.tmp) <- unique(bfs.ann.can.info$canyon)
  colnames(beta.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.vgty.can[[i]] <- beta.can.tmp
  beta.vgty.null.can[[i]] <- beta.can.null.tmp
}
names(beta.vgty.can) <- guild.list
names(beta.vgty.null.can) <- guild.list

# calculate temporal beta diversity values at transect scale for butterfly vagility guilds
guild.list <- unique(bfs.guilds$vagility)
beta.vgty.tran <- vector("list", length = length(guild.list))
beta.vgty.null.tran <- vector("list", length = length(guild.list))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.tran.data[, match(bfs.guilds$sppcode[which(bfs.guilds$vagility ==
                                                                  guild.list[i])],
                                       colnames(bfs.ann.tran.data))]
  beta.tran.tmp <- matrix(NA, ncol = 3, nrow = length(unique(bfs.ann.tran.info$transect)))
  beta.tran.null.tmp <- vector("list", length = length(unique(bfs.ann.tran.info$transect)))
  for (j in 1:length(unique(bfs.ann.tran.info$transect))) {
    data.sub <- bfs.sub[which(bfs.ann.tran.info$transect == unique(bfs.ann.tran.info$transect)[j]), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.tran.tmp[j, ] <- beta.set$out.beta
      beta.tran.null.tmp[[j]] <- beta.set$out.null
    } else {
      beta.tran.tmp[j, ] <- rep(NA, ncol(beta.tran.tmp))
      beta.tran.null.tmp[[j]] <- NA
    }
  }
  rownames(beta.tran.tmp) <- unique(bfs.ann.tran.info$transect)
  colnames(beta.tran.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.vgty.tran[[i]] <- beta.tran.tmp
  beta.vgty.null.tran[[i]] <- beta.tran.null.tmp
}
names(beta.vgty.tran) <- guild.list
names(beta.vgty.null.tran) <- guild.list

# calculate spatial beta diversity values at canyon scale for butterfly overwintering guilds
guild.list <- c("pupae", "adults", "larvae", "eggs")
beta.spat.owint.can <- vector("list", length = length(guild.list))
beta.spat.owint.null.can <- vector("list", length = length(guild.list))
vals.list <- expand.grid(unique(bfs.cov.can.all$range),
                         unique(bfs.cov.can.all$year))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.can.data[, match(bfs.guilds$sppcode[which(bfs.guilds$overwintering.stage ==
                                                                 guild.list[i])],
                                      colnames(bfs.ann.can.data))]
  beta.spat.can.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.can.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bfs.sub[which((bfs.cov.can.all$range == vals.list[j, 1]) &
                                (bfs.cov.can.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.can.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.can.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.can.tmp[j, ] <- rep(NA, ncol(beta.spat.can.tmp))
      beta.spat.null.can.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.can.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.owint.can[[i]] <- beta.spat.can.tmp
  beta.spat.owint.null.can[[i]] <- beta.spat.null.can.tmp
}
names(beta.spat.owint.can) <- guild.list
names(beta.spat.owint.null.can) <- guild.list

# calculate spatial beta diversity values at transect scale for butterfly overwintering guilds
guild.list <- c("pupae", "adults", "larvae", "eggs")
beta.spat.owint.tran <- vector("list", length = length(guild.list))
beta.spat.owint.null.tran <- vector("list", length = length(guild.list))
vals.list <- expand.grid(unique(bfs.cov.all$canyon),
                         unique(bfs.cov.all$year))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.tran.data[, match(bfs.guilds$sppcode[which(bfs.guilds$overwintering.stage ==
                                                                  guild.list[i])],
                                       colnames(bfs.ann.tran.data))]
  beta.spat.pt.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.pt.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bfs.sub[which((bfs.cov.all$canyon == vals.list[j, 1]) &
                                (bfs.cov.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.pt.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.pt.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.pt.tmp[j, ] <- rep(NA, ncol(beta.spat.pt.tmp))
      beta.spat.null.pt.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.pt.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.owint.tran[[i]] <- beta.spat.pt.tmp
  beta.spat.owint.null.tran[[i]] <- beta.spat.null.pt.tmp
}
names(beta.spat.owint.tran) <- guild.list
names(beta.spat.owint.null.tran) <- guild.list

# calculate spatial beta diversity values at canyon scale for butterfly vagility guilds
guild.list <- unique(bfs.guilds$vagility)
beta.spat.vgty.can <- vector("list", length = length(guild.list))
beta.spat.vgty.null.can <- vector("list", length = length(guild.list))
vals.list <- expand.grid(unique(bfs.cov.can.all$range),
                         unique(bfs.cov.can.all$year))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.can.data[, match(bfs.guilds$sppcode[which(bfs.guilds$vagility ==
                                                                 guild.list[i])],
                                      colnames(bfs.ann.can.data))]
  beta.spat.can.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.can.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bfs.sub[which((bfs.cov.can.all$range == vals.list[j, 1]) &
                                (bfs.cov.can.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.can.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.can.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.can.tmp[j, ] <- rep(NA, ncol(beta.spat.can.tmp))
      beta.spat.null.can.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.can.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.can.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.vgty.can[[i]] <- beta.spat.can.tmp
  beta.spat.vgty.null.can[[i]] <- beta.spat.null.can.tmp
}
names(beta.spat.vgty.can) <- guild.list
names(beta.spat.vgty.null.can) <- guild.list

# calculate spatial beta diversity values at transect scale for butterfly vagility guilds
guild.list <- unique(bfs.guilds$vagility)
beta.spat.vgty.tran <- vector("list", length = length(guild.list))
beta.spat.vgty.null.tran <- vector("list", length = length(guild.list))
vals.list <- expand.grid(unique(bfs.cov.all$canyon),
                         unique(bfs.cov.all$year))
for (i in seq(along = guild.list)) {
  bfs.sub <- bfs.ann.tran.data[, match(bfs.guilds$sppcode[which(bfs.guilds$vagility ==
                                                                  guild.list[i])],
                                       colnames(bfs.ann.tran.data))]
  beta.spat.pt.tmp <- matrix(NA, ncol = 3, nrow = nrow(vals.list))
  beta.spat.null.pt.tmp <- vector("list", length = nrow(vals.list))
  for (j in 1:nrow(vals.list)) {
    data.sub <- bfs.sub[which((bfs.cov.all$canyon == vals.list[j, 1]) &
                                (bfs.cov.all$year == vals.list[j, 2])), ]
    data.sub <- ifelse(data.sub > 0, 1, 0)
    if (length(data.sub) > ncol(bfs.sub)) {
      beta.set <- beta.fun(data.sub, n.sub = 2, n.sample = 100, n.perm = 999)
      beta.spat.pt.tmp[j, ] <- beta.set$out.beta
      beta.spat.null.pt.tmp[[j]] <- beta.set$out.null
    } else {
      beta.spat.pt.tmp[j, ] <- rep(NA, ncol(beta.spat.pt.tmp))
      beta.spat.null.pt.tmp[[j]] <- NA
    }
  }
  rownames(beta.spat.pt.tmp) <- paste(vals.list[, 1], vals.list[, 2], sep = "_")
  colnames(beta.spat.pt.tmp) <- c("beta_turnover", "beta_nestedness", "beta_total")
  beta.spat.vgty.tran[[i]] <- beta.spat.pt.tmp
  beta.spat.vgty.null.tran[[i]] <- beta.spat.null.pt.tmp
}
names(beta.spat.vgty.tran) <- guild.list
names(beta.spat.vgty.null.tran) <- guild.list

rm("beta.can.null.tmp", "beta.can.tmp", "beta.set", "beta.spat.can.tmp",
   "beta.pt.null.tmp", "beta.pt.tmp",
   "beta.spat.null.can.tmp", "beta.spat.null.pt.tmp", "beta.spat.pt.tmp",
   "beta.tran.null.tmp", "beta.tran.tmp", "bfs.sub", "data.sub", "bird.sub",
   "guild.list", "i", "j", "vals.list")
