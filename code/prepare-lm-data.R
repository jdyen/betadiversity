# beta_mat, preds, coords
# bird pt (turnover + nestedness)
beta.set <- beta.all.pt[[1]]
coords <- matrix(NA, nrow = nrow(beta.set), ncol = 2)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
range <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
canyon <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
match.rows <- match(rownames(beta.set), cent.pt.coords$point)
region[which(!is.na(match.rows)), 1] <- 1
range[, 1] <- bird.cov.all$mtnrange[match(rownames(beta.set), bird.cov.all$point)]
canyon[, 1] <- bird.cov.all$canyon[match(rownames(beta.set), bird.cov.all$point)]
coords[which(!is.na(match.rows)), 1] <- as.numeric(cent.pt.coords[match.rows[which(!is.na(match.rows))], 2])
coords[which(!is.na(match.rows)), 2] <- as.numeric(cent.pt.coords[match.rows[which(!is.na(match.rows))], 3])
match.rows <- match(rownames(beta.set), west.pt.coords$point)
region[which(!is.na(match.rows)), 1] <- 2
coords[which(!is.na(match.rows)), 1] <- as.numeric(west.pt.coords[match.rows[which(!is.na(match.rows))], 2])
coords[which(!is.na(match.rows)), 2] <- as.numeric(west.pt.coords[match.rows[which(!is.na(match.rows))], 3])
range <- as.matrix(as.integer(as.factor(range)))
canyon <- as.matrix(as.integer(as.factor(canyon)))
rownames(coords) <- rownames(beta.set)
rownames(region) <- rownames(beta.set)
rownames(range) <- rownames(beta.set)
rownames(canyon) <- rownames(beta.set)

beta.pred.mean <- matrix(NA, nrow = nrow(beta.set), ncol = (ncol(bird.cov.all) - 5))
beta.pred.sd <- matrix(NA, nrow = nrow(beta.set), ncol = ncol(beta.pred.mean))
for (i in 1:ncol(beta.pred.mean)) {
  beta.pred.mean[, i] <- tapply(as.numeric(bird.cov.all[, (i + 5)]),
                                bird.cov.all$point,
                                mean,
                                na.rm = TRUE)
  beta.pred.sd[, i] <- tapply(as.numeric(bird.cov.all[, (i + 5)]),
                              bird.cov.all$point,
                              sd,
                              na.rm = TRUE)
}
colnames(beta.pred.mean) <- paste0(colnames(bird.cov.all)[6:ncol(bird.cov.all)], '_mean')
colnames(beta.pred.sd) <- paste0(colnames(bird.cov.all)[6:ncol(bird.cov.all)], '_sd')
beta.pred.mean <- beta.pred.mean[, -c(6, 9, 10, 12)]
beta.pred.sd <- beta.pred.sd[, -c(1:10, 12)]
beta.pred <- cbind(beta.pred.mean, beta.pred.sd)

bird.temp.pt <- list(beta = beta.set,
                     pred = beta.pred,
                     region = region,
                     range = range,
                     canyon = canyon,
                     coords = coords)

# bird can (turnover + nestedness)
beta.set <- beta.all.can[[1]]
coords <- matrix(NA, nrow = nrow(beta.set), ncol = 2)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
range <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
match.rows <- match(rownames(beta.set), cent.can.coords$canyon)
region[which(!is.na(match.rows)), 1] <- 1
range[, 1] <- bird.cov.can.all$mtnrange[match(rownames(beta.set), bird.cov.can.all$canyon)]
coords[which(!is.na(match.rows)), 1] <- as.numeric(cent.can.coords[match.rows[which(!is.na(match.rows))], 3])
coords[which(!is.na(match.rows)), 2] <- as.numeric(cent.can.coords[match.rows[which(!is.na(match.rows))], 4])
match.rows <- match(rownames(beta.set), west.can.coords$canyon)
region[which(!is.na(match.rows)), 1] <- 2
coords[which(!is.na(match.rows)), 1] <- as.numeric(west.can.coords[match.rows[which(!is.na(match.rows))], 3])
coords[which(!is.na(match.rows)), 2] <- as.numeric(west.can.coords[match.rows[which(!is.na(match.rows))], 4])
range <- as.matrix(as.integer(as.factor(range)))
rownames(coords) <- rownames(beta.set)
rownames(region) <- rownames(beta.set)
rownames(range) <- rownames(beta.set)

beta.pred.mean <- matrix(NA, nrow = nrow(beta.set), ncol = (ncol(bird.cov.can.all) - 5))
beta.pred.sd <- matrix(NA, nrow = nrow(beta.set), ncol = ncol(beta.pred.mean))
for (i in 1:ncol(beta.pred.mean)) {
  beta.pred.mean[, i] <- tapply(as.numeric(bird.cov.can.all[, (i + 5)]),
                                bird.cov.can.all$canyon,
                                mean,
                                na.rm = TRUE)
  beta.pred.sd[, i] <- tapply(as.numeric(bird.cov.can.all[, (i + 5)]),
                              bird.cov.can.all$canyon,
                              sd,
                              na.rm = TRUE)
}
colnames(beta.pred.mean) <- paste0(colnames(bird.cov.can.all)[6:ncol(bird.cov.can.all)], '_mean')
colnames(beta.pred.sd) <- paste0(colnames(bird.cov.can.all)[6:ncol(bird.cov.can.all)], '_sd')
beta.pred.mean <- beta.pred.mean[, -c(6, 9, 10, 12)]
beta.pred.sd <- beta.pred.sd[, -c(1:7, 9:11)]
beta.pred <- cbind(beta.pred.mean, beta.pred.sd)
beta.pred <- beta.pred[, -c(9)]

bird.temp.can <- list(beta = beta.set,
                      pred = beta.pred,
                      region = region,
                      range = range,
                      coords = coords)

# bfs tran (turnover + nestedness)
beta.set <- beta.all.pt[[2]]
coords <- matrix(NA, nrow = nrow(beta.set), ncol = 2)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
range <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
canyon <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
match.rows <- match(rownames(beta.set), cent.bfs.coords$transect)
region[which(!is.na(match.rows)), 1] <- 1
range[, 1] <- bfs.cov.all$range[match(rownames(beta.set), bfs.cov.all$transect)]
canyon[, 1] <- bfs.cov.all$canyon[match(rownames(beta.set), bfs.cov.all$transect)]
coords[which(!is.na(match.rows)), 1] <- as.numeric(cent.bfs.coords[match.rows[which(!is.na(match.rows))], 4])
coords[which(!is.na(match.rows)), 2] <- as.numeric(cent.bfs.coords[match.rows[which(!is.na(match.rows))], 5])
match.rows <- match(rownames(beta.set), west.bfs.coords$transect)
region[which(!is.na(match.rows)), 1] <- 2
coords[which(!is.na(match.rows)), 1] <- as.numeric(west.bfs.coords[match.rows[which(!is.na(match.rows))], 4])
coords[which(!is.na(match.rows)), 2] <- as.numeric(west.bfs.coords[match.rows[which(!is.na(match.rows))], 5])
range <- as.matrix(as.integer(as.factor(range)))
canyon <- as.matrix(as.integer(as.factor(canyon)))
rownames(coords) <- rownames(beta.set)
rownames(region) <- rownames(beta.set)
rownames(range) <- rownames(beta.set)
rownames(canyon) <- rownames(beta.set)

beta.pred.mean <- matrix(NA, nrow = nrow(beta.set), ncol = (ncol(bfs.cov.all) - 5))
beta.pred.sd <- matrix(NA, nrow = nrow(beta.set), ncol = ncol(beta.pred.mean))
for (i in 1:ncol(beta.pred.mean)) {
  beta.pred.mean[, i] <- tapply(as.numeric(bfs.cov.all[, (i + 5)]),
                                bfs.cov.all$transect,
                                mean,
                                na.rm = TRUE)
  beta.pred.sd[, i] <- tapply(as.numeric(bfs.cov.all[, (i + 5)]),
                              bfs.cov.all$transect,
                              sd,
                              na.rm = TRUE)
}
colnames(beta.pred.mean) <- paste0(colnames(bfs.cov.all)[6:ncol(bfs.cov.all)], '_mean')
colnames(beta.pred.sd) <- paste0(colnames(bfs.cov.all)[6:ncol(bfs.cov.all)], '_sd')
beta.pred.mean <- beta.pred.mean[, -5]
beta.pred.sd <- beta.pred.sd[, -c(1:3, 5, 7)]
beta.pred <- cbind(beta.pred.mean, beta.pred.sd)

bfs.temp.tran <- list(beta = beta.set,
                      pred = beta.pred,
                      region = region,
                      range = range,
                      canyon = canyon,
                      coords = coords)

# bfs can (turnover + nestedness)
beta.set <- beta.all.can[[2]]
coords <- matrix(NA, nrow = nrow(beta.set), ncol = 2)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
range <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
match.rows <- match(rownames(beta.set), cent.bfs.coords$canyon)
region[which(!is.na(match.rows)), 1] <- 1
range[, 1] <- bfs.cov.all$range[match(rownames(beta.set), bfs.cov.all$canyon)]
coords[which(!is.na(match.rows)), 1] <- as.numeric(cent.bfs.coords[match.rows[which(!is.na(match.rows))], 4])
coords[which(!is.na(match.rows)), 2] <- as.numeric(cent.bfs.coords[match.rows[which(!is.na(match.rows))], 5])
match.rows <- match(rownames(beta.set), west.bfs.coords$canyon)
region[which(!is.na(match.rows)), 1] <- 2
coords[which(!is.na(match.rows)), 1] <- as.numeric(west.bfs.coords[match.rows[which(!is.na(match.rows))], 4])
coords[which(!is.na(match.rows)), 2] <- as.numeric(west.bfs.coords[match.rows[which(!is.na(match.rows))], 5])
range <- as.matrix(as.integer(as.factor(range)))
rownames(coords) <- rownames(beta.set)
rownames(region) <- rownames(beta.set)
rownames(range) <- rownames(beta.set)

beta.pred.mean <- matrix(NA, nrow = nrow(beta.set), ncol = (ncol(bfs.cov.can.all) - 5))
beta.pred.sd <- matrix(NA, nrow = nrow(beta.set), ncol = ncol(beta.pred.mean))
pred.col.id <- c(5:13, 15:ncol(bfs.cov.can.all))
for (i in 1:ncol(beta.pred.mean)) {
  beta.pred.mean[, i] <- tapply(as.numeric(bfs.cov.can.all[, pred.col.id[i]]),
                                bfs.cov.can.all$canyon,
                                mean,
                                na.rm = TRUE)
  beta.pred.sd[, i] <- tapply(as.numeric(bfs.cov.can.all[, pred.col.id[i]]),
                              bfs.cov.can.all$canyon,
                              sd,
                              na.rm = TRUE)
}
colnames(beta.pred.mean) <- paste0(colnames(bfs.cov.can.all)[pred.col.id], '_mean')
colnames(beta.pred.sd) <- paste0(colnames(bfs.cov.can.all)[pred.col.id], '_sd')
beta.pred.mean <- beta.pred.mean[, -c(6, 8, 11)]
beta.pred.sd <- beta.pred.sd[, -c(1:11)]
beta.pred <- cbind(beta.pred.mean, beta.pred.sd)
colnames(beta.pred)[ncol(beta.pred)] <- 'mean_mintemp_sd'

bfs.temp.can <- list(beta = beta.set,
                     pred = beta.pred,
                     region = region,
                     range = range,
                     coords = coords)

# bird pt spat (turnover + nestedness)
beta.set <- beta.spat.pt[[1]]
coords <- matrix(NA, nrow = nrow(beta.set), ncol = 2)
range <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
year <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
beta.row.names <- sapply(strsplit(rownames(beta.set), '_'),
                         function(x) x[1])
beta.row.names <- gsub(' ', '', beta.row.names)
range[, 1] <- bird.cov.can.all$mtnrange[match(beta.row.names, bird.cov.can.all$canyon)]
year[, 1] <- bird.cov.can.all$year[match(beta.row.names, bird.cov.can.all$canyon)]
match.rows <- match(beta.row.names, gsub(' ', '', cent.can.coords$canyon))
region[which(!is.na(match.rows)), 1] <- 1
coords[which(!is.na(match.rows)), 1] <- as.numeric(cent.can.coords[match.rows[which(!is.na(match.rows))], 3])
coords[which(!is.na(match.rows)), 2] <- as.numeric(cent.can.coords[match.rows[which(!is.na(match.rows))], 4])
match.rows <- match(beta.row.names, gsub(' ', '', west.can.coords$canyon))
region[which(!is.na(match.rows)), 1] <- 2
coords[which(!is.na(match.rows)), 1] <- as.numeric(west.can.coords[match.rows[which(!is.na(match.rows))], 3])
coords[which(!is.na(match.rows)), 2] <- as.numeric(west.can.coords[match.rows[which(!is.na(match.rows))], 4])
rownames(coords) <- rownames(beta.set)
range <- as.integer(as.factor(range))
year <- as.integer(as.factor(year))

beta.pred <- matrix(NA, nrow = nrow(beta.set), ncol = (ncol(bird.cov.can.all) - 5))
year.list <- sapply(strsplit(rownames(beta.set), '_'),
                    function(x) x[2])
for (i in 1:nrow(beta.set)) {
    row.id <- which((bird.cov.can.all$year == year.list[i]) &
                      (gsub(' ', '', bird.cov.can.all$canyon) == beta.row.names[i]))
    beta.pred[i, ] <- as.numeric(bird.cov.can.all[row.id, (6:ncol(bird.cov.can.all))])
}
colnames(beta.pred) <- colnames(bird.cov.can.all)[6:ncol(bird.cov.can.all)]
beta.pred <- beta.pred[, -c(6, 9, 10, 12)]

bird.spat.pt <- list(beta = beta.set,
                     pred = beta.pred,
                     year = year,
                     range = range,
                     region = region,
                     coords = coords)

# bfs tran spat (turnover + nestedness)
beta.set <- beta.spat.pt[[2]]
coords <- matrix(NA, nrow = nrow(beta.set), ncol = 2)
range <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
year <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
beta.row.names <- sapply(strsplit(rownames(beta.set), '_'),
                         function(x) x[1])
beta.row.names <- gsub(' ', '', beta.row.names)
range[, 1] <- bfs.cov.can.all$range[match(beta.row.names, bfs.cov.can.all$canyon)]
year[, 1] <- bfs.cov.can.all$year[match(beta.row.names, bfs.cov.can.all$canyon)]
match.rows <- match(beta.row.names, gsub(' ', '', cent.can.coords$canyon))
region[which(!is.na(match.rows)), 1] <- 1
coords[which(!is.na(match.rows)), 1] <- as.numeric(cent.can.coords[match.rows[which(!is.na(match.rows))], 3])
coords[which(!is.na(match.rows)), 2] <- as.numeric(cent.can.coords[match.rows[which(!is.na(match.rows))], 4])
match.rows <- match(beta.row.names, gsub(' ', '', west.can.coords$canyon))
region[which(!is.na(match.rows)), 1] <- 2
coords[which(!is.na(match.rows)), 1] <- as.numeric(west.can.coords[match.rows[which(!is.na(match.rows))], 3])
coords[which(!is.na(match.rows)), 2] <- as.numeric(west.can.coords[match.rows[which(!is.na(match.rows))], 4])
rownames(coords) <- rownames(beta.set)
range <- as.integer(as.factor(range))
year <- as.integer(as.factor(year))

beta.pred <- matrix(NA, nrow = nrow(beta.set), ncol = (ncol(bfs.cov.can.all) - 5))
year.list <- sapply(strsplit(rownames(beta.set), '_'),
                    function(x) x[2])
for (i in 1:nrow(beta.set)) {
  row.id <- which((bfs.cov.can.all$year == year.list[i]) &
                    (gsub(' ', '', bfs.cov.can.all$canyon) == beta.row.names[i]))
  beta.pred[i, ] <- as.numeric(bfs.cov.can.all[row.id, c(5:13, 15:ncol(bfs.cov.can.all))])
}
colnames(beta.pred) <- colnames(bfs.cov.can.all)[c(5:13, 15:ncol(bfs.cov.can.all))]
beta.pred <- beta.pred[, -c(6, 9, 11)]

bfs.spat.tran <- list(beta = beta.set,
                      pred = beta.pred,
                      year = year,
                      range = range,
                      region = region,
                      coords = coords)

# bird can spat (turnover + nestedness)
beta.set <- beta.spat.can[[1]]
coords <- NULL
range <- NULL
beta.pred <- NULL
year <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
beta.row.names <- sapply(strsplit(rownames(beta.set), '_'),
                         function(x) x[1])
beta.row.names <- gsub(' ', '', beta.row.names)
year[, 1] <- bird.cov.can.all$year[match(beta.row.names, bird.cov.can.all$mtnrange)]
match.rows <- match(beta.row.names, gsub(' ', '', cent.can.coords$range))
region[which(!is.na(match.rows)), 1] <- 1
match.rows <- match(beta.row.names, gsub(' ', '', west.can.coords$range))
region[which(!is.na(match.rows)), 1] <- 2
year <- as.integer(as.factor(year))

bird.spat.can <- list(beta = beta.set,
                      pred = beta.pred,
                      year = year,
                      range = range,
                      region = region,
                      coords = coords)

# bfs can spat (turnover + nestedness)
beta.set <- beta.spat.can[[2]]
coords <- NULL
range <- NULL
beta.pred <- NULL
year <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
region <- matrix(NA, nrow = nrow(beta.set), ncol = 1)
beta.row.names <- sapply(strsplit(rownames(beta.set), '_'),
                         function(x) x[1])
beta.row.names <- gsub(' ', '', beta.row.names)
year[, 1] <- bfs.cov.can.all$year[match(beta.row.names, bfs.cov.can.all$range)]
match.rows <- match(beta.row.names, gsub(' ', '', cent.can.coords$range))
region[which(!is.na(match.rows)), 1] <- 1
match.rows <- match(beta.row.names, gsub(' ', '', west.can.coords$range))
region[which(!is.na(match.rows)), 1] <- 2
year <- as.integer(as.factor(year))

bfs.spat.can <- list(beta = beta.set,
                     pred = beta.pred,
                     year = year,
                     range = range,
                     region = region,
                     coords = coords)
