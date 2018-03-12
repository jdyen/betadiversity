# calculate total species richness
sum.fun <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(sum(ifelse(x > 0, 1, 0), na.rm = TRUE))
  }
}

# function to calculate null beta diversity for each sample
null.fun <- function(non.null.data) {
  null.data <- matrix(0, nrow = nrow(non.null.data), ncol = ncol(non.null.data))
  for (i in 1:nrow(non.null.data)) {
    null.set <- sample(1:ncol(non.null.data), size = sum(non.null.data[i, ]), replace = FALSE)
    null.data[i, null.set] <- 1
  }
  null.data
}

# function to calculate spatial weights for CAR error term
wt.fun <- function(x, l) {
  exp(-((x / (0.4 * l)) ** 2))
}

# function to subsample surveys used to calculate beta diversity values
beta.fun <- function(sad.data, n.sub = 3, n.sample = 100, n.perm = 999) {
  if (n.sub == 1) {
    stop("Must have more than one site to calculate beta diversity", call. = FALSE)
  }
  beta.null <- matrix(NA, nrow = n.perm, ncol = 3)
  beta.tmp <- matrix(NA, nrow = n.sample, ncol = 3)
  for (i in 1:n.sample) {
    data.sub <- sad.data[sample(1:nrow(sad.data), size = n.sub, replace = FALSE), ]
    if (!is.matrix(data.sub)) {
      data.sub <- matrix(data.sub, ncol = 1)
    }
    beta.tmp[i, ] <- unlist(beta.multi(data.sub, index.family = "sorensen"))
  }
  for (i in 1:n.perm) {
    data.sub2 <- sad.data[sample(1:nrow(sad.data), size = n.sub, replace = FALSE), ]
    if (!is.matrix(data.sub2)) {
      data.sub2 <- matrix(data.sub2, ncol = 1)
    }
    data.null <- null.fun(data.sub2)
    beta.null[i, ] <- unlist(beta.multi(data.null, index.family = "sorensen"))
  }
  out.beta <- apply(beta.tmp, 2, mean, na.rm = TRUE)
  out.null <- beta.null
  return(list(out.beta = out.beta, out.null = out.null))
}
