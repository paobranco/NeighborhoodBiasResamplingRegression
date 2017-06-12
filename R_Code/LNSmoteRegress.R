## ===================================================
## Smote for regression with Neighborhood bias
# ---------------------------------------------------
LNRFNFSmoteRegress <- function(form, dat, method = method, npts = npts, 
                           control.pts=control.pts, thr.rel = 0.5,
                           C.perc = "balance", k = 5, repl = FALSE,
                           dist = "Euclidean", p = 2)
  
  # INPUTS:
  # form    a model formula
  # dat     the original training set (with the unbalanced distribution)
  # method  parameter previously defined for the relevance function
  # npts    parameter previously defined for the relevance function
  # control.pts parameter previously defined for the relevance function
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc  is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the 
  #         threshold are increased by this percentage. The under sampling
  #         percentage means that the normal cases (cases below the threshold)
  #         are under-sampled by this percentage. Alternatively it may be
  #         "balance" or "extreme", cases where the sampling percentages
  #         are automatically estimated.
  # k       is the number of neighbors to consider as the pool from where
  #         the new synthetic examples are generated
  # repl    is it allowed to perform sampling with replacement
  # dist    is the distance measure to be used (defaults to "Euclidean")
  # p       is a parameter used when a p-norm is computed
{
 
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  if (tgt < ncol(dat)) {
    orig.order <- colnames(dat)
    cols <- 1:ncol(dat)
    cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
    dat <- dat[, cols]
  }
  if (is.na(thr.rel)) {
    stop("Future work!")
  }
  
  y <- dat[, ncol(dat)]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)

  pc <- list()
  pc$method <- method
  pc$npts <- npts
  pc$control.pts <- control.pts
  
  temp <- y.relev <- phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. 
         Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. 
         Please, redefine your relevance function!")
  }

  bumps <- c()
  for (i in 1:(length(y) - 1)) { 
    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
   }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
    obs.ind <- as.list(rep(NA, nbump))
    last <- 1
    for (i in 1:length(bumps)) {
      obs.ind[[i]] <- s.y[last:bumps[i]]
      last <- bumps[i] + 1
    }
    obs.ind[[nbump]] <- s.y[last:length(s.y)]

  newdata <- data.frame()
  
  if (is.list(C.perc)) {
    if (length(C.perc) != nbump){
      stop("The percentages provided must be the same length as the number
           of bumps!")
    }
  } else if (C.perc == "balance") {
    # estimate the percentages of over/under sampling
    B <- round(nrow(dat)/nbump, 0)
    C.perc <- B/sapply(obs.ind, length)        
  } else if (C.perc == "extreme") {
    B <- round(nrow(dat)/nbump, 0)
    rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
    obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
    C.perc <- round(obj/sapply(obs.ind, length), 1)
  }
  
  for (i in 1:nbump) {
    if (C.perc[[i]] == 1) {
      newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] > 1) {
      newExs <- Smote.exsRegressFr(dat, names(obs.ind[[i]]),
                                 ncol(dat),
                                 C.perc[[i]],
                                 k,
                                 dist,
                                 p,
                                 Fr=TRUE)
      # add original rare examples and synthetic generated examples
      newdata <- rbind(newdata, newExs, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] < 1) {
      
      kNNs <- neighbours(tgt, dat, dist, p, k)
      ind <- names(obs.ind[[i]])
      indpos <- match(ind, rownames(dat))
      # determine nr of examples out of the partition in the KNN
      # the higher the r the higher the number of neighbors that are 
      # not in the same partition
      r <- c()
      for(ex in indpos){
          r <- c(r, length(which(!(kNNs[ex,] %in% (indpos))))/k)
      }
      
      # to avoid having error of "too few positive probabilities"
      new.r <- r+0.1
      
      # case vector r is always zero!
      if (any(is.na(new.r))){
        new.r <- NULL
      }
      # this reinforces the frontier
      sel.maj <- sample(indpos,
                        as.integer(C.perc[[i]] * length(obs.ind[[i]])),
                        replace = repl,
                        prob=new.r)
      newdata <- rbind(newdata, dat[sel.maj, ])
      
    }
  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[, cols]
    dat <- dat[, cols]
  }
  
  newdata
}


LNRNFSmoteRegress <- function(form, dat, method = method, npts = npts, 
                               control.pts=control.pts, thr.rel = 0.5,
                               C.perc = "balance", k = 5, repl = FALSE,
                               dist = "Euclidean", p = 2)
  
  # INPUTS:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # method  parameter previously defined for the relevance function
  # npts    parameter previously defined for the relevance function
  # control.pts parameter previously defined for the relevance function
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the 
  #         threshold are increased by this percentage. The under sampling
  #         percentage means that the normal cases (cases below the threshold)
  #         are under-sampled by this percentage. Alternatively it may be
  #         "balance" or "extreme", cases where the sampling percentages
  #         are automatically estimated.
  # k       is the number of neighbors to consider as the pool from where
  #         the new synthetic examples are generated
  # repl    is it allowed to perform sampling with replacement
  # dist    is the distance measure to be used (defaults to "Euclidean")
  # p       is a parameter used when a p-norm is computed
  {
  
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  if (tgt < ncol(dat)) {
    orig.order <- colnames(dat)
    cols <- 1:ncol(dat)
    cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
    dat <- dat[, cols]
  }
  if (is.na(thr.rel)) {
    stop("Future work!")
  }

  y <- dat[, ncol(dat)]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  pc <- list()
  pc$method <- method
  pc$npts <- npts
  pc$control.pts <- control.pts
  
  temp <- y.relev <- phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. 
         Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. 
         Please, redefine your relevance function!")
  }
  
  bumps <- c()
  for (i in 1:(length(y) - 1)) { 
    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]
  
  newdata <- data.frame()
  
  if (is.list(C.perc)) {
    if (length(C.perc) != nbump){
      stop("The percentages provided must be the same length as the number
           of bumps!")
    }
    } else if (C.perc == "balance") {
      # estimate the percentages of over/under sampling
      B <- round(nrow(dat)/nbump, 0)
      C.perc <- B/sapply(obs.ind, length)        
  } else if (C.perc == "extreme") {
    B <- round(nrow(dat)/nbump, 0)
    rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
    obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
    C.perc <- round(obj/sapply(obs.ind, length), 1)
  }
  
  for (i in 1:nbump) {
    if (C.perc[[i]] == 1) {
      newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] > 1) {
      newExs <- Smote.exsRegressFr(dat, names(obs.ind[[i]]),
                                   ncol(dat),
                                   C.perc[[i]],
                                   k,
                                   dist,
                                   p,
                                   Fr=FALSE)
      # add original rare examples and synthetic generated examples
      newdata <- rbind(newdata, newExs, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] < 1) {
      
      kNNs <- neighbours(tgt, dat, dist, p, k)
      ind <- names(obs.ind[[i]])
      indpos <- match(ind, rownames(dat))
      # determine nr of examples out of the partition in the KNN
      # the higher the r the higher the number of neighbors that are 
      # not in the same partition
      r <- c()
      for(ex in indpos){
        r <- c(r, length(which(!(kNNs[ex,] %in% (indpos))))/k)
      }

      # to avoid error of "too few positive probabilities"
      new.r <- r+0.1      

      # case vector r is always zero!
      if (any(is.na(new.r))){
        new.r <- NULL
      }
      # this reinforces the frontier
      sel.maj <- sample(indpos,
                        as.integer(C.perc[[i]] * length(obs.ind[[i]])),
                        replace = repl,
                        prob=new.r)
      newdata <- rbind(newdata, dat[sel.maj, ])
      
    }
  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[, cols]
    dat <- dat[, cols]
  }
  
  newdata
  }



LNRNSmoteRegress <- function(form, dat, method = method, npts = npts, 
                               control.pts=control.pts, thr.rel = 0.5,
                               C.perc = "balance", k = 5, repl = FALSE,
                               dist = "Euclidean", p = 2)
  
  # INPUTS:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # method  parameter previously defined for the relevance function
  # npts    parameter previously defined for the relevance function
  # control.pts parameter previously defined for the relevance function
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the 
  #         threshold are increased by this percentage. The under sampling
  #         percentage means that the normal cases (cases below the threshold)
  #         are under-sampled by this percentage. Alternatively it may be
  #         "balance" or "extreme", cases where the sampling percentages
  #         are automatically estimated.
  # k       is the number of neighbors to consider as the pool from where
  #         the new synthetic examples are generated
  # repl    is it allowed to perform sampling with replacement
  # dist    is the distance measure to be used (defaults to "Euclidean")
  # p       is a parameter used when a p-norm is computed
  {
  
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  if (tgt < ncol(dat)) {
    orig.order <- colnames(dat)
    cols <- 1:ncol(dat)
    cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
    dat <- dat[, cols]
  }
  if (is.na(thr.rel)) {
    stop("Future work!")
  }
  
  y <- dat[, ncol(dat)]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  pc <- list()
  pc$method <- method
  pc$npts <- npts
  pc$control.pts <- control.pts
  
  temp <- y.relev <- phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. 
         Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. 
         Please, redefine your relevance function!")
  }
  
  bumps <- c()
  for (i in 1:(length(y) - 1)) { 
    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]
  
  newdata <- data.frame()
  
  if (is.list(C.perc)) {
    if (length(C.perc) != nbump){
      stop("The percentages provided must be the same length as the number
           of bumps!")
    }
    } else if (C.perc == "balance") {
      # estimate the percentages of over/under sampling
      B <- round(nrow(dat)/nbump, 0)
      C.perc <- B/sapply(obs.ind, length)        
  } else if (C.perc == "extreme") {
    B <- round(nrow(dat)/nbump, 0)
    rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
    obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
    C.perc <- round(obj/sapply(obs.ind, length), 1)
  }
  
  for (i in 1:nbump) {
    if (C.perc[[i]] == 1) {
      newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] > 1) {
      newExs <- Smote.exsRegressFr(dat, names(obs.ind[[i]]),
                                   ncol(dat),
                                   C.perc[[i]],
                                   k,
                                   dist,
                                   p,
                                   Fr=FALSE)
      # add original rare examples and synthetic generated examples
      newdata <- rbind(newdata, newExs, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] < 1) {
      
      kNNs <- neighbours(tgt, dat, dist, p, k)
      ind <- names(obs.ind[[i]])
      indpos <- match(ind, rownames(dat))
      # determine nr of examples out of the partition in the KNN
      # the higher the r the higher the number of neighbors that are 
      # in the same partition
      r <- c()
      for(ex in indpos){
        r <- c(r, length(which((kNNs[ex,] %in% (indpos))))/k)
      }
      
      # to avoid error of "too few positive probabilities"
      new.r <- r+0.1      

      # case vector r is always zero!
      if (any(is.na(new.r))){
        new.r <- NULL
      }
      # this reinforces the far from frontier cases
      sel.maj <- sample(indpos,
                        as.integer(C.perc[[i]] * length(obs.ind[[i]])),
                        replace = repl,
                        prob=new.r)
      newdata <- rbind(newdata, dat[sel.maj, ])
      
    }
  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[, cols]
    dat <- dat[, cols]
  }
  
  newdata
  }


LNRFNSmoteRegress <- function(form, dat, method = method, npts = npts, 
                               control.pts=control.pts, thr.rel = 0.5,
                               C.perc = "balance", k = 5, repl = FALSE,
                               dist = "Euclidean", p = 2)
  
  # INPUTS:
  # form    a model formula
  # dat    the original training set (with the unbalanced distribution)
  # method  parameter previously defined for the relevance function
  # npts    parameter previously defined for the relevance function
  # control.pts parameter previously defined for the relevance function
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the percentage of under- or/and 
  #         over-sampling to apply to each "class" obtained with the threshold.
  #         The over-sampling percentage means that the examples above the 
  #         threshold are increased by this percentage. The under sampling
  #         percentage means that the normal cases (cases below the threshold)
  #         are under-sampled by this percentage. Alternatively it may be
  #         "balance" or "extreme", cases where the sampling percentages
  #         are automatically estimated.
  # k       is the number of neighbors to consider as the pool from where
  #         the new synthetic examples are generated
  # repl    is it allowed to perform sampling with replacement
  # dist    is the distance measure to be used (defaults to "Euclidean")
  # p       is a parameter used when a p-norm is computed
  {
  
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }
  
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  if (tgt < ncol(dat)) {
    orig.order <- colnames(dat)
    cols <- 1:ncol(dat)
    cols[c(tgt, ncol(dat))] <- cols[c(ncol(dat), tgt)]
    dat <- dat[, cols]
  }
  if (is.na(thr.rel)) {
    stop("Future work!")
  }

  y <- dat[, ncol(dat)]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  pc <- list()
  pc$method <- method
  pc$npts <- npts
  pc$control.pts <- control.pts
  
  temp <- y.relev <- phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. 
         Please, redefine your relevance function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. 
         Please, redefine your relevance function!")
  }
  
  bumps <- c()
  for (i in 1:(length(y) - 1)) { 
    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) || 
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }
  }
  nbump <- length(bumps) + 1 # number of different "classes"
  
  # collect the indexes in each "class"
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]
  
  newdata <- data.frame()
  
  if (is.list(C.perc)) {
    if (length(C.perc) != nbump){
      stop("The percentages provided must be the same length as the number
           of bumps!")
    }
    } else if (C.perc == "balance") {
      # estimate the percentages of over/under sampling
      B <- round(nrow(dat)/nbump, 0)
      C.perc <- B/sapply(obs.ind, length)        
  } else if (C.perc == "extreme") {
    B <- round(nrow(dat)/nbump, 0)
    rescale <- nbump * B/sum(B^2/sapply(obs.ind, length))
    obj <- round((B^2/sapply(obs.ind, length)) * rescale, 2)
    C.perc <- round(obj/sapply(obs.ind, length), 1)
  }
  
  for (i in 1:nbump) {
    if (C.perc[[i]] == 1) {
      newdata <- rbind(newdata, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] > 1) {
      newExs <- Smote.exsRegressFr(dat, names(obs.ind[[i]]),
                                   ncol(dat),
                                   C.perc[[i]],
                                   k,
                                   dist,
                                   p,
                                   Fr=TRUE)
      # add original rare examples and synthetic generated examples
      newdata <- rbind(newdata, newExs, dat[names(obs.ind[[i]]), ])
    } else if (C.perc[[i]] < 1) {
      
      kNNs <- neighbours(tgt, dat, dist, p, k)
      ind <- names(obs.ind[[i]])
      indpos <- match(ind, rownames(dat))
      # determine nr of examples out of the partition in the KNN
      # the higher the r the higher the number of neighbors that are 
      # not in the same partition
      r <- c()
      for(ex in indpos){
        r <- c(r, length(which((kNNs[ex,] %in% (indpos))))/k)
      }
      
      # to avoid error of "too few positive probabilities"
      new.r <- r+0.1      
      
      # case vector r is always zero!
      if (any(is.na(new.r))){
        new.r <- NULL
      }
      # this reinforces the far from frontier cases
      sel.maj <- sample(indpos,
                        as.integer(C.perc[[i]] * length(obs.ind[[i]])),
                        replace = repl,
                        prob=new.r)
      newdata <- rbind(newdata, dat[sel.maj, ])
      
    }
  }
  
  if (tgt < ncol(dat)) {
    newdata <- newdata[, cols]
    dat <- dat[, cols]
  }
  
  newdata
  }


# ===================================================
# Obtain a set of smoted examples for a set of rare cases.
#
# ---------------------------------------------------
Smote.exsRegressFr <- function(orig, ind, tgt, N, k, dist, p, Fr)
  # INPUTS:
  # orig the original data set with rare and normal cases
  # ind are character indexes of the rare cases (the minority "class" cases)
  # tgt the column nr of the target variable
  # N is the percentage of over-sampling to carry out;
  # and k is the number of nearest neighours
  # dist is the distance function used for the neighours computation
  # p is an integer used when a "p-norm" distance is selected
  # Fr logical indicating the behaviour that should be applied. 
  # If TRUE the frontier is reinforced, else the safe points far away
  # from the frontier are reinforced
  # OUTPUTS:
  # The result of the function is a (N-1)*nrow(dat) set of generate
  # examples with rare values on the target
{
  indpos <- match(ind, rownames(orig))
  dat <- orig[indpos,]
  # check for constant features and remove them, if any
  # add the constant value of those features in the returned synthetic examples
  
  ConstFeat <- which(apply(dat, 2, function(col){length(unique(col)) == 1}))
  
  tgtOld <- tgt
  
  if(length(ConstFeat)){
    badds <- dat
    ConstRes <- dat[1,ConstFeat]
    dat <- dat[,apply(dat, 2, function(col) { length(unique(col)) > 1 })]
    tgt <- ncol(dat)
  }
  
  nomatr <- c()
  T <- matrix(nrow = dim(dat)[1], ncol = dim(dat)[2])
  for (col in seq.int(dim(T)[2])){
    if (class(dat[, col]) %in% c('factor', 'character')) {
      T[, col] <- as.integer(dat[, col])
      nomatr <- c(nomatr, col)
    } else {
      T[, col] <- dat[, col]
    }
  }
  nC <- dim(T)[2]
  nT <- dim(T)[1]
  

  # ranges evaluated only inside examples to oversample
  ranges <- rep(1, nC)
  if (length(nomatr)) {
    for (x in (1:nC)[-c(nomatr)]) {
      ranges[x] <- max(T[, x]) - min(T[, x])
    }
  } else {
    for(x in (1:nC)) {
      ranges[x] <- max(T[, x]) - min(T[, x])
    }
  }
  

  # use tgtOld because orig data still has all constFeat columns
  kNNs <- neighbours(tgtOld, orig, dist, p, k)

  # for each example in the partition under consideration
  # determine nr of examples out of the partition in the KNN
  r <- c()
  if(Fr){
    for(ex in indpos){
      r <- c(r, length(which(!(kNNs[ex,] %in% (indpos))))/k)
    }
  } else {
    for(ex in indpos){
      r <- c(r, length(which((kNNs[ex,] %in% (indpos))))/k)
    }
    
  }
    
  # normalize
  new.r <- r/sum(r)
  
  # case vector r is always zero!
  if (any(is.na(new.r))){
    r <- rep(1, nrow(dat))
    new.r <- r/sum(r)
  }
      
  nexs <- as.integer(N - 1) # nr of examples to generate for each rare case
  extra <- as.integer(nT * (N - 1 - nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  newM <- matrix(nrow = nexs * nT + extra, ncol = nC)    # the new cases

  NrSynth <- nexs*nT+extra


  # for each ex generate g_i new examples
  g <- round(new.r*NrSynth, 0)
    
  # correct g for matching NrSynth if necessary
  if(sum(g) > NrSynth){
    rem <- sum(g) - NrSynth
    if(all(g != 0)){
      sp <- sample(1:nrow(dat), rem)
      g[sp] <- g[sp] - 1
    } else {
      g.zero <- which(g == 0)
      sp <- sample(setdiff(1:nrow(dat), g.zero), rem)
      g[sp] <- g[sp] - 1
    }
  }
  if(sum(g)< NrSynth){
    add <- NrSynth - sum(g)
    sp <- sample(1:nrow(dat), add)
    g[sp] <- g[sp]+1
  }
  
  if (nrow(dat) > k){
    kNNsLoc <- neighbours(tgt, dat, dist, p, k)
  } else if (nrow(dat)>=2){
    kNNsLoc <- neighbours(tgt, dat, dist, p, (nrow(dat)-1))
    k <- nrow(dat)-1
  } else if(nrow(dat)<2){
    stop("Unable to determine neighbors using only one example!")
  }

  count <- 1
  for(i in 1:nrow(dat)){
    if(g[i] != 0){
      for(n in 1:g[i]){
          # select randomly one of the k NNs
        neig <- sample(1:k, 1)
        # the attribute values of the generated case
        difs <- T[kNNsLoc[i, neig], -tgt] - T[i, -tgt]
        newM[count, -tgt] <- T[i, -tgt] + runif(1) * difs
        for (a in nomatr) {
          # nominal attributes are randomly selected among the existing
          # values of seed and the selected neighbour 
          newM[count, a] <- c(T[kNNsLoc[i, neig], a],
                              T[i, a])[1 + round(runif(1), 0)]
        }
        # now the target value (weighted (by inverse distance) average)
        d1 <- d2 <- 0
        for (x in (1:nC)[-c(nomatr, tgt)]) {
          d1 <- abs(T[i, x] - newM[count, x])/ranges[x]
          d2 <- abs(T[kNNsLoc[i, neig], x] - newM[count, x])/ranges[x]
        }
        if (length(nomatr)) {
          d1 <- d1 + sum(T[i, nomatr] != newM[count, nomatr])
          d2 <- d2 + sum(T[kNNsLoc[i, neig], nomatr] != newM[count, nomatr])
        }
        # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
        if (d1 == d2) {
          newM[count, tgt] <- (T[i, tgt] + T[kNNsLoc[i, neig], tgt])/2
        } else {
          newM[count, tgt] <- (d2 * T[i, tgt] + 
                               d1 * T[kNNsLoc[i, neig], tgt])/(d1 + d2)
        }
        count <- count + 1
      }
    }
  }
  
  newCases <- data.frame(newM)
  for (a in nomatr) {
    newCases[, a] <- factor(newCases[, a],
                            levels = 1:nlevels(dat[, a]),
                            labels = levels(dat[, a]))
  }
  
  if(length(ConstFeat)){ # add constant features that were removed in the beginning
    
    newCases <- cbind(newCases, 
                      as.data.frame(lapply(ConstRes,
                                           function(x){rep(x, nrow(newCases))})))
    colnames(newCases) <- c(colnames(dat), names(ConstFeat))
    newCases <- newCases[colnames(badds)]
    
  } else {
    colnames(newCases) <- colnames(dat)
  }
  
  newCases
  
}

