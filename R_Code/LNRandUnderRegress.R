# ===================================================
# Performs a random under-sampling strategy for regression problems
# with a bias towards the frontier or towards the safe cases.
# ---------------------------------------------------
LNRandUnderRegress <- function(form, dat, rel = "auto", thr.rel = 0.5,
                             C.perc = "balance", repl = FALSE, Fr=TRUE, dist="HEOM", p=2, k=5)
  # Args:
  # form a model formula
  # dat the original training set (with the unbalanced distribution)
  # rel is relevance determined automatically (default) 
  #     or provided by the user
  # thr.rel is the relevance threshold above which a case is considered
  #         as belonging to the rare "class"
  # C.perc is a list containing the under-sampling percentage/s to apply 
  #        to all/each "class" obtained with the relevance threshold. This
  #        percentage represents the percentage of examples that is maintained
  #        in each "class". Examples are randomly removed in each "class".
  #        Moreover, different percentages may be provided for each "class". 
  #        Alternatively, it may be "balance" or "extreme", cases 
  #        where the under-sampling percentages are automatically estimated.
  # repl is it allowed to perform sampling with replacement.
  # fr  logical. Is the frontier reinforced (TRUE) or are the safe cases(FALSE).
  # dist  character specifying the distance metric to use when computing neighbors
  # p     numeric required when a p-norm is computed
  # k    numeric. Nr. of neighbors to compute 
  # Returns: a data frame modified by the random under-sampling strategy
  #           with neighborhood bias

{
  if(is.list(C.perc) & any(unlist(C.perc)>1)){
    stop("The under-sampling percentages provided in parameter C.perc
         can not be higher than 1!")
  }
  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))
  
  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)
  s.y <- sort(y)
  
  if (is.matrix(rel)) { 
    pc <- phi.control(y, method = "range", control.pts = rel)
  } else if (is.list(rel)) { 
    pc <- rel
  } else if (rel == "auto") {
    pc <- phi.control(y, method = "extremes")
  } else {
    stop("future work!")
  }
  
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
  
  
  imp <- sapply(obs.ind, function(x) mean(phi(x, pc)))
  
  und <- which(imp < thr.rel)
  ove <- which(imp >= thr.rel)
  
  newdata <- NULL
  for (j in 1:length(ove)) {
    newdata <- rbind(newdata, dat[names(obs.ind[[ove[j]]]), ])
    # start with the examples from the minority "classes"
  }
  
  # set the under-sampling percentages
  if (is.list(C.perc)) {
    if (length(und) > 1 & length(C.perc) == 1) { 
      # the same under-sampling percentage is applied to all the "classes"
      C.perc <- rep(C.perc[1], length(und))
    } else if (length(und) > length(C.perc) & length(C.perc) > 1) {
      stop("The number of under-sampling percentages must be equal 
           to the number of bumps below the threshold defined!")      
    }else if (length(und) < length(C.perc)) {
      stop("the number of under-sampling percentages must be at most
           the number of bumps below the threshold defined!")
    }
    } else if (C.perc == "balance") {
      B <- sum(sapply(obs.ind[ove], length))
      obj <- B/length(und)
      C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 5))
  } else if (C.perc == "extreme") {
    Bove <- sum(sapply(obs.ind[ove], length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length), 5))
  }
  
  for (j in 1:length(und)) {
    kNNs <- neighbours(tgt, dat, dist, p, k)
    ind <- names(obs.ind[[und[j]]])
    indpos <- match(ind, rownames(dat))
    if(Fr){ # reinforce the frontier cases
      
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
      # this reinforces the frontier cases
      sel <- sample(indpos,
                    as.integer(C.perc[[j]] * length(obs.ind[[und[j]]])),
                    replace = repl,
                    prob=new.r)
      newdata <- rbind(newdata, dat[sel, ])
      
    } else { # reinforce the far from frontier cases (the safe cases)

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
      # this reinforces the far from frontier cases
      sel <- sample(indpos,
                    as.integer(C.perc[[j]] * length(obs.ind[[und[j]]])),
                    replace = repl,
                    prob=new.r)
      newdata <- rbind(newdata, dat[sel, ])
      
    }
  }
  newdata
  }

