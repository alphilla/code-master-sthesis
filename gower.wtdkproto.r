rm(list= ls())
ls() #character(0)
library(StatMatch) #gower.dist
library(questionr) #wtd.table

gower.wtdkproto <- function (x, ...)
  UseMethod("gowerwtdkproto")

gowerwtdkproto.default<- function(x, k, iter.max = 100, nstart = 1, na.rm = TRUE, keep.data = TRUE, verbose = TRUE, w, ...){

  ############## initial error checks ##############

  if(!is.data.frame(x)) stop("x should be a data frame!")
  if(ncol(x) < 2) stop("For clustering x should contain at least two variables!")
  if(iter.max < 1 | nstart < 1) stop("iter.max and nstart must not be specified < 1!")
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  if(!anynum) stop("\n No numeric variables in x! Try using kmodes() from package klaR...\n\n")
  if(!anyfact) stop("\n No factor variables in x! Try using kmeans()...\n\n")

  ############## treatment of missings ###############

  # NAcount <- apply(x, 2, function(z) sum(is.na(z)))
  # if(verbose){
  #   cat("# NAs in variables:\n")
  #   print(NAcount)
  # }
  # if(any(NAcount == nrow(x))) stop(paste("Variable(s) have only NAs please remove them:",names(NAcount)[NAcount == nrow(x)],"!"))
  # if(na.rm) {
  #   miss <- apply(x, 1, function(z) any(is.na(z)))
  #   if(verbose){
  #     cat(sum(miss), "observation(s) with NAs.\n")
  #     if(sum(miss) > 0) message("Observations with NAs are removed.\n")
  #     cat("\n")
  #   }
  #   x <- x[!miss,]
  # }

  ############ remove missings ################

  # if(!na.rm){
  #   allNAs <- apply(x,1,function(z) all(is.na(z)))
  #   if(sum(allNAs) > 0){
  #     if(verbose) cat(sum(allNAs), "observation(s) where all variables NA.\n")
  #     warning("No meaningful cluster assignment possible for observations where all variables NA.\n")
  #     if(verbose) cat("\n")
  #
  #   }
  # }

  if(nrow(x) == 1) stop("Only one observation clustering not meaningful.")

  k_input <- k # store input k for nstart > 1 as clusters can be merged

  ############# initialize prototypes ##################

  if(!is.data.frame(k)){
    if (length(k) == 1){
      if(as.integer(k) != k){k <- as.integer(k); warning(paste("k has been set to", k,"!"))}
      if(sum(complete.cases(x)) < k) stop("Data frame has less complete observations than clusters!")
      ids <- sample(row.names(x[complete.cases(x),]), k)
      if(nrow(x) < k) stop("Data frame has less observations than clusters!")
      ids <- sample(nrow(x), k)
      protos <- x[ids,]
    }
    if (length(k) > 1){
      if(nrow(x) < length(k)) stop("Data frame has less observations than clusters!")
      ids <- k
      k <- length(ids)
      if(length(unique(ids)) != length(ids)) stop("If k is specified as a vector it should contain different indices!")
      if(any(ids<1)|any(ids>nrow(x))) stop("If k is specified as a vector all elements must be valid indices of x!")
      #check for integer
      protos <- x[ids,]
      if(any(!complete.cases(protos))) stop("Choose initial prototypes without missing values!")
    }
    rm(ids)
  }
  if(is.data.frame(k)){
    if(nrow(x) < nrow(k)) stop("Data frame has less observations than clusters!")
    if(length(names(k)) != length(names(x))) stop("k and x have different numbers of columns!")
    if(any(names(k) != names(x))) stop("k and x have different column names!")
    if(anynum) {if( any(sapply(k, is.numeric) != numvars)) stop("Numeric variables of k and x do not match!")}
    if(anyfact) {if( any(sapply(k, is.factor) != catvars)) stop("Factor variables of k and x do not match!")}
    protos <- k
    if(any(!complete.cases(protos))) stop("Prototypes with missing values. Choose initial prototypes without missing values!")
    k <- nrow(protos)
  }
  if(k < 1) stop("Number of clusters k must not be smaller than 1!")



  ############## initialize clusters ################

  clusters  <- numeric(nrow(x))
  tot.dists <- NULL
  moved   <- NULL
  iter <- 1

  ## check for any equal prototypes and reduce cluster number in case of occurence
  if(k > 1){
    keep.protos <- rep(TRUE,k)
    for(l in 1:(k-1)){
      for(m in (l+1):k){
        d1 <- as.numeric(gower.dist(protos[l,],protos[m,])) # euclidean for numerics

        if(d1 == 0) keep.protos[m] <- FALSE
      }
    }
    if(!all(keep.protos)){
      protos <- protos[keep.protos,]
      k <- sum(keep.protos)
      if(verbose) message("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")
    }
  }

  # special case only one cluster
  if(k == 1){clusters <- rep(1, nrow(x)); size  <- table(clusters); iter <- iter.max} # REM: named vector size is needed later...

  # start iterations for standard case (i.e. k > 1)
  while(iter < iter.max){

    ########## compute distances #################
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      dists[,i]<- as.matrix(gower.dist(x,protos[i,]))
    }

    #############  assign clusters ##############

    old.clusters  <- clusters
    # clusters      <- apply(dists, 1, function(z) which.min(z))
    clusters      <- apply(dists, 1, function(z) {a <- which(z == min(z)); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
    size          <- table(clusters)
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
    # prevent from empty classes
    #tot.within    <- numeric(k)
    #totw.list     <- by(min.dists, clusters, sum)
    #tot.within[names(totw.list)] <- as.numeric(totw.list)

    ######## ...check for empty clusters and eventually reduce number of prototypes ###########

    if (length(size) < k){
      k <- length(size)
      protos <- protos[1:length(size),]
      if(verbose) cat("Empty clusters occur. Cluster number reduced to:", k, "\n\n")
    }

    ####### trace #######

    tot.dists <- c(tot.dists, sum(tot.within))
    moved <- c(moved, sum(clusters != old.clusters))

    ####### weighted mean for numerical variables #######

    wtd.mean<-function(z){
        weighted_columns <- weighted.mean(z, w=w)
      return(weighted_columns)
    }

    ####### weighted frequencies for categorical variables #######
    #This function returns the most frequent factor
    wtd.freq<-function(z){
      weighted_mode<- levels(z)[which.max(wtd.table(z, weights=w))]
      return(weighted_mode)
    }

    ####### compute new prototypes #######

    remids <- as.integer(names(size))
    for(i in remids){

      protos[which(remids == i), numvars]<- sapply(x[clusters==i, numvars, drop = FALSE], wtd.mean)
      # protos[which(remids == i), numvars] <- sapply(x[clusters==i, numvars, drop = FALSE], mean, na.rm = TRUE)

      protos[which(remids == i), catvars] <-sapply(x[clusters==i, catvars, drop = FALSE], wtd.freq)
      # protos[which(remids == i), catvars] <- sapply(x[clusters==i, catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
    }

    if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}

    # check for any equal prototypes and reduce cluster number in case of occurence
    if(iter == (iter.max-1)){ # REM: for last iteration equal prototypes are allowed. otherwise less prototypes than assigned clusters.
      keep.protos <- rep(TRUE,k)
      for(l in 1:(k-1)){
        for(m in (l+1):k){
          d1 <- as.numeric(gower.dist(protos[l,],protos[m,]))
          if(d1 == 0) keep.protos[m] <- FALSE
        }
      }
      if(!all(keep.protos)){
        protos <- protos[keep.protos,]
        k <- sum(keep.protos)
        if(verbose) cat("Equal prototypes merged. Cluster number reduced to:", k, "\n\n")
      }
    }
    ######## add stopping rules  ############

    if(moved[length(moved)] ==  0) break

    if(k == 1){clusters <- rep(1, length(clusters)); size <- table(clusters); iter <- iter.max; break}

    #cat("iter", iter, "moved", moved[length(moved)], "tot.dists",tot.dists[length(tot.dists)],"\n" )
    iter <- iter+1
  }


  ########### Final update of prototypes and dists ###########

  if(iter == iter.max){ # otherwise there have been no moves anymore and prototypes correspond to cluster assignments
    # compute new prototypes
    remids <- as.integer(names(size))
    for(i in remids){

      protos[which(remids == i), numvars]<- sapply(x[clusters==i, numvars, drop = FALSE], wtd.mean)
      protos[which(remids == i), catvars] <-sapply(x[clusters==i, catvars, drop = FALSE], wtd.freq)
      # protos[which(remids == i), numvars] <- sapply(x[clusters==i, numvars, drop = FALSE], mean, na.rm = TRUE)
      # protos[which(remids == i), catvars] <- sapply(x[clusters==i, catvars, drop = FALSE], function(z) levels(z)[which.max(table(z))])
    }

    ########## compute distances #############
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      dists[,i]<- as.matrix(gower.dist(x,protos[i,]))
    }

    size          <- table(clusters)
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
  }


  names(clusters) <- row.names(dists) <- row.names(x)
  rownames(protos) <- NULL
  # create result:
  res <- list(cluster = clusters,
              centers = protos,
              size = size,
              withinss = within,
              tot.withinss = tot.within,
              dists = dists,
              iter = iter,
              trace = list(tot.dists = tot.dists, moved = moved))

  # loop: if nstart > 1:
  if(nstart > 1)
    for(j in 2:nstart){
      res.new <- gower.kproto(x=x, k=k_input, iter.max = iter.max, nstart=1, verbose=verbose, na.rm = na.rm)
      if(res.new$tot.withinss < res$tot.withinss) res <- res.new
    }

  if(keep.data) res$data = x
  class(res) <- "gowerwtdkproto"
  return(res)
}


summary.kproto <- function(object, data = NULL, pct.dig = 3, ...){
  if(class(object) != "gowerwtdkproto") stop("object must be of class gowerwtdkproto!")

  if(is.null(data)){
    data <- object$data
    cluster <- object$cluster
  }
  if(!is.null(data)) cluster <- predict(object, data)$cluster

  numvars <- sapply(data, is.numeric)
  #anynum <- any(numvars)
  catvars <- sapply(data, is.factor)
  #anyfact <- any(catvars)

  res <- NULL
  for(i in 1:ncol(data)){
    cat(names(data)[i],"\n")
    if(numvars[i]){
      resi <- by(data[,i], cluster, summary, ...)
      res[[i]] <- matrix(unlist(resi), nrow = length(unique(cluster)), byrow=TRUE)
      colnames(res[[i]]) <- names(resi[[1]])
      rownames(res[[i]]) <- sort(unique(cluster))
    }
    if(catvars[i])  res[[i]] <- round(prop.table(table(cluster, data[,i]),1), digits = pct.dig)
    print(res[[i]])
    cat("\n-----------------------------------------------------------------\n")
  }
  names(res) <- names(data)

  #return(res)
  invisible(res)
}

print.kproto <- function(x, ...){
  cat("Numeric predictors:", sum(sapply(x$centers, is.numeric)), "\n")
  cat("Categorical predictors:", sum(sapply(x$centers, is.factor)), "\n")
  cat("Number of Clusters:", length(x$size), "\n")
  cat("Cluster sizes:", x$size, "\n")
  cat("Within cluster error:", x$withinss, "\n\n")

  cat("Cluster prototypes:\n")
  print(x$centers)
}
