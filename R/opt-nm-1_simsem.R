#' Search for the optimal missing pattern with only one missing measured variable. An internal function for forward assembly.
#' 
#'\code{opt1.simsem} is an internal function that runs simulations using \code{lavaan} and \code{simsem}. It returns the optimal missing pattern that only contains one missing measured variable. This is the first step of forward assembly.
#'
#' @inheritParams balance.miss.l
#' 
#' @return An object containing the information of the optimal missing
#'    data pattern containing only one missing observed variable. The
#'    optimal pattern is the one that yields highest statistical power for
#'    testing the focal parameters, compared to other patterns with only
#'    one missing observed variable.
#' 
#' @seealso \code{\link{simPM}} which is a warpper function for this
#'    function.
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export opt1.simsem
#' @examples


opt1.simsem <- function(
  popModel,
  analyzeModel,
  NAMES,    # a vactor containing the measured variable names
  distal.var,
  n,
  nreps,
  seed,
  Time,      # total planned time points
  k,         # original number of measured variables at each time point
  Time.complete,      # time points with completed data
  costmx,             # a vector containing the cost of the remaining measured variables acorss time points
  pc,                 # proportion of completers
  pd,                 # proportion of droppers
  focal.param,        # user identified focal parameters, in a matrix form. User needs to identify the rows in the readModel parameter object that are of focal interest
  complete.var=NULL
#  multicore=T
  ) { # a list of the variables that need to be complete

  num.miss <- 1   #this function only deals with one missing slot

  VNAMES <- NAMES

  # calculated info based on user supplied info
  future.k <- (Time - Time.complete) * k   # data points not yet completed, also maximum possible # of missing
  ms.range <- c((Time.complete * k + 1):(Time * k))  # The available time slots to plant missingness

  ms.combn <- combn(ms.range, num.miss) #all possible combinations of missing slots given num.miss
  all.pattern <- matrix(0, nrow = choose(future.k, num.miss), ncol = length(VNAMES)) # place holder, all possible patterns with a certain ms

  # update the missing patterns with 1s
  for (q in seq_len(nrow(all.pattern))) {
    all.pattern[q, ms.combn[, q]] = 1
  }                                         # add missingness
  completers <- rep(0, ncol(all.pattern))       # completers pattern
  dropper <- c(rep(0, Time.complete * k), rep(1, future.k))   #droppers pattern

  #### if user specify certain variables to be complete

  if (length(complete.var)==1) {
    # find which column these variables are
    complete.cols <- which(VNAMES %in% complete.var)
    # whether to keep the rows
    keep <- all.pattern[, complete.cols]==0

    all.pattern <- all.pattern[keep, , drop=FALSE]
    # take out these designs columns out of the ms.combn
    ms.combn <- ms.combn[, keep, drop=FALSE]
  }

  if (length(complete.var)>1) {
    complete.cols <- which(VNAMES %in% complete.var)
    keep <- rowSums(all.pattern[, complete.cols]==0)==length(complete.var)

    temp.pattern <- all.pattern[keep, , drop = FALSE]
    if (is.null(dim(temp.pattern))==F) {
      all.pattern <- temp.pattern
    }
    if (is.null(dim(temp.pattern))==T) {
      all.pattern <- t(as.matrix(temp.pattern))
    }

    temp.combs <- ms.combn[ , keep, drop = FALSE]
    if (is.null(dim(temp.combs))==F) {
      ms.combn <- temp.combs
    }
    if (is.null(dim(temp.combs))==T) {
      ms.combn <- as.matrix(temp.combs)
    }
  }

  ### storage bins for simulation results
  convergence.rate <- rep(NA, nrow(all.pattern))   #convergence rate
  weakest.param.name <- rep(NA, nrow(all.pattern))
  weakest.para.power <- rep(NA, nrow(all.pattern))
  cost.design <- rep(NA, nrow(all.pattern)) # cost of each design
  miss.num <- rep(num.miss, nrow(all.pattern))
  miss.name <- matrix(NA, nrow(all.pattern), num.miss)
  sim.seq <- rep(NA, nrow(all.pattern))
  miss.loc <- matrix(NA, nrow(all.pattern), num.miss)
  sim.out <- vector("list", nrow(all.pattern))   # simulation output storage


  for (i in seq_len(nrow(all.pattern))) {
    # because num.miss=1, no previously selected pattern applicable
    if(pd!=0) {
      patmx <- rbind(all.pattern[i, ], completers, dropper)   # missing patterns
    }
    if(pd==0) {
      patmx <- rbind(all.pattern[i, ], completers)
    }

    #FNAME=paste0("missing-",num.miss,"-sim.seq-",VNAMES[ms.combn[,i]])  #file name
    #FNAME=paste0("missing-",num.miss,"-sim.seq-",i)  #file name

    ### distal variables
    if (is.null(distal.var)==F) {
      dis.pat <- matrix(0, nrow = nrow(patmx), ncol = length(distal.var))
      patmx <- cbind(patmx, dis.pat)
      VNAMES <- c(VNAMES, distal.var)
    }

    # pattern probs

    if (pd==0) {
      p.probs <- c(rep(round((1-pc-pd)/(nrow(patmx)-1),6), nrow(patmx)-1), pc)
    }
    if (pd!=0) {
      p.probs <- c(rep(round((1-pc-pd)/(nrow(patmx)-2),6), nrow(patmx)-2), pc, pd)
    }

    # ns in each pattern
    pn <- rep(0, length(p.probs))
    for (pp in 1:(length(pn)-1)) {
      pn[pp] <- floor(n*p.probs[pp])
    }
    pn[length(pn)] <- n-sum(pn[1:(length(pn)-1)])

    logical.mx <- matrix(0, nrow = n, ncol = ncol(patmx))
    logical.mx[1:pn[1],] <- patmx[rep(1,pn[1]), ]

    if (length(pn)>=3) {
    for (pi in 2:(length(pn)-1)) {
      logical.mx[(sum(pn[1:(pi-1)])+1):sum(pn[1:pi]),] <- patmx[rep(pi,pn[pi]), ]
    }
    }

    # rest are completers

    # need to make sure the order of the variable is consistent
    colnames(logical.Mx) <- colnames(logical.mx) <- VNAMES
    logical.Matrix[[d]] <- logical.Mx
    
    get_data <- simsem::sim(nRep = 1, 
                            model = analyzeModel, 
                            generate = popModel, 
                            n = 2, 
                            dataOnly = T)
    
    # important! Need the variable to be ordered the same as the generated data!
    logical.Mx.reorder <- logical.Mx[, colnames(get_data[[1]])]
    
    misstemplate <- miss(logical = logical.Mx.reorder, m=0)
    output <- simsem::sim(nreps, 
                          n = n, 
                          model = analyzeModel, 
                          generate = popModel, 
                          miss = misstemplate,
                          seed = seed + i
                          #, multicore=multicore
                          )
    sim.out[[i]] <- output    #save the simulation output

    sim.param <- summaryParam(output)
    name.param <- rownames(summaryParam(output))
    converged <- output@converged == 0

    if (sum(converged)==0) {
      convergence.rate[i] <- 0
      weakest.param.name[i] <- "NA"
      weakest.para.power[i] <- 0
    }
    if (sum(converged)>0) {
      f.param <- sim.param[name.param %in% focal.param, ]
      weakest.f.param <- f.param[f.param$`Power (Not equal 0)`==min(f.param$`Power (Not equal 0)`), ]
      if (nrow(weakest.f.param)>1) {
        weakest.f.param <- weakest.f.param[1, ]      ########## may need to be changed later
      }
      convergence.rate[i] <- sum(converged)/nreps    #converged number of simulations
      weakest.param.name[i] <- rownames(weakest.f.param)
      weakest.para.power[i] <- weakest.f.param$`Power (Not equal 0)`
    }

      if (pd==0) {
        cost.design[i] <- sum(c((1-pc)*n,pc*n)*((1-patmx[,ms.range])%*%costmx))   #patmx depends on i
      }
      if (pd!=0) {
        cost.design[i] <- sum(c((1-pc-pd)*n,pc*n,pd*n)*((1-patmx[, ms.range]) %*% costmx))
      }

      miss.name[i] <- VNAMES[ms.combn[,i]]
      sim.seq[i] <- i  # location as specified in the miss.combn matrix
      miss.loc[i,] <- ms.combn[,i]
    }



  ### combine the results
  sim.results.out <- cbind.data.frame(convergence.rate,   #convergence rate
                                   weakest.param.name,
                                   weakest.para.power,
                                   cost.design, # cost of each design
                                   miss.num,
                                   miss.name,
                                   sim.seq,
                                   miss.loc)

  opt.design.1 <- sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]), ]

  if (nrow(opt.design.1)==1) {
    opt.design <- opt.design.1
  }

  if (nrow(opt.design.1)>1) {
    n.min.cost <- nrow(opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ])
    if (n.min.cost==1) {
      opt.design <- opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ]
    } else {
      opt.min.cost <- opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ]

      ### only applies to num.miss=1
      opt.design <- opt.min.cost[opt.min.cost[,"miss.loc"]==max(opt.min.cost[,"miss.loc"]), ]
    }
  }

  op <- opt.design[, "sim.seq"]

  opt.output <- sim.out[[op]]

  if (pd==0) {
    opt.pattern <- rbind(all.pattern[op, ], completers)
    opt.probs <- c(rep(round((1-pc)/(nrow(patmx)-1),6), nrow(patmx)-1), pc)
  }
  if (pd!=0) {
    opt.pattern <- rbind(all.pattern[op, ], completers, dropper)
    opt.probs <- c(rep(round((1-pc-pd)/(nrow(patmx)-2), 6), nrow(patmx)-2), pc, pd)
  }
  colnames(opt.pattern) <- NAMES

  # ns in each pattern
  pn <- rep(0,length(opt.probs))
  for (pp in 1:(length(pn)-1)){
    pn[pp] <- floor(n*opt.probs[pp])
  }
  pn[length(pn)] <- n-sum(pn[1:(length(pn)-1)])

  misc <- list(time = Time, k = k, focal.param = focal.param)
  
  re.ob <- list("results" = sim.results.out,
             "opt.design" = opt.design,
             "opt.pattern" = opt.pattern,
             "opt.probs" = opt.probs,
             "opt.ns" = pn,
             "design.order" = op,
             "opt.output" = opt.output,
             "misc" = misc)

  class(re.ob) <- append(class(re.ob),"simpm")
  
  return(re.ob)

}


