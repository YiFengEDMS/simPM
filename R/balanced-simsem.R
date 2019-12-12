
#' Search for the optimal balanced item-level planned missing designs.
#' 
#'  \code{balance.miss.l} runs simulations using lavaan and simsem. It returns the search results for optimal balanced item-level PM designs. 
#'  
#' @param popModel The data generation model (population model) specified #'    using lavaan script.
#' @param analyzeModel The analysis model, specified using lavaan script. #'    The analysis model can be different from the population model. 
#' @param NAMES A character vector containing the names of the observed
#'    variables. The variable names must be ordered chronologically, by
#'    the time (wave) they are measured.
#' @param Time Numeric. The total number of time points (or total number
#'    of waves of data collection).
#' @param Time.complete Numeric. Number of waves of data collection that
#'    have been completed before the funding cut occurs.
#' @param k Numeric. The number of observed variables collected at each
#'    wave.
#' @param pc Numeric. Proportion of completers: the proportion of subjects
#'    that will participate in all of the following waves of data
#'    collection and provide complete data. This must be greater than 0.
#' @param pd Numeric. The proportion of subjects that will not participate
#'    in any of the following waves of data collection (i.e., people who
#'    will drop from the longitudinal study). This value can be 0.
#' @param costmx  A numeric vector containing the unit cost of each
#'    observed variable that is yet to be measured (post the funding cut).
#'    The cost is assumed to be constant across subjects, but it is
#'    allowed to vary across variables and across waves.
#' @param n The total sample size as initially planned.
#' @param nreps Number of replications for Monte Carlo simulations.
#' @param focal.param Character vector. The parameters of focal interest.
#'    The focal parameters should be specified in the format of 
#'    \code{lavaan} script.
#' @param complete.var Char vector. Specify the name(s) of the 
#'    variable(s) if there are any variable(s) that need to have complete
#'    data collected from all the participating subjects.
#' @param eval.budget Logical scalar (TRUE or FALSE), indicating whether
#'    there is any budget constraint. If the user wishes to search for PM
#'    designs under the budget limit, they need to specify the amount of
#'    the remaining available budget that can be used for future data
#'    collection.
#' @param rm.budget Numeric. The amount of remaining budget avaialbe for
#'    future data collection. User must supply a value for this argument
#'    if \code{eval.budget = TRUE}. 
#' @param distal.var Char vector. Specify the name(s) of the distal
#'    variables. User needs to specify this argument if there are any 
#'    time-independent distal variables included in the model that are not
#'    subject to planned missingness.
#' @param seed The random seed for random number generation.
#' 
#' @return An object containing the information of the optimal balanced
#'    item-level missing design. The optimal design is the one that yields
#'    highest power for testing the focal parameters, compared to other
#'    plausible candidate PM designs.
#' 
#' @seealso \code{\link{simPM}} which is a warpper function for this
#'    function.
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export balance.miss.l
#' @examples


balance.miss.l <- function(
  popModel,
  analyzeModel,
  NAMES,
  Time,
  Time.complete,
  k,
  pc,
  pd,
  costmx,
  n,
  nreps,
  focal.param,
  complete.var=NULL,
  eval.budget=T,
  rm.budget=NULL,
  distal.var=NULL,
  seed=1234
 # multicore=T
  ) {

  n.miss.point <- 1:(Time*k - Time.complete * k - 1)
  ms.range <- c((Time.complete * k + 1):(Time * k))

  # storage bins
  designs <- list()
  probs <- list()
  pat.n <- list()
  #num.miss.wave=c()
  cost.design <- c()                   # cost of each design


  rs <- 1

  for (i in n.miss.point) {             # loop over different # of missing waves (designs)

    mpoint <- combn(ms.range, i)           # possible combinations of missing points when n.miss.point=i
    pattern <- matrix(0, nrow = ncol(mpoint), ncol = k * Time)  #pattern matrix

    for (j in seq_len(nrow(pattern))) {
      for (m in seq_len(nrow(mpoint)))
        pattern[j, mpoint[m, j]]=1   #put missing in pattern matrix
    }

    completers <- rep(0, ncol(pattern))       # completers pattern
    dropper <- c(rep(0, Time.complete * k),rep(1, (Time - Time.complete) * k))   #droppers pattern

    # when user needs to have complete data on certain variables

    if (is.null(complete.var)==F) {

      if (length(complete.var)==1) {
        # find which column these variables are
        complete.cols <- which(VNAMES %in% complete.var)
        # whether to keep the rows
        keep <- pattern[, complete.cols]==0

        pattern <- pattern[keep, , drop=FALSE]
      }

      if (length(complete.var) > 1) {
        complete.cols <- which(VNAMES %in% complete.var)
        keep <- rowSums(pattern[ , complete.cols]==0)==length(complete.var)

        temp.pattern <- pattern[keep, , drop=FALSE]
        if (is.null(dim(temp.pattern))==F) {
          pattern <- temp.pattern
        }
        if (is.null(dim(temp.pattern))==T) {
          pattern <- t(as.matrix(temp.pattern))
        }

      }

    }

    if(pd!=0) {      # if there are droppers
      patmx <- rbind(pattern, completers, dropper)   # missing patterns for Mplus later
    }
    if(pd==0) {      # if there are no droppers
      patmx <- rbind(pattern,completers)
    }

    designs[[rs]]=patmx

    #### pattern probs
    # p.probs

    if (pd==0) {
      p.probs <- c(rep(round((1 - pc - pd) / (nrow(patmx) - 1), 6), nrow(patmx) - 1), pc)
    }
    if (pd!=0) {
      p.probs <- c(rep(round((1 - pc - pd) / (nrow(patmx) - 2), 6), nrow(patmx) - 2), pc, pd)
    }

    probs[[rs]] <- p.probs


    # ns in each pattern
    pn <- rep(0,length(p.probs))
    for (pp in 1:(length(pn)-1)){
      pn[pp] <- floor(n * p.probs[pp])
    }
    pn[length(pn)] <- n-sum(pn[1:(length(pn)-1)])

    pat.n[[rs]] <- pn

    ### cost of each design
    cost.design[rs] <- sum((1-patmx[,ms.range])%*%costmx*pn)

    rs <- rs+1
  }
  #### evaluate cost
  # only simulate for those designs that are below the budget limit

  if (eval.budget==T) {
    if (sum(cost.design>rm.budget)==length(cost.design)) {
      stop ("all wave missing designs cost more than the avaiable remaing budget. Try other designs.")
    }

    designs2 <- designs[cost.design<=rm.budget] #select the designs that are below the budget limit
    probs2 <- probs[cost.design <= rm.budget]
    miss.point <- n.miss.point[cost.design <= rm.budget]
    cost.design2 <- cost.design[cost.design <= rm.budget]
    pat.n2 <- pat.n[cost.design <= rm.budget]
  }

  if (eval.budget==F) {
    designs2 <- designs
    probs2 <- probs
    miss.point <- n.miss.point
    cost.design2 <- cost.design
    pat.n2 <- pat.n
  }

  convergence.rate <- c()   #convergence rate
  weakest.param.name <- c()
  weakest.para.power <- c()
  template <- list()
  logical.Matrix <- list()
  sim.out <- list()

  for (d in seq_len(length(designs2))) {

    patmx <- designs2[[d]]
    p.probs <- probs2[[d]]
    pn <- pat.n2[[d]]
    VNAMES <- NAMES

    if (prod(pn>=1)==0) {
      convergence.rate[d] <- NA   #convergence rate
      weakest.param.name[d] <- NA
      weakest.para.power[d] <- -99

      template[d] <- NA
      logical.Matrix[d] <- NA
      sim.out[d] <- NA
    }
    if (prod(pn >= 1) != 0) {

    ###distal variables
    if (is.null(distal.var)==F) {
      dis.pat <- matrix(0,nrow=nrow(patmx),ncol=length(distal.var))
      patmx <- cbind(patmx,dis.pat)
      VNAMES <- c(VNAMES,distal.var)
    }

    #  FNAME=paste0("missing-waves-",miss.waves[d])  #file name


    logical.mx <- matrix(0,nrow=n, ncol=ncol(patmx))
    #  data <- data[sample(1:nrow(data)),]   #shuffle the rows
    logical.mx[1:pn[1], ] <- patmx[rep(1,pn[1]), ]

    if (length(pn) >= 3) {
    for (pi in 2:(length(pn)-1)) {
      logical.mx[(sum(pn[1:(pi-1)])+1):sum(pn[1:pi]), ] <- patmx[rep(pi, pn[pi]), ]
    }
    }

    # rest are completers

    logical.Mx <- logical.mx==1
    
    # need to make sure the order of the variable is consistent
    colnames(logical.Mx) <- colnames(logical.mx) <- VNAMES
    logical.Matrix[[d]] <- logical.Mx
    
    get_data <- simsem::sim(nRep = 1, 
                            model = analyzeModel, 
                            generate = popModel, 
                            n = 2, 
                            dataOnly = T)
    

    #  miss.model=miss(logical=logical.Mx)
    #  data.miss=impose(miss.model,data)
    
    # important! Need the variable to be ordered the same as the generated data!
    logical.Mx.reorder <- logical.Mx[, colnames(get_data[[1]])]
    
    misstemplate <- miss(logical = logical.Mx.reorder, m=0)
    output <- simsem::sim(nreps, 
                          n = n, 
                          model = analyzeModel, 
                          generate=popModel,
                          miss=misstemplate, 
                          seed = seed + d
                          #,multicore=multicore
                          )

    template[[d]] <- misstemplate
    sim.out[[d]] <- output

    sim.param <- summaryParam(output)   #turn off
    name.param <- rownames(summaryParam(output))   #turn off
    converged <- output@converged == 0

    if (sum(converged)==0) {
      convergence.rate[d] <- 0
      weakest.param.name[d] <- "NA"
      weakest.para.power[d] <- 0
    }
    if (sum(converged)>0) {
      f.param <- sim.param[name.param%in%focal.param,]
      weakest.f.param <- f.param[f.param$`Power (Not equal 0)`==min(f.param$`Power (Not equal 0)`),]
      if (nrow(weakest.f.param)>1) {
        weakest.f.param <- weakest.f.param[1,]      ########## may need to be changed later
      }
      convergence.rate[d] <- sum(converged)/nreps    #converged number of simulations
      weakest.param.name[d] <- rownames(weakest.f.param)
      weakest.para.power[d] <- weakest.f.param$`Power (Not equal 0)`
    }
  }
}


  sim.results.out <- cbind.data.frame(convergence.rate,   #convergence rate
                                   weakest.param.name,
                                   weakest.para.power,
                                   "cost.design"=cost.design2, # cost of each design
                                   miss.point)

  opt.design.1 <- sim.results.out[sim.results.out[, "weakest.para.power"]==max(sim.results.out[, "weakest.para.power"]), ]

  if (nrow(opt.design.1)==1) {
    opt.design <- opt.design.1
  }
  if (nrow(opt.design.1)>1) {
    opt.design <- opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ]
  }

  op <- which(miss.point==opt.design$miss.point)    #which design is chosen
  opt.pattern <- designs2[[op]]
  colnames(opt.pattern) <- NAMES
  opt.probs <- probs2[[op]]
  opt.patns <- pat.n2[[op]]

  opt.template <- template[[op]]
  opt.logical <- logical.Matrix[[op]]
  opt.output <- sim.out[[op]]
  
  misc <- list(time=Time,k=k,focal.param=focal.param)

  re.ob <- list(
    "results"=sim.results.out,
    "opt.design"=opt.design,
    "opt.pattern"=opt.pattern,
    "opt.probs"=opt.probs,
    "opt.ns"=opt.patns,
    "n.miss.point"=opt.design$miss.point,
    "opt.template"=opt.template,
    'opt.logical'=opt.logical,
    'opt.output'=opt.output,
    "misc"=misc)
  
  class(re.ob) <- append(class(re.ob),"simpm")
  
  return(re.ob)
}
