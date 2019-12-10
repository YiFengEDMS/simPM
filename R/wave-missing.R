#' Search for the optimal wave-level PM designs.
#'
#' \code{wave.miss} runs simulations using M\emph{plus}. It returns the search results for optimal wave-level PM designs. 
#'
#' @inheritParams balance.miss
#'
#' @param complete.wave Numeric vector. Specify which wave(s) that the
#'    user wish to have complete data collected from all the participants.
#'    
#' @return An object containing the information of the optimal
#'    wave-level missing design. The optimal design is the one that yields
#'    highest power for testing the focal parameters, compared to other
#'    plausible candidate PM designs.
#'    
#' @seealso \code{\link{simPM}} which is a warpper function for this
#'    function.
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export wave.miss
#' @examples



wave.miss <- function(
  VNAMES,
  distal.var=NULL,
  n,
  nreps,
  seed,
  Time,
  k,
  Time.complete,
  costmx,
  pc,
  pd,
  design0.out,
  focal.param,
  #max.mk,                 # maximum number of missing slots allowed
  eval.budget=T,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
  rm.budget=NULL,          # remaining available budget
  complete.wave=NULL
) {

  # number of avaialbe waves for PM
  NAMES <- VNAMES
  n.miss.waves <- 1:(Time-Time.complete-1)  # possible number of waves missing
  ms.range <- c((Time.complete*k+1):(Time*k))

  # storage bins
  designs <- vector(mode = "list", length = n.miss.waves)
  probs <- vector(mode = "list", length = n.miss.waves)
  #num.miss.wave=c()
  cost.design <- rep(NA, n.miss.waves) # cost of each design

  rs <- 1

  for (i in n.miss.waves) {               # loop over different # of missing waves (designs)

    if (Time.complete==0) {
      mwave <- combn(c(1:Time), i)
    }
    if (Time.complete>0) {
      mwave <- combn(c(1:Time)[-c(1:Time.complete)], i)
    }

    pattern <- matrix(0, nrow = ncol(mwave), ncol = k*Time)  #pattern matrix

    for (j in seq_len(nrow(pattern))) {
      for (m in seq_len(nrow(mwave)))
        pattern[j, ((mwave[m,j]-1)*k+1):(mwave[m,j]*k)] <- 1   #put missing in pattern matrix
    }

    completers <- rep(0,ncol(pattern))       # completers pattern
    dropper <- c(rep(0,Time.complete*k), rep(1, (Time-Time.complete)*k))   #droppers pattern

    ### users may wish to specify certain waves to have complete data

    if (is.null(complete.wave)==F) {

      if (i==1) {

        evalpattern <- matrix(mwave %in% c(complete.wave), nrow = nrow(mwave), byrow = F)
        # keep the patterns
        keep <- pattern[evalpattern==F, , drop = FALSE]

        if (is.null(dim(keep))==F){
          pattern <- keep
        }
        if (is.null(dim(keep))==T){    # it may become a non-matrix object, if only one row
          pattern <- t(as.matrix(keep))
        }
      }

      if (i>1) {
        evalpattern <- matrix(mwave %in% c(complete.wave), nrow = nrow(mwave), byrow = F)
        # keep the patterns
        keep <- pattern[colSums(evalpattern)==0, , drop = FALSE]

        if (is.null(dim(keep))==F) {
          pattern <- keep
        }
        if (is.null(dim(keep))==T) {    # it may become a non-matrix object, if only one row
          pattern <- t(as.matrix(keep))
        }
      }
    }

    ### design matrix

    if (pd!=0) {      # if there are droppers
      patmx <- rbind(pattern, completers, dropper)   # missing patterns for Mplus later
    }
    if (pd==0) {      # if there are no droppers
      patmx <- rbind(pattern, completers)
    }

    designs[[rs]] <- patmx

    #### pattern probs
    # p.probs

    if (pd==0) {
      p.probs <- c(rep(round((1-pc-pd)/(nrow(patmx)-1),6), nrow(patmx)-1), pc)
    }
    if (pd!=0) {
      p.probs <- c(rep(round((1-pc-pd)/(nrow(patmx)-2),6), nrow(patmx)-2), pc, pd)
    }

    probs[[rs]] <- p.probs

    # cost of each design
    cost.design[rs] <- sum((1-patmx[,ms.range]) %*% costmx*p.probs*n)

    rs <- rs+1
  }

  #### evaluate cost
  # only simulate for those designs that are below the budget limit

  if (eval.budget==T) {
    if (sum(cost.design>rm.budget)==length(cost.design)) {
      stop ("All wave missing designs cost more than the avaiable remaing budget. Try other designs.")
    }

    designs2 <- designs[cost.design<=rm.budget] #select the designs that are below the budget limit
    probs2 <- probs[cost.design<=rm.budget]
    miss.waves <- n.miss.waves[cost.design<=rm.budget]
    cost.design2 <- cost.design[cost.design<=rm.budget]
  }

  if (eval.budget==F) {
    designs2 <- designs
    probs2 <- probs
    miss.waves <- n.miss.waves
    cost.design2 <- cost.design
  }

  convergence.rate <- c()   #convergence rate
  weakest.param.DV <- c()
  weakest.param.IV <- c()
  weakest.para.power <- c()

  for (d in seq_len(length(designs2))) {

    patmx <- designs2[[d]]
    p.probs <- probs2[[d]]
    VNAMES <- NAMES

    ###distal variables
    if (is.null(distal.var)==F) {
      dis.pat <- matrix(0, nrow=nrow(patmx), ncol=length(distal.var))
      patmx <- cbind(patmx,dis.pat)
      VNAMES <- c(VNAMES,distal.var)
    }

    FNAME <- paste0("missing-waves-", miss.waves[d])  #file name

    #### need to fit in mplus format

    max.col <- 9        # the break-up point for the PATMISS matrix to meet the 90 character limitation of MPLUS
    r.PAT <- ceiling(length(VNAMES)/max.col)  #number of blocks

    if (r.PAT==1) {
      PATMISS <- rep(NA, nrow(patmx))
      for (j in 1:(nrow(patmx) - 1)) {
        Pat <- paste0(VNAMES, "(", patmx[j,], ")", collapse = " ")
        PATMISS[j] <- paste(Pat, "|", sep=" ")
      }
      lastPatLine <- paste0(VNAMES, "(", patmx[nrow(patmx),], ")", collapse = " ") #last line
      PATMISS[nrow(patmx)] <- paste(lastPatLine , ";")
      VNAMES.inp=c(paste(VNAMES,collapse=" "),";")
    }

    if (r.PAT==2) {
      PATMISS <- rep(NA, nrow(patmx)*2)
      patmx.1 <- patmx[,1:max.col]
      patmx.2 <- patmx[,(max.col+1):ncol(patmx)]

      VNAMES.1 <- VNAMES[1:max.col]
      VNAMES.2 <- VNAMES[(max.col+1):length(VNAMES)]

      for (j in 1:(nrow(patmx) - 1)) {
        Pat.1 <- paste0(VNAMES.1, "(", patmx.1[j,], ")", collapse = " ")
        Pat.2 <- paste0(VNAMES.2, "(", patmx.2[j,], ")", collapse = " ")
        PATMISS[2*j-1] <- Pat.1
        PATMISS[2*j] <- paste(Pat.2, "|", sep=" ")
      }
      lastPatLine.1 <- paste0(VNAMES.1, "(", patmx.1[nrow(patmx),], ")", collapse = " ") #last line
      lastPatLine.2 <- paste0(VNAMES.2,"(",patmx.2[nrow(patmx),], ")", collapse = " ")

      PATMISS[nrow(patmx)*2-1] <- lastPatLine.1
      PATMISS[nrow(patmx)*2] <- paste(lastPatLine.2 , ";")

      VNAMES.inp <- c(paste(VNAMES.1,collapse=" "),paste(VNAMES.2,collapse=" "),";")
    }

    if (r.PAT>2) {
      pat.list <- vector(mode = "list", length = r.PAT)
      name.list <- vector(mode = "list", length = r.PAT)
      name.list2 <- vector(mode = "list", length = r.PAT)

      PATMISS <- rep(NA, nrow(patmx)*r.PAT)

      pat.list[[1]] <- patmx[,1:max.col]
      name.list[[1]] <- VNAMES[1:max.col]

      for (rp in 2:(r.PAT-1)) {
        pat.list[[rp]] <- patmx[,((rp-1)*max.col+1):(rp*max.col)]
        name.list[[rp]] <- VNAMES[((rp-1)*max.col+1):(rp*max.col)]
      }
      pat.list[[r.PAT]] <- patmx[,((r.PAT-1)*max.col+1):length(VNAMES)]
      name.list[[r.PAT]] <- VNAMES[((r.PAT-1)*max.col+1):length(VNAMES)]

      PATMISS.j <- list()

      for (j in 1:(nrow(patmx) - 1)) {
        Pat.j <- list()  #storage

        for (rp in 1:(r.PAT-1)) {
          Pat.j[[rp]] <- paste0(name.list[[rp]], "(", pat.list[[rp]][j,], ")", collapse = " ")
        }
        Pat.j[[r.PAT]] <- paste(paste0(name.list[[r.PAT]], "(", pat.list[[r.PAT]][j,], ")", collapse = " "),"|",sep=" ")
        PATMISS.j[[j]] <- unlist(Pat.j)
      }

      # last line
      j <- nrow(patmx)
      for (rp in 1:(r.PAT-1)){
      Pat.j[[rp]] <- paste0(name.list[[rp]], "(", pat.list[[rp]][j,], ")", collapse = " ")
      }
      Pat.j[[r.PAT]] <- paste(paste0(name.list[[r.PAT]], "(", pat.list[[r.PAT]][j,], ")", collapse = " "), ";", sep=" ")
      PATMISS.j[[j]] <- unlist(Pat.j)

      PATMISS <- unlist(PATMISS.j)

      # variable names

      for (l in 1:r.PAT) {
        name.list2[[l]] <- paste0(name.list[[l]], collapse=" ")
      }
      VNAMES.inp <- c(unlist(name.list2), ";")
    }

    # pattern probs
    r.probs <- ceiling(nrow(patmx)/max.col)


    if (r.probs==1) {
      A <- paste(p.probs[1:(length(p.probs)-1)], "|", collapse = "")
      B <- paste(p.probs[length(p.probs)], ";")
      PATPROBS <- paste(A,B)
    }

    if (r.probs==2) {

      if ((max.col*(r.probs-1))<(length(p.probs)-1)) {
        A <- paste(p.probs[1:max.col], "|", collapse = "")
        B <- paste(p.probs[(max.col+1):(length(p.probs)-1)],"|", collapse="")
        C <- paste(p.probs[length(p.probs)],";")
        PATPROBS <- c(A,paste(B,C))
      }

      if ((max.col*(r.probs-1))==(length(p.probs)-1)) {
        A <- paste(p.probs[1:max.col], "|", collapse = "")
        #B=paste(p.probs[(max.col+1):(length(p.probs)-1)],"|", collapse="")
        C <- paste(p.probs[length(p.probs)], ";")
        PATPROBS <- c(A,paste(B,C))
      }
    }

    if (r.probs>2) {

      prob.list <- list()
      for (rp in 1:(r.probs-1)) {
        prob.list[[rp]] <- paste(p.probs[((rp-1)*max.col+1):(rp*max.col)],"|", collapse="")
      }

      if ((max.col*(r.probs-1))==(length(p.probs)-1)) {
        #B=paste(p.probs[((r.probs-1)*max.col+1):(length(p.probs)-1)],"|", collapse="")
        C <- paste(p.probs[length(p.probs)], ";")
        prob.list[[r.probs]] <- C
        PATPROBS <- unlist(prob.list)
      } else if ((max.col*(r.probs-1)) < (length(p.probs)-1)) {
        B <- paste(p.probs[((r.probs-1)*max.col+1):(length(p.probs)-1)],"|", collapse="")
        C <- paste(p.probs[length(p.probs)],";")
        prob.list[[r.probs]] <- paste(B,C)
        PATPROBS <- unlist(prob.list)
      }
    }


    scriptMplus <- c(
      paste0("TITLE: ", FNAME, ";"),
      "MONTECARLO: ",
      #paste0("NAMES ARE ", paste(VNAMES, collapse = " "), ";"),
      "NAMES ARE ",
      VNAMES.inp,
      paste0("NOBSERVATIONS = ", n, ";"),
      paste0("NREPS = ", nreps, ";"),
      paste0("SEED = ", seed, ";"),
      "PATMISS =",
      PATMISS,
      "PATPROBS =",
      PATPROBS,
      "MODEL POPULATION: ",
      design0.out$input$model.population,
      "MODEL: ",
      design0.out$input$model,
      "OUTPUT: TECH9;")

    fileConn <- file(paste0(FNAME, ".inp"))
    writeLines(scriptMplus, fileConn)  #write the input file into the folder
    close(fileConn)

    wd.dir <- getwd()
    MplusAutomation::runModels(target=paste0(wd.dir,"/",FNAME,".inp"),replaceOutfile = F)                        #run the input file in Mplus

    filename <- paste0(paste0(FNAME, ".out"))  #output file name
    design.out <- MplusAutomation::readModels(filename)


    temp <- design.out$parameters$unstandardized
    if (is.null(temp)==T) {
      convergence.rate[d] <- 0
      weakest.param.DV[d] <- "NA"
      weakest.param.IV[d] <- "NA"
      weakest.para.power[d] <- 0
    }
    if (is.null(temp)==F) {
      f.param <- temp[temp$paramHeader%in%focal.param$paramHeader&temp$param%in%focal.param$param,]
      weakest.f.param <- f.param[f.param$pct_sig_coef==min(f.param$pct_sig_coef),]
      if (nrow(weakest.f.param)>1) {
        weakest.f.param <- weakest.f.param[1, ]      ########## may need to be changed later
      }
      convergence.rate[d] <- design.out$summaries$ChiSqM_NumComputations/nreps    #converged number of simulations
      weakest.param.DV[d] <- weakest.f.param[, "paramHeader"]
      weakest.param.IV[d] <- weakest.f.param[, "param"]
      weakest.para.power[d] <- weakest.f.param[, "pct_sig_coef"]
    }
  }

  VNAMES <- NAMES

  sim.results.out <- cbind.data.frame(convergence.rate,   #convergence rate
                                   weakest.param.DV,
                                   weakest.param.IV,
                                   weakest.para.power,
                                   "cost.design"=cost.design2, # cost of each design
                                   miss.waves)

  opt.design.1 <- sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]),]

  if (nrow(opt.design.1)==1) {
    opt.design <- opt.design.1
  }
  if (nrow(opt.design.1)>1) {
    opt.design <- opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ]
  }

  op <- which(miss.waves==opt.design$miss.waves)    #which design is chosen
  opt.pattern <- designs2[[op]]
  colnames(opt.pattern) <- VNAMES
  opt.probs <- probs2[[op]]

  misc <- list(time = Time, k = k, focal.param = focal.param)
  
  re.ob=list("results" = sim.results.out,
             "opt.design" = opt.design,
             "opt.pattern" = opt.pattern,
             "opt.probs" = opt.probs,
             "n.miss.waves" = opt.design$miss.waves,
             "misc" = misc)
  
  class(re.ob) <- append(class(re.ob),"simpm")
  
  return(re.ob)

}
