#' Search for the optimal item-level PM design via forward assembly.
#' 
#'  \code{forward.opt} runs simulations using M\emph{plus}. It returns the search results for optimal item-level PM designs via forward assembly. 
#'  
#' @inheritParams balance.miss
#' @param max.mk Specify the maximum number of unique missing data
#'   patterns in the selected design. Only applicable if forward assembly
#'   is used. 
#'
#' @return An object containing the information of the optimal 
#'    item-level missing design. The optimal design is the one that yields
#'    highest power for testing the focal parameters, compared to other
#'    plausible candidate PM designs.
#' @seealso \code{\link{simPM}} which is a warpper function for this
#'    function.
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export forward.opt
#' @examples

forward.opt <- function(
  VNAMES,
  distal.var,
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
  max.mk,                 # maximum number of missing slots allowed
  eval.budget=F,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
  rm.budget=NULL,          # remaining available budget
  complete.var=NULL)
{
  NAMES <- VNAMES
  previous <- list()
  if (max.mk==1) {          # when num.miss=1
    output <- opt.nm.1(VNAMES,
                    distal.var,
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
                    complete.var)

    ### if user wish to evaluate the budget
    if (eval.budget==T) {
      if ((1 - prod(output$results$cost.design > rm.budget))==0) {
        return(output)
        warning("All designs would cost more than the avaiable remaing budget. Consider increasing max.mk.")
      }
    }
    return(output)
  } else if (max.mk >= 2) {     # when max number of missing slots is greater than 1

    future.k <- (Time - Time.complete) * k   # data points not yet completed, also maximum possible # of missing
    ms.range <- c((Time.complete*k + 1):(Time * k))  # The available time slots to plant missingness

    design.order <- rep(NA, max.mk)    # keep record of the order of the selected designs

    #num.miss must be less than future.k
    if (max.mk >= future.k) {
      stop ("number of missing points exceeds the maximum possible value.")
    }

    output1 <- opt.nm.1(VNAMES,   #get the best design with one missing slot
                     distal.var,
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
                     complete.var)
    
    previous[[1]] <- output1$opt.pattern
    design.order[1] <- output1$design.order

    for (num.miss in 2:max.mk) {

      if (eval.budget == T) {            # evaluate the cost and the budget
        if (num.miss == max.mk) {        # cannot evaluate the cost without knowing the previous selected patterns

          ms.combn <- combn(ms.range, num.miss) #all possible combinations of missing slots given num.miss
          all.pattern <- matrix(0, nrow = choose(future.k, num.miss), ncol = length(VNAMES)) # place holder, all possible patterns with a certain ms

          for (q in seq_len(nrow(all.pattern))) {
            all.pattern[q, ms.combn[, q]] <- 1
          }                                         # add missingness

          if (length(complete.var)==1) {
            # find which column these variables are
            complete.cols <- which(VNAMES %in% complete.var)
            # whether to keep the rows
            keep <- all.pattern[, complete.cols]==0

            temp.pattern <- all.pattern[keep, , drop=FALSE]
            if (is.null(dim(temp.pattern))==F) {
              all.pattern <- temp.pattern
            }
            if (is.null(dim(temp.pattern))==T) {    # it may become a non-matrix object, if only one row
              all.pattern <- t(as.matrix(temp.pattern))
            }
          }

          if (length(complete.var)>1) {
            complete.cols <- which(VNAMES %in% complete.var)
            keep <- rowSums(all.pattern[, complete.cols]==0)==length(complete.var)

            temp.pattern <- all.pattern[keep, , drop=FALSE]
            if (is.null(dim(temp.pattern))==F) {
              all.pattern <- temp.pattern
            }
            if (is.null(dim(temp.pattern))==T) {
              all.pattern <- t(as.matrix(temp.pattern))       # need to be a matrix object even if it's only one row
            }
          }


          # the optimal pattern obtained in previous round of selection
          opt1.pattern <- previous[[num.miss-1]]

          # unit cost of each pattern
          if (nrow(all.pattern)>1) {
            unit.cost <- rowSums((1 - all.pattern[, ms.range]) * costmx) * ((1 - pc - pd) * n / num.miss)
          }
          if (nrow(all.pattern)==1) {
            temp <- t(as.matrix(1 - all.pattern[, ms.range]))
            unit.cost <- rowSums(temp * costmx) * ((1 - pc - pd) * n / num.miss)
          }

          # cost of previous chosen design
          if (pd==0) {
            opt1.cost <- sum(rowSums((1 - opt1.pattern[, ms.range]) * costmx) * c(rep((1 - pc) * n / num.miss, num.miss - 1), pc * n))
          }
          if (pd!=0) {
            opt1.cost <- sum(rowSums((1 - opt1.pattern[, ms.range]) * costmx) * c(rep((1 - pc - pd) * n / num.miss, num.miss - 1), pc * n, pd * n))
          }

          sum.cost <- unit.cost + opt1.cost
          if ((1 - prod(sum.cost > rm.budget))==0) {
            stop("All designs cost more than the avaiable remaing budget, please allow a a larger max.mk. To turn off the budget evaluation, set eval.budget=F and rm.budget=NULL")
          }
        }
      }


      ms.combn <- combn(ms.range, num.miss) #all possible combinations of missing slots given num.miss
      all.pattern <- matrix(0, nrow = choose(future.k, num.miss), ncol = length(VNAMES)) # place holder, all possible patterns with a certain ms

      # update the missing patterns with 1s
      for (q in seq_len(nrow(all.pattern))) {
        all.pattern[q, ms.combn[, q]] <- 1
      }                                         # add missingness

      # the optimal pattern obtained in previous round of selection
      opt1.pattern <- previous[[num.miss-1]]

      if (length(complete.var)==1) {
        # find which column these variables are
        complete.cols <- which(VNAMES %in% complete.var)
        # whether to keep the rows
        keep <- all.pattern[, complete.cols]==0

        temp.pattern <- all.pattern[keep, , drop=FALSE]
        if (is.null(dim(temp.pattern))==F) {
          all.pattern <- temp.pattern
        }
        if (is.null(dim(temp.pattern))==T) {
          all.pattern <- t(as.matrix(temp.pattern))
        }


        # take out these designs columns out of the ms.combn

        temp.combs <- ms.combn[ , keep, drop=FALSE]

        if(is.null(dim(temp.combs))==F) {
          ms.combn <- temp.combs
        }
        if(is.null(dim(temp.combs))==T) {        #may become a non-matrix object if only one column
          ms.combn <- as.matrix(temp.combs)
        }
      }


      if (length(complete.var)>1) {
        complete.cols <- which(VNAMES %in% complete.var)
        keep <- rowSums(all.pattern[, complete.cols]==0)==length(complete.var)
        temp.pattern <- all.pattern[keep, , drop=FALSE]

        if (is.null(dim(temp.pattern))==F) {
          all.pattern <- temp.pattern
        }
        if (is.null(dim(temp.pattern))==T) {
          all.pattern <- t(as.matrix(temp.pattern))
        }


        # take out these designs columns out of the ms.combn

        temp.combs <- ms.combn[, keep, drop=FALSE]
        if (is.null(dim(temp.combs))==F) {
          ms.combn <- temp.combs
        }
        if (is.null(dim(temp.combs))==T) {
          ms.combn <- as.matrix(temp.combs)
        }
      }

      # storage bins for Mplus simulation results
      convergence.rate <- rep(NA, ncol(ms.combn))   #convergence rate
      weakest.param.DV <- rep(NA, ncol(ms.combn))
      weakest.param.IV <- rep(NA, ncol(ms.combn))
      weakest.para.power <- rep(NA, ncol(ms.combn))
      cost.design <- rep(NA, ncol(ms.combn)) # cost of each design
      miss.num <- rep(num.miss, ncol(ms.combn))
      miss.name <- matrix(NA, ncol(ms.combn), num.miss)
      sim.seq <- rep(NA, ncol(ms.combn))
      miss.loc <- matrix(NA, ncol(ms.combn), num.miss)


      # generate Mplus input files
      for (i in seq_len(nrow(all.pattern))) {

        ###distal variables
        if (is.null(distal.var)==F) {
          dis.pat <- rep(0,length(distal.var))
          patmx <- rbind(c(all.pattern[i, ], dis.pat), cbind(opt1.pattern, t(replicate(nrow(opt1.pattern), dis.pat))))
          VNAMES <- c(VNAMES, distal.var)
        }
        if (is.null(distal.var)==T) {
          patmx <- rbind(all.pattern[i, ], opt1.pattern)   # missing patterns
        }

        FNAME <- paste0("missing-", num.miss, "-sim.seq-", i)  #file name

        max.col <- 9        # the break-up point for the PATMISS matrix to meet the 90 character limitation of MPLUS
        r.PAT <- ceiling(length(VNAMES) / max.col)  #number of blocks

        if (r.PAT==1) {
          PATMISS <- rep(NA, nrow(patmx))
          for (j in 1:(nrow(patmx) - 1)) {
            Pat <- paste0(VNAMES, "(", patmx[j,], ")", collapse = " ")
            PATMISS[j] <- paste(Pat, "|", sep=" ")
          }
          lastPatLine <- paste0(VNAMES, "(", patmx[nrow(patmx),], ")", collapse = " ") #last line
          PATMISS[nrow(patmx)] <- paste(lastPatLine , ";")
          VNAMES.inp <- c(paste(VNAMES,collapse=" "),";")
        }

        if (r.PAT==2) {
          PATMISS <- rep(NA, nrow(patmx)*2)
          patmx.1 <- patmx[, 1:max.col]
          patmx.2 <- patmx[, (max.col+1):ncol(patmx)]

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

          PATMISS[nrow(patmx) * 2 - 1] <- lastPatLine.1
          PATMISS[nrow(patmx) * 2] <- paste(lastPatLine.2, ";")

          VNAMES.inp <- c(paste(VNAMES.1,collapse=" "), paste(VNAMES.2, collapse=" "), ";")
        }

        if (r.PAT > 2) {
          pat.list <- list()
          name.list <- list()
          name.list2 <- list()

          PATMISS <- rep(NA, nrow(patmx) * r.PAT)

          pat.list[[1]] <- patmx[, 1:max.col]
          name.list[[1]] <- VNAMES[1:max.col]

          for (rp in 2:(r.PAT-1)) {
            pat.list[[rp]] <- patmx[, ((rp - 1) * max.col + 1):(rp * max.col)]
            name.list[[rp]] <- VNAMES[((rp - 1) * max.col + 1):(rp * max.col)]
          }
          pat.list[[r.PAT]] <- patmx[, ((r.PAT - 1) * max.col + 1):length(VNAMES)]
          name.list[[r.PAT]] <- VNAMES[((r.PAT - 1) * max.col + 1):length(VNAMES)]

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
          for (rp in 1:(r.PAT-1)) {
            Pat.j[[rp]] <- paste0(name.list[[rp]], "(", pat.list[[rp]][j,], ")", collapse = " ")
          }
          Pat.j[[r.PAT]] <- paste(paste0(name.list[[r.PAT]], "(", pat.list[[r.PAT]][j,], ")", collapse = " "),";",sep=" ")
          PATMISS.j[[j]] <- unlist(Pat.j)

          PATMISS <- unlist(PATMISS.j)

          # variable names

          for (l in 1:r.PAT) {
            name.list2[[l]] <- paste0(name.list[[l]], collapse=" ")
          }
          VNAMES.inp <- c(unlist(name.list2), ";")
        }

        # pattern probs
        r.probs <- ceiling(nrow(patmx) / max.col)

        if (pd==0) {
          p.probs <- c(rep(round((1 - pc - pd) / (nrow(patmx) - 1), 6), nrow(patmx) - 1), pc)
        }
        if (pd!=0) {
          p.probs <- c(rep(round((1 - pc - pd) / (nrow(patmx) - 2), 6), nrow(patmx) - 2), pc, pd)
        }

        if (r.probs==1) {
          A <- paste(p.probs[1:(length(p.probs)-1)], "|", collapse = "")
          B <- paste(p.probs[length(p.probs)],";")
          PATPROBS <- paste(A,B)
        }

        if (r.probs==2) {

          if ((max.col*(r.probs-1))<(length(p.probs)-1)) {
            A <- paste(p.probs[1:max.col], "|", collapse = "")
            B <- paste(p.probs[(max.col+1):(length(p.probs)-1)], "|", collapse="")
            C <- paste(p.probs[length(p.probs)], ";")
            PATPROBS <- c(A, paste(B,C))
          }

          if ((max.col*(r.probs-1))==(length(p.probs)-1)) {
            A <- paste(p.probs[1:max.col], "|", collapse = "")
            #B=paste(p.probs[(max.col+1):(length(p.probs)-1)],"|", collapse="")
            C <- paste(p.probs[length(p.probs)],";")
            PATPROBS <- c(A,paste(B,C))
          }
        }

        if (r.probs>2) {

          prob.list <- list()
          for (rp in 1:(r.probs-1)) {
            prob.list[[rp]] <- paste(p.probs[((rp-1)*max.col+1):(rp*max.col)],"|",collapse="")
          }

          if ((max.col*(r.probs-1))==(length(p.probs)-1)) {
            #B=paste(p.probs[((r.probs-1)*max.col+1):(length(p.probs)-1)],"|", collapse="")
            C <- paste(p.probs[length(p.probs)], ";")
            prob.list[[r.probs]] <- C
            PATPROBS <- unlist(prob.list)
          } else if ((max.col*(r.probs-1))<(length(p.probs)-1)) {
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
          convergence.rate[i] <- 0
          weakest.param.DV[i] <- "NA"
          weakest.param.IV[i] <- "NA"
          weakest.para.power[i] <- 0
        }
        if (is.null(temp)==F) {
          f.param <- temp[temp$paramHeader %in% focal.param$paramHeader&temp$param %in% focal.param$param, ]
          weakest.f.param <- f.param[f.param$pct_sig_coef==min(f.param$pct_sig_coef), ]
          if (nrow(weakest.f.param)>1) {
            weakest.f.param <- weakest.f.param[1, ]      ########## may need to be changed later
          }

          convergence.rate[i] <- design.out$summaries$ChiSqM_NumComputations/nreps    #converged number of simulations
          weakest.param.DV[i] <- weakest.f.param[,"paramHeader"]
          weakest.param.IV[i] <- weakest.f.param[,"param"]
          weakest.para.power[i] <- weakest.f.param[,"pct_sig_coef"]

        }

        if (pd==0) {
          cost.design[i] <- sum(c(rep(n*(1-pc)/num.miss,num.miss), pc*n)*((1-patmx[,ms.range]) %*% costmx))   #patmx depends on i
        }
        if (pd!=0) {
          cost.design[i] <- sum(c(rep(n*(1-pc-pd)/num.miss,num.miss),pc*n,pd*n)*((1-patmx[,ms.range]) %*% costmx))   #patmx depends on i
        }

        miss.name[i,] <- VNAMES[ms.combn[,i]]
        sim.seq[i] <- i  # location as specified in the miss.combn matrix
        miss.loc[i,] <- ms.combn[,i]
        VNAMES <- NAMES
      }


      colnames(miss.name) <- paste0("miss.var", 1:num.miss)
      colnames(miss.loc) <- paste0("miss.loc", 1:num.miss)

      sim.results.out <- cbind.data.frame(convergence.rate,   #convergence rate
                                       weakest.param.DV,
                                       weakest.param.IV,
                                       weakest.para.power,
                                       cost.design, # cost of each design
                                       miss.num,
                                       miss.name,
                                       sim.seq,
                                       miss.loc)

      if (eval.budget==F) {
      opt.design.1 <- sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[, "weakest.para.power"]), ]
      }

      if (eval.budget==T) {
        sim.results.out <- sim.results.out[sim.results.out$cost.design <= rm.budget, ]
        opt.design.1 <- sim.results.out[sim.results.out[, "weakest.para.power"]==max(sim.results.out[, "weakest.para.power"]), ]
        }

      if (nrow(opt.design.1)==1) {
        opt.design <- opt.design.1
      }

      if (nrow(opt.design.1)>1) {
        n.min.cost <- nrow(opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ])
        if (n.min.cost==1) {
          opt.design <- opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ]
        } else {
          opt.min.cost <- opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design), ]
          sum.loc <- rowSums(opt.min.cost[, colnames(miss.loc)])
          opt.design <- opt.min.cost[sum.loc==max(sum.loc), ]
        }
      }

      op <- opt.design[, "sim.seq"]

      opt.pattern <- rbind(all.pattern[op, ], opt1.pattern)
      colnames(opt.pattern) <- VNAMES

      if (pd==0) {
        opt.probs <- c(rep(round((1-pc)/(nrow(patmx)-1), 6), num.miss), pc)
      }
      if (pd!=0) {
        opt.probs <- c(rep(round((1-pc-pd)/(nrow(patmx)-2), 6), num.miss), pc, pd)
      }

      design.order[num.miss] <- op
      
      misc <- list(time = Time, k = k, focal.param = focal.param, max.mk = max.mk)

      re.ob <- list(
        "results"=sim.results.out,
        "opt.design"=opt.design,
        "opt.pattern"=opt.pattern,
        "opt.probs"=opt.probs,
        "design.order"=design.order,
        "misc"=misc)

      previous[[num.miss]]=re.ob$opt.pattern   # update the previous pattern list

    }
    return(re.ob)
  }
}



