#' Select the optimal pattern with only one missing indicator using Mplus

#' @param design0.out Mplus output file which contains the a priori power analysis/sample size planning (simulation) results for this
#' specific model assuming a complete data design. Theoretically, such analysis was supposed to be conducted before the study began.
#' @param VNAMES The names of the observed variables, ordered by the time they are measured
#' @param Time The total number of time points (waves of data collection)
#' @param Time.complete Number of waves of data collection that have been completed before the funding cut occured
#' @param k The number of observed variables collected at each wave
#' @param pc The proportion of subjects that will participate in all of the following waves of data collection and provide complete data (must be greater than 0)
#' @param pd The proportion of subjects that will not participate in any of the following waves of data collection (i.e., attritors). This can be 0.
#' @param costmx  The vector containing the unit cost of each observed variable which has no data collected yet. They are constant across subjects, but they can vary across variables and across time.
#' @param n The total sample size as initially planned
#' @param nreps Number of replications for Monte Carlo simulation for each possible design
#' @param focal.param The parameters of focal interest. If engine="l", the focal parameters should be specified using
#' the lavaan script. If engine="m", the focal parameters should be specified based on the Mplus output file design0.out.
#' @param complete.var Specify if there are any variables that need to have complete data collected across all the participating subjects
#' @param distal.var Any distal variables included in the model that would have complete data
#' @param seed Random seet for simulation
#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export opt.nm.1
#' @examples


opt.nm.1=function(
  VNAMES,    # a vactor containing the measured variable names
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
  design0.out,        # an readModel object, from the user supplied Mplus output file
  focal.param,        # user identified focal parameters, in a matrix form. User needs to identify the rows in the readModel parameter object that are of focal interest
  complete.var=NULL){ # a list of the variables that need to be complete

  num.miss=1   #this function only deals with one missing slot
  NAMES=VNAMES

  # calculated info based on user supplied info
  future.k=(Time-Time.complete)*k   # data points not yet completed, also maximum possible # of missing
  ms.range=c((Time.complete*k+1):(Time*k))  # The available time slots to plant missingness

  ms.combn=combn(ms.range,num.miss) #all possible combinations of missing slots given num.miss
  all.pattern=matrix(0,nrow=choose(future.k,num.miss),ncol=length(VNAMES)) # place holder, all possible patterns with a certain ms

  # update the missing patterns with 1s
  for (q in 1:nrow(all.pattern)){
    all.pattern[q,ms.combn[,q]]=1
  }                                         # add missingness
  completers=rep(0,ncol(all.pattern))       # completers pattern
  dropper=c(rep(0,Time.complete*k),rep(1,future.k))   #droppers pattern

  #### if user specify certain variables to be complete

  if (length(complete.var)==1){
    # find which column these variables are
    complete.cols=which(VNAMES%in%complete.var)
    # whether to keep the rows
    keep=all.pattern[,complete.cols]==0

    all.pattern=all.pattern[keep,]
    # take out these designs columns out of the ms.combn
    ms.combn=t(as.matrix(ms.combn[,keep]))
  }

  if (length(complete.var)>1){
    complete.cols=which(VNAMES%in%complete.var)
    keep=rowSums(all.pattern[,complete.cols]==0)==length(complete.var)

    temp.pattern=all.pattern[keep,]
    if (is.null(dim(temp.pattern))==F){
      all.pattern=temp.pattern
    }
    if (is.null(dim(temp.pattern))==T){
      all.pattern=t(as.matrix(temp.pattern))
    }

    temp.combs=ms.combn[,keep]
    if(is.null(dim(temp.combs))==F){
      ms.combn=temp.combs
    }
    if(is.null(dim(temp.combs))==T){
      ms.combn=as.matrix(temp.combs)
    }
  }

  ### storage bins for Mplus simulation results
  convergence.rate=rep(NA,nrow(all.pattern))   #convergence rate
  weakest.param.DV=rep(NA,nrow(all.pattern))
  weakest.param.IV=rep(NA,nrow(all.pattern))
  weakest.para.power=rep(NA,nrow(all.pattern))
  cost.design=rep(NA,nrow(all.pattern)) # cost of each design
  miss.num=rep(num.miss,nrow(all.pattern))
  miss.name=matrix(NA,nrow(all.pattern),num.miss)
  sim.seq=rep(NA,nrow(all.pattern))
  miss.loc=matrix(NA,nrow(all.pattern),num.miss)

  ### generate Mplus input files
  for (i in 1:nrow(all.pattern)){
    # because num.miss=1, no previously selected pattern applicable
    if(pd!=0){
      patmx=rbind(all.pattern[i,],completers,dropper)   # missing patterns
    }
    if(pd==0){
      patmx=rbind(all.pattern[i,],completers)
    }

    #FNAME=paste0("missing-",num.miss,"-sim.seq-",VNAMES[ms.combn[,i]])  #file name
    FNAME=paste0("missing-",num.miss,"-sim.seq-",i)  #file name


    ### distal variables
    if (is.null(distal.var)==F){
      dis.pat=matrix(0,nrow=nrow(patmx),ncol=length(distal.var))
      patmx=cbind(patmx,dis.pat)
      VNAMES=c(VNAMES,distal.var)
    }

    ### get the patterns into the Mplus format

    max.col=9        # the break-up point for the PATMISS matrix to meet the 90 character limitation of MPLUS
    r.PAT=ceiling(length(VNAMES)/max.col)  #number of blocks

    if (r.PAT==1){
      PATMISS <- rep(NA, nrow(patmx))
      for (j in 1:(nrow(patmx) - 1)){
        Pat <- paste0(VNAMES, "(", patmx[j,], ")", collapse = " ")
        PATMISS[j] <- paste(Pat, "|", sep=" ")
      }
      lastPatLine <- paste0(VNAMES, "(", patmx[nrow(patmx),], ")", collapse = " ") #last line
      PATMISS[nrow(patmx)] <- paste(lastPatLine , ";")
      VNAMES.inp=c(paste(VNAMES,collapse=" "),";")
    }

    if (r.PAT==2){
      PATMISS <- rep(NA, nrow(patmx)*2)
      patmx.1=patmx[,1:max.col]
      patmx.2=patmx[,(max.col+1):ncol(patmx)]

      VNAMES.1=VNAMES[1:max.col]
      VNAMES.2=VNAMES[(max.col+1):length(VNAMES)]

      for (j in 1:(nrow(patmx) - 1)){
        Pat.1 <- paste0(VNAMES.1, "(", patmx.1[j,], ")", collapse = " ")
        Pat.2 <- paste0(VNAMES.2, "(", patmx.2[j,], ")", collapse = " ")
        PATMISS[2*j-1] <- Pat.1
        PATMISS[2*j] <- paste(Pat.2, "|", sep=" ")
      }
      lastPatLine.1 <- paste0(VNAMES.1, "(", patmx.1[nrow(patmx),], ")", collapse = " ") #last line
      lastPatLine.2 <- paste0(VNAMES.2,"(",patmx.2[nrow(patmx),], ")", collapse = " ")

      PATMISS[nrow(patmx)*2-1] <- lastPatLine.1
      PATMISS[nrow(patmx)*2] <- paste(lastPatLine.2 , ";")

      VNAMES.inp=c(paste(VNAMES.1,collapse=" "),paste(VNAMES.2,collapse=" "),";")
    }

    if (r.PAT>2){
      pat.list=list()
      name.list=list()
      name.list2=list()

      PATMISS <- rep(NA, nrow(patmx)*r.PAT)

      pat.list[[1]]=patmx[,1:max.col]
      name.list[[1]]=VNAMES[1:max.col]

      for (rp in 2:(r.PAT-1)){
        pat.list[[rp]]=patmx[,((rp-1)*max.col+1):(rp*max.col)]
        name.list[[rp]]=VNAMES[((rp-1)*max.col+1):(rp*max.col)]
      }
      pat.list[[r.PAT]]=patmx[,((r.PAT-1)*max.col+1):length(VNAMES)]
      name.list[[r.PAT]]=VNAMES[((r.PAT-1)*max.col+1):length(VNAMES)]

      PATMISS.j=list()

      for (j in 1:(nrow(patmx) - 1)){
        Pat.j=list()  #storage

        for (rp in 1:(r.PAT-1)){
          Pat.j[[rp]] <- paste0(name.list[[rp]], "(", pat.list[[rp]][j,], ")", collapse = " ")
        }
        Pat.j[[r.PAT]] <- paste(paste0(name.list[[r.PAT]], "(", pat.list[[r.PAT]][j,], ")", collapse = " "),"|",sep=" ")
        PATMISS.j[[j]]=unlist(Pat.j)
      }

      # last line
      j=nrow(patmx)
      for (rp in 1:(r.PAT-1)){
        Pat.j[[rp]] <- paste0(name.list[[rp]], "(", pat.list[[rp]][j,], ")", collapse = " ")
      }
      Pat.j[[r.PAT]] <- paste(paste0(name.list[[r.PAT]], "(", pat.list[[r.PAT]][j,], ")", collapse = " "),";",sep=" ")
      PATMISS.j[[j]]=unlist(Pat.j)

      PATMISS=unlist(PATMISS.j)

      # variable names

      for (l in 1:r.PAT){
        name.list2[[l]]=paste0(name.list[[l]],collapse=" ")
      }
      VNAMES.inp=c(unlist(name.list2),";")
    }

    # pattern probs
    r.probs=ceiling(nrow(patmx)/max.col)

    if (pd==0){
      p.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-1),6), nrow(patmx)-1),pc)
    }
    if (pd!=0){
      p.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-2),6), nrow(patmx)-2),pc,pd)
    }

    if (r.probs==1){
      A=paste(p.probs[1:(length(p.probs)-1)], "|", collapse = "")
      B=paste(p.probs[length(p.probs)],";")
      PATPROBS=paste(A,B)
    }

    if (r.probs==2){

      if ((max.col*(r.probs-1))<(length(p.probs)-1)){
        A=paste(p.probs[1:max.col], "|", collapse = "")
        B=paste(p.probs[(max.col+1):(length(p.probs)-1)],"|", collapse="")
        C=paste(p.probs[length(p.probs)],";")
        PATPROBS=c(A,paste(B,C))
      }

      if ((max.col*(r.probs-1))==(length(p.probs)-1)){
        A=paste(p.probs[1:max.col], "|", collapse = "")
        #B=paste(p.probs[(max.col+1):(length(p.probs)-1)],"|", collapse="")
        C=paste(p.probs[length(p.probs)],";")
        PATPROBS=c(A,paste(B,C))
      }
    }

    if (r.probs>2){

      prob.list=list()
      for (rp in 1:(r.probs-1)){
        prob.list[[rp]]=paste(p.probs[((rp-1)*max.col+1):(rp*max.col)],"|",collapse="")
      }

      if ((max.col*(r.probs-1))==(length(p.probs)-1)){
        #B=paste(p.probs[((r.probs-1)*max.col+1):(length(p.probs)-1)],"|", collapse="")
        C=paste(p.probs[length(p.probs)],";")
        prob.list[[r.probs]]=C
        PATPROBS=unlist(prob.list)
      }else if ((max.col*(r.probs-1))<(length(p.probs)-1)){
        B=paste(p.probs[((r.probs-1)*max.col+1):(length(p.probs)-1)],"|", collapse="")
        C=paste(p.probs[length(p.probs)],";")
        prob.list[[r.probs]]=paste(B,C)
        PATPROBS=unlist(prob.list)
      }
    }

    ### write the input scripts
    scriptMplus <- c(
      paste0("TITLE: ", FNAME, ";"),
      "MONTECARLO: ",
      "NAMES ARE ",
      VNAMES.inp,
      #paste0("NAMES ARE ", paste(VNAMES, collapse = " "), ";"),
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

    wd.dir=getwd()
    MplusAutomation::runModels(target=paste0(wd.dir,"/",FNAME,".inp"),replaceOutfile = F)                        #run the input file in Mplus

    filename=paste0(paste0(FNAME, ".out"))  #output file name
    design.out=MplusAutomation::readModels(filename)

    temp=design.out$parameters$unstandardized
    if (is.null(temp)==T){
      convergence.rate[i]=0
      weakest.param.DV[i]="NA"
      weakest.param.IV[i]="NA"
      weakest.para.power[i]=0
    }
    if (is.null(temp)==F){
      ### focal parameter
      f.param=temp[temp$paramHeader%in%focal.param$paramHeader&temp$param%in%focal.param$param,]

      weakest.f.param=f.param[f.param$pct_sig_coef==min(f.param$pct_sig_coef),]
      if (nrow(weakest.f.param)>1){
        weakest.f.param=weakest.f.param[1,]      ########## may need to be changed later
      }
      convergence.rate[i]=design.out$summaries$ChiSqM_NumComputations/nreps    #converged number of simulations
      weakest.param.DV[i]=weakest.f.param[,"paramHeader"]
      weakest.param.IV[i]=weakest.f.param[,"param"]
      weakest.para.power[i]=weakest.f.param[,"pct_sig_coef"]

      if (pd==0){
        cost.design[i]=sum(c((1-pc)*n,pc*n)*((1-patmx[,ms.range])%*%costmx))   #patmx depends on i
      }
      if (pd!=0){
        cost.design[i]=sum(c((1-pc-pd)*n,pc*n,pd*n)*((1-patmx[,ms.range])%*%costmx))
      }

      miss.name[i]=VNAMES[ms.combn[,i]]
      sim.seq[i]=i  # location as specified in the miss.combn matrix
      miss.loc[i,]=ms.combn[,i]
    }
    VNAMES=NAMES
  }


  ### combine the results
  sim.results.out=cbind.data.frame(convergence.rate,   #convergence rate
                                   weakest.param.DV,
                                   weakest.param.IV,
                                   weakest.para.power,
                                   cost.design, # cost of each design
                                   miss.num,
                                   miss.name,
                                   sim.seq,
                                   miss.loc)

  opt.design.1=sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]),]

  if (nrow(opt.design.1)==1){
    opt.design=opt.design.1
  }

  if (nrow(opt.design.1)>1){
    n.min.cost=nrow(opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),])
    if (n.min.cost==1){
      opt.design=opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),]
    }else{
      opt.min.cost=opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),]

      ### only applies to num.miss=1
      opt.design=opt.min.cost[opt.min.cost[,"miss.loc"]==max(opt.min.cost[,"miss.loc"]),]
    }
  }

  op=opt.design[,"sim.seq"]

  if (pd==0){
    opt.pattern=rbind(all.pattern[op,],completers)
    opt.probs=c(rep(round((1-pc)/(nrow(patmx)-1),6),nrow(patmx)-1),pc)
    #opt.probs=c(A,B)
  }
  if (pd!=0){
    opt.pattern=rbind(all.pattern[op,],completers,dropper)
    opt.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-2),6),nrow(patmx)-2),pc,pd)
  }
  colnames(opt.pattern)=VNAMES


  re.ob=list("results"=sim.results.out,"opt.design"=opt.design,"opt.pattern"=opt.pattern,"opt.probs"=opt.probs,"design.order"=op)

  return(re.ob)

}


