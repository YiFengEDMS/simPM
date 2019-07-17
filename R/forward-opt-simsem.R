#' Function to search for the optimal PM design with forward selection (lavaan-based)

#' @param popModel The data generation model (population model) specified using lavaan script
#' @param analyzeModel The analysis model,specified using lavaan script. The analysis model can be different from the population model.
#' @param NAMES A vector containing the names of the observed variables. The variable names must be ordered chronologically, by the time (wave) they are measured.
#' @param Time The total number of time points (or waves of data collection).
#' @param Time.complete Number of waves of data collection that have been completed before the funding cut occurs.
#' @param k The number of observed variables collected at each wave.
#' @param pc Proportion of completers. The proportion of subjects that will participate in all of the following waves of data collection and provide complete data. This must be greater than 0.
#' @param pd The proportion of subjects that will not participate in any of the following waves of data collection (i.e., drop from the longitudinal study). This can be 0.
#' @param costmx  A vector containing the unit cost of each observed variable that is yet to be measured (post the funding cut). The cost is assumed to be constant across subjects, but it is allowed to vary across variables and across waves.
#' @param n The total sample size as initially planned.
#' @param nreps Number of replications for Monte Carlo simulations.
#' @param focal.param The parameters of focal interest. The focal parameters should be specified #' using the lavaan script.
#' @param complete.var Specify the names of the variable(s) if there are any variable(s) that need to have complete data collected across all the participating subjects.
#' @param max.mk Specify the maximum number of unique missing data patterns in the selected design. Only applicable if forward selection is used. 
#' @param eval.budget Logical, indicating whether there is any budget constraint. If the user wishes to search for PM designs under the budget limit, they need to specify the amount of the remaining available budget that can be used for future data collection.
#' @param rm.budget The amount of remaining budget avaialbe for future data collection.
#' @param distal.var Specify the names of the variables, if there are any time-independent distal variables included in the model that are not subject planned missingness.
#' @param seed Random seet for simulation
#' @return An object containing the information of the optimal item-level PM design, with highest power for testing the
#' focal parameters, compared with other candidate PM designs
#' 
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export forward.opt.simsem
#' @examples


forward.opt.simsem=function(
  popModel,
  analyzeModel,
  NAMES,
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
  focal.param,
  max.mk,                 # maximum number of missing slots allowed
  eval.budget=F,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
  rm.budget=NULL,          # remaining available budget
  complete.var=NULL
  #multicore=T
  )
{
  VNAMES=NAMES
  previous=list()

  if (max.mk>n-(pc+pd)*n){
    stop("Consider using a smaller value for max.mk")
  }

  if (max.mk==1){          # when num.miss=1
    output=opt1.simsem(popModel=popModel,
                    analyzeModel=analyzeModel,
                    NAMES=VNAMES,
                    distal.var=distal.var,
                    n=n,
                    nreps=nreps,
                    seed=seed,
                    Time=Time,
                    k=k,
                    Time.complete=Time.complete,
                    costmx=costmx,
                    pc=pc,
                    pd=pd,
                    focal.param=focal.param,
                    complete.var=complete.var
                    #multicore=multicore
                    )

    ### if user wish to evaluate the budget
    if (eval.budget==T){
      if ((1-prod(output$results$cost.design>rm.budget))==0){
        return(output)
        warning("all designs would cost more than the avaiable remaing budget. Consider increasing max.mk.")
      }
    }
    return(output)
  }else if (max.mk>=2){     # when max number of missing slots is greater than 1

    future.k=(Time-Time.complete)*k   # data points not yet completed, also maximum possible # of missing
    ms.range=c((Time.complete*k+1):(Time*k))  # The available time slots to plant missingness

    design.order=rep(NA,max.mk)    # keep record of the order of the selected designs

    #num.miss must be less than future.k
    if (max.mk>=future.k){
      stop ("number of missing points exceeds the maximum possible value.")
    }

    output1=opt1.simsem(popModel=popModel,
                        analyzeModel=analyzeModel,
                        NAMES=VNAMES,
                        distal.var=distal.var,
                        n=n,
                        nreps=nreps/5,
                        seed=seed,
                        Time=Time,
                        k=k,
                        Time.complete=Time.complete,
                        costmx=costmx,
                        pc=pc,
                        pd=pd,
                        focal.param=focal.param,
                        complete.var=complete.var
                        #multicore=multicore
                        )
    previous[[1]]=output1$opt.pattern
    design.order[1]=output1$design.order

    for (num.miss in 2:max.mk){

      if (eval.budget==T){            # evaluate the cost and the budget
        if (num.miss==max.mk){        # cannot evaluate the cost without knowing the previous selected patterns

          ms.combn=combn(ms.range,num.miss) #all possible combinations of missing slots given num.miss
          all.pattern=matrix(0,nrow=choose(future.k,num.miss),ncol=length(VNAMES)) # place holder, all possible patterns with a certain ms

          for (q in 1:nrow(all.pattern)){
            all.pattern[q,ms.combn[,q]]=1
          }                                         # add missingness

          if (length(complete.var)==1){
            # find which column these variables are
            complete.cols=which(VNAMES%in%complete.var)
            # whether to keep the rows
            keep=all.pattern[,complete.cols]==0

            temp.pattern=all.pattern[keep,]
            if (is.null(dim(temp.pattern))==F){
              all.pattern=temp.pattern
            }
            if (is.null(dim(temp.pattern))==T){    # it may become a non-matrix object, if only one row
              all.pattern=t(as.matrix(temp.pattern))
            }
          }

          if (length(complete.var)>1){
            complete.cols=which(VNAMES%in%complete.var)
            keep=rowSums(all.pattern[,complete.cols]==0)==length(complete.var)

            temp.pattern=all.pattern[keep,]
            if (is.null(dim(temp.pattern))==F){
              all.pattern=temp.pattern
            }
            if (is.null(dim(temp.pattern))==T){
              all.pattern=t(as.matrix(temp.pattern))       # need to be a matrix object even if it's only one row
            }
          }


          # the optimal pattern obtained in previous round of selection
          opt1.pattern=previous[[num.miss-1]]

          # unit cost of each pattern
          if (nrow(all.pattern)>1){
            unit.cost=rowSums((1-all.pattern[,ms.range])*costmx)*((1-pc-pd)*n/num.miss)
          }
          if (nrow(all.pattern)==1){
            temp=t(as.matrix(1-all.pattern[,ms.range]))
            unit.cost=rowSums(temp*costmx)*((1-pc-pd)*n/num.miss)
          }

          # cost of previous chosen design
          if (pd==0){
            opt1.cost=sum(rowSums((1-opt1.pattern[,ms.range])*costmx)*c(rep((1-pc)*n/num.miss,num.miss-1),pc*n))
          }
          if (pd!=0){
            opt1.cost=sum(rowSums((1-opt1.pattern[,ms.range])*costmx)*c(rep((1-pc-pd)*n/num.miss,num.miss-1),pc*n,pd*n))
          }

          sum.cost=unit.cost+opt1.cost
          if ((1-prod(sum.cost>rm.budget))==0){
            stop("all designs cost more than the avaiable remaing budget, please allow a a larger max.mk. To turn off the budget evaluation, set eval.budget=F and rm.budget=NULL")
          }
        }
      }


      ms.combn=combn(ms.range,num.miss) #all possible combinations of missing slots given num.miss
      all.pattern=matrix(0,nrow=choose(future.k,num.miss),ncol=length(VNAMES)) # place holder, all possible patterns with a certain ms

      # update the missing patterns with 1s
      for (q in 1:nrow(all.pattern)){
        all.pattern[q,ms.combn[,q]]=1
      }                                         # add missingness

      # the optimal pattern obtained in previous round of selection
      opt1.pattern=previous[[num.miss-1]]

      if (length(complete.var)==1){
        # find which column these variables are
        complete.cols=which(VNAMES%in%complete.var)
        # whether to keep the rows
        keep=all.pattern[,complete.cols]==0

        temp.pattern=all.pattern[keep,]
        if (is.null(dim(temp.pattern))==F){
          all.pattern=temp.pattern
        }
        if (is.null(dim(temp.pattern))==T){
          all.pattern=t(as.matrix(temp.pattern))
        }


        # take out these designs columns out of the ms.combn

        temp.combs=ms.combn[,keep]

        if(is.null(dim(temp.combs))==F){
          ms.combn=temp.combs
        }
        if(is.null(dim(temp.combs))==T){        #may become a non-matrix object if only one column
          ms.combn=as.matrix(temp.combs)
        }
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


        # take out these designs columns out of the ms.combn

        temp.combs=ms.combn[,keep]
        if(is.null(dim(temp.combs))==F){
          ms.combn=temp.combs
        }
        if(is.null(dim(temp.combs))==T){
          ms.combn=as.matrix(temp.combs)
        }
      }

      # storage bins for simulation results

      ct=ncol(ms.combn)

      convergence.rate=rep(NA,ct)   #convergence rate
      weakest.param.name=rep(NA,ct)
      weakest.para.power=rep(NA,ct)
      cost.design=rep(NA,ct) # cost of each design
      miss.num=rep(num.miss,ct)
      miss.name=matrix(NA,ct,num.miss)
      sim.seq=rep(NA,ct)
      miss.loc=matrix(NA,ct,num.miss)
      sim.out=vector("list",ct)   # simulation output storage


      for (i in 1:nrow(all.pattern)){

        VNAMES=NAMES
        ###distal variables
        if (is.null(distal.var)==F){
          dis.pat=rep(0,length(distal.var))
          patmx=rbind(c(all.pattern[i,],dis.pat),cbind(opt1.pattern,t(replicate(nrow(opt1.pattern),dis.pat))))
          VNAMES=c(VNAMES,distal.var)
        }
        if (is.null(distal.var)==T){
          patmx=rbind(all.pattern[i,],opt1.pattern)   # missing patterns
        }

#        FNAME=paste0("missing-",num.miss,"-sim.seq-",i)  #file name


        # pattern probs

        if (pd==0){
          p.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-1),6), nrow(patmx)-1),pc)
        }
        if (pd!=0){
          p.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-2),6), nrow(patmx)-2),pc,pd)
        }

        # ns in each pattern
        pn=rep(0,length(p.probs))
        for (pp in 1:(length(pn)-1)){
          pn[pp]=floor(n*p.probs[pp])
        }
        pn[length(pn)]=n-sum(pn[1:(length(pn)-1)])

        logical.mx=matrix(0,nrow=n,ncol=ncol(patmx))
        logical.mx[1:pn[1],]=patmx[rep(1,pn[1]),]

        if (length(pn)>=3){
          for (pi in 2:(length(pn)-1)){
            logical.mx[(sum(pn[1:(pi-1)])+1):sum(pn[1:pi]),]=patmx[rep(pi,pn[pi]),]
          }
        }
        # rest are completers

        logical.Mx=logical.mx==1

        misstemplate <- miss(logical=logical.Mx, m=0)

        if (num.miss==max.mk){
        output <- simsem::sim(nRep = nreps, n=n, model=analyzeModel, generate=popModel,miss=misstemplate,
                              #multicore=multicore,
                              seed=seed)
        }else if (num.miss<max.mk){
        output <- simsem::sim(nRep=nreps/5, n=n, model=analyzeModel, generate=popModel,miss=misstemplate,
                              #multicore=multicore,
                              seed=seed)
        }

        sim.out[[i]] <- output

        sim.param=summaryParam(output)
        name.param=rownames(summaryParam(output))
        converged <- output@converged == 0

        if (sum(converged)==0){
          convergence.rate[i]=0
          weakest.param.name[i]="NA"
          weakest.para.power[i]=0
        }
        if (sum(converged)>0){
          f.param=sim.param[name.param%in%focal.param,]
          weakest.f.param=f.param[f.param$`Power (Not equal 0)`==min(f.param$`Power (Not equal 0)`),]
          if (nrow(weakest.f.param)>1){
            weakest.f.param=weakest.f.param[1,]      ########## may need to be changed later
          }
          convergence.rate[i]=sum(converged)/nreps    #converged number of simulations
          weakest.param.name[i]=rownames(weakest.f.param)
          weakest.para.power[i]=weakest.f.param$`Power (Not equal 0)`
        }

        if (pd==0){
          cost.design[i]=sum(c(rep(n*(1-pc)/num.miss,num.miss),pc*n)*((1-patmx[,ms.range])%*%costmx))   #patmx depends on i
        }
        if (pd!=0){
          cost.design[i]=sum(c(rep(n*(1-pc-pd)/num.miss,num.miss),pc*n,pd*n)*((1-patmx[,ms.range])%*%costmx))   #patmx depends on i
        }

        miss.name[i,]=VNAMES[ms.combn[,i]]
        sim.seq[i]=i  # location as specified in the miss.combn matrix
        miss.loc[i,]=ms.combn[,i]
      }


      colnames(miss.name)=paste0("miss.var",1:num.miss)
      colnames(miss.loc)=paste0("miss.loc",1:num.miss)

      sim.results.out=cbind.data.frame(convergence.rate,   #convergence rate
                                       weakest.param.name,
                                       weakest.para.power,
                                       cost.design, # cost of each design
                                       miss.num,
                                       miss.name,
                                       sim.seq,
                                       miss.loc)
      if(num.miss<max.mk){
        opt.design.1=sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]),]
      } else if (num.miss==max.mk){

      if (eval.budget==F){
      opt.design.1=sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]),]
      }

      if (eval.budget==T){
        sim.results.out=sim.results.out[sim.results.out$cost.design<=rm.budget,]
        opt.design.1=sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]),]
      }
      }

      if (nrow(opt.design.1)==1){
        opt.design=opt.design.1
      }

      if (nrow(opt.design.1)>1){
        n.min.cost=nrow(opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),])
        if (n.min.cost==1){
          opt.design=opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),]
        }else{
          opt.min.cost=opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),]
          sum.loc=rowSums(opt.min.cost[,colnames(miss.loc)])
          opt.design=opt.min.cost[sum.loc==max(sum.loc),]
        }
      }

      op=opt.design[,"sim.seq"]

      opt.output=sim.out[[op]]     #save the simulation output

      opt.pattern=rbind(all.pattern[op,],opt1.pattern)
      colnames(opt.pattern)=NAMES

      if (pd==0){
        opt.probs=c(rep(round((1-pc)/(nrow(patmx)-1),6),num.miss),pc)
      }
      if (pd!=0){
        opt.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-2),6),num.miss),pc,pd)
      }

      # ns in each pattern
      pn=rep(0,length(opt.probs))
      for (pp in 1:(length(pn)-1)){
        pn[pp]=floor(n*opt.probs[pp])
      }
      pn[length(pn)]=n-sum(pn[1:(length(pn)-1)])


      design.order[num.miss]=op
      
      misc=list(time=Time,k=k,focal.param=focal.param)
      
      re.ob=list("results"=sim.results.out,"opt.design"=opt.design,"opt.pattern"=opt.pattern,"opt.probs"=opt.probs,
                 "opt.ns"=pn,
                 "design.order"=design.order,
                 "opt.output"=opt.output,
                 "misc"=misc)
      previous[[num.miss]]=re.ob$opt.pattern   # update the previous pattern list
    }
    return(re.ob)
  }
}



