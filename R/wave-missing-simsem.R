

#' Searching for wave missing designs, using simsem/lavaan

#' @param popModel The data generation model (population model) specified using lavaan script
#' @param analyzeModel The analysis model, which can be different from the population model, using lavaan script
#' @param NAMES The names of the observed variables, ordered by the time they are measured
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
#' @param complete.wave Specify if there are any waves that need to have complete data collected across all the participating subjects
#' @param eval.budget Budget constraints. If the researcher wishes to search for designs under the budget limit, you to provide the remaining available budget that can be used for data collection
#' @param rm.budget The amount of remaining budget avaialbe for data collection
#' @param distal.var Any distal variables included in the model that would have complete data
#' @param seed Random seet for simulation
#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export wave.miss.l
#' @examples

wave.miss.l=function(popModel,
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
                     complete.wave=NULL,
                     eval.budget=T,
                     rm.budget=NULL,
                     distal.var=NULL,
                     seed=1234
#                     multicore=T
){

n.miss.waves=1:(Time-Time.complete-1)  # possible number of waves missing
ms.range=c((Time.complete*k+1):(Time*k))

# storage bins
designs=list()
probs=list()
pat.n=list()
#num.miss.wave=c()
cost.design=c() # cost of each design


rs=1

for (i in n.miss.waves){               # loop over different # of missing waves (designs)

  if (Time.complete==0){
    mwave=combn(c(1:Time),i)
  }
  if (Time.complete>0){
    mwave=combn(c(1:Time)[-c(1:Time.complete)],i)
  }

  pattern=matrix(0,nrow=ncol(mwave),ncol=k*Time)  #pattern matrix

  for (j in 1:nrow(pattern)){
    for (m in 1:nrow(mwave))
      pattern[j,((mwave[m,j]-1)*k+1):(mwave[m,j]*k)]=1   #put missing in pattern matrix
  }

  completers=rep(0,ncol(pattern))       # completers pattern
  dropper=c(rep(0,Time.complete*k),rep(1,(Time-Time.complete)*k))   #droppers pattern

  ### users may wish to specify certain waves to have complete data

  if(is.null(complete.wave)==F){

    if (i==1){

      evalpattern=matrix(mwave%in%c(complete.wave),nrow=nrow(mwave),byrow=F)
      # keep the patterns
      keep=pattern[evalpattern==F,]

      if (is.null(dim(keep))==F){
        pattern=keep
      }
      if (is.null(dim(keep))==T){    # it may become a non-matrix object, if only one row
        pattern=t(as.matrix(keep))
      }
    }

    if (i>1){
      evalpattern=matrix(mwave%in%c(complete.wave),nrow=nrow(mwave),byrow=F)
      # keep the patterns
      keep=pattern[colSums(evalpattern)==0,]

      if (is.null(dim(keep))==F){
        pattern=keep
      }
      if (is.null(dim(keep))==T){    # it may become a non-matrix object, if only one row
        pattern=t(as.matrix(keep))
      }
    }
  }

  ### design matrix

  if(pd!=0){      # if there are droppers
    patmx=rbind(pattern,completers,dropper)   # missing patterns for Mplus later
  }
  if(pd==0){      # if there are no droppers
    patmx=rbind(pattern,completers)
  }

  designs[[rs]]=patmx

  #### pattern probs
  # p.probs

  if (pd==0){
    p.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-1),6), nrow(patmx)-1),pc)
  }
  if (pd!=0){
    p.probs=c(rep(round((1-pc-pd)/(nrow(patmx)-2),6), nrow(patmx)-2),pc,pd)
  }

  probs[[rs]]=p.probs

  # cost of each design
  cost.design[rs]=sum((1-patmx[,ms.range])%*%costmx*p.probs*n)

  # ns in each pattern
  pn=rep(0,length(p.probs))
  for (pp in 1:(length(pn)-1)){
    pn[pp]=floor(n*p.probs[pp])
  }
  pn[length(pn)]=n-sum(pn[1:(length(pn)-1)])

  pat.n[[rs]]=pn

  # cost of each design
  cost.design[rs]=sum((1-patmx[,ms.range])%*%costmx*pn)


  rs=rs+1
}

#### evaluate cost
# only simulate for those designs that are below the budget limit

if (eval.budget==T){
  if (sum(cost.design>rm.budget)==length(cost.design)){
    stop ("all wave missing designs cost more than the avaiable remaing budget. Try other designs.")
  }

  designs2=designs[cost.design<=rm.budget] #select the designs that are below the budget limit
  probs2=probs[cost.design<=rm.budget]
  miss.waves=n.miss.waves[cost.design<=rm.budget]
  cost.design2=cost.design[cost.design<=rm.budget]
  pat.n2=pat.n[cost.design<=rm.budget]
}

if (eval.budget==F){
  designs2=designs
  probs2=probs
  miss.waves=n.miss.waves
  cost.design2=cost.design
  pat.n2=pat.n
}

convergence.rate=c()   #convergence rate
weakest.param.name=c()
weakest.para.power=c()

template=list()
logical.Matrix=list()
sim.out=list()

for (d in 1:length(designs2)){

  patmx=designs2[[d]]
  p.probs=probs2[[d]]
  pn=pat.n2[[d]]
  VNAMES=NAMES

if (prod(pn>=1)!=0){
  ###distal variables
  if (is.null(distal.var)==F){
    dis.pat=matrix(0,nrow=nrow(patmx),ncol=length(distal.var))
    patmx=cbind(patmx,dis.pat)
    VNAMES=c(VNAMES,distal.var)
  }

#  FNAME=paste0("missing-waves-",miss.waves[d])  #file name


  logical.mx=matrix(0,nrow=n,ncol=ncol(patmx))
  #  data=data[sample(1:nrow(data)),]   #shuffle the rows
  logical.mx[1:pn[1],]=patmx[rep(1,pn[1]),]

  for (pi in 2:(length(pn)-1)){
    logical.mx[(sum(pn[1:(pi-1)])+1):sum(pn[1:pi]),]=patmx[rep(pi,pn[pi]),]
  }

  # rest are completers

  logical.Mx=logical.mx==1
  logical.Matrix[[d]]=logical.Mx

  #  miss.model=miss(logical=logical.Mx)
  #  data.miss=impose(miss.model,data)

  misstemplate <- miss(logical=logical.Mx, m=0)
  output <- simsem::sim(nreps, n=n, model=analyzeModel, generate=popModel,miss=misstemplate,
#                        multicore=multicore,
seed=seed+d)

  template[[d]]=misstemplate
  sim.out[[d]]=output

  sim.param=summaryParam(output)
  name.param=rownames(summaryParam(output))
  converged <- output@converged == 0

  if (sum(converged)==0){
    convergence.rate[d]=0
    weakest.param.name[d]="NA"
    weakest.para.power[d]=0
  }
  if (sum(converged)>0){
    f.param=sim.param[name.param%in%focal.param,]
    weakest.f.param=f.param[f.param$`Power (Not equal 0)`==min(f.param$`Power (Not equal 0)`),]
    if (nrow(weakest.f.param)>1){
      weakest.f.param=weakest.f.param[1,]      ########## may need to be changed later
    }
    convergence.rate[d]=sum(converged)/nreps    #converged number of simulations
    weakest.param.name[d]=rownames(weakest.f.param)
    weakest.para.power[d]=weakest.f.param$`Power (Not equal 0)`
  }
}

if (prod(pn>=1)==0){
  convergence.rate[d]=NA   #convergence rate
  weakest.param.name[d]=NA
  weakest.para.power[d]=-99

  template[d]=NA
  logical.Matrix[d]=NA
  sim.out[d]=NA
}

}


sim.results.out=cbind.data.frame(convergence.rate,   #convergence rate
                                 weakest.param.name,
                                 weakest.para.power,
                                 "cost.design"=cost.design2, # cost of each design
                                 miss.waves)

opt.design.1=sim.results.out[sim.results.out[,"weakest.para.power"]==max(sim.results.out[,"weakest.para.power"]),]

if (nrow(opt.design.1)==1){
  opt.design=opt.design.1
}
if (nrow(opt.design.1)>1){
  opt.design=opt.design.1[opt.design.1$cost.design==min(opt.design.1$cost.design),]
}

op=which(miss.waves==opt.design$miss.waves)    #which design is chosen
opt.pattern=designs2[[op]]
colnames(opt.pattern)=NAMES
opt.probs=probs2[[op]]
opt.patns=pat.n2[[op]]

opt.template=template[[op]]
opt.logical=logical.Matrix[[op]]
opt.output=sim.out[[op]]

re.ob=list("results"=sim.results.out,"opt.design"=opt.design,"opt.pattern"=opt.pattern,
           "opt.probs"=opt.probs,"opt.ns"=opt.patns,"n.miss.waves"=opt.design$miss.waves,
           "opt.template"=opt.template,'opt.logical'=opt.logical,'opt.output'=opt.output)

return(re.ob)
}


