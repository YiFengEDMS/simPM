#' The wrapper function for all the searching strategies


#' @param popModel The data generation model (population model) specified using lavaan script
#' @param analyzeModel The analysis model, which can be different from the population model, using lavaan script
#' @param design0.out Mplus output file which contains the a priori power analysis/sample size planning (simulation) results for this
#' specific model assuming a complete data design. Theoretically, such analysis was supposed to be conducted before the study began.
#' @param VarNAMES The names of the observed variables, ordered by the time they are measured
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
#' @param complete.var Specify if there are any variables that need to have complete data collected across all the participating subjects
#' @param max.mk If using forward selection, specify the maximum number of unique missing data patterns in the design
#' @param eval.budget Budget constraints. If the researcher wishes to search for designs under the budget limit, you to provide the remaining available budget that can be used for data collection
#' @param rm.budget The amount of remaining budget avaialbe for data collection
#' @param distal.var Any distal variables included in the model that would have complete data
#' @param seed Random seet for simulation
#' @param engine Specify the whether the simulation should be completed using lavaan/simsem (engine="l") or Mplus (engine="m")
#' @param methods Specify which searching method should be used ("wave","indicator","forward")
#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' @seealso \code{\link{simPM}} which wraps this function
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export simPM
#' @examples


simPM=function(
  popModel,               # lavaan and simsem specific
  analyzeModel,           # lavaan and simsem specific
  design0.out=NULL,            # Mplus specific
  VarNAMES,
  Time,
  Time.complete,
  k,
  pc,
  pd,
  costmx,
  n,
  nreps,
  focal.param,
  complete.wave=NULL,      #wave missing specific
  complete.var=NULL,       #balanced missing specific
  max.mk=NULL,                  #forward selection specific
  eval.budget=T,
  rm.budget=NULL,
  distal.var=NULL,
  seed=1234,
  engine="l",              #which engine to use, Mplus or lavaan
  methods="wave"
#  multicore=T
  ){         #which searching method to use, wave missing, balanced missing, or forward selection


  if (engine=="l"&methods=="wave"){
    output.w=wave.miss.l(popModel=popModel,
                         analyzeModel=analyzeModel,
                         NAMES=VarNAMES,
                         Time=Time,
                         Time.complete=Time.complete,
                         k=k,
                         pc=pc,
                         pd=pd,
                         costmx=costmx,
                         n=n,
                         nreps=nreps,
                         focal.param=focal.param,
                         complete.wave=complete.wave,
                         eval.budget=eval.budget,
                         rm.budget=rm.budget,
                         distal.var=distal.var,
                         seed=seed
                         #multicore=multicore
                         )
  }

  else if (engine=="l"&methods=="item"){
    output.w=balance.miss.l(
      popModel=popModel,
      analyzeModel=analyzeModel,
      NAMES=VarNAMES,
      Time=Time,
      Time.complete=Time.complete,
      k=k,
      pc=pc,
      pd=pd,
      costmx=costmx,
      n=n,
      nreps=nreps,
      focal.param=focal.param,
      complete.var=complete.var,
      eval.budget=eval.budget,
      rm.budget=rm.budget,
      distal.var=distal.var,
      seed=seed
      #multicore=multicore
      )
  }

  else if (engine=="l"&methods=="forward"){
    output.w=forward.opt.simsem(
      popModel=popModel,
      analyzeModel=analyzeModel,
      NAMES=VarNAMES,
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
      max.mk=max.mk,                 # maximum number of missing slots allowed
      eval.budget=eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget=rm.budget,          # remaining available budget
      complete.var=complete.var
      #multicore=multicore
      )
  }

  else if (engine=="m"&methods=="wave"){
    output.w=wave.miss(
      VNAMES=VarNAMES,
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
      design0.out=design0.out,
      focal.param=focal.param,
      eval.budget=eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget=rm.budget,          # remaining available budget
      complete.wave=complete.wave)
  }

  else if (engine=="m"&methods=="item"){
    output.w=balance.miss(
      VNAMES=VarNames,
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
      design0.out=design0.out,
      focal.param=focal.param,
      eval.budget=eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget=rm.budget,          # remaining available budget
      complete.var=complete.var)
  }

  else if (engine=="m"&methods=="forward"){
    output.w=forward.opt(
      VNAMES=VarNAMES,
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
      design0.out=design0.out,
      focal.param=focal.param,
      max.mk=max.mk,                 # maximum number of missing slots allowed
      eval.budget=eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget=rm.budget,          # remaining available budget
      complete.var=complete.var)
  }


  return(output.w)
}

#summary.opt(re.ob)
