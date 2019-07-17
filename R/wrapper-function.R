#' The wrapper function for all the searching strategies


#' @param popModel The data generation model (population model) specified using lavaan script
#' @param analyzeModel The analysis model,specified using lavaan script. The analysis model can be different from the population model. 
#' @param design0.out Mplus output file which contains the a priori power analysis/sample size planning (simulation) results for this
#' specific model assuming a complete data design. Theoretically, such analysis was supposed to be conducted before the study began.
#' @param VarNAMES A vector containing the names of the observed variables. The variable names must be ordered chronologically, by the time (wave) they are measured.
#' @param Time The total number of time points (or waves of data collection).
#' @param Time.complete Number of waves of data collection that have been completed before the funding cut occurs.
#' @param k The number of observed variables collected at each wave.
#' @param pc Proportion of completers. The proportion of subjects that will participate in all of the following waves of data collection and provide complete data. This must be greater than 0.
#' @param pd The proportion of subjects that will not participate in any of the following waves of data collection (i.e., drop from the longitudinal study). This can be 0.
#' @param costmx  A vector containing the unit cost of each observed variable that is yet to be measured (post the funding cut). The cost is assumed to be constant across subjects, but it is allowed to vary across variables and across waves.
#' @param n The total sample size as initially planned.
#' @param nreps Number of replications for Monte Carlo simulations.
#' @param focal.param The parameters of focal interest. If engine="l", the focal parameters should be specified using
#' the lavaan script. If engine="m", the focal parameters should be specified based on the Mplus output file design0.out.
#' @param complete.wave Specify the wave(s) if there are any waves that need to have complete data collected across all the participants.
#' @param complete.var Specify the names of the variable(s) if there are any variable(s) that need to have complete data collected across all the participating subjects.
#' @param max.mk Specify the maximum number of unique missing data patterns in the selected design. Only applicable if forward selection is used. 
#' @param eval.budget Logical, indicating whether there is any budget constraint. If the user wishes to search for PM designs under the budget limit, they need to specify the amount of the remaining available budget that can be used for future data collection.
#' @param rm.budget The amount of remaining budget avaialbe for future data collection.
#' @param distal.var Specify the names of the variables, if there are any time-independent distal variables included in the model that are not subject planned missingness.
#' @param seed seed for random number generation.
#' @param engine Specify the whether the simulations should be conducted using lavaan/simsem (engine="l") or Mplus (engine="m").
#' @param methods Specify which searching strategy should be used ("wave","indicator","forward").
#' 
#' @return An object containing the information of the optimal PM design, with highest power for testing the
#' focal parameters, compared with other PM designs
#' 
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
