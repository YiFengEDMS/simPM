#' Search for an optimal PM design.
#' 
#' \code{simPM} runs Monte Carlo simulations and returns the search results for optimal PM designs. This is a wrapper function for all the available searching methods. 
#' 
#' @inheritParams balance.miss.l
#' @param design0.out An object returned by \code{\link[MplusAutomation]{readModels}}. 
#' To obtain this object, the user need to have a Mplus 
#'    output file which contains the \emph{a priori} power analysis 
#'    results for this specific model assuming a complete data design 
#'    (i.e., simulation-based power analysis for sample size planning). 
#'    In principle, \emph{a priori} power analysis is supposed to be 
#'    conducted before the study began.
#' @param VarNAMES A character vector containing the names of the observed
#'    variables. The variable names must be ordered chronologically, by
#'    the time (wave) they are measured.
#' @param focal.param Character Vector. Specify the parameters of focal
#'    interest. If engine="l", the focal parameters should be specified 
#'    following the format of \code{lavaan} script. If engine="m", the
#'    focal parameters should be specified in the specific format based on
#'    the Mplus output object \code{design0.out}.
#' @param complete.wave Numeric vector. Specify which wave(s) that the
#'    user wish to have complete data collected from all the participants.
#'    Only applicable for wave-level PM designs.
#' @param max.mk Specify the maximum number of unique missing data
#'    patterns in the selected design. Only applicable if forward assembly
#'    is used. 
#' @param engine Specify the whether the simulations should be conducted
#'    using \code{lavaan}/\code{simsem} (engine="l") or M\emph{plus}
#'   (engine="m").
#' @param methods Specify which searching strategy should be used
#'   (wave-level PM designs only: methods = "wave"; balanced item-level
#'   PM designs only: methods = "indicator"; item-level PM designs via
#'   forward assembly: methods = "forward").
#' 
#' @return An object containing the information of the optimal PM design.
#'    The optimal design is the one that yields highest statistical power
#'    for testing the focal parameters, compared to other plausible
#'    candidate PM designs.
#' 
#' @import MplusAutomation
#' @import simsem
#' @import lavaan
#' @export simPM
#' @examples


simPM <- function(
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
  ) {         #which searching method to use, wave missing, balanced missing, or forward selection


  if (engine=="l"&methods=="wave") {
    output.w <- wave.miss.l(popModel = popModel,
                         analyzeModel = analyzeModel,
                         NAMES = VarNAMES,
                         Time = Time,
                         Time.complete = Time.complete,
                         k = k,
                         pc = pc,
                         pd = pd,
                         costmx = costmx,
                         n = n,
                         nreps = nreps,
                         focal.param = focal.param,
                         complete.wave = complete.wave,
                         eval.budget = eval.budget,
                         rm.budget = rm.budget,
                         distal.var = distal.var,
                         seed = seed
                         #multicore=multicore
                         )
  } else if (engine=="l" & methods=="item") {
    output.w <- balance.miss.l(
      popModel = popModel,
      analyzeModel = analyzeModel,
      NAMES = VarNAMES,
      Time = Time,
      Time.complete = Time.complete,
      k = k,
      pc = pc,
      pd = pd,
      costmx = costmx,
      n = n,
      nreps = nreps,
      focal.param = focal.param,
      complete.var = complete.var,
      eval.budget = eval.budget,
      rm.budget = rm.budget,
      distal.var = distal.var,
      seed = seed
      #multicore=multicore
      )
  }

  else if (engine=="l" & methods=="forward") {
    output.w <- forward.opt.simsem(
      popModel = popModel,
      analyzeModel = analyzeModel,
      NAMES = VarNAMES,
      distal.var = distal.var,
      n = n,
      nreps = nreps,
      seed = seed,
      Time = Time,
      k = k,
      Time.complete = Time.complete,
      costmx = costmx,
      pc = pc,
      pd = pd,
      focal.param = focal.param,
      max.mk = max.mk,                 # maximum number of missing slots allowed
      eval.budget = eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget = rm.budget,          # remaining available budget
      complete.var = complete.var
      #multicore=multicore
      )
  } else if (engine=="m" & methods=="wave") {
    output.w = wave.miss(
      VNAMES = VarNAMES,
      distal.var = distal.var,
      n = n,
      nreps = nreps,
      seed = seed,
      Time = Time,
      k = k,
      Time.complete = Time.complete,
      costmx = costmx,
      pc = pc,
      pd = pd,
      design0.out = design0.out,
      focal.param = focal.param,
      eval.budget = eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget = rm.budget,          # remaining available budget
      complete.wave = complete.wave)
    
  } else if (engine=="m"&methods=="item") {
    output.w = balance.miss(
      VNAMES = VarNames,
      distal.var = distal.var,
      n = n,
      nreps = nreps,
      seed = seed,
      Time = Time,
      k = k,
      Time.complete = Time.complete,
      costmx = costmx,
      pc = pc,
      pd = pd,
      design0.out = design0.out,
      focal.param = focal.param,
      eval.budget = eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget = rm.budget,          # remaining available budget
      complete.var = complete.var)
  
    } else if (engine=="m" & methods=="forward") {
    output.w = forward.opt(
      VNAMES = VarNAMES,
      distal.var = distal.var,
      n = n,
      nreps = nreps,
      seed = seed,
      Time = Time,
      k = k,
      Time.complete = Time.complete,
      costmx = costmx,
      pc = pc,
      pd = pd,
      design0.out = design0.out,
      focal.param = focal.param,
      max.mk = max.mk,                 # maximum number of missing slots allowed
      eval.budget = eval.budget,          # logical, whether the user would like to evaluate the budget constraints. If =T, the function will stop with a warning if all possible patterns would exceed the avaialbe remaining budget.
      rm.budget = rm.budget,          # remaining available budget
      complete.var = complete.var)
  }

  return(output.w)
}


